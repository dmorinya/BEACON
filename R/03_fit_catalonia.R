helpers_path <- c("helpers.R", "R/helpers.R", "../R/helpers.R")
helpers_path <- helpers_path[file.exists(helpers_path)][1]
if (is.na(helpers_path)) stop("helpers.R not found.")
source(helpers_path)

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(posterior)
})

dir.create("results", recursive = TRUE, showWarnings = FALSE)

input <- "data/derived/catalonia_beacon_analysis.csv"
if (!file.exists(input)) stop("Run R/02_download_prepare_catalonia.R first.")
d <- fread(input)
d[, week_date := as.Date(week_date)]

max_areas <- as.integer(Sys.getenv("BEACON_MAX_AREAS", "60"))
area_size <- d[, .(population = median(population, na.rm = TRUE), n = .N),
               by = abs_code][order(-population)]
keep <- head(area_size$abs_code, max_areas)
d <- d[abs_code %in% keep & week_in_season >= 1 & week_in_season <= 52]
d <- d[
  is.finite(y) & is.finite(population) & population > 0 &
    is.finite(immunisation_intensity) & is.finite(treated) &
    is.finite(week_in_season) & is.finite(ses_z)
]
if (nrow(d) < 500) {
  stop("Too few complete records after filtering. Inspect immunisation-",
       "intensity and age-group preparation.")
}

d[, area_id := abs_code]
d[, season_id := season]
d[, ses := ses_z]
d[, immunisation_intensity := pmax(immunisation_intensity, 0)]

# ---------------------------------------------------------------------------
# Overlap first. If the identifying contrast is degenerate, the effect estimate
# is uninformative regardless of how well the sampler behaves.
# ---------------------------------------------------------------------------
overlap <- overlap_diagnostic(d, area_col = "abs_code", season_col = "season")
fwrite(overlap, "results/catalonia_overlap_diagnostic.csv")
message("Overlap in scaled immunisation intensity by programme season:")
print(overlap)

# ---------------------------------------------------------------------------
# Stan data and model
# ---------------------------------------------------------------------------
stan_H <- as.integer(Sys.getenv("BEACON_STAN_H", "4"))
stan_harmonic_decay <- as.numeric(Sys.getenv("BEACON_HARMONIC_DECAY", "1"))
prior_sd_beta <- as.numeric(Sys.getenv("BEACON_PRIOR_SD_BETA", "0.5"))

built <- build_beacon_stan_data(
  d,
  H = stan_H,
  harmonic_decay = stan_harmonic_decay,
  intensity_quadratic = as.integer(Sys.getenv("BEACON_INTENSITY_QUADRATIC", "0")),
  # Primary analysis: neutral prior. Centring on log(0.6) would put the prior
  # mean on the hypothesis under test.
  prior_mean_beta_intensity = 0,
  prior_sd_beta_intensity = prior_sd_beta,
  prior_sd_ascertainment = 0,
  compute_log_lik = 1L,
  compute_y_rep = 1L,
  dynamic_gformula = 1L
)
stan_data <- built$stan_data
# build_beacon_stan_data() reorders rows into contiguous series blocks; use its
# frame from here on so every downstream join matches the Stan row order.
d <- built$data

stan_file <- Sys.getenv("BEACON_STAN_FILE", "")
if (!nzchar(stan_file)) {
  candidates <- c("Stan/beacon_affine_v2.stan", "../Stan/beacon_affine_v2.stan",
                  "beacon_affine_v2.stan")
  stan_file <- candidates[file.exists(candidates)][1]
  if (is.na(stan_file)) stop("beacon_affine_v2.stan not found.")
}
stan_file <- normalizePath(stan_file, mustWork = TRUE)
assert_beacon_v2(stan_file)

message("Stan source: ", stan_file)
message("Stan source MD5: ", unname(tools::md5sum(stan_file)))
message("Series blocks: ", stan_data$M, "; weeks per season: ", stan_data$W)
message("Lag transform: (log1p(y) - ", round(stan_data$lag_center, 4), ") / ",
        round(stan_data$lag_scale, 4))

force_recompile <- tolower(Sys.getenv("BEACON_FORCE_RECOMPILE", "true")) %in%
  c("1", "true", "yes")
mod <- cmdstan_model(stan_file, force_recompile = force_recompile)

chains <- as.integer(Sys.getenv("BEACON_CHAINS", "4"))
parallel_chains <- min(
  chains,
  as.integer(Sys.getenv("BEACON_CORES", as.character(60)))
)

set.seed(17072026)

run_beacon <- function(stan_data, label, seed = 17072026) {
  init_list <- make_beacon_inits_v2(stan_data, chains, seed + 991L)

  init_list <- lapply(init_list, function(chain_inits) {
    chain_inits$reciprocal_phi <- 0.5
    return(chain_inits)
  })
  mod$sample(
    data = stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = as.integer(Sys.getenv("BEACON_WARMUP", "1000")),
    iter_sampling = as.integer(Sys.getenv("BEACON_SAMPLING", "1000")),
    init = init_list,
    adapt_delta = as.numeric(Sys.getenv("BEACON_ADAPT_DELTA", "0.97")),
    max_treedepth = as.integer(Sys.getenv("BEACON_MAX_TREEDEPTH", "13")),
    refresh = 100
  )
}

# ---------------------------------------------------------------------------
# Prior predictive check. Run before the posterior: if the prior implies weekly
# incidence orders of magnitude away from the observed range, fix that first.
# ---------------------------------------------------------------------------
if (tolower(Sys.getenv("BEACON_PRIOR_CHECK", "true")) %in% c("1", "true", "yes")) {
  prior_data <- stan_data
  prior_data$prior_only <- 1L
  prior_data$compute_y_rep <- 1L
  prior_data$compute_log_lik <- 0L
  prior_data$dynamic_gformula <- 0L

  message("Prior predictive check ...")
  prior_fit <- try(
    mod$sample(
      data = prior_data,
      seed = 17072026,
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      init = make_beacon_inits_v2(prior_data, 2, 4242L),
      refresh = 0
    ),
    silent = TRUE
  )
  if (inherits(prior_fit, "try-error")) {
    warning("Prior predictive run failed; inspect before trusting the posterior.")
  } else {
    yrep <- as_draws_matrix(prior_fit$draws("y_rep"))
    prior_totals <- rowSums(yrep)
    observed_total <- sum(stan_data$y)
    fwrite(
      data.table(
        quantity = c("observed_total", "prior_predictive_median",
                     "prior_predictive_lower", "prior_predictive_upper"),
        value = c(
          observed_total,
          median(prior_totals),
          unname(quantile(prior_totals, 0.025)),
          unname(quantile(prior_totals, 0.975))
        )
      ),
      "results/catalonia_prior_predictive_check.csv"
    )
    if (observed_total < quantile(prior_totals, 0.025) ||
        observed_total > quantile(prior_totals, 0.975)) {
      warning(
        "The observed total lies outside the central 95% of the prior ",
        "predictive distribution. Revisit the prior on alpha and sigma_f ",
        "before interpreting the posterior."
      )
    }
    rm(prior_fit, yrep)
  }
}

# ---------------------------------------------------------------------------
# Primary fit
# ---------------------------------------------------------------------------
fit <- run_beacon(stan_data, "primary")
fit$save_object("results/catalonia_beacon_fit.rds")

parameter_variables <- c(
  "beta_intensity", "beta_intensity_cos", "beta_intensity_sin",
  "beta_lag", "beta_ses",
  "sigma_area", "sigma_season", "sigma_series",
  "sigma_shift", "sigma_speed", "sigma_f", "phi"
)
pars <- fit$summary(parameter_variables)
fwrite(as.data.table(pars), "results/catalonia_parameter_summary.csv")

scalar_draws <- as_draws_matrix(fit$draws(c(
  "beta_intensity", "beta_intensity_cos", "beta_intensity_sin",
  "rr_reference",
  "prevented_total", "prevented_fraction",
  "prevented_total_dynamic", "prevented_fraction_dynamic",
  "sigma_shift", "sigma_speed", "sigma_area", "sigma_season", "sigma_series",
  "phi"
)))

beta <- scalar_draws[, "beta_intensity"]
modulation <- sqrt(
  scalar_draws[, "beta_intensity_cos"]^2 +
    scalar_draws[, "beta_intensity_sin"]^2
)

summ <- function(x) {
  q <- q_summary(x)
  data.table(mean = mean(x, na.rm = TRUE), lower = q[1], median = q[2],
             upper = q[3])
}

metric_rows <- list(
  cbind(metric = "intensity_log_rate_ratio_cycle_average", summ(beta)),
  cbind(metric = "rate_ratio_at_reference_intensity",
        summ(scalar_draws[, "rr_reference"])),
  cbind(metric = "effect_modulation_amplitude", summ(modulation)),
  cbind(metric = "prevented_cases", summ(scalar_draws[, "prevented_total"])),
  cbind(metric = "prevented_fraction",
        summ(scalar_draws[, "prevented_fraction"])),
  cbind(metric = "prevented_cases_dynamic",
        summ(scalar_draws[, "prevented_total_dynamic"])),
  cbind(metric = "prevented_fraction_dynamic",
        summ(scalar_draws[, "prevented_fraction_dynamic"])),
  cbind(metric = "epidemic_shift_sd", summ(scalar_draws[, "sigma_shift"])),
  cbind(metric = "epidemic_speed_sd", summ(scalar_draws[, "sigma_speed"])),
  cbind(metric = "area_sd", summ(scalar_draws[, "sigma_area"])),
  cbind(metric = "season_sd", summ(scalar_draws[, "sigma_season"])),
  cbind(metric = "area_season_sd", summ(scalar_draws[, "sigma_series"])),
  cbind(metric = "dispersion", summ(scalar_draws[, "phi"]))
)

real_summary <- rbindlist(metric_rows, fill = TRUE)
real_summary[, status := "fitted_open_catalonia_data"]
fwrite(real_summary, "results/catalonia_posterior_summary.csv")

# How much did the data move the intervention coefficient away from its prior?
contraction <- prior_posterior_contraction(beta, prior_sd = prior_sd_beta)
fwrite(contraction, "results/catalonia_prior_posterior_contraction.csv")
message("Prior-posterior contraction for beta_intensity: ",
        round(contraction$contraction, 3))
if (contraction$contraction < 0.25) {
  warning(
    "The posterior SD for beta_intensity is within 25% of its prior SD. The ",
    "data carry little information about the intervention effect; report this ",
    "as weak identification rather than as a null result."
  )
}

# ---------------------------------------------------------------------------
# Per-season estimands
# ---------------------------------------------------------------------------
# build_beacon_stan_data() indexes seasons by match(season_id, sort(unique(.))),
# so season index j corresponds to the j-th sorted season label.
season_levels <- sort(unique(d$season_id))
season_map <- data.table(season_index = seq_along(season_levels),
                         season = season_levels)

season_estimand <- function(varname, metric) {
  m <- as_draws_matrix(fit$draws(varname))
  rbindlist(lapply(seq_len(ncol(m)), function(j) {
    cbind(metric = metric, season_index = j, summ(m[, j]))
  }))
}

season_summary <- rbindlist(list(
  season_estimand("prevented_by_season", "prevented_cases"),
  season_estimand("peak_attenuation", "peak_attenuation"),
  season_estimand("peak_shift_weeks", "peak_shift_weeks"),
  season_estimand("epidemic_peak_week", "epidemic_peak_week")
), fill = TRUE)
season_summary <- merge(season_summary, season_map, by = "season_index",
                        all.x = TRUE)
fwrite(season_summary, "results/catalonia_season_estimands.csv")

# ---------------------------------------------------------------------------
# Weekly national trajectory, from the season-by-week aggregates.
# mu_fit_sw and mu_cf_sw are S x W matrices summed over areas inside Stan.
# ---------------------------------------------------------------------------
mu_fit <- as_draws_matrix(fit$draws("mu_fit_sw"))
mu_cf <- as_draws_matrix(fit$draws("mu_cf_sw"))

# cmdstanr names matrix columns "mu_fit_sw[s,w]"; parse the indices rather than
# assuming a storage order.
parse_index <- function(nms) {
  inside <- sub("^[^\\[]+\\[(.*)\\]$", "\\1", nms)
  parts <- do.call(rbind, strsplit(inside, ",", fixed = TRUE))
  data.table(col = seq_along(nms),
             season_index = as.integer(parts[, 1]),
             week_in_season = as.integer(parts[, 2]))
}
idx <- parse_index(colnames(mu_fit))
idx <- merge(idx, season_map, by = "season_index", all.x = TRUE)

observed_weekly <- d[, .(observed = sum(y), population = sum(population),
                         treated = max(treated),
                         week_date = min(week_date)),
                     by = .(season = season_id, week_in_season)]

traj <- rbindlist(lapply(seq_len(nrow(idx)), function(k) {
  j <- idx$col[k]
  fitted_week <- mu_fit[, j]
  counterfactual_week <- mu_cf[, j]
  data.table(
    season = idx$season[k],
    week_in_season = idx$week_in_season[k],
    fitted_mean = mean(fitted_week),
    fitted_lower = unname(quantile(fitted_week, .025)),
    fitted_upper = unname(quantile(fitted_week, .975)),
    counterfactual_mean = mean(counterfactual_week),
    counterfactual_lower = unname(quantile(counterfactual_week, .025)),
    counterfactual_upper = unname(quantile(counterfactual_week, .975))
  )
}))

traj <- merge(traj, observed_weekly, by = c("season", "week_in_season"),
              all.x = TRUE)
traj <- traj[!is.na(week_date)]
setorder(traj, week_date)
fwrite(traj, "results/catalonia_weekly_trajectory.csv")

# ---------------------------------------------------------------------------
# Declared sensitivity analyses.
#
#   1. Protective prior on the intervention coefficient, N(log 0.6, 0.5^2).
#   2. Differential-ascertainment bias analysis. lambda is collinear with
#      beta_intensity in the likelihood, so the data cannot separate them; the
#      point is that prior uncertainty about ascertainment is added to the
#      posterior for the effect rather than assumed away.
#   3. Two harmonics with no decay, reproducing the v1 epidemic curve.
# ---------------------------------------------------------------------------
run_sensitivity <- tolower(Sys.getenv("BEACON_SENSITIVITY", "true")) %in%
  c("1", "true", "yes")

if (run_sensitivity) {
  sens_specs <- list(
    list(label = "protective_prior",
         mutate = function(sd) { sd$prior_mean_beta_intensity <- log(0.6); sd }),
    list(label = "ascertainment_sd_0.10",
         mutate = function(sd) { sd$prior_sd_ascertainment <- 0.10; sd }),
    list(label = "two_harmonics",
         mutate = function(sd) { sd$H <- 2L; sd$harmonic_decay <- 0; sd })
  )

  sens_rows <- list()
  for (spec in sens_specs) {
    message("Sensitivity analysis: ", spec$label)
    sd_i <- spec$mutate(stan_data)
    sd_i$compute_log_lik <- 0L
    sd_i$compute_y_rep <- 0L

    f_i <- try(run_beacon(sd_i, spec$label, seed = 17072026), silent = TRUE)
    if (inherits(f_i, "try-error")) {
      warning("Sensitivity fit failed: ", spec$label)
      next
    }
    dr <- as_draws_matrix(f_i$draws(c("rr_reference", "prevented_fraction",
                                      "prevented_fraction_dynamic")))
    sens_rows[[length(sens_rows) + 1L]] <- rbindlist(list(
      cbind(analysis = spec$label, metric = "rate_ratio_at_reference_intensity",
            summ(dr[, "rr_reference"])),
      cbind(analysis = spec$label, metric = "prevented_fraction",
            summ(dr[, "prevented_fraction"])),
      cbind(analysis = spec$label, metric = "prevented_fraction_dynamic",
            summ(dr[, "prevented_fraction_dynamic"]))
    ))
    rm(f_i)
  }

  if (length(sens_rows)) {
    sens <- rbindlist(sens_rows, fill = TRUE)
    primary <- rbindlist(list(
      cbind(analysis = "primary", metric = "rate_ratio_at_reference_intensity",
            summ(scalar_draws[, "rr_reference"])),
      cbind(analysis = "primary", metric = "prevented_fraction",
            summ(scalar_draws[, "prevented_fraction"])),
      cbind(analysis = "primary", metric = "prevented_fraction_dynamic",
            summ(scalar_draws[, "prevented_fraction_dynamic"]))
    ))
    fwrite(rbindlist(list(primary, sens), fill = TRUE),
           "results/catalonia_sensitivity_summary.csv")
  }
}

# Diagnostics
diag <- fit$summary()[, c("variable", "rhat", "ess_bulk", "ess_tail")]
fwrite(as.data.table(diag), "results/catalonia_diagnostics.csv")

sampler_diag <- fit$diagnostic_summary()
diag_value <- function(x, f = sum) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(f(x, na.rm = TRUE))
}
fwrite(
  data.table(
    quantity = c("divergences", "max_treedepth_hits", "min_ebfmi"),
    value = c(
      diag_value(sampler_diag$num_divergent),
      diag_value(sampler_diag$num_max_treedepth),
      diag_value(sampler_diag$ebfmi, min)
    )
  ),
  "results/catalonia_sampler_diagnostics.csv"
)

message("Empirical posterior summaries written to results/.")

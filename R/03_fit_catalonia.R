source("helpers.R")

library(data.table)
library(cmdstanr)
library(posterior)

dir.create("results", recursive = TRUE, showWarnings = FALSE)

env_bool <- function(name, default) {
  tolower(Sys.getenv(name, if (default) "true" else "false")) %in%
    c("1", "true", "yes")
}

input <- "data/derived/catalonia_beacon_analysis.csv"
if (!file.exists(input)) stop("Run R/02_download_prepare_catalonia.R first.")
d_all <- fread(input)
d_all[, week_date := as.Date(week_date)]

## ---------------------------------------------------------------------------
## Analysis population
## ---------------------------------------------------------------------------
max_areas <- as.integer(Sys.getenv("BEACON_MAX_AREAS", "60"))
area_size <- d_all[, .(population = median(population, na.rm = TRUE), n = .N),
                   by = abs_code][order(-population)]
keep <- head(area_size$abs_code, max_areas)
d_all <- d_all[abs_code %in% keep & week_in_season >= 1 & week_in_season <= 52]
d_all <- d_all[
  is.finite(y) & is.finite(population) & population > 0 &
    is.finite(immunisation_intensity) & is.finite(treated) &
    is.finite(week_in_season) & is.finite(ses_z)
]

## Seasons with fewer than 52 distinct weeks in the extract are incomplete
## and are not used in the model.
season_weeks <- d_all[, .(n_weeks = uniqueN(week_in_season)), by = season]
incomplete_seasons <- season_weeks[n_weeks < 52, season]
if (length(incomplete_seasons) &&
    !env_bool("BEACON_KEEP_INCOMPLETE_SEASONS", FALSE)) {
  message("Dropping incomplete seasons: ", paste(incomplete_seasons, collapse = ", "))
  d_all <- d_all[!season %in% incomplete_seasons]
}

## The 2020-2021 season is excluded from the primary fit and
## included in a sensitivity fit.
excluded_seasons <- trimws(strsplit(
  Sys.getenv("BEACON_EXCLUDE_SEASONS", "2020-2021"), ",", fixed = TRUE)[[1]])
excluded_seasons <- excluded_seasons[nzchar(excluded_seasons)]

prepare_model_data <- function(d_in, exclude) {
  d <- copy(d_in)
  if (length(exclude)) d <- d[!season %in% exclude]
  if (nrow(d) < 500) {
    stop("Too few complete records after filtering. Inspect immunisation-",
         "activity and age-group preparation.")
  }
  d[, area_id := abs_code]
  d[, season_id := season]
  d[, ses := ses_z]
  d[, immunisation_intensity := pmax(immunisation_intensity, 0)]
  d
}

d <- prepare_model_data(d_all, excluded_seasons)
fwrite(
  data.table(
    item = c("incomplete_seasons_dropped", "seasons_excluded", "seasons_used"),
    value = c(paste(incomplete_seasons, collapse = ";"),
              paste(excluded_seasons, collapse = ";"),
              paste(sort(unique(d$season)), collapse = ";"))
  ),
  "results/catalonia_model_seasons.csv"
)

## ---------------------------------------------------------------------------
## Overlap first. If the cross-sectional contrast is degenerate, the effect is
## informed mainly by the pre-post component (see the unpooled-season fit).
## ---------------------------------------------------------------------------
overlap <- overlap_diagnostic(d, area_col = "abs_code", season_col = "season")
fwrite(overlap, "results/catalonia_overlap_diagnostic.csv")
message("Overlap in scaled programme activity by programme season:")
print(overlap)

## ---------------------------------------------------------------------------
## Stan data and model
## ---------------------------------------------------------------------------
stan_H <- as.integer(Sys.getenv("BEACON_STAN_H", "4"))
stan_harmonic_decay <- as.numeric(Sys.getenv("BEACON_HARMONIC_DECAY", "1"))
prior_sd_beta <- as.numeric(Sys.getenv("BEACON_PRIOR_SD_BETA", "0.5"))
ascertainment_sd_sensitivity <- as.numeric(
  Sys.getenv("BEACON_ASCERTAINMENT_SD", "0.10"))
unpooled_season_sd <- as.numeric(Sys.getenv("BEACON_UNPOOLED_SEASON_SD", "5"))

build_model_data <- function(d_in, compute_log_lik = 1L, compute_y_rep = 1L) {
  build_beacon_stan_data(
    d_in,
    H = stan_H,
    harmonic_decay = stan_harmonic_decay,
    intensity_quadratic = as.integer(Sys.getenv("BEACON_INTENSITY_QUADRATIC", "0")),
    prior_mean_beta_intensity = 0,
    prior_sd_beta_intensity = prior_sd_beta,
    prior_sd_ascertainment = 0,
    compute_log_lik = compute_log_lik,
    compute_y_rep = compute_y_rep,
    dynamic_gformula = 1L
  )
}

built <- build_model_data(d)
stan_data <- built$stan_data
d <- built$data

stan_file <- Sys.getenv("BEACON_STAN_FILE", "")
stan_file <- normalizePath("Stan/beacon_affine_v2.stan", mustWork = TRUE)
assert_beacon_v2(stan_file)
stan_data$unpooled_season_sd <- unpooled_season_sd

message("Stan source: ", stan_file)
message("Stan source MD5: ", unname(tools::md5sum(stan_file)))
message("Series blocks: ", stan_data$M, "; weeks per season: ", stan_data$W)
message("Lag transform: (log1p(y) - ", round(stan_data$lag_center, 4), ") / ",
        round(stan_data$lag_scale, 4))

force_recompile <- env_bool("BEACON_FORCE_RECOMPILE", TRUE)
mod <- cmdstan_model(stan_file, force_recompile = force_recompile)

chains <- as.integer(Sys.getenv("BEACON_CHAINS", "4"))
parallel_chains <- min(
  chains,
  as.integer(Sys.getenv("BEACON_CORES", as.character(60)))
)
adapt_delta <- as.numeric(Sys.getenv("BEACON_ADAPT_DELTA", "0.97"))

set.seed(17072026)

run_beacon <- function(stan_data, label, seed = 17072026) {
  message("Fitting: ", label)
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
    adapt_delta = adapt_delta,
    max_treedepth = as.integer(Sys.getenv("BEACON_MAX_TREEDEPTH", "13")),
    refresh = 100
  )
}

summ <- function(x) {
  q <- q_summary(x)
  data.table(mean = mean(x, na.rm = TRUE), lower = q[1], median = q[2],
             upper = q[3])
}

sampler_diagnostics <- function(fit) {
  sd_ <- fit$diagnostic_summary()
  val <- function(x, f = sum) {
    if (is.null(x) || !length(x) || all(is.na(x))) return(NA_real_)
    as.numeric(f(x, na.rm = TRUE))
  }
  data.table(
    divergences = val(sd_$num_divergent),
    max_treedepth_hits = val(sd_$num_max_treedepth),
    min_ebfmi = val(sd_$ebfmi, min)
  )
}

summarise_fit <- function(fit, stan_data, d_fit, label, prior_sd = prior_sd_beta) {
  dr <- as_draws_matrix(fit$draws(c(
    "beta_intensity", "beta_intensity_cos", "beta_intensity_sin", "beta_lag",
    "rr_reference", "prevented_total", "prevented_fraction",
    "prevented_total_dynamic", "prevented_fraction_dynamic",
    "sigma_shift", "sigma_speed", "sigma_season", "sigma_series", "phi"
  )))
  modulation <- sqrt(dr[, "beta_intensity_cos"]^2 + dr[, "beta_intensity_sin"]^2)

  centiles <- activity_centiles_by_season(
    data.table(area = stan_data$area, season = stan_data$season,
               treated = stan_data$treated,
               immunisation_intensity = stan_data$immunisation_intensity))
  ws <- within_support_pf(fit, stan_data, centiles)

  contraction <- prior_posterior_contraction(dr[, "beta_intensity"],
                                             prior_sd = prior_sd)

  rows <- rbindlist(list(
    cbind(metric = "log_rate_ratio_cycle_average", summ(dr[, "beta_intensity"])),
    cbind(metric = "rate_ratio_at_reference_activity", summ(dr[, "rr_reference"])),
    cbind(metric = "effect_modulation_amplitude", summ(modulation)),
    cbind(metric = "prevented_cases", summ(dr[, "prevented_total"])),
    cbind(metric = "prevented_fraction", summ(dr[, "prevented_fraction"])),
    cbind(metric = "prevented_cases_dynamic", summ(dr[, "prevented_total_dynamic"])),
    cbind(metric = "prevented_fraction_dynamic",
          summ(dr[, "prevented_fraction_dynamic"])),
    cbind(metric = "prevented_fraction_within_support_p90_vs_p10", ws$overall),
    cbind(metric = "beta_lag", summ(dr[, "beta_lag"])),
    cbind(metric = "exp_beta_lag", summ(exp(dr[, "beta_lag"]))),
    cbind(metric = "epidemic_shift_sd", summ(dr[, "sigma_shift"])),
    cbind(metric = "epidemic_speed_sd", summ(dr[, "sigma_speed"])),
    cbind(metric = "season_sd", summ(dr[, "sigma_season"])),
    cbind(metric = "area_season_sd", summ(dr[, "sigma_series"])),
    cbind(metric = "dispersion", summ(dr[, "phi"])),
    data.table(metric = "beta_prior_posterior_contraction",
               mean = contraction$contraction, lower = NA_real_,
               median = NA_real_, upper = NA_real_)
  ), fill = TRUE)

  season_levels <- sort(unique(d_fit$season_id))
  ws_by_season <- ws$by_season
  ws_by_season[, season := season_levels[season_index]]

  diag <- sampler_diagnostics(fit)
  list(
    metrics = cbind(analysis = label, rows),
    within_support_by_season = cbind(analysis = label, ws_by_season),
    diagnostics = cbind(analysis = label, diag)
  )
}

if (env_bool("BEACON_PRIOR_CHECK", TRUE)) {
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

fit <- run_beacon(stan_data, "primary")
fit$save_object("results/catalonia_beacon_fit.rds")

season_levels <- sort(unique(d$season_id))
season_map <- data.table(season_index = seq_along(season_levels),
                         season = season_levels)

parameter_variables <- c(
  "beta_intensity", "beta_intensity_cos", "beta_intensity_sin",
  "beta_lag", "beta_ses",
  "sigma_area", "sigma_season", "sigma_series",
  "sigma_shift", "sigma_speed", "sigma_f", "phi"
)
pars <- fit$summary(parameter_variables)
fwrite(as.data.table(pars), "results/catalonia_parameter_summary.csv")

primary <- summarise_fit(fit, stan_data, d, "primary")
real_summary <- copy(primary$metrics)
real_summary[, status := "fitted_open_catalonia_data"]
fwrite(real_summary, "results/catalonia_posterior_summary.csv")
fwrite(primary$within_support_by_season,
       "results/catalonia_within_support_pf_by_season.csv")

contraction <- real_summary[metric == "beta_prior_posterior_contraction", mean]
message("Prior-posterior contraction for beta_intensity: ", round(contraction, 3))
if (contraction < 0.25) {
  warning(
    "The posterior SD for beta_intensity is within 25% of its prior SD. The ",
    "data carry little information about the intervention effect; report this ",
    "as weak identification rather than as a null result."
  )
}

clock_draws <- as_draws_matrix(fit$draws(c("shift", "log_speed")))
season_clocks <- rbindlist(lapply(seq_len(stan_data$S), function(s) {
  sh <- clock_draws[, paste0("shift[", s, "]")]
  sp <- exp(clock_draws[, paste0("log_speed[", s, "]")])
  rbind(
    cbind(season_index = s, parameter = "shift_weeks", summ(sh)),
    cbind(season_index = s, parameter = "speed_ratio", summ(sp))
  )
}))
season_clocks <- merge(season_clocks, season_map, by = "season_index")
season_clocks[, programme := season_index %in%
                which(season_programme_indicator(stan_data) == 1L)]
setorder(season_clocks, parameter, season_index)
fwrite(season_clocks, "results/catalonia_season_clocks.csv")

timing_draws <- as.data.table(as_draws_df(fit$draws(c("sigma_shift", "sigma_speed"))))
timing_draws <- timing_draws[, .(sigma_shift, sigma_speed)]
fwrite(timing_draws, "results/catalonia_timing_sd_draws.csv")
fwrite(
  data.table(parameter = c("sigma_shift", "sigma_speed"),
             prior = "half-normal",
             prior_sd = c(stan_data$prior_sd_shift, stan_data$prior_sd_speed)),
  "results/catalonia_timing_priors.csv"
)

season_estimand <- function(fit, varname, metric) {
  m <- as_draws_matrix(fit$draws(varname))
  rbindlist(lapply(seq_len(ncol(m)), function(j) {
    cbind(metric = metric, season_index = j, summ(m[, j]))
  }))
}

season_summary <- rbindlist(list(
  season_estimand(fit, "prevented_by_season", "prevented_cases"),
  season_estimand(fit, "peak_attenuation", "peak_attenuation"),
  season_estimand(fit, "peak_shift_weeks", "peak_shift_weeks"),
  season_estimand(fit, "epidemic_peak_week", "epidemic_peak_week")
), fill = TRUE)
season_summary <- merge(season_summary, season_map, by = "season_index",
                        all.x = TRUE)
fwrite(season_summary, "results/catalonia_season_estimands.csv")

mu_fit <- as_draws_matrix(fit$draws("mu_fit_sw"))
mu_cf <- as_draws_matrix(fit$draws("mu_cf_sw"))

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
rm(mu_fit, mu_cf)

run_ppc <- env_bool("BEACON_PPC", TRUE)
n_ppc <- as.integer(Sys.getenv("BEACON_PPC_DRAWS", "400"))
if (run_ppc) {
  ppc_primary <- ppc_by_season(fit, stan_data, season_map, n_ppc)
  fwrite(cbind(analysis = "primary", ppc_primary$weekly),
         "results/catalonia_ppc_weekly.csv")
  ppc_all <- list(cbind(analysis = "primary", ppc_primary$season))
}

all_metrics <- list(primary$metrics)
all_ws <- list(primary$within_support_by_season)
all_diag <- list(primary$diagnostics)

calendar_data <- beacon_clock_off(stan_data)
calendar_data$compute_y_rep <- as.integer(run_ppc)
fit_cal <- run_beacon(calendar_data, "BEACON-calendar")
if (run_ppc) {
  ppc_cal <- ppc_by_season(fit_cal, calendar_data, season_map, n_ppc)
  fwrite(cbind(analysis = "calendar", ppc_cal$weekly),
         "results/catalonia_ppc_weekly_calendar.csv")
  ppc_all[[2]] <- cbind(analysis = "calendar", ppc_cal$season)
  ppc_season_all <- rbindlist(ppc_all, fill = TRUE)
  ppc_season_all[, abs_peak_week_error := abs(rep_peak_week_median -
                                                observed_peak_week)]
  fwrite(ppc_season_all, "results/catalonia_ppc_season.csv")
  fwrite(ppc_season_all[, .(mean_abs_peak_week_error = mean(abs_peak_week_error),
                            max_abs_peak_week_error = max(abs_peak_week_error),
                            seasons_peak_outside_90pct = sum(
                              observed_peak_week < rep_peak_week_lower |
                                observed_peak_week > rep_peak_week_upper)),
                        by = analysis],
         "results/catalonia_ppc_peak_timing_comparison.csv")
}
cal <- summarise_fit(fit_cal, calendar_data, d, "calendar")
all_metrics[[length(all_metrics) + 1L]] <- cal$metrics
all_ws[[length(all_ws) + 1L]] <- cal$within_support_by_season
all_diag[[length(all_diag) + 1L]] <- cal$diagnostics

cal_season <- rbindlist(list(
  season_estimand(fit_cal, "peak_attenuation", "peak_attenuation"),
  season_estimand(fit_cal, "peak_shift_weeks", "peak_shift_weeks"),
  season_estimand(fit_cal, "prevented_by_season", "prevented_cases")
), fill = TRUE)
cal_season <- merge(cal_season, season_map, by = "season_index", all.x = TRUE)
fwrite(rbindlist(list(
  cbind(analysis = "primary", season_summary[metric != "epidemic_peak_week"]),
  cbind(analysis = "calendar", cal_season)
), fill = TRUE), "results/catalonia_clock_vs_calendar_season_estimands.csv")

if (env_bool("BEACON_LOSO", FALSE)) {
loso <- try({
  loo_primary <- grouped_loso(fit, stan_data)
  loo_calendar <- grouped_loso(fit_cal, calendar_data)
  cmp <- loo::loo_compare(list(BEACON = loo_primary,
                               `BEACON-calendar` = loo_calendar))
  k_primary <- loo::pareto_k_values(loo_primary)
  k_calendar <- loo::pareto_k_values(loo_calendar)
  list(
    compare = data.table(model = rownames(cmp),
                         elpd_diff = cmp[, "elpd_diff"],
                         se_diff = cmp[, "se_diff"],
                         elpd_loso = cmp[, "elpd_loo"],
                         se_elpd_loso = cmp[, "se_elpd_loo"]),
    pareto_k = data.table(season_index = seq_len(stan_data$S),
                          season = season_levels,
                          k_beacon = k_primary, k_calendar = k_calendar)
  )
}, silent = TRUE)
if (inherits(loso, "try-error")) {
  warning("Leave-one-season-out comparison failed: ", as.character(loso))
} else {
  fwrite(loso$compare, "results/catalonia_loso_comparison.csv")
  fwrite(loso$pareto_k, "results/catalonia_loso_pareto_k.csv")
  if (any(unlist(loso$pareto_k[, .(k_beacon, k_calendar)]) > 0.7)) {
    warning(
      "Pareto k > 0.7 for at least one held-out season: the PSIS estimate is ",
      "unreliable for that season. Report it as such, alongside the ",
      "posterior predictive comparison of peak timing."
    )
  }
}
}
rm(fit_cal)

unpooled_data <- stan_data
unpooled_data$compute_y_rep <- 0L
unpooled_data$compute_log_lik <- 0L
unpooled_data$unpool_programme_seasons <- 1L
fit_unp <- run_beacon(unpooled_data, "unpooled programme seasons")
unp <- summarise_fit(fit_unp, unpooled_data, d, "unpooled_programme_seasons")
all_metrics[[length(all_metrics) + 1L]] <- unp$metrics
all_ws[[length(all_ws) + 1L]] <- unp$within_support_by_season
all_diag[[length(all_diag) + 1L]] <- unp$diagnostics
rm(fit_unp)

if (env_bool("BEACON_SENSITIVITY", TRUE)) {
  sens_specs <- list(
    list(label = "protective_prior", prior_sd = prior_sd_beta,
         mutate = function(sd) { sd$prior_mean_beta_intensity <- log(0.6); sd }),
    list(label = paste0("ascertainment_sd_", ascertainment_sd_sensitivity),
         prior_sd = prior_sd_beta,
         mutate = function(sd) {
           sd$prior_sd_ascertainment <- ascertainment_sd_sensitivity; sd }),
    list(label = "two_harmonics", prior_sd = prior_sd_beta,
         mutate = function(sd) { sd$H <- 2L; sd$harmonic_decay <- 0; sd })
  )

  for (spec in sens_specs) {
    sd_i <- spec$mutate(stan_data)
    sd_i$compute_log_lik <- 0L
    sd_i$compute_y_rep <- 0L
    f_i <- try(run_beacon(sd_i, spec$label, seed = 17072026), silent = TRUE)
    if (inherits(f_i, "try-error")) {
      warning("Sensitivity fit failed: ", spec$label)
      next
    }
    s_i <- summarise_fit(f_i, sd_i, d, spec$label, prior_sd = spec$prior_sd)
    all_metrics[[length(all_metrics) + 1L]] <- s_i$metrics
    all_ws[[length(all_ws) + 1L]] <- s_i$within_support_by_season
    all_diag[[length(all_diag) + 1L]] <- s_i$diagnostics
    rm(f_i)
  }

  if (length(excluded_seasons)) {
    d_incl <- prepare_model_data(d_all, character())
    built_incl <- build_model_data(d_incl, compute_log_lik = 0L, compute_y_rep = 0L)
    built_incl$stan_data$unpooled_season_sd <- unpooled_season_sd
    f_i <- try(run_beacon(built_incl$stan_data, "including excluded seasons"),
               silent = TRUE)
    if (inherits(f_i, "try-error")) {
      warning("Sensitivity fit failed: including excluded seasons")
    } else {
      s_i <- summarise_fit(f_i, built_incl$stan_data, built_incl$data,
                           paste0("including_", paste(excluded_seasons, collapse = "_")))
      all_metrics[[length(all_metrics) + 1L]] <- s_i$metrics
      all_ws[[length(all_ws) + 1L]] <- s_i$within_support_by_season
      all_diag[[length(all_diag) + 1L]] <- s_i$diagnostics
      rm(f_i)
    }
  }
}

if (env_bool("BEACON_PLACEBO", FALSE)) {
  placebo_seasons <- trimws(strsplit(
    Sys.getenv("BEACON_PLACEBO_SEASONS", "2016-2017,2017-2018,2018-2019"),
    ",", fixed = TRUE)[[1]])
  d_pl <- prepare_model_data(d_all, excluded_seasons)
  d_pl <- d_pl[treated == 0]
  if (all(placebo_seasons %in% d_pl$season)) {
    d_pl[, treated := as.integer(season %in% placebo_seasons)]
    d_pl[, immunisation_intensity := as.numeric(treated)]
    built_pl <- build_model_data(d_pl, compute_log_lik = 0L, compute_y_rep = 0L)
    built_pl$stan_data$unpooled_season_sd <- unpooled_season_sd
    f_i <- try(run_beacon(built_pl$stan_data, "in-time placebo"), silent = TRUE)
    if (inherits(f_i, "try-error")) {
      warning("Placebo fit failed.")
    } else {
      pl_label <- paste0("placebo_", paste(placebo_seasons, collapse = "_"))
      s_i <- summarise_fit(f_i, built_pl$stan_data, built_pl$data, pl_label)
      all_metrics[[length(all_metrics) + 1L]] <- s_i$metrics
      all_diag[[length(all_diag) + 1L]] <- s_i$diagnostics
      pl_levels <- sort(unique(built_pl$data$season_id))
      pl_est <- rbindlist(list(
        season_estimand(f_i, "peak_attenuation", "peak_attenuation"),
        season_estimand(f_i, "peak_shift_weeks", "peak_shift_weeks")
      ), fill = TRUE)
      pl_est[, season := pl_levels[season_index]]
      fwrite(pl_est[season %in% placebo_seasons],
             "results/catalonia_placebo_season_estimands.csv")
      rm(f_i)
    }
  } else {
    warning("Placebo seasons not all available among pre-programme seasons.")
  }
}

fwrite(rbindlist(all_metrics, fill = TRUE),
       "results/catalonia_sensitivity_summary.csv")
fwrite(rbindlist(all_ws, fill = TRUE),
       "results/catalonia_within_support_pf_all_analyses.csv")
fwrite(rbindlist(all_diag, fill = TRUE),
       "results/catalonia_sampler_diagnostics_all_analyses.csv")

diag <- fit$summary()[, c("variable", "rhat", "ess_bulk", "ess_tail")]
fwrite(as.data.table(diag), "results/catalonia_diagnostics.csv")
fwrite(
  melt(primary$diagnostics, id.vars = "analysis",
       variable.name = "quantity")[, .(quantity, value)],
  "results/catalonia_sampler_diagnostics.csv"
)
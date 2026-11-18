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
invisible(
  data.table(
    item = c("incomplete_seasons_dropped", "seasons_excluded", "seasons_used"),
    value = c(paste(incomplete_seasons, collapse = ";"),
              paste(excluded_seasons, collapse = ";"),
              paste(sort(unique(d$season)), collapse = ";"))
  )
)

## ---------------------------------------------------------------------------
## Overlap first. If the cross-sectional contrast is degenerate, the effect is
## informed mainly by the pre-post component (see the unpooled-season fit).
## ---------------------------------------------------------------------------
overlap <- overlap_diagnostic(d, area_col = "abs_code", season_col = "season")
# (overlap diagnostic not rewritten by the placebo script)
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
# adapt_delta is the target acceptance rate. The manuscript states 0.97.
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

season_estimand <- function(fit, varname, metric) {
  m <- as_draws_matrix(fit$draws(varname))
  rbindlist(lapply(seq_len(ncol(m)), function(j) {
    cbind(metric = metric, season_index = j, summ(m[, j]))
  }))
}

all_metrics <- list()
all_diag <- list()

if (env_bool("BEACON_PLACEBO", TRUE)) {
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

if (length(all_metrics)) {
  fwrite(rbindlist(all_metrics, fill = TRUE), "results/catalonia_placebo_summary.csv")
  fwrite(rbindlist(all_diag, fill = TRUE),
         "results/catalonia_placebo_sampler_diagnostics.csv")
} 
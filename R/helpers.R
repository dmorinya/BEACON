clean_names_base <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("(^_+|_+$)", "", x)
}

find_col <- function(nms, patterns, required = TRUE, label = paste(patterns, collapse = "/")) {
  nms_clean <- clean_names_base(nms)
  hit <- unique(unlist(lapply(patterns, function(p) grep(p, nms_clean, perl = TRUE))))
  if (!length(hit)) {
    if (required) stop("Could not identify column '", label, "'. Available columns: ", paste(nms, collapse = ", "))
    return(NA_character_)
  }
  nms[hit[1]]
}

parse_date_flex <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.Date(x))
  bad <- is.na(out)
  if (any(bad)) out[bad] <- suppressWarnings(as.Date(x[bad], format = "%d/%m/%Y"))
  bad <- is.na(out)
  if (any(bad)) out[bad] <- suppressWarnings(as.Date(substr(x[bad], 1, 10), format = "%Y-%m-%d"))
  out
}

epidemic_season <- function(date, start_month = 7L) {
  y <- as.integer(format(date, "%Y"))
  m <- as.integer(format(date, "%m"))
  start <- ifelse(m >= start_month, y, y - 1L)
  paste0(start, "-", start + 1L)
}

safe_download <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading ", url)
    tryCatch(
      utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
      error = function(e) stop("Download failed. Download the CSV manually from the open-data portal and save it as: ", dest, "\nOriginal error: ", conditionMessage(e))
    )
  }
  invisible(dest)
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

q_summary <- function(x, prob = c(.025, .5, .975)) {
  stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE)
}

## ===========================================================================
## BEACON v2 helpers
##
## Added for beacon_affine_v2.stan. The v1 data block did not carry an
## area-season index, a week pivot, or the contiguous series blocks that the
## dynamic g-formula needs, so the Stan data list must be rebuilt rather than
## extended.
## ===========================================================================

BEACON_V2_MARKERS <- c("sum_to_zero_helmert", "week_pivot", "series_start",
                       "reciprocal_phi", "beta_intensity_cos", "prior_only")

## Guard replacing the v1 checks for inv_phi and the bounded-phi declaration,
## which no longer exist in v2.
assert_beacon_v2 <- function(stan_file) {
  if (!file.exists(stan_file)) stop("Stan model not found: ", stan_file)
  code <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")

  if (grepl("tanh(shift_raw", code, fixed = TRUE)) {
    stop("This is the v1 model (soft tanh clock). Use beacon_affine_v2.stan.")
  }
  missing <- BEACON_V2_MARKERS[
    !vapply(BEACON_V2_MARKERS, grepl, logical(1), x = code, fixed = TRUE)
  ]
  if (length(missing)) {
    stop("Stan file does not look like BEACON v2; missing: ",
         paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## Build the v2 Stan data list.
##
## dat must contain: area_id, season_id, week_in_season, y, population,
## immunisation_intensity, treated, ses. The lagged covariate is rebuilt here
## rather than taken from the caller, because the recursion inside Stan has to
## use exactly the same centring and scaling constants.
##
## Returns the Stan list AND the reordered data frame, so downstream joins use
## the same row order.
## ---------------------------------------------------------------------------
build_beacon_stan_data <- function(dat,
                                   H = 4L,
                                   harmonic_decay = 1,
                                   intensity_quadratic = 0L,
                                   prior_mean_beta_intensity = 0,
                                   prior_sd_beta_intensity = 0.5,
                                   prior_sd_shift = 3,
                                   prior_sd_speed = 0.15,
                                   prior_sd_ascertainment = 0,
                                   prior_only = 0L,
                                   compute_log_lik = 0L,
                                   compute_y_rep = 0L,
                                   dynamic_gformula = 1L,
                                   lag_init = 0) {

  required <- c("area_id", "season_id", "week_in_season", "y", "population",
                "immunisation_intensity", "treated", "ses")
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols)) {
    stop("build_beacon_stan_data(): missing columns: ",
         paste(missing_cols, collapse = ", "))
  }

  dat <- data.table::as.data.table(data.table::copy(dat))
  dat[, area := match(area_id, sort(unique(area_id)))]
  dat[, season := match(season_id, sort(unique(season_id)))]

  ## Contiguous, time-ordered blocks, one per area-season.
  data.table::setorder(dat, area, season, week_in_season)
  dat[, series := .GRP, by = .(area, season)]

  if (any(diff(dat$series) < 0)) {
    stop("series blocks are not contiguous after sorting")
  }
  if (dat[, any(diff(week_in_season) <= 0), by = series][, any(V1)]) {
    stop("weeks are not strictly increasing within at least one series")
  }

  ## Lag WITHIN series. Lagging by area alone, as the v1 preparation did,
  ## carries the last week of June into the first week of July of the next
  ## season and breaks the recursion.
  dat[, lag_raw := data.table::shift(log1p(y)), by = series]
  lag_center <- mean(dat$lag_raw, na.rm = TRUE)
  lag_scale <- stats::sd(dat$lag_raw, na.rm = TRUE)
  if (!is.finite(lag_scale) || lag_scale <= 0) lag_scale <- 1
  dat[, lag_log_y := (lag_raw - lag_center) / lag_scale]
  dat[!is.finite(lag_log_y), lag_log_y := lag_init]

  starts <- dat[, .I[1], by = series]$V1
  lens <- dat[, .N, by = series]$N

  stan_data <- list(
    N = nrow(dat),
    A = max(dat$area),
    S = max(dat$season),
    M = max(dat$series),
    W = max(dat$week_in_season),

    area = as.integer(dat$area),
    season = as.integer(dat$season),
    series = as.integer(dat$series),
    week_in_season = as.integer(dat$week_in_season),

    series_start = as.integer(starts),
    series_len = as.integer(lens),

    y = as.integer(round(dat$y)),
    log_population = log(dat$population),
    week = as.numeric(dat$week_in_season),
    ## Pivot in the middle of the season rather than at week zero; see [C11].
    week_pivot = mean(as.numeric(dat$week_in_season)),

    immunisation_intensity = as.numeric(dat$immunisation_intensity),
    treated = as.numeric(dat$treated),
    ses = as.numeric(dat$ses),
    lag_log_y = as.numeric(dat$lag_log_y),

    lag_center = lag_center,
    lag_scale = lag_scale,
    lag_init = lag_init,

    H = as.integer(H),
    harmonic_decay = harmonic_decay,
    intensity_quadratic = as.integer(intensity_quadratic),

    prior_mean_beta_intensity = prior_mean_beta_intensity,
    prior_sd_beta_intensity = prior_sd_beta_intensity,
    prior_sd_shift = prior_sd_shift,
    prior_sd_speed = prior_sd_speed,
    prior_sd_ascertainment = prior_sd_ascertainment,

    prior_only = as.integer(prior_only),
    compute_log_lik = as.integer(compute_log_lik),
    compute_y_rep = as.integer(compute_y_rep),
    dynamic_gformula = as.integer(dynamic_gformula)
  )

  finite_inputs <- c("log_population", "week", "immunisation_intensity",
                     "treated", "ses", "lag_log_y")
  bad <- vapply(finite_inputs,
                function(nm) sum(!is.finite(stan_data[[nm]])), integer(1))
  if (any(bad > 0L)) {
    stop("Non-finite values in Stan data: ",
         paste(names(bad)[bad > 0L], bad[bad > 0L], sep = "=", collapse = ", "))
  }
  if (any(stan_data$immunisation_intensity < 0)) {
    stop("Immunisation intensity must be non-negative.")
  }
  if (!all(stan_data$treated %in% c(0, 1))) {
    stop("treated must contain only 0 and 1.")
  }

  list(stan_data = stan_data, data = dat)
}

## ---------------------------------------------------------------------------
## Dispersed initial values on the v2 parameter names.
##
## These are deliberately NOT centred on the data-generating values: a coverage
## study initialised at the truth is not a neutral test. Only alpha is anchored,
## from the crude event rate, because the offset makes its scale predictable.
## ---------------------------------------------------------------------------
make_beacon_inits_v2 <- function(stan_data, chains, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  crude_rate <- (sum(stan_data$y) + 0.5) /
    (sum(exp(stan_data$log_population)) + 0.5)
  alpha_init <- log(crude_rate)
  S <- stan_data$S
  H <- stan_data$H

  lapply(seq_len(chains), function(chain_id) {
    list(
      alpha = alpha_init + rnorm(1, 0, 0.25),
      z_area = rnorm(stan_data$A, 0, 0.5),
      z_season = rnorm(S, 0, 0.5),
      z_series = rnorm(stan_data$M, 0, 0.5),
      sigma_area = abs(rnorm(1, 0.25, 0.10)),
      sigma_season = abs(rnorm(1, 0.20, 0.10)),
      sigma_series = abs(rnorm(1, 0.15, 0.08)),
      ## length S - 1 under the exact zero-sum parameterisation
      shift_raw = rnorm(S - 1L, 0, 0.5),
      log_speed_raw = rnorm(S - 1L, 0, 0.5),
      sigma_shift = abs(rnorm(1, 1.0, 0.5)),
      sigma_speed = abs(rnorm(1, 0.06, 0.03)),
      a_raw = rnorm(H, 0, 0.5),
      b_raw = rnorm(H, 0, 0.5),
      sigma_f = abs(rnorm(1, 0.6, 0.2)),
      beta_lag = rnorm(1, 0, 0.10),
      beta_ses = rnorm(1, 0, 0.10),
      beta_intensity = rnorm(1, 0, 0.25),
      beta_intensity_cos = rnorm(1, 0, 0.10),
      beta_intensity_sin = rnorm(1, 0, 0.10),
      beta_intensity_sq = rnorm(1, 0, 0.10),
      lambda_raw = rnorm(1, 0, 0.5),
      reciprocal_phi = abs(rnorm(1, 0.10, 0.05)) + 1e-3
    )
  })
}

## ---------------------------------------------------------------------------
## Overlap / positivity diagnostic.
##
## Run this BEFORE interpreting any effect estimate. beta_intensity is
## identified by between-area variation in intensity within programme seasons,
## because the season effects absorb anything programme-wide. If almost every
## area sits at the season median there is little such variation and the
## posterior will revert towards its prior.
## ---------------------------------------------------------------------------
overlap_diagnostic <- function(dat,
                               area_col = "abs_code",
                               season_col = "season") {
  d <- data.table::as.data.table(dat)
  d <- d[treated == 1]
  d <- unique(d[, c(area_col, season_col, "immunisation_intensity"),
                with = FALSE])
  data.table::setnames(d, c(area_col, season_col), c("area_key", "season_id"))

  d[is.finite(immunisation_intensity),
    .(
      n_abs = .N,
      median_intensity = stats::median(immunisation_intensity),
      sd_intensity = stats::sd(immunisation_intensity),
      iqr_intensity = stats::IQR(immunisation_intensity),
      p10 = unname(stats::quantile(immunisation_intensity, 0.10)),
      p90 = unname(stats::quantile(immunisation_intensity, 0.90)),
      prop_within_5pct = mean(
        abs(immunisation_intensity /
              stats::median(immunisation_intensity) - 1) < 0.05
      )
    ),
    by = season_id][order(season_id)]
}

## Ratio of posterior to prior SD for the intervention coefficient. Values near
## 1 mean the data carry almost no information about the effect.
prior_posterior_contraction <- function(draws, prior_sd) {
  post_sd <- stats::sd(draws)
  data.table::data.table(
    prior_sd = prior_sd,
    posterior_sd = post_sd,
    contraction = 1 - post_sd / prior_sd
  )
}

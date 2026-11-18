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

BEACON_V2_MARKERS <- c("sum_to_zero_helmert", "week_pivot", "series_start",
                       "reciprocal_phi", "beta_intensity_cos", "prior_only",
                       "use_clock", "unpool_programme_seasons",
                       "tau_grid[g] = -26.0",
                       "fourier_grid", "sw_index")

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
    stop("Stan file is not the revised (v3) BEACON program; missing: ",
         paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

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
    dynamic_gformula = as.integer(dynamic_gformula),

    use_clock = 1L,
    unpool_programme_seasons = 0L,
    season_programme = as.integer(dat[, max(treated), by = season][order(season)]$V1),
    unpooled_season_sd = 5
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
      min_intensity = min(immunisation_intensity),
      max_intensity = max(immunisation_intensity),
      prop_within_5pct = mean(
        abs(immunisation_intensity /
              stats::median(immunisation_intensity) - 1) < 0.05
      )
    ),
    by = season_id][order(season_id)]
}

prior_posterior_contraction <- function(draws, prior_sd) {
  post_sd <- stats::sd(draws)
  data.table::data.table(
    prior_sd = prior_sd,
    posterior_sd = post_sd,
    contraction = 1 - post_sd / prior_sd
  )
}

beacon_clock_off <- function(stan_data) {
  stan_data$use_clock <- 0L
  stan_data
}

stan_supports <- function(stan_file, marker) {
  code <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  grepl(marker, code, fixed = TRUE)
}

fit_has_variable <- function(fit, variable) {
  vars <- tryCatch(fit$metadata()$stan_variables, error = function(e) character())
  variable %in% vars
}

season_programme_indicator <- function(stan_data) {
  v <- tapply(stan_data$treated, stan_data$season, max)
  as.integer(v[as.character(seq_len(stan_data$S))])
}

activity_centiles_by_season <- function(d, probs = c(0.10, 0.90)) {
  dd <- unique(data.table::as.data.table(d)[treated == 1,
                                             .(area, season, immunisation_intensity)])
  dd[, .(lo = unname(stats::quantile(immunisation_intensity, probs[1])),
         hi = unname(stats::quantile(immunisation_intensity, probs[2])),
         n_area = .N), by = season][order(season)]
}

within_support_pf <- function(fit, stan_data, centiles, n_draws = 1000L,
                              seed = 1L) {
  quadratic <- isTRUE(stan_data$intensity_quadratic == 1L)
  vars <- c("beta_intensity", "beta_intensity_cos", "beta_intensity_sin",
            "shift", "log_speed", "mu_cf_sw",
            if (quadratic) "beta_intensity_sq")
  dr <- posterior::as_draws_matrix(fit$draws(vars))
  dr <- matrix(as.numeric(dr), nrow = nrow(dr),
               dimnames = list(NULL, posterior::variables(dr)))
  set.seed(seed)
  keep <- sort(sample.int(nrow(dr), min(n_draws, nrow(dr))))
  dr <- dr[keep, , drop = FALSE]

  W <- stan_data$W
  pivot <- stan_data$week_pivot
  weeks <- seq_len(W)
  out_season <- list()
  num_total <- numeric(nrow(dr))
  den_total <- numeric(nrow(dr))

  for (k in seq_len(nrow(centiles))) {
    s <- centiles$season[k]
    mu_cols <- paste0("mu_cf_sw[", s, ",", weeks, "]")
    mu <- dr[, mu_cols, drop = FALSE]
    shift <- dr[, paste0("shift[", s, "]")]
    speed <- exp(dr[, paste0("log_speed[", s, "]")])
    tau <- outer(-shift, weeks - pivot, "+") * speed      # draws x W
    m <- dr[, "beta_intensity"] +
      dr[, "beta_intensity_cos"] * cos(2 * pi * tau / 52) +
      dr[, "beta_intensity_sin"] * sin(2 * pi * tau / 52)
    q_hi <- if (quadratic) dr[, "beta_intensity_sq"] * centiles$hi[k]^2 else 0
    q_lo <- if (quadratic) dr[, "beta_intensity_sq"] * centiles$lo[k]^2 else 0
    num <- rowSums(mu * exp(centiles$hi[k] * m + q_hi))
    den <- rowSums(mu * exp(centiles$lo[k] * m + q_lo))
    num_total <- num_total + num
    den_total <- den_total + den
    pf <- 1 - num / den
    out_season[[k]] <- data.table::data.table(
      season_index = s, activity_lo = centiles$lo[k],
      activity_hi = centiles$hi[k],
      mean = mean(pf), lower = unname(stats::quantile(pf, .025)),
      median = unname(stats::quantile(pf, .5)),
      upper = unname(stats::quantile(pf, .975))
    )
  }
  pf_all <- 1 - num_total / den_total
  list(
    by_season = data.table::rbindlist(out_season),
    overall = data.table::data.table(
      mean = mean(pf_all), lower = unname(stats::quantile(pf_all, .025)),
      median = unname(stats::quantile(pf_all, .5)),
      upper = unname(stats::quantile(pf_all, .975))
    )
  )
}

grouped_loso <- function(fit, stan_data) {
  if (!requireNamespace("loo", quietly = TRUE)) stop("Package 'loo' needed.")
  ll <- fit$draws("log_lik", format = "draws_array")   # iter x chain x N
  n_iter <- dim(ll)[1]; n_chain <- dim(ll)[2]
  season <- stan_data$season
  S <- stan_data$S
  grouped <- array(NA_real_, dim = c(n_iter, n_chain, S))
  for (s in seq_len(S)) {
    idx <- which(season == s)
    grouped[, , s] <- apply(ll[, , idx, drop = FALSE], c(1, 2), sum)
  }
  scaled <- grouped
  for (s in seq_len(S)) scaled[, , s] <- exp(grouped[, , s] - max(grouped[, , s]))
  r_eff <- loo::relative_eff(scaled)
  loo::loo(grouped, r_eff = r_eff)
}

half_normal_density <- function(x, sd) {
  ifelse(x < 0, 0, 2 * stats::dnorm(x, 0, sd))
}

pc_by_season <- function(fit, stan_data, season_map, n_ppc = 400L,
                          seed = 20260101) {
  yrep <- posterior::as_draws_matrix(fit$draws("y_rep"))
  set.seed(seed)
  yrep <- yrep[sort(sample.int(nrow(yrep), min(n_ppc, nrow(yrep)))), , drop = FALSE]
  group_key <- paste(stan_data$season, stan_data$week_in_season, sep = "_")
  agg <- rowsum(t(unclass(yrep)), group = group_key)   
  rm(yrep)
  key_dt <- data.table::data.table(group_key = rownames(agg))
  key_dt[, c("season_index", "week_in_season") :=
           data.table::tstrsplit(group_key, "_", type.convert = TRUE)]
  obs <- data.table::data.table(season_index = stan_data$season,
                                week_in_season = stan_data$week_in_season,
                                y = stan_data$y)[, .(observed = sum(y)),
                                                 by = .(season_index, week_in_season)]

  weekly <- cbind(key_dt[, .(season_index, week_in_season)], data.table::data.table(
    rep_median = apply(agg, 1, stats::median),
    rep_lower = apply(agg, 1, stats::quantile, 0.05),
    rep_upper = apply(agg, 1, stats::quantile, 0.95)
  ))
  weekly <- merge(weekly, obs, by = c("season_index", "week_in_season"))
  weekly <- merge(weekly, season_map, by = "season_index")
  data.table::setorder(weekly, season_index, week_in_season)

  season <- data.table::rbindlist(lapply(seq_len(stan_data$S), function(s) {
    rows <- which(key_dt$season_index == s)
    wk <- key_dt$week_in_season[rows]
    ord <- order(wk)
    m <- agg[rows[ord], , drop = FALSE]
    wk <- wk[ord]
    o <- obs[season_index == s][order(week_in_season)]
    rep_peak_week <- wk[apply(m, 2, which.max)]
    rep_peak_height <- apply(m, 2, max)
    rep_total <- colSums(m)
    data.table::data.table(
      season_index = s,
      observed_total = sum(o$observed),
      rep_total_median = stats::median(rep_total),
      ppp_total = mean(rep_total >= sum(o$observed)),
      observed_peak_week = o$week_in_season[which.max(o$observed)],
      rep_peak_week_median = stats::median(rep_peak_week),
      rep_peak_week_lower = unname(stats::quantile(rep_peak_week, 0.05)),
      rep_peak_week_upper = unname(stats::quantile(rep_peak_week, 0.95)),
      observed_peak_height = max(o$observed),
      rep_peak_height_median = stats::median(rep_peak_height),
      ppp_peak_height = mean(rep_peak_height >= max(o$observed))
    )
  }))
  season <- merge(season, season_map, by = "season_index")
  list(weekly = weekly, season = season)
}
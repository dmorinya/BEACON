library(data.table)
library(MASS)
library(mgcv)
library(future.apply)

source("helpers.R")

dir.create("results", recursive = TRUE, showWarnings = FALSE)

get_env_int <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value < minimum) value <- as.integer(default)
  value
}

get_env_num <- function(name, default, minimum = -Inf, maximum = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value < minimum || value > maximum) value <- default
  value
}

get_env_bool <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
  value %in% c("true", "1", "yes", "y")
}

base_seed <- get_env_int("BEACON_SIM_SEED", 17072026L)
n_rep_fast <- get_env_int("BEACON_FAST_NREP", 500L)
n_rep_stan <- get_env_int("BEACON_STAN_NREP", 200L)
n_fast_cores <- get_env_int(
  "BEACON_FAST_CORES",
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
)
skip_stan <- get_env_bool("BEACON_SKIP_STAN", FALSE)
skip_fast <- get_env_bool("BEACON_SKIP_FAST", FALSE)
collect_only <- get_env_bool("BEACON_STAN_COLLECT_ONLY", FALSE)
# Number of Bayesian fits run at the same time within this job.
stan_workers <- get_env_int("BEACON_STAN_WORKERS", 1L)
# Optional repetition range for array jobs, e.g. "1-50".
stan_reps_env <- Sys.getenv("BEACON_STAN_REPS", "")
write_summaries <- !nzchar(stan_reps_env) || collect_only

# Four chains by default. Rank-normalised R-hat is a weak diagnostic with two.
stan_chains <- get_env_int("BEACON_STAN_CHAINS", 4L)
# With several fits in parallel, run the chains of each fit sequentially by
# default, so that the job uses about BEACON_STAN_WORKERS cores.
stan_parallel_chains <- min(
  stan_chains,
  get_env_int("BEACON_STAN_PARALLEL_CHAINS",
              if (stan_workers > 1L) 1L else stan_chains)
)
stan_warmup <- get_env_int("BEACON_STAN_WARMUP", 300L)
stan_sampling <- get_env_int("BEACON_STAN_SAMPLING", 300L)
# adapt_delta is the target acceptance rate (the manuscript previously called
# it "step size"). Raised from 0.90 to 0.95 to reduce divergences [C15].
stan_adapt_delta <- get_env_num("BEACON_STAN_ADAPT_DELTA", 0.95, 0.5, 0.999)
stan_max_treedepth <- get_env_int("BEACON_STAN_MAX_TREEDEPTH", 11L)
stan_refresh <- get_env_int("BEACON_STAN_REFRESH", 0L, minimum = 0L)
force_recompile <- get_env_bool("BEACON_FORCE_RECOMPILE", FALSE)
resume_stan <- get_env_bool("BEACON_STAN_RESUME", TRUE)

# Long-chain robustness subset: the first n repetitions of each scenario are
# re-fitted (BEACON only) with a full-length schedule and identical seeds.
n_rep_stan_long <- get_env_int("BEACON_STAN_LONG_NREP", 20L, minimum = 0L)
stan_warmup_long <- get_env_int("BEACON_STAN_LONG_WARMUP", 1000L)
stan_sampling_long <- get_env_int("BEACON_STAN_LONG_SAMPLING", 1000L)
stan_adapt_delta_long <- get_env_num("BEACON_STAN_LONG_ADAPT_DELTA",
                                     0.95, 0.5, 0.999)

# Fitted-model configuration, kept here so the manuscript and the simulation
# cannot drift apart.
stan_H <- get_env_int("BEACON_STAN_H", 4L)
stan_harmonic_decay <- get_env_num("BEACON_STAN_HARMONIC_DECAY", 1.0, 0)
stan_prior_mean_beta <- get_env_num("BEACON_PRIOR_MEAN_BETA", 0.0)
stan_prior_sd_beta <- get_env_num("BEACON_PRIOR_SD_BETA", 0.5, 1e-6)
bias_analysis_sd <- get_env_num("BEACON_BIAS_ANALYSIS_SD", 0.15, 0)

quick_mode <- get_env_bool("BEACON_QUICK", FALSE)
if (quick_mode) {
  n_rep_fast <- min(n_rep_fast, 20L)
  n_rep_stan <- min(n_rep_stan, 2L)
  n_rep_stan_long <- 0L
  stan_warmup <- min(stan_warmup, 100L)
  stan_sampling <- min(stan_sampling, 100L)
}

## ---------------------------------------------------------------------------
## Definition for the scenarios.
## ---------------------------------------------------------------------------
scenarios <- data.table(
  scenario = c(
    "aligned",
    "phase_shift",
    "phase_and_ascertainment",
    "phase_and_confounding",
    "phase_dependent_effect",
    "limited_overlap"
  ),
  scenario_label = c(
    "1. Aligned",
    "2. Phase shift",
    "3. Phase + ascertainment",
    "4. Phase + confounding",
    "5. Phase-dependent effect",
    "6. Limited overlap"
  ),
  shift_sd = c(0, 2.5, 2.5, 2.5, 2.5, 2.5),
  speed_sd = c(0, 0.08, 0.08, 0.08, 0.08, 0.08),
  ascertainment_change = c(0, 0, -0.15, 0, 0, 0),
  confounding = c(0, 0, 0, 0.60, 0, 0),
  # Amplitude of the phase-locked modulation (largest reduction at the peak).
  modulation_amplitude = c(0, 0, 0, 0, 0.25, 0),
  # Activity A*_r = exp{ses_coef Z_r + c U_r + eps_r}, eps_r ~ N(0, noise_sd^2).
  activity_ses_coef = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.03),
  activity_noise_sd = c(0.18, 0.18, 0.18, 0.18, 0.18, 0.04),
  # Extra delay (weeks) of the programme seasons before centring.
  programme_delay = c(0, 0, 0, 0, 0, 3)
)
scenario_levels <- scenarios$scenario_label

method_levels <- c("ITS", "GAM", "Aligned NB", "BEACON-calendar", "BEACON",
                   "BEACON + bias analysis")

requested_scenarios <- trimws(
  strsplit(
    Sys.getenv("BEACON_STAN_SCENARIOS", paste(scenarios$scenario, collapse = ",")),
    ",",
    fixed = TRUE
  )[[1]]
)
requested_scenarios <- intersect(requested_scenarios, scenarios$scenario)
if (!length(requested_scenarios)) {
  stop("BEACON_STAN_SCENARIOS did not match any configured scenario.")
}

simulation_constants <- list(
  n_area = 10L,
  n_season = 6L,
  n_week = 52L,
  intervention_season = 5L,
  # Same pivot as the fitted model: mean week of the analysis window.
  week_pivot = (1 + 52) / 2,
  alpha = -6.35,
  sigma_area = 0.25,
  sigma_season = 0.15,
  beta_lag = 0.10,
  beta_ses = 0.15,
  beta_intensity = log(0.55),
  phi = 12,
  fourier = c(0.75, -0.15, 0.30, -0.10),
  # Area-season disturbance w_rs
  sigma_area_season = 0.15
)

## Baseline curve on the model's latent scale. Written as f(tau + t0) so that,
## with delta = kappa = 0, the calendar epidemic is identical to the first
## submission's.
dgp_f <- function(tau, K = simulation_constants) {
  x <- tau + K$week_pivot
  K$fourier[1] * sin(2 * pi * x / 52) +
    K$fourier[2] * cos(2 * pi * x / 52) +
    K$fourier[3] * sin(4 * pi * x / 52) +
    K$fourier[4] * cos(4 * pi * x / 52)
}

## Latent peak of f, used to lock the modulation to the epidemic peak.
tau_grid <- seq(-simulation_constants$week_pivot, 52 - simulation_constants$week_pivot,
                by = 0.01)
tau_star <- tau_grid[which.max(dgp_f(tau_grid))]

## True modulation coefficients for a scenario:
## m(tau) = beta0 - amp cos{2 pi (tau - tau*)/52}
##        = beta0 + beta_cos cos(2 pi tau/52) + beta_sin sin(2 pi tau/52).
true_modulation <- function(amp) {
  c(beta_cos = -amp * cos(2 * pi * tau_star / 52),
    beta_sin = -amp * sin(2 * pi * tau_star / 52))
}

## Clip, then re-centre, until the vector both sums to zero and is in bounds.
## (The first submission centred and then clipped, which could break the zero
## sum that the fitted model imposes.)
bounded_center <- function(x, lower, upper, max_iter = 100L) {
  for (i in seq_len(max_iter)) {
    x <- pmin(pmax(x - mean(x), lower), upper)
    if (abs(mean(x)) < 1e-10) break
  }
  x - mean(x)
}

simulate_one <- function(seed, sc) {
  set.seed(seed)

  K <- simulation_constants
  A <- K$n_area
  S <- K$n_season
  W <- K$n_week
  mod <- true_modulation(sc$modulation_amplitude)

  population_area <- pmax(
    300L,
    round(exp(rnorm(A, log(950), 0.22)))
  )
  ses_area <- as.numeric(scale(rnorm(A)))
  latent_access <- rnorm(A)

  area_re <- rnorm(A, 0, K$sigma_area)
  season_re <- rnorm(S, 0, K$sigma_season)
  # w_rs: one draw per area-season, constant across the 52 weeks of a block,
  # i.e. at exactly the level at which programme activity varies.
  area_season_re <- matrix(
    rnorm(A * S, 0, K$sigma_area_season),
    nrow = A, ncol = S
  )

  programme_season <- seq_len(S) >= K$intervention_season

  shift <- if (sc$shift_sd > 0) {
    raw <- rnorm(S, 0, sc$shift_sd) + sc$programme_delay * programme_season
    bounded_center(raw, -6, 6)
  } else {
    rep(0, S)
  }

  log_speed <- if (sc$speed_sd > 0) {
    bounded_center(rnorm(S, 0, sc$speed_sd), -0.25, 0.25)
  } else {
    rep(0, S)
  }

  intensity_area_raw <- exp(
    sc$activity_ses_coef * ses_area +
      sc$confounding * latent_access +
      rnorm(A, 0, sc$activity_noise_sd)
  )
  # Scale by the median so that activity = 1 represents typical programme
  # activity and activity = 0 represents no programme.
  intensity_area <- intensity_area_raw / median(intensity_area_raw)

  out <- vector("list", A * S)
  pos <- 1L

  for (a in seq_len(A)) {
    for (s in seq_len(S)) {
      previous_y <- 0
      treated_s <- as.integer(s >= K$intervention_season)
      intensity_s <- if (treated_s == 1L) intensity_area[a] else 0

      rows <- vector("list", W)

      for (w in seq_len(W)) {
        tau <- (w - K$week_pivot - shift[s]) * exp(log_speed[s])
        epidemic <- dgp_f(tau)

        # A fixed, interpretable transform is used during generation; the
        # fitted model re-standardises empirically.
        lag_log_y <- if (w == 1L) 0 else (log1p(previous_y) - 1.0) / 1.25

        intervention <- treated_s * intensity_s *
          (
            K$beta_intensity +
              mod[["beta_cos"]] * cos(2 * pi * tau / 52) +
              mod[["beta_sin"]] * sin(2 * pi * tau / 52)
          )

        latent_risk_confounding <-
          0.22 * sc$confounding * latent_access[a] * treated_s
        differential_ascertainment <-
          sc$ascertainment_change * treated_s * intensity_s

        # True outcome scale (no ascertainment), with and without programme.
        eta_true <-
          log(population_area[a]) +
          K$alpha +
          area_re[a] +
          season_re[s] +
          epidemic +
          K$beta_lag * lag_log_y +
          K$beta_ses * ses_area[a] +
          area_season_re[a, s] +
          latent_risk_confounding +
          intervention

        mu <- exp(eta_true + differential_ascertainment)
        y <- rnbinom(1L, mu = mu, size = K$phi)

        rows[[w]] <- data.table(
          y = y,
          population = population_area[a],
          area = a,
          season = s,
          week = w,
          treated = treated_s,
          immunisation_intensity = intensity_s,
          ses = ses_area[a],
          lag_log_y = lag_log_y,
          tau_true = tau,
          shift_true = shift[s],
          log_speed_true = log_speed[s],
          # Controlled contrast at the observed lag, outcome scale.
          mu_true_programme = exp(eta_true),
          mu_true_zero = exp(eta_true - intervention)
        )

        previous_y <- y
      }

      out[[pos]] <- rbindlist(rows)
      pos <- pos + 1L
    }
  }

  dat <- rbindlist(out)
  setattr(dat, "true_pf", dat[treated == 1L,
                              1 - sum(mu_true_programme) / sum(mu_true_zero)])
  dat
}

true_values <- function(dat, sc) {
  mod <- true_modulation(sc$modulation_amplitude)
  list(
    beta0 = simulation_constants$beta_intensity,
    prevented_fraction = attr(dat, "true_pf"),
    modulation_amplitude = sc$modulation_amplitude,
    beta_cos = unname(mod["beta_cos"]),
    beta_sin = unname(mod["beta_sin"])
  )
}

## ---------------------------------------------------------------------------
## Comparators
## ---------------------------------------------------------------------------

extract_estimates <- function(fit, dat, coefficient_pattern) {
  na4 <- c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, se = NA_real_)
  empty <- list(beta0 = na4, prevented_fraction = na4)
  if (inherits(fit, "try-error")) return(empty)

  coefficient_name <- grep(coefficient_pattern, names(stats::coef(fit)),
                           value = TRUE)[1]
  if (is.na(coefficient_name)) return(empty)

  estimate <- unname(stats::coef(fit)[coefficient_name])
  se <- sqrt(stats::vcov(fit)[coefficient_name, coefficient_name])
  lo <- estimate - 1.96 * se
  hi <- estimate + 1.96 * se

  mu_hat <- as.numeric(stats::predict(fit, newdata = dat, type = "response"))
  prog <- dat$treated == 1L
  a <- dat$immunisation_intensity[prog]
  mu0 <- mu_hat[prog] * exp(-a * estimate)
  pf_at <- function(b) 1 - sum(mu0 * exp(a * b)) / sum(mu0)
  pf_limits <- sort(c(pf_at(lo), pf_at(hi)))

  list(
    beta0 = c(estimate = estimate, lower = lo, upper = hi, se = se),
    prevented_fraction = c(estimate = pf_at(estimate), lower = pf_limits[1],
                           upper = pf_limits[2],
                           se = (pf_limits[2] - pf_limits[1]) / 3.92)
  )
}

fit_its <- function(dat) {
  fit <- try(
    MASS::glm.nb(
      y ~ offset(log(population)) +
        splines::ns(week, df = 6) +
        factor(season) +
        factor(area) +
        lag_log_y +
        ses +
        treated:immunisation_intensity,
      data = dat,
      control = glm.control(maxit = 80)
    ),
    silent = TRUE
  )

  extract_estimates(fit, dat, "treated:immunisation_intensity")
}

fit_gam <- function(dat) {
  dd <- copy(dat)
  # Random-effect smooths need factors; the first submission passed the
  # integer codes, which mgcv treats as a single linear "random slope".
  dd[, area_f := factor(area)]
  dd[, season_f := factor(season)]
  fit <- try(
    mgcv::gam(
      y ~ offset(log(population)) +
        s(week, k = 12, bs = "cc") +
        s(area_f, bs = "re") +
        s(season_f, bs = "re") +
        lag_log_y +
        ses +
        treated:immunisation_intensity,
      family = mgcv::nb(),
      data = dd,
      method = "REML",
      knots = list(week = c(0.5, 52.5))
    ),
    silent = TRUE
  )

  extract_estimates(fit, dd, "treated:immunisation_intensity")
}

estimate_epidemic_clock <- function(dat) {
  t0 <- simulation_constants$week_pivot
  pre <- dat[treated == 0L]

  template <- pre[
    ,
    .(incidence = sum(y) / sum(population)),
    by = week
  ][order(week)]

  shift_grid <- seq(-5, 5, by = 1)
  speed_grid <- c(0.90, 0.95, 1.00, 1.05, 1.10)

  season_clock <- data.table(
    season = sort(unique(dat$season)),
    shift_hat = 0,
    log_speed_hat = 0
  )

  for (ss in season_clock$season) {
    dd <- dat[
      season == ss,
      .(incidence = sum(y) / sum(population)),
      by = week
    ][order(week)]

    grid <- CJ(shift = shift_grid, speed = speed_grid)

    grid[, score := vapply(seq_len(.N), function(i) {
      tau_candidate <- t0 + (dd$week - t0 - grid$shift[i]) * grid$speed[i]
      template_prediction <- approx(
        x = template$week,
        y = log1p(template$incidence),
        xout = tau_candidate,
        rule = 2
      )$y

      suppressWarnings(
        cor(
          log1p(dd$incidence),
          template_prediction,
          use = "complete.obs"
        )
      )
    }, numeric(1))]

    grid[!is.finite(score), score := -Inf]
    best <- grid[which.max(score)]

    season_clock[
      season == ss,
      `:=`(
        shift_hat = best$shift,
        log_speed_hat = log(best$speed)
      )
    ]
  }

  season_clock
}

fit_aligned_nb <- function(dat) {
  t0 <- simulation_constants$week_pivot
  clock <- estimate_epidemic_clock(dat)
  dd <- merge(dat, clock, by = "season", all.x = TRUE, sort = FALSE)
  dd[, tau_hat := t0 + (week - t0 - shift_hat) * exp(log_speed_hat)]

  fit <- try(
    MASS::glm.nb(
      y ~ offset(log(population)) +
        sin(2 * pi * tau_hat / 52) +
        cos(2 * pi * tau_hat / 52) +
        sin(4 * pi * tau_hat / 52) +
        cos(4 * pi * tau_hat / 52) +
        factor(season) +
        factor(area) +
        lag_log_y +
        ses +
        treated:immunisation_intensity,
      data = dd,
      control = glm.control(maxit = 80)
    ),
    silent = TRUE
  )

  extract_estimates(fit, dd, "treated:immunisation_intensity")
}

estimates_to_rows <- function(estimates, truth, base) {
  rbindlist(lapply(names(estimates), function(estimand) {
    z <- estimates[[estimand]]
    data.table(
      base,
      estimand = estimand,
      true = truth[[estimand]],
      estimate = unname(z["estimate"]),
      lower = unname(z["lower"]),
      upper = unname(z["upper"]),
      se = unname(z["se"])
    )
  }))
}

run_fast_rep <- function(i, sc) {
  scenario_index <- match(sc$scenario, scenarios$scenario)
  dat <- simulate_one(
    base_seed + 100000L * scenario_index + i,
    sc
  )
  truth <- true_values(dat, sc)

  methods <- list(
    ITS = fit_its,
    GAM = fit_gam,
    `Aligned NB` = fit_aligned_nb
  )

  rbindlist(lapply(names(methods), function(method_name) {
    est <- methods[[method_name]](dat)
    estimates_to_rows(
      est, truth,
      data.table(
        rep = i,
        scenario = sc$scenario,
        method = method_name,
        component = "large_scale_comparison",
        interval_type = "Wald",
        schedule = NA_character_
      )
    )
  }))
}

n_fast_cores <- min(n_fast_cores, get_env_int("BEACON_FAST_MAX_CORES", 8L))
if (skip_fast || collect_only) n_fast_cores <- 1L
if (!(skip_fast || collect_only)) message(
  "Running large-scale comparison: ", n_fast_cores, " worker(s), ",
  n_rep_fast, " repetitions x ", nrow(scenarios), " scenarios."
)

single_thread <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  data.table::setDTthreads(1L)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  invisible(NULL)
}

run_fast_rep_safe <- function(i, sc) {
  single_thread()
  tryCatch(run_fast_rep(i, sc), error = function(e) data.table(
    rep = i, scenario = sc$scenario, method = NA_character_,
    component = "large_scale_comparison", error = conditionMessage(e)))
}

if (n_fast_cores > 1L) {
  future::plan(future::multisession, workers = n_fast_cores)
} else {
  single_thread()
  future::plan(future::sequential)
}

fast_chunks <- vector("list", nrow(scenarios))
for (j in seq_len(nrow(scenarios))) {
  sc <- scenarios[j]
  chk <- sprintf("results/simulation_fast_raw_v3_%s.csv", sc$scenario)
  if (skip_fast || collect_only) {
    if (file.exists(chk)) fast_chunks[[j]] <- fread(chk)
    next
  }
  if (resume_stan && file.exists(chk)) {
    fast_chunks[[j]] <- fread(chk)
    if (uniqueN(fast_chunks[[j]]$rep) >= n_rep_fast) {
      message("  ", sc$scenario_label, ": loaded from checkpoint.")
      next
    }
  }
  t0 <- Sys.time()
  message("  ", sc$scenario_label, ": started ", format(t0, "%H:%M:%S"))
  ## Blocks of 50 repetitions, so progress is reported while it runs.
  blocks <- split(seq_len(n_rep_fast), ceiling(seq_len(n_rep_fast) / 50))
  res <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    res[[b]] <- rbindlist(
      future_lapply(blocks[[b]], run_fast_rep_safe, sc = sc,
                    future.seed = TRUE),
      fill = TRUE
    )
    message(sprintf("    %d/%d repetitions, %.1f min elapsed",
                    max(blocks[[b]]), n_rep_fast,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  fast_chunks[[j]] <- rbindlist(res, fill = TRUE)
  fwrite(fast_chunks[[j]], chk)
}
future::plan(future::sequential)

fast_raw <- rbindlist(fast_chunks, fill = TRUE)
if (!nrow(fast_raw)) fast_raw <- data.table(estimate = numeric(), lower = numeric(),
                                            upper = numeric(), true = numeric())
if ("error" %in% names(fast_raw) && any(!is.na(fast_raw$error))) {
  warning(sum(!is.na(fast_raw$error)), " comparator repetitions failed; see ",
          "the 'error' column of results/simulation_fast_raw_v3.csv.")
  fast_raw <- fast_raw[is.na(error)]
}

if (write_summaries && nrow(fast_raw)) fwrite(fast_raw, "results/simulation_fast_raw_v3.csv")

## ---------------------------------------------------------------------------
## Bayesian tier
## ---------------------------------------------------------------------------

## Which Bayesian methods are fitted in which scenarios.
stan_methods <- list(
  list(method = "BEACON", clock = TRUE, ascertainment_sd = 0,
       scenarios = scenarios$scenario, long = TRUE),
  list(method = "BEACON-calendar", clock = FALSE, ascertainment_sd = 0,
       scenarios = scenarios$scenario, long = FALSE),
  list(method = "BEACON + bias analysis", clock = TRUE,
       ascertainment_sd = bias_analysis_sd,
       scenarios = "phase_and_ascertainment", long = FALSE)
)
requested_methods <- trimws(strsplit(
  Sys.getenv("BEACON_STAN_METHODS",
             paste(vapply(stan_methods, `[[`, "", "method"), collapse = ",")),
  ",", fixed = TRUE)[[1]])
stan_methods <- Filter(function(m) m$method %in% requested_methods, stan_methods)

stan_raw <- data.table()
stan_fit_raw <- data.table()
stan_long <- data.table()

if (!skip_stan) {

  required_stan_packages <- c("cmdstanr", "posterior")
  missing_stan_packages <- required_stan_packages[
    !vapply(required_stan_packages, requireNamespace, quietly = TRUE,
            FUN.VALUE = logical(1))
  ]
  if (length(missing_stan_packages)) {
    stop("Install the following packages before running the BEACON ",
         "validation: ", paste(missing_stan_packages, collapse = ", "))
  }

  stan_file <- Sys.getenv("BEACON_STAN_FILE", "")
  if (!nzchar(stan_file)) {
    candidates <- c(
      "../Stan/beacon_affine_v2.stan",
      "Stan/beacon_affine_v2.stan",
      "beacon_affine_v2.stan"
    )
    stan_file <- candidates[file.exists(candidates)][1]
    if (is.na(stan_file)) {
      stop("beacon_affine_v2.stan not found; set BEACON_STAN_FILE.")
    }
  }

  assert_beacon_v2(stan_file)

  cmdstan_path_ok <- tryCatch(
    {
      path <- cmdstanr::cmdstan_path()
      nzchar(path) && dir.exists(path)
    },
    error = function(e) FALSE
  )
  if (!cmdstan_path_ok) {
    stop(
      "CmdStan is not configured. Run cmdstanr::install_cmdstan() or ",
      "cmdstanr::set_cmdstan_path()."
    )
  }

  message("Compiling BEACON model once: ", normalizePath(stan_file))
  beacon_model <- cmdstanr::cmdstan_model(
    stan_file,
    force_recompile = force_recompile,
    quiet = TRUE
  )
 
  beacon_exe <- beacon_model$exe_file()
  cmdstan_dir <- cmdstanr::cmdstan_path()
  get_beacon_model <- function() {
    if (!exists(".beacon_model", envir = globalenv())) {
      cmdstanr::set_cmdstan_path(cmdstan_dir)
      assign(".beacon_model", cmdstanr::cmdstan_model(exe_file = beacon_exe),
             envir = globalenv())
    }
    get(".beacon_model", envir = globalenv())
  }

  make_stan_data <- function(dat, clock = TRUE, ascertainment_sd = 0) {
    dat <- copy(dat)
    dat[, area_id := match(area, sort(unique(area)))]
    dat[, season_id := match(season, sort(unique(season)))]
    dat[, week_in_season := as.integer(week)]

    sd <- build_beacon_stan_data(
      dat,
      H = stan_H,
      harmonic_decay = stan_harmonic_decay,
      prior_mean_beta_intensity = stan_prior_mean_beta,
      prior_sd_beta_intensity = stan_prior_sd_beta,
      prior_sd_ascertainment = ascertainment_sd,
      compute_log_lik = 0L,
      compute_y_rep = 0L,
      dynamic_gformula = 0L
    )$stan_data
    if (!isTRUE(all.equal(sd$week_pivot, simulation_constants$week_pivot))) {
      stop("Fitted pivot differs from the generating pivot.")
    }
    if (!clock) sd <- beacon_clock_off(sd)
    sd
  }

  safe_diagnostic_value <- function(diagnostics, candidate_names) {
    candidate <- intersect(candidate_names, names(diagnostics))[1]
    if (is.na(candidate)) return(NA_real_)
    sum(diagnostics[[candidate]], na.rm = TRUE)
  }

  interval_row <- function(x) {
    c(estimate = mean(x), lower = unname(quantile(x, 0.025)),
      upper = unname(quantile(x, 0.975)), se = sd(x))
  }

  rhat_variables <- c("beta_intensity", "beta_intensity_cos",
                      "beta_intensity_sin", "sigma_area", "sigma_season",
                      "sigma_series", "sigma_shift", "sigma_speed", "phi")

  fit_beacon_stan <- function(dat, seed, spec,
                              warmup = stan_warmup,
                              sampling = stan_sampling,
                              adapt_delta = stan_adapt_delta) {
    stan_data <- make_stan_data(dat, spec$clock, spec$ascertainment_sd)
    init_list <- make_beacon_inits_v2(stan_data, stan_chains, seed + 900000L)
    n_transitions <- stan_chains * sampling

    fit <- try(
      get_beacon_model()$sample(
        data = stan_data,
        seed = seed,
        chains = stan_chains,
        parallel_chains = stan_parallel_chains,
        init = init_list,
        iter_warmup = warmup,
        iter_sampling = sampling,
        adapt_delta = adapt_delta,
        max_treedepth = stan_max_treedepth,
        refresh = stan_refresh,
        save_warmup = FALSE,
        show_messages = FALSE,
        output_dir = tempdir()
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(list(estimates = NULL, diagnostics = data.table(
        n_transitions = n_transitions, error = as.character(fit))))
    }

    output_files <- try(fit$output_files(), silent = TRUE)
    if (!inherits(output_files, "try-error")) {
      on.exit(unlink(output_files, force = TRUE), add = TRUE)
    }

    has_pf <- fit_has_variable(fit, "prevented_fraction")
    variables <- c("beta_intensity", "beta_intensity_cos", "beta_intensity_sin",
                   "sigma_series", "phi", "shift", "log_speed",
                   if (has_pf) "prevented_fraction")
    draws <- posterior::as_draws_df(fit$draws(variables = variables))

    modulation_draws <- sqrt(
      draws[["beta_intensity_cos"]]^2 + draws[["beta_intensity_sin"]]^2
    )

    na4 <- c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, se = NA_real_)
    estimates <- list(
      beta0 = interval_row(draws[["beta_intensity"]]),
      modulation_amplitude = interval_row(modulation_draws),
      prevented_fraction = if (has_pf) interval_row(draws[["prevented_fraction"]]) else na4
    )

    summ <- fit$summary(variables = intersect(
      rhat_variables, fit$metadata()$stan_variables))
    rhat_beta <- summ$rhat[summ$variable == "beta_intensity"]
    ess_beta <- summ$ess_bulk[summ$variable == "beta_intensity"]
    ## In the calendar ablation the clock scales are sampled from their
    ## priors only (use_clock = 0); exclude them from its R-hat summary.
    rhat_set <- if (spec$clock) summ$rhat else
      summ$rhat[!summ$variable %in% c("sigma_shift", "sigma_speed")]

    shift_names <- grep("^shift\\[", names(draws), value = TRUE)
    speed_names <- grep("^log_speed\\[", names(draws), value = TRUE)
    draws_plain <- as.data.frame(draws)
    shift_estimate <- colMeans(as.matrix(draws_plain[, shift_names, drop = FALSE]))
    speed_estimate <- colMeans(as.matrix(draws_plain[, speed_names, drop = FALSE]))
    shift_true <- dat[, shift_true[1], by = season][order(season)]$V1
    speed_true <- dat[, log_speed_true[1], by = season][order(season)]$V1

    diagnostics <- try(fit$diagnostic_summary(), silent = TRUE)
    if (inherits(diagnostics, "try-error")) {
      divergences <- NA_integer_
      max_depth <- NA_integer_
    } else {
      divergences <- as.integer(safe_diagnostic_value(
        diagnostics, c("num_divergent", "num_divergences", "divergent__")))
      max_depth <- as.integer(safe_diagnostic_value(
        diagnostics, c("num_max_treedepth", "num_max_treedepth_hits")))
    }

    list(
      estimates = estimates,
      diagnostics = data.table(
        rhat_beta = if (length(rhat_beta)) rhat_beta else NA_real_,
        rhat_max = max(rhat_set, na.rm = TRUE),
        ess_bulk = if (length(ess_beta)) ess_beta else NA_real_,
        divergences = divergences,
        n_transitions = n_transitions,
        max_treedepth_hits = max_depth,
        shift_rmse = sqrt(mean((shift_estimate - shift_true)^2)),
        speed_rmse = sqrt(mean((speed_estimate - speed_true)^2)),
        phi_estimate = mean(draws[["phi"]]),
        beta_cos_estimate = mean(draws[["beta_intensity_cos"]]),
        beta_sin_estimate = mean(draws[["beta_intensity_sin"]]),
        sigma_series_estimate = mean(draws[["sigma_series"]]),
        error = NA_character_
      )
    )
  }

  fit_dir <- "results/beacon_fits"
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  fit_file <- function(schedule, method, scenario_name, i) {
    file.path(fit_dir, sprintf("%s__%s__%s__%03d.rds", schedule,
                               gsub("[^A-Za-z0-9]+", "-", method),
                               scenario_name, i))
  }

  run_and_store <- function(scenario_name, spec, i, schedule, ...) {
    sc <- scenarios[scenario == scenario_name]
    scenario_index <- match(scenario_name, scenarios$scenario)
    # Same seed for every Bayesian method, so they see the same datasets.
    seed_i <- base_seed + 5000000L + 100000L * scenario_index + i
    dat <- simulate_one(seed_i, sc)
    truth <- true_values(dat, sc)
    z <- fit_beacon_stan(dat, seed_i, spec, ...)

    base <- data.table(rep = i, scenario = scenario_name, method = spec$method,
                       component = "bayesian_validation",
                       interval_type = "credible", schedule = schedule)
    est_rows <- if (is.null(z$estimates)) data.table() else
      estimates_to_rows(z$estimates, truth, base)
    fit_row <- cbind(base, z$diagnostics,
                     true_beta_cos = truth$beta_cos,
                     true_beta_sin = truth$beta_sin)
    list(est = est_rows, fit = fit_row)
  }

  legacy <- list(short = c("results/simulation_beacon_raw_v3.csv",
                           "results/simulation_beacon_fits_v3.csv"),
                 long = c("results/simulation_beacon_long_raw_v3.csv",
                          "results/simulation_beacon_long_fits_v3.csv"))
  if (resume_stan && !collect_only && !length(list.files(fit_dir))) {
    for (sched in names(legacy)) {
      if (!all(file.exists(legacy[[sched]]))) next
      est_old <- fread(legacy[[sched]][1])
      fit_old <- fread(legacy[[sched]][2])
      for (k in seq_len(nrow(fit_old))) {
        r <- fit_old[k]
        saveRDS(list(est = est_old[scenario == r$scenario & method == r$method &
                                     rep == r$rep],
                     fit = r),
                fit_file(sched, r$method, r$scenario, r$rep))
      }
      message("Imported ", nrow(fit_old), " ", sched, " fits from ",
              legacy[[sched]][2])
    }
  }

  ## Task list: every (schedule, method, scenario, rep) still to be fitted.
  rep_range <- seq_len(n_rep_stan)
  if (nzchar(stan_reps_env)) {
    lims <- as.integer(strsplit(stan_reps_env, "-", fixed = TRUE)[[1]])
    rep_range <- intersect(seq(lims[1], lims[length(lims)]), rep_range)
  }
  tasks <- rbindlist(c(
    lapply(stan_methods, function(spec) CJ(
      schedule = "short", method = spec$method,
      scenario = intersect(requested_scenarios, spec$scenarios),
      rep = rep_range)),
    lapply(Filter(function(m) isTRUE(m$long), stan_methods), function(spec) CJ(
      schedule = "long", method = spec$method, scenario = requested_scenarios,
      rep = intersect(rep_range, seq_len(min(n_rep_stan_long, n_rep_stan)))))
  ), fill = TRUE)
  if (nrow(tasks)) {
    tasks[, file := fit_file(schedule, method, scenario, rep)]
    tasks <- tasks[!file.exists(file)]
  }
  if (collect_only) tasks <- tasks[0]

  run_task <- function(k) tryCatch(run_task_unsafe(k), error = function(e) {
    message("Task ", k, " failed: ", conditionMessage(e))
    invisible(NULL)
  })
  run_task_unsafe <- function(k) {
    single_thread()
    t <- tasks[k]
    spec <- Filter(function(m) m$method == t$method, stan_methods)[[1]]
    res <- if (t$schedule == "long") {
      run_and_store(t$scenario, spec, t$rep, "long",
                    warmup = stan_warmup_long, sampling = stan_sampling_long,
                    adapt_delta = stan_adapt_delta_long)
    } else {
      run_and_store(t$scenario, spec, t$rep, "short")
    }
    tmp <- paste0(t$file, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, t$file)       
    invisible(t$file)
  }

  message("Running BEACON validation: ", nrow(tasks), " fits to run with ",
          stan_workers, " fit(s) at a time; ", stan_chains, " chains, ",
          stan_warmup, " + ", stan_sampling, " iterations (long: ",
          stan_warmup_long, " + ", stan_sampling_long, "), adapt_delta = ",
          stan_adapt_delta, ".")
  if (nrow(tasks)) {
    if (stan_workers > 1L) {
      future::plan(future::multisession, workers = stan_workers)
    } else {
      future::plan(future::sequential)
    }
    t0 <- Sys.time()
    blocks <- split(seq_len(nrow(tasks)),
                    ceiling(seq_len(nrow(tasks)) / max(10L, stan_workers)))
    for (b in seq_along(blocks)) {
      future_lapply(blocks[[b]], run_task, future.seed = TRUE)
      message(sprintf("  %d/%d fits done, %.1f h elapsed", max(blocks[[b]]),
                      nrow(tasks),
                      as.numeric(difftime(Sys.time(), t0, units = "hours"))))
    }
    future::plan(future::sequential)
  }

  ## Collect every finished fit.
  collected <- lapply(list.files(fit_dir, pattern = "[.]rds$", full.names = TRUE),
                      readRDS)
  all_est <- rbindlist(lapply(collected, `[[`, "est"), fill = TRUE, use.names = TRUE)
  all_fit <- rbindlist(lapply(collected, `[[`, "fit"), fill = TRUE, use.names = TRUE)
  if (nrow(all_fit)) {
    stan_raw <- all_est[schedule == "short"]
    stan_fit_raw <- all_fit[schedule == "short"]
    stan_long <- all_est[schedule == "long"]
    stan_long_fit <- all_fit[schedule == "long"]
  } else {
    stan_long_fit <- data.table()
  }
  if (write_summaries && nrow(all_fit)) {
    fwrite(stan_raw, "results/simulation_beacon_raw_v3.csv")
    fwrite(stan_fit_raw, "results/simulation_beacon_fits_v3.csv")
    if (nrow(stan_long)) {
      fwrite(stan_long, "results/simulation_beacon_long_raw_v3.csv")
      fwrite(stan_long_fit, "results/simulation_beacon_long_fits_v3.csv")
    }
  }

  long_specs <- Filter(function(m) isTRUE(m$long), stan_methods)
  if (write_summaries && n_rep_stan_long > 0L && length(long_specs)) {
    spec <- long_specs[[1]]
    if (nrow(stan_long)) {
      short <- stan_raw[method == spec$method & estimand == "beta0"]
      long <- stan_long[estimand == "beta0"]
      matched <- merge(
        short[, .(scenario, rep, estimate_short = estimate,
                  width_short = upper - lower,
                  covered_short = lower <= true & upper >= true)],
        long[, .(scenario, rep, estimate_long = estimate,
                 width_long = upper - lower,
                 covered_long = lower <= true & upper >= true)],
        by = c("scenario", "rep")
      )
      schedule_check <- matched[
        is.finite(estimate_short) & is.finite(estimate_long),
        .(
          n_pairs = .N,
          mean_abs_difference = mean(abs(estimate_short - estimate_long)),
          max_abs_difference = max(abs(estimate_short - estimate_long)),
          mean_width_short = mean(width_short),
          mean_width_long = mean(width_long),
          coverage_short = mean(covered_short),
          coverage_long = mean(covered_long)
        ),
        by = scenario
      ]
      long_diag <- stan_long_fit[, .(
        long_total_divergences = sum(divergences, na.rm = TRUE),
        long_pct_fits_divergent = 100 * mean(divergences > 0, na.rm = TRUE),
        long_max_rhat = max(rhat_max, na.rm = TRUE)
      ), by = scenario]
      schedule_check <- merge(schedule_check, long_diag, by = "scenario",
                              all.x = TRUE)
      schedule_check[, scenario_label := scenarios$scenario_label[
        match(scenario, scenarios$scenario)]]
      fwrite(schedule_check, "results/simulation_schedule_check_v3.csv")
    }
  }
}

## ---------------------------------------------------------------------------
## Performance summaries
## ---------------------------------------------------------------------------
if (!write_summaries) {
  message("Array job finished (BEACON_STAN_REPS = ", stan_reps_env, "); ",
          "run once with BEACON_STAN_COLLECT_ONLY=true to write the summaries.")
} else {
raw <- rbindlist(list(fast_raw, stan_raw), fill = TRUE, use.names = TRUE)
raw[, covered := is.finite(lower) & is.finite(upper) &
      lower <= true & upper >= true]
raw[, interval_width := upper - lower]
raw[, scenario_label := factor(
  scenarios$scenario_label[match(scenario, scenarios$scenario)],
  levels = scenario_levels)]
raw[, method := factor(method, levels = method_levels)]

n_requested <- function(method) ifelse(
  method %in% c("ITS", "GAM", "Aligned NB"), n_rep_fast, n_rep_stan)

performance <- raw[
  is.finite(estimate) & is.finite(lower) & is.finite(upper),
  {
    errors <- estimate - true
    empirical_coverage <- mean(covered)
    list(
      n_success = .N,
      bias = mean(errors),
      bias_mcse = if (.N > 1L) sd(errors) / sqrt(.N) else NA_real_,
      rmse = sqrt(mean(errors^2)),
      coverage = empirical_coverage,
      coverage_mcse = sqrt(empirical_coverage * (1 - empirical_coverage) / .N),
      mean_interval_width = mean(interval_width),
      mean_estimate = mean(estimate),
      mean_true = mean(true)
    )
  },
  by = .(scenario, scenario_label, method, estimand, component, interval_type)
]
performance[, n_rep_requested := n_requested(as.character(method))]
performance[, n_failed := n_rep_requested - n_success]
performance[estimand == "modulation_amplitude" & mean_true == 0,
            `:=`(coverage = NA_real_, coverage_mcse = NA_real_)]
setorder(performance, scenario_label, estimand, method)
fwrite(performance, "results/simulation_performance_v3.csv")

if (nrow(stan_fit_raw)) {
  stan_fit_raw[, scenario_label := factor(
    scenarios$scenario_label[match(scenario, scenarios$scenario)],
    levels = scenario_levels)]
  beta_rows <- stan_raw[estimand == "beta0" & schedule == "short",
                        .(scenario, method, rep, estimate, true)]
  recovery_input <- merge(stan_fit_raw[schedule == "short"], beta_rows,
                          by = c("scenario", "method", "rep"), all.x = TRUE)
  beacon_parameter_recovery <- recovery_input[
    ,
    .(
      n_fits = .N,
      n_success = sum(is.finite(estimate)),
      beta_bias = mean(estimate - true, na.rm = TRUE),
      beta_rmse = sqrt(mean((estimate - true)^2, na.rm = TRUE)),
      beta_cos_bias = mean(beta_cos_estimate - true_beta_cos, na.rm = TRUE),
      beta_sin_bias = mean(beta_sin_estimate - true_beta_sin, na.rm = TRUE),
      phi_bias = mean(phi_estimate - simulation_constants$phi, na.rm = TRUE),
      sigma_area_season_bias = mean(
        sigma_series_estimate - simulation_constants$sigma_area_season,
        na.rm = TRUE),
      shift_rmse = mean(shift_rmse, na.rm = TRUE),
      speed_rmse = mean(speed_rmse, na.rm = TRUE),
      mean_rhat_max = mean(rhat_max, na.rm = TRUE),
      max_rhat = max(rhat_max, na.rm = TRUE),
      pct_fits_rhat_gt_1.01 = 100 * mean(rhat_max > 1.01, na.rm = TRUE),
      median_ess_bulk = median(ess_bulk, na.rm = TRUE),
      pct_fits_divergent = 100 * mean(divergences > 0, na.rm = TRUE),
      total_divergences = sum(divergences, na.rm = TRUE),
      total_transitions = sum(n_transitions, na.rm = TRUE),
      pct_divergent_transitions = 100 * sum(divergences, na.rm = TRUE) /
        sum(n_transitions, na.rm = TRUE),
      total_max_treedepth_hits = sum(max_treedepth_hits, na.rm = TRUE),
      n_excluded = 0L
    ),
    by = .(scenario, scenario_label, method)
  ]
  setorder(beacon_parameter_recovery, scenario_label, method)
  fwrite(beacon_parameter_recovery,
         "results/simulation_beacon_parameter_recovery_v3.csv")
}

example <- simulate_one(base_seed + 9917L, scenarios[scenario == "phase_shift"])
example_summary <- example[
  ,
  .(
    observed = sum(y),
    population = sum(population),
    immunisation_intensity = mean(immunisation_intensity),
    true_shift = shift_true[1],
    true_log_speed = log_speed_true[1]
  ),
  by = .(season, week, treated)
]
fwrite(example_summary, "results/simulation_example_v3.csv")

settings <- data.table(
  setting = c(
    "fast_repetitions_per_scenario", "stan_repetitions_per_scenario",
    "stan_chains", "stan_warmup", "stan_sampling", "stan_adapt_delta",
    "stan_max_treedepth", "stan_long_repetitions_per_scenario",
    "stan_long_warmup", "stan_long_sampling", "stan_long_adapt_delta",
    "fitted_harmonics", "fitted_harmonic_decay", "prior_mean_beta_intensity",
    "prior_sd_beta_intensity", "bias_analysis_sd",
    "week_pivot", "tau_star", "sigma_area_season_true", "areas", "seasons",
    "weeks_per_season"
  ),
  value = c(
    n_rep_fast, n_rep_stan, stan_chains, stan_warmup, stan_sampling,
    stan_adapt_delta, stan_max_treedepth, n_rep_stan_long, stan_warmup_long,
    stan_sampling_long, stan_adapt_delta_long, stan_H, stan_harmonic_decay,
    stan_prior_mean_beta, stan_prior_sd_beta, bias_analysis_sd,
    simulation_constants$week_pivot, tau_star,
    simulation_constants$sigma_area_season, simulation_constants$n_area,
    simulation_constants$n_season, simulation_constants$n_week
  )
)
fwrite(settings, "results/simulation_settings_v3.csv")
fwrite(scenarios, "results/simulation_scenarios_v3.csv")
} 

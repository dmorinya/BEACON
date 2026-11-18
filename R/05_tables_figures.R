library(data.table)
library(ggplot2)

dir.create("results", recursive = TRUE, showWarnings = FALSE)
fig_dir <- Sys.getenv("BEACON_FIG_DIR", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- c("ITS", "GAM", "Aligned NB", "BEACON-calendar", "BEACON",
                   "BEACON + bias analysis")
fmt <- function(x, d = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "--")
fmt_pct <- function(x, d = 1) ifelse(is.finite(x), formatC(100 * x, format = "f", digits = d), "--")
write_rows <- function(lines, path) {
  writeLines(lines, path)
}

theme_beacon <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0))

## ---------------------------------------------------------------------------
## Simulation
## ---------------------------------------------------------------------------
perf_path <- "results/simulation_performance_v3.csv"
if (file.exists(perf_path)) {
  perf <- fread(perf_path)
  scen <- fread("results/simulation_scenarios_v3.csv")
  perf[, scenario_label := factor(scenario_label, levels = scen$scenario_label)]
  perf[, method := factor(method, levels = method_levels)]
  setorder(perf, scenario_label, method)

  bias_plot <- function(est, xlab, file, cap) {
    p <- perf[estimand == est & is.finite(bias)]
    if (!nrow(p)) return(invisible(NULL))
    p[, off := abs(bias) > cap]
    p[, x := pmax(pmin(bias, cap), -cap)]
    p[, lo := pmax(bias - 2 * bias_mcse, -cap)]
    p[, hi := pmin(bias + 2 * bias_mcse, cap)]
    g <- ggplot(p, aes(y = method)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      geom_errorbarh(data = p[off == FALSE], aes(xmin = lo, xmax = hi), height = 0) +
      geom_point(data = p[off == FALSE], aes(x = x), size = 1.8) +
      geom_point(data = p[off == TRUE], aes(x = x), shape = 17, size = 2.2) +
      geom_text(data = p[off == TRUE], aes(x = x, label = sprintf("%.2f", bias)),
                hjust = ifelse(p[off == TRUE]$bias > 0, 1.3, -0.3), size = 2.6) +
      scale_y_discrete(limits = rev(intersect(method_levels, unique(as.character(p$method))))) +
      coord_cartesian(xlim = c(-cap, cap)) +
      facet_wrap(~ scenario_label, ncol = 2) +
      labs(x = xlab, y = NULL) +
      theme_beacon
    ggsave(file.path(fig_dir, file), g, width = 7, height = 7.5)
  }
  bias_plot("beta0", "Bias in the log rate ratio", "fig-simulation-bias-1.pdf",
            cap = as.numeric(Sys.getenv("BEACON_FIG_BIAS_CAP", "0.30")))
  bias_plot("prevented_fraction", "Bias in the prevented fraction",
            "fig-simulation-pf-bias-1.pdf",
            cap = as.numeric(Sys.getenv("BEACON_FIG_PF_BIAS_CAP", "0.20")))

  ## Table 2 body: scenario, method, fits, bias (MCSE), RMSE, coverage (MCSE), width
  t2 <- perf[estimand == "beta0"]
  lines <- character()
  for (sl in levels(t2$scenario_label)) {
    block <- t2[scenario_label == sl]
    if (!nrow(block)) next
    for (k in seq_len(nrow(block))) {
      b <- block[k]
      lines <- c(lines, sprintf(
        "%s & %s & %d/%d & %s (%s) & %s & %s (%s) & %s\\\\",
        if (k == 1) sl else "", as.character(b$method), b$n_success,
        b$n_rep_requested, fmt(b$bias), fmt(b$bias_mcse), fmt(b$rmse),
        fmt_pct(b$coverage), fmt_pct(b$coverage_mcse), fmt(b$mean_interval_width)))
    }
    lines <- c(lines, "\\addlinespace")
  }
  write_rows(head(lines, -1), "results/tex_table_simulation_rows.tex")

  ## Prevented fraction and modulation amplitude
  t_pf <- perf[estimand %in% c("prevented_fraction", "modulation_amplitude")]
  fwrite(t_pf, "results/simulation_pf_modulation_summary_v3.csv")
}

rec_path <- "results/simulation_beacon_parameter_recovery_v3.csv"
if (file.exists(rec_path)) {
  rec <- fread(rec_path)
  scen <- fread("results/simulation_scenarios_v3.csv")
  rec[, scenario_label := factor(scenario_label, levels = scen$scenario_label)]
  rec[, method := factor(method, levels = method_levels)]
  setorder(rec, scenario_label, method)
  r3 <- rec[method == "BEACON"]
  lines <- sprintf(
    "%s & %d & %s (%s) & %s & %s & %s / %s & %s / %s & %s & %.0f & %s & %d / %d\\\\",
    as.character(r3$scenario_label), r3$n_fits, fmt(r3$beta_bias), fmt(r3$beta_rmse),
    fmt(r3$beta_cos_bias), fmt(r3$phi_bias, 2), fmt(r3$shift_rmse, 2),
    fmt(r3$speed_rmse, 3), fmt(r3$mean_rhat_max), fmt(r3$max_rhat),
    fmt(r3$pct_fits_rhat_gt_1.01, 1), r3$median_ess_bulk,
    fmt(r3$pct_fits_divergent, 1), as.integer(r3$total_divergences),
    as.integer(r3$total_max_treedepth_hits))
  write_rows(lines, "results/tex_table_recovery_rows.tex")
}

sched_path <- "results/simulation_schedule_check_v3.csv"
if (file.exists(sched_path)) {
  sc <- fread(sched_path)
  scen <- fread("results/simulation_scenarios_v3.csv")
  sc[, scenario_label := factor(scenario_label, levels = scen$scenario_label)]
  setorder(sc, scenario_label)
  lines <- sprintf(
    "%s & %d & %s & %s & %s / %s & %s / %s & %d & %s\\\\",
    as.character(sc$scenario_label), sc$n_pairs, fmt(sc$mean_abs_difference),
    fmt(sc$max_abs_difference), fmt(sc$mean_width_short), fmt(sc$mean_width_long),
    fmt_pct(sc$coverage_short), fmt_pct(sc$coverage_long),
    as.integer(sc$long_total_divergences), fmt(sc$long_max_rhat))
  write_rows(lines, "results/tex_table_s1_schedule_rows.tex")
}

## ---------------------------------------------------------------------------
## Catalan data
## ---------------------------------------------------------------------------
act_path <- "results/catalonia_immunisation_intensity_by_abs.csv"
if (file.exists(act_path)) {
  act <- fread(act_path)[treated == 1 & is.finite(immunisation_intensity)]
  analysed_only <- tolower(Sys.getenv("BEACON_FIG_ANALYSED_ABS_ONLY", "false")) %in%
    c("1", "true", "yes")
  if (analysed_only && file.exists("data/derived/catalonia_beacon_analysis.csv")) {
    an <- fread("data/derived/catalonia_beacon_analysis.csv")
    sz <- an[, .(p = median(population, na.rm = TRUE)), by = abs_code][order(-p)]
    act <- act[abs_code %in% head(sz$abs_code,
                                  as.integer(Sys.getenv("BEACON_MAX_AREAS", "60")))]
  }
  g <- ggplot(act, aes(x = season, y = immunisation_intensity)) +
    geom_boxplot(outlier.size = 0.6, width = 0.5) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    labs(x = "Programme season", y = "Scaled RSV immunisation activity") +
    theme_beacon
  ggsave(file.path(fig_dir, "fig-activity-1.pdf"), g, width = 6, height = 4)
}

## Figure 2: weekly incidence in programme seasons against the distribution
## of complete pre-programme seasons (10th-90th percentile ribbon and median),
## excluding the incomplete 2011-2012 and the pandemic 2020-2021 seasons.
prof_path <- "results/catalonia_descriptive_weekly.csv"
if (file.exists(prof_path)) {
  d <- fread(prof_path)[week_in_season >= 1 & week_in_season <= 52]
  excl <- trimws(strsplit(Sys.getenv("BEACON_PROFILE_EXCLUDE",
                                     "2011-2012,2020-2021"), ",")[[1]])
  pre <- d[treated == 0 & !season %in% excl]
  band <- pre[, .(lo = quantile(incidence_per_100k, 0.1),
                  hi = quantile(incidence_per_100k, 0.9),
                  med = median(incidence_per_100k)), by = week_in_season]
  prog <- d[treated == 1]
  g <- ggplot() +
    geom_ribbon(data = band, aes(week_in_season, ymin = lo, ymax = hi), fill = "grey80") +
    geom_line(data = band, aes(week_in_season, med), colour = "grey40", linetype = 2) +
    geom_line(data = prog, aes(week_in_season, incidence_per_100k, colour = season),
              linewidth = 0.7) +
    scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                        name = "Programme season") +
    labs(x = "Week of season (from 1 July)", y = "Weekly incidence per 100,000 infants") +
    theme_beacon + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "fig-epidemic-profile-1.pdf"), g, width = 6.5, height = 4)
}

clock_path <- "results/catalonia_season_clocks.csv"
if (file.exists(clock_path)) {
  cl <- fread(clock_path)
  cl[, parameter := factor(parameter, levels = c("shift_weeks", "speed_ratio"),
                           labels = c("Shift (weeks)", "Speed relative to average"))]
  ref <- data.table(parameter = levels(cl$parameter), ref = c(0, 1))
  ref[, parameter := factor(parameter, levels = levels(cl$parameter))]
  g <- ggplot(cl, aes(x = season, y = mean, colour = programme)) +
    geom_hline(data = ref, aes(yintercept = ref), linetype = 2, colour = "grey50") +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.3) +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "firebrick"),
                        labels = c("Pre-programme", "Programme"), name = NULL) +
    facet_wrap(~ parameter, ncol = 1, scales = "free_y") +
    labs(x = "Season", y = "Posterior mean and 95% CrI") +
    theme_beacon +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ggsave(file.path(fig_dir, "fig-season-clocks-1.pdf"), g, width = 6.5, height = 5.5)
}

td_path <- "results/catalonia_timing_sd_draws.csv"
if (file.exists(td_path)) {
  td <- fread(td_path)
  pri <- fread("results/catalonia_timing_priors.csv")
  long <- melt(td, measure.vars = c("sigma_shift", "sigma_speed"),
               variable.name = "parameter", value.name = "value")
  dens_prior <- rbindlist(lapply(seq_len(nrow(pri)), function(i) {
    par <- pri$parameter[i]
    upper <- max(4 * pri$prior_sd[i], max(long[parameter == par, value]))
    x <- seq(0, upper, length.out = 400)
    data.table(parameter = par, value = x,
               density = 2 * dnorm(x, 0, pri$prior_sd[i]))
  }))
  labs_map <- c(sigma_shift = "sigma[delta]~(weeks)", sigma_speed = "sigma[kappa]")
  long[, parameter := factor(labs_map[as.character(parameter)], levels = labs_map)]
  dens_prior[, parameter := factor(labs_map[parameter], levels = labs_map)]
  g <- ggplot(long, aes(x = value)) +
    geom_density(aes(linetype = "Posterior"), fill = "grey80", colour = "black") +
    geom_line(data = dens_prior, aes(y = density, linetype = "Prior")) +
    scale_linetype_manual(values = c(Posterior = 1, Prior = 2), name = NULL) +
    facet_wrap(~ parameter, scales = "free", labeller = label_parsed) +
    labs(x = NULL, y = "Density") +
    theme_beacon + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "fig-timing-prior-posterior-1.pdf"), g, width = 6.5, height = 3.5)
}

ppc_path <- "results/catalonia_ppc_weekly.csv"
if (file.exists(ppc_path)) {
  pp <- fread(ppc_path)
  cal_path <- "results/catalonia_ppc_weekly_calendar.csv"
  if (file.exists(cal_path)) pp <- rbind(pp, fread(cal_path), fill = TRUE)
  pp[, model := factor(ifelse(analysis == "calendar", "BEACON-calendar", "BEACON"),
                       levels = c("BEACON", "BEACON-calendar"))]
  obs <- unique(pp[, .(season, week_in_season, observed)])
  g <- ggplot(pp, aes(x = week_in_season)) +
    geom_ribbon(aes(ymin = rep_lower, ymax = rep_upper, fill = model), alpha = 0.35) +
    geom_line(aes(y = rep_median, colour = model)) +
    geom_point(data = obs, aes(y = observed), size = 0.4) +
    scale_fill_manual(values = c(BEACON = "grey40", `BEACON-calendar` = "firebrick"),
                      name = NULL) +
    scale_colour_manual(values = c(BEACON = "grey20", `BEACON-calendar` = "firebrick"),
                        name = NULL) +
    facet_wrap(~ season, ncol = 3, scales = "free_y") +
    labs(x = "Week of season (from 1 July)",
         y = "Weekly bronchiolitis diagnoses (analysed ABS)") +
    theme_beacon + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "fig-ppc-1.pdf"), g, width = 7, height = 8.5)
}

sens_path <- "results/catalonia_sensitivity_summary.csv"
if (file.exists(sens_path)) {
  se <- fread(sens_path)
  pl_path <- "results/catalonia_placebo_summary.csv"
  if (file.exists(pl_path)) se <- rbind(se, fread(pl_path), fill = TRUE)
  wanted <- c("rate_ratio_at_reference_activity", "effect_modulation_amplitude",
              "prevented_fraction", "prevented_fraction_dynamic",
              "prevented_fraction_within_support_p90_vs_p10",
              "beta_prior_posterior_contraction")
  w <- dcast(se[metric %in% wanted], analysis ~ metric,
             value.var = c("median", "lower", "upper", "mean"))
  ## Row order and labels of Supplementary Table S2; the placebo row last.
  labels <- c(primary = "Primary analysis",
              calendar = "Clock switched off (BEACON-calendar)",
              unpooled_programme_seasons = "Unpooled programme-season effects",
              protective_prior = "Protective prior, $\\beta_0\\sim N(\\log 0.6,0.5^2)$",
              `ascertainment_sd_0.1` = "Ascertainment bias analysis, $\\sigma_\\lambda=0.10$",
              two_harmonics = "Two harmonics, no decay",
              `including_2020-2021` = "Including 2020--2021")
  w[, ord := match(analysis, names(labels))]
  w[is.na(ord), ord := 100L]
  setorder(w, ord, analysis)
  w[, label := ifelse(analysis %in% names(labels), labels[analysis],
                      ifelse(grepl("^placebo_", analysis),
                             "In-time placebo (2016--2017 to 2018--2019)",
                             gsub("_", " ", analysis)))]
  ci <- function(m, l, u, d = 2, pct = FALSE) {
    k <- if (pct) 100 else 1
    sprintf("%s (%s, %s)", fmt(k * m, d), fmt(k * l, d), fmt(k * u, d))
  }
  pf9010 <- ifelse(grepl("^placebo_", w$analysis), "--",
                   ci(w$median_prevented_fraction_within_support_p90_vs_p10,
                      w$lower_prevented_fraction_within_support_p90_vs_p10,
                      w$upper_prevented_fraction_within_support_p90_vs_p10, 1, TRUE))
  lines <- sprintf(
    "%s & %s & %s & %s & %s & %s & %s\\\\",
    w$label,
    ci(w$median_rate_ratio_at_reference_activity, w$lower_rate_ratio_at_reference_activity,
       w$upper_rate_ratio_at_reference_activity),
    ci(w$median_effect_modulation_amplitude, w$lower_effect_modulation_amplitude,
       w$upper_effect_modulation_amplitude),
    ci(w$median_prevented_fraction, w$lower_prevented_fraction,
       w$upper_prevented_fraction, 1, TRUE),
    ci(w$median_prevented_fraction_dynamic, w$lower_prevented_fraction_dynamic,
       w$upper_prevented_fraction_dynamic, 1, TRUE),
    pf9010,
    fmt(w$mean_beta_prior_posterior_contraction, 2))
  write_rows(lines, "results/tex_table_s2_sensitivity_rows.tex")
}
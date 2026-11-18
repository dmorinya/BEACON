library(data.table)
library(posterior)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
fig_dir <- Sys.getenv("BEACON_FIG_DIR", "figures")
dir.create(fig_dir, showWarnings = FALSE)

fit_path <- "results/catalonia_beacon_fit.rds"
desc_path <- "results/catalonia_descriptive_weekly.csv"
if (!file.exists(fit_path)) stop("Run 03_fit_catalonia.R first (", fit_path, ").")
if (!file.exists(desc_path)) stop("Run 02_download_prepare_catalonia.R first.")

frac_main <- as.numeric(Sys.getenv("BEACON_WINDOW_FRACTION", "0.5"))
frac_broad <- as.numeric(Sys.getenv("BEACON_WINDOW_FRACTION_BROAD", "0.25"))
excl <- trimws(strsplit(Sys.getenv("BEACON_PROFILE_EXCLUDE", "2011-2012,2020-2021"),
                        ",")[[1]])

## ---- epidemic window from the pre-programme seasons -----------------------
d <- fread(desc_path)[week_in_season >= 1 & week_in_season <= 52]
pre <- d[treated == 0 & !season %in% excl]
avg <- pre[, .(incidence = mean(incidence_per_100k)), by = week_in_season][
  order(week_in_season)]
peak <- max(avg$incidence)

window_for <- function(frac) {
  wpk <- avg[which.max(incidence), week_in_season]
  above <- avg$incidence >= frac * peak
  lo <- wpk; hi <- wpk
  while (lo > 1 && above[lo - 1]) lo <- lo - 1
  while (hi < 52 && above[hi + 1]) hi <- hi + 1
  data.table(fraction_of_peak = frac, first_week = lo, last_week = hi,
             share_of_cases = sum(avg[week_in_season %between% c(lo, hi), incidence]) /
               sum(avg$incidence))
}
windows <- rbind(window_for(frac_main), window_for(frac_broad))
windows[, window := c("epidemic", "broad")]
fwrite(windows, "results/catalonia_epidemic_window.csv")
print(windows)

## ---- rate ratio profile from the posterior draws ---------------------------
fit <- readRDS(fit_path)
dr <- as_draws_df(fit$draws(c("beta_intensity", "beta_intensity_cos",
                              "beta_intensity_sin")))
pivot <- (1 + 52) / 2
weeks <- 1:52
tau <- weeks - pivot
lrr <- outer(dr$beta_intensity, rep(1, 52)) +
  outer(dr$beta_intensity_cos, cos(2 * pi * tau / 52)) +
  outer(dr$beta_intensity_sin, sin(2 * pi * tau / 52))
rr <- exp(lrr)
profile <- data.table(
  week_in_season = weeks,
  rr_median = apply(rr, 2, median),
  rr_lower = apply(rr, 2, quantile, 0.025),
  rr_upper = apply(rr, 2, quantile, 0.975),
  in_epidemic_window = weeks %between% c(windows[window == "epidemic", first_week],
                                         windows[window == "epidemic", last_week]),
  in_broad_window = weeks %between% c(windows[window == "broad", first_week],
                                      windows[window == "broad", last_week])
)
fwrite(profile, "results/catalonia_effect_profile.csv")

summ_window <- function(w) {
  idx <- which(weeks %between% c(windows[window == w, first_week],
                                 windows[window == w, last_week]))
  wmin <- idx[which.min(profile$rr_median[idx])]
  wmax <- idx[which.max(profile$rr_median[idx])]
  data.table(
    window = w,
    weeks = paste0(min(idx), "-", max(idx)),
    rr_min_week = wmin, rr_min = profile$rr_median[wmin],
    rr_min_lower = profile$rr_lower[wmin], rr_min_upper = profile$rr_upper[wmin],
    rr_max_week = wmax, rr_max = profile$rr_median[wmax],
    rr_max_lower = profile$rr_lower[wmax], rr_max_upper = profile$rr_upper[wmax],
    last_week_rr_below_one = suppressWarnings(max(idx[profile$rr_median[idx] < 1]))
  )
}
wsum <- rbind(summ_window("epidemic"), summ_window("broad"))
fwrite(wsum, "results/catalonia_effect_window_summary.csv")
print(wsum)

## ---- figure: profile over the broad window, epidemic window shaded --------
ew <- windows[window == "epidemic"]
pb <- profile[in_broad_window == TRUE]
g <- ggplot(pb, aes(week_in_season)) +
  annotate("rect", xmin = ew$first_week - 0.5, xmax = ew$last_week + 0.5,
           ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = rr_lower, ymax = rr_upper), fill = "grey70", alpha = 0.6) +
  geom_line(aes(y = rr_median)) +
  labs(x = "Week of season (from 1 July), average-timed season",
       y = "Rate ratio at reference activity") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig-effect-profile-1.pdf"), g, width = 6, height = 3.5)

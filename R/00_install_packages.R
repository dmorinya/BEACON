required <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "MASS", "mgcv",
  "splines", "future", "future.apply", "cmdstanr", "posterior", "loo",
  "readr", "stringr", "lubridate", "knitr", "kableExtra", "scales", "curl"
)

installed <- rownames(installed.packages())
missing <- setdiff(required, installed)
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  ok <- tryCatch({
    v <- cmdstanr::cmdstan_version()
    if (utils::compareVersion(v, "2.32.0") < 0) {
      warning("CmdStan ", v, " is too old for beacon_affine_v2.stan. ",
              "Run cmdstanr::install_cmdstan() to update.")
    } else {
      message("CmdStan ", v, " detected.")
    }
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    message("CmdStan is not configured. Run cmdstanr::install_cmdstan().")
  }
}

message("Packages available.")

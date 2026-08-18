helpers_path <- c("helpers.R", "R/helpers.R", "../R/helpers.R")
helpers_path <- helpers_path[file.exists(helpers_path)][1]
if (is.na(helpers_path)) {
  stop("helpers.R not found; it now supplies overlap_diagnostic().")
}
source(helpers_path)

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(curl)
  library(jsonlite)
})

# -----------------------------------------------------------------------------
# Fallback helpers: used only when R/helpers.R is not available.
# -----------------------------------------------------------------------------
if (!exists("clean_names_base", mode = "function")) {
  clean_names_base <- function(x) {
    out <- iconv(x, from = "", to = "ASCII//TRANSLIT")
    out[is.na(out)] <- x[is.na(out)]
    out <- tolower(out)
    out <- gsub("[^a-z0-9]+", "_", out)
    out <- gsub("^_+|_+$", "", out)
    make.unique(out, sep = "_")
  }
}

if (!exists("find_col", mode = "function")) {
  find_col <- function(nms, patterns, required = TRUE, label = "column") {
    for (pattern in patterns) {
      hit <- grep(pattern, nms, ignore.case = TRUE, perl = TRUE, value = TRUE)
      if (length(hit)) return(hit[1L])
    }
    if (required) {
      stop(
        "Could not identify column '", label, "'. Available columns: ",
        paste(nms, collapse = ", ")
      )
    }
    NA_character_
  }
}

if (!exists("parse_date_flex", mode = "function")) {
  parse_date_flex <- function(x) {
    x <- as.character(x)
    out <- suppressWarnings(as.Date(x))
    missing <- is.na(out) & nzchar(x)
    if (any(missing)) {
      parsed <- suppressWarnings(lubridate::parse_date_time(
        x[missing],
        orders = c("Ymd HMS", "Ymd HM", "Ymd", "dmY HMS", "dmY", "mdY"),
        quiet = TRUE,
        tz = "UTC"
      ))
      out[missing] <- as.Date(parsed)
    }
    out
  }
}

if (!exists("epidemic_season", mode = "function")) {
  epidemic_season <- function(x) {
    x <- as.Date(x)
    y <- lubridate::year(x)
    start_y <- ifelse(lubridate::month(x) >= 7L, y, y - 1L)
    paste0(start_y, "-", start_y + 1L)
  }
}

if (!exists("weighted_mean_safe", mode = "function")) {
  weighted_mean_safe <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(x[ok], w[ok])
  }
}

is_true_env <- function(x, default = "false") {
  tolower(trimws(Sys.getenv(x, default))) %in% c("1", "true", "yes", "y")
}

first_year_in_text <- function(x) {
  x <- as.character(x)
  m <- regexpr("20[0-9]{2}", x, perl = TRUE)
  len <- attr(m, "match.length")
  ans <- rep(NA_integer_, length(x))
  ok <- m > 0L & len == 4L
  ans[ok] <- as.integer(substr(x[ok], m[ok], m[ok] + len[ok] - 1L))
  ans
}

# -----------------------------------------------------------------------------
# Socrata/SODA helpers.
# -----------------------------------------------------------------------------
socrata_domain <- Sys.getenv(
  "BEACON_SOCRATA_DOMAIN",
  "analisi.transparenciacatalunya.cat"
)
sivic_id <- Sys.getenv("BEACON_SIVIC_DATASET_ID", "fa7i-d8gc")
immun_id <- Sys.getenv("BEACON_IMMUN_DATASET_ID", "hxw8-ertk")

sivic_file <- "data/raw/sivic_primary_care_bronchiolitis_infants.csv"
immun_file <- "data/raw/sivic_rsv_immunisation_campaigns.csv"

age_pattern <- Sys.getenv(
  "BEACON_INFANT_AGE_PATTERN",
  "(^0$|0 anys|menors? d.?1|0[ -]?11|< ?1)"
)

start_date <- as.Date(Sys.getenv("BEACON_START_DATE", "2011-01-01"))
requested_end_date <- as.Date(Sys.getenv("BEACON_END_DATE", as.character(Sys.Date())))
chunk_size <- as.integer(Sys.getenv("BEACON_SOCRATA_CHUNK_SIZE", "25000"))
force_download <- is_true_env("BEACON_FORCE_DOWNLOAD", "false")

if (!is.finite(chunk_size) || chunk_size < 1000L) chunk_size <- 25000L
if (is.na(start_date) || is.na(requested_end_date) || requested_end_date < start_date) {
  stop("Invalid BEACON_START_DATE or BEACON_END_DATE.")
}

# Retain complete July--June epidemic seasons only. For example, on any date
# after 30 June 2026 the most recent complete season is 2025--2026; the newly
# started 2026--2027 season is excluded automatically.
default_last_complete_start <- lubridate::year(requested_end_date) - 1L
candidate_end <- as.Date(sprintf("%d-06-30", default_last_complete_start + 1L))
if (requested_end_date < candidate_end) {
  default_last_complete_start <- default_last_complete_start - 1L
}
last_complete_season_start <- as.integer(Sys.getenv(
  "BEACON_LAST_COMPLETE_SEASON_START",
  as.character(default_last_complete_start)
))
if (!is.finite(last_complete_season_start)) {
  stop("BEACON_LAST_COMPLETE_SEASON_START must be an integer year.")
}
last_complete_season <- paste0(
  last_complete_season_start, "-", last_complete_season_start + 1L
)
complete_season_end <- as.Date(sprintf(
  "%d-06-30", last_complete_season_start + 1L
))
end_date <- min(requested_end_date, complete_season_end)
message(
  "Last complete epidemic season retained: ", last_complete_season,
  " (data through ", end_date, ")."
)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

encode_query <- function(params) {
  # Some optional Socrata parameters, such as $where, are deliberately NULL.
  # Filter them before coercion: as.character(NULL) has length zero and causes
  # vapply(..., character(1)) to fail.
  if (!length(params)) return("")

  keep <- vapply(
    params,
    function(x) {
      !is.null(x) && length(x) > 0L && !all(is.na(x)) &&
        nzchar(paste(as.character(x), collapse = ","))
    },
    logical(1)
  )
  params <- params[keep]
  if (!length(params)) return("")

  paste(
    vapply(seq_along(params), function(i) {
      nm <- names(params)[i]
      value <- paste(as.character(params[[i]]), collapse = ",")
      paste0(
        utils::URLencode(nm, reserved = TRUE), "=",
        utils::URLencode(value, reserved = TRUE)
      )
    }, character(1)),
    collapse = "&"
  )
}

soda_url <- function(domain, dataset_id, params = list(), format = c("csv", "json")) {
  format <- match.arg(format)
  base <- sprintf("https://%s/resource/%s.%s", domain, dataset_id, format)
  query <- encode_query(params)
  if (nzchar(query)) paste0(base, "?", query) else base
}

fetch_to_file <- function(url, extension, timeout = 300L) {
  tmp <- tempfile(fileext = extension)
  handle <- curl::new_handle(
    useragent = "BEACON-Catalonia/1.2",
    followlocation = TRUE,
    connecttimeout = 30,
    timeout = timeout,
    failonerror = FALSE
  )
  response <- curl::curl_fetch_disk(url, tmp, handle = handle)
  if (response$status_code < 200L || response$status_code >= 300L) {
    detail <- tryCatch(
      paste(readLines(tmp, n = 12L, warn = FALSE), collapse = " "),
      error = function(e) ""
    )
    unlink(tmp)
    stop(
      "HTTP ", response$status_code, " returned by the Socrata API.",
      if (nzchar(detail)) paste0("\n", substr(detail, 1L, 1200L)) else "",
      "\nURL: ", url
    )
  }
  tmp
}

fetch_csv_dt <- function(url, timeout = 300L) {
  tmp <- fetch_to_file(url, ".csv", timeout = timeout)
  on.exit(unlink(tmp), add = TRUE)

  first_line <- readLines(tmp, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (!length(first_line) || grepl("^\\s*(<html|<!doctype|\\{\\s*\"error\")", first_line, ignore.case = TRUE)) {
    stop("The Socrata endpoint did not return a CSV file. URL: ", url)
  }

  data.table::fread(
    tmp,
    encoding = "UTF-8",
    showProgress = FALSE,
    fill = FALSE
  )
}

fetch_json_dt <- function(url, timeout = 300L) {
  tmp <- fetch_to_file(url, ".json", timeout = timeout)
  on.exit(unlink(tmp), add = TRUE)

  txt <- paste(readLines(tmp, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!nzchar(trimws(txt)) || grepl("^\\s*<", txt)) {
    stop("The Socrata endpoint did not return JSON. URL: ", url)
  }

  obj <- tryCatch(
    jsonlite::fromJSON(txt, simplifyDataFrame = TRUE),
    error = function(e) stop("Invalid JSON returned by Socrata: ", conditionMessage(e), "\nURL: ", url)
  )
  if (is.list(obj) && !is.data.frame(obj) && !is.null(obj$message)) {
    stop("Socrata API error: ", obj$message, "\nURL: ", url)
  }
  if (is.null(obj) || (is.list(obj) && !length(obj))) return(data.table())
  data.table::as.data.table(obj)
}

fetch_sample <- function(dataset_id) {
  fetch_json_dt(soda_url(
    socrata_domain,
    dataset_id,
    list("$limit" = 1),
    format = "json"
  ))
}

fetch_distinct_values <- function(dataset_id, field, where = NULL, limit = 5000L) {
  ans <- fetch_json_dt(soda_url(
    socrata_domain,
    dataset_id,
    list(
      "$select" = paste0(field, " as value"),
      "$where" = where,
      "$group" = field,
      "$order" = field,
      "$limit" = limit
    ),
    format = "json"
  ))
  if (!nrow(ans)) return(character())
  out <- trimws(as.character(ans[[1L]]))
  out <- gsub('^"+|"+$', '', out)
  unique(out[nzchar(out) & !is.na(out)])
}

split_env_values <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(character())
  out <- trimws(unlist(strsplit(x, "\\s*(?:\\|\\||;|\\n)\\s*", perl = TRUE)))
  unique(out[nzchar(out)])
}

extract_age_numbers <- function(x) {
  hits <- regmatches(x, gregexpr("[0-9]+(?:\\.[0-9]+)?", x, perl = TRUE))[[1L]]
  if (!length(hits) || identical(hits, character(0))) return(numeric())
  suppressWarnings(as.numeric(hits))
}

is_infant_age_level <- function(x, fallback_pattern = NULL) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub('^"+|"+$', '', z)
  if (!nzchar(z) || grepl("no disponible|desconegut|unknown|missing", z)) return(FALSE)

  textual <- grepl(
    "(<\\s*1|menor(?:s)?\\s+(?:de|d')?\\s*1|menys\\s+de\\s+1|0\\s*(?:any|anys|year|years|mes|mesos|month|months))",
    z,
    perl = TRUE
  )
  exact_zero <- grepl("^\\s*0(?:\\.0+)?\\s*$", z, perl = TRUE)
  nums <- extract_age_numbers(z)
  interval_zero_one <- length(nums) >= 2L && is.finite(nums[1L]) && is.finite(nums[2L]) &&
    abs(nums[1L]) < 1e-10 && nums[2L] > 0 && nums[2L] <= 1.01

  fallback <- FALSE
  if (!is.null(fallback_pattern) && nzchar(fallback_pattern)) {
    fallback <- grepl(fallback_pattern, z, ignore.case = TRUE, perl = TRUE)
  }
  isTRUE(textual || exact_zero || interval_zero_one || fallback)
}

choose_infant_levels <- function(levels, label, fallback_pattern = NULL) {
  levels <- unique(trimws(as.character(levels)))
  levels <- levels[nzchar(levels) & !is.na(levels)]

  override <- split_env_values(Sys.getenv("BEACON_INFANT_AGE_VALUES", ""))
  if (length(override)) {
    missing <- setdiff(override, levels)
    if (length(missing)) {
      stop(
        "BEACON_INFANT_AGE_VALUES contains values not present in ", label, ": ",
        paste(dQuote(missing), collapse = ", "),
        ". Available values are: ", paste(dQuote(levels), collapse = ", ")
      )
    }
    return(override)
  }

  keep <- vapply(levels, is_infant_age_level, logical(1), fallback_pattern = fallback_pattern)
  selected <- levels[keep]
  if (!length(selected)) {
    stop(
      "No under-one age interval could be identified in ", label, ". ",
      "Available values are: ", paste(dQuote(levels), collapse = ", "),
      ". Set BEACON_INFANT_AGE_VALUES to the exact value shown by the API, ",
      "for example Sys.setenv(BEACON_INFANT_AGE_VALUES = '[0,1)')."
    )
  }
  selected
}

soql_string <- function(x) {
  paste0("'", gsub("'", "''", as.character(x), fixed = TRUE), "'")
}

soql_in <- function(field, values) {
  paste0(field, " in (", paste(soql_string(values), collapse = ","), ")")
}

filtered_count <- function(dataset_id, where = NULL) {
  ans <- fetch_json_dt(soda_url(
    socrata_domain,
    dataset_id,
    list(
      "$select" = "count(*) as n",
      "$where" = where,
      "$limit" = 1
    ),
    format = "json"
  ))
  if (!nrow(ans)) return(0L)
  as.integer(ans[[1L]][1L])
}

# Downloads a filtered query with a known finite row count. The explicit count
# prevents an endless pagination loop if the endpoint ignores an offset.
download_filtered_query <- function(dataset_id, dest, select, where = NULL,
                                    order = NULL, chunk_size = 25000L) {
  total <- filtered_count(dataset_id, where)
  message("Rows selected from ", dataset_id, ": ", format(total, big.mark = ","))

  if (total == 0L) {
    stop("The filtered Socrata query returned no records for dataset ", dataset_id, ".")
  }

  if (file.exists(dest)) unlink(dest)
  offset <- 0L
  first <- TRUE
  page <- 0L
  max_pages <- ceiling(total / chunk_size) + 2L

  while (offset < total) {
    page <- page + 1L
    if (page > max_pages) {
      stop("Pagination safety limit reached for dataset ", dataset_id, ".")
    }

    n_request <- min(chunk_size, total - offset)
    url <- soda_url(
      socrata_domain,
      dataset_id,
      list(
        "$select" = select,
        "$where" = where,
        "$order" = order,
        "$limit" = n_request,
        "$offset" = offset
      )
    )

    chunk <- fetch_csv_dt(url, timeout = 600L)
    if (!nrow(chunk)) {
      stop(
        "The Socrata endpoint returned an empty page after ", offset,
        " of ", total, " expected rows."
      )
    }

    data.table::fwrite(
      chunk,
      dest,
      append = !first,
      col.names = first
    )
    first <- FALSE
    offset <- offset + nrow(chunk)

    message(
      "  ", format(offset, big.mark = ","), " / ",
      format(total, big.mark = ","), " rows"
    )
  }

  if (!file.exists(dest) || file.info(dest)$size < 20L) {
    stop("The filtered download did not create a valid file: ", dest)
  }
  invisible(dest)
}

# Downloads SIVIC in separate calendar-year requests. This avoids large offsets
# and makes retries local to one year rather than restarting the whole dataset.
download_sivic_target <- function(dest) {
  sample <- fetch_sample(sivic_id)
  api_names <- names(sample)
  clean <- clean_names_base(api_names)
  lookup <- setNames(api_names, clean)

  c_date_clean <- find_col(clean, c("^data$", "^date$", "data_diagnostic", "dia"), label = "date")
  c_abs_clean <- find_col(clean, c("^abs$", "area_basica", "codi_abs"), label = "ABS")
  c_diag_clean <- find_col(clean, c("diagnostic", "diagnosi", "malaltia"), label = "diagnosis")
  c_age_clean <- find_col(clean, c("grup.*edat", "edat.*grup", "^edat$"), label = "age group")
  c_cases_clean <- find_col(clean, c("^recompte$", "^casos$", "nombre.*cas", "num.*cas"), label = "cases")
  c_pop_clean <- find_col(clean, c("poblacio", "population"), label = "population")
  c_ses_clean <- find_col(clean, c("socio", "index.*socio"), required = FALSE, label = "socioeconomic index")
  c_sex_clean <- find_col(clean, c("^sexe$", "^sex$"), required = FALSE, label = "sex")

  c_date <- unname(lookup[c_date_clean])
  c_abs <- unname(lookup[c_abs_clean])
  c_diag <- unname(lookup[c_diag_clean])
  c_age <- unname(lookup[c_age_clean])
  c_cases <- unname(lookup[c_cases_clean])
  c_pop <- unname(lookup[c_pop_clean])
  c_ses <- if (!is.na(c_ses_clean)) unname(lookup[c_ses_clean]) else NA_character_
  c_sex <- if (!is.na(c_sex_clean)) unname(lookup[c_sex_clean]) else NA_character_

  diagnosis_where <- paste0("lower(", c_diag, ") like '%bronquiol%'")

  # Distinct categorical values are requested as JSON. The SIVIC age labels
  # are interval strings containing commas (for example [0,1)), which made the
  # previous CSV discovery query split one value into several columns.
  age_levels <- fetch_distinct_values(
    dataset_id = sivic_id,
    field = c_age,
    where = diagnosis_where
  )
  infant_levels <- choose_infant_levels(
    age_levels,
    label = "SIVIC bronchiolitis age groups",
    fallback_pattern = age_pattern
  )
  message("Infant age groups selected: ", paste(dQuote(infant_levels), collapse = ", "))

  # c_age is deliberately omitted from the bulk CSV. It has already been
  # filtered server-side and its interval label may contain an unquoted comma.
  fields <- c(c_date, c_abs, c_diag, c_cases, c_pop, c_ses, c_sex)
  fields <- unique(fields[!is.na(fields) & nzchar(fields)])
  select <- paste(fields, collapse = ",")
  order <- paste(c_date, c_abs, sep = ",")

  if (file.exists(dest)) unlink(dest)
  first_file <- TRUE
  years <- seq.int(year(start_date), year(end_date))

  for (yy in years) {
    lower <- max(start_date, as.Date(sprintf("%d-01-01", yy)))
    upper <- min(end_date + 1L, as.Date(sprintf("%d-01-01", yy + 1L)))
    where <- paste(
      diagnosis_where,
      soql_in(c_age, infant_levels),
      sprintf("%s >= '%sT00:00:00.000'", c_date, lower),
      sprintf("%s < '%sT00:00:00.000'", c_date, upper),
      sep = " and "
    )

    total <- filtered_count(sivic_id, where)
    message("SIVIC ", yy, ": ", format(total, big.mark = ","), " selected rows")
    if (total == 0L) next

    year_file <- tempfile(fileext = ".csv")
    on.exit(unlink(year_file), add = TRUE)
    download_filtered_query(
      dataset_id = sivic_id,
      dest = year_file,
      select = select,
      where = where,
      order = order,
      chunk_size = chunk_size
    )

    year_dt <- fread(year_file, encoding = "UTF-8", showProgress = FALSE)
    fwrite(year_dt, dest, append = !first_file, col.names = first_file)
    first_file <- FALSE
    unlink(year_file)
  }

  if (first_file) stop("No SIVIC records were downloaded for the requested period.")
  invisible(dest)
}

download_immun_target <- function(dest) {
  sample <- fetch_sample(immun_id)
  api_names <- names(sample)
  clean <- clean_names_base(api_names)
  lookup <- setNames(api_names, clean)

  c_campaign_clean <- find_col(clean, c("campanya", "campaign", "vacuna", "immun"), label = "campaign")
  c_abs_clean <- find_col(clean, c("^abs$", "area_basica", "codi_abs"), label = "ABS")
  c_count_clean <- find_col(clean, c("^recompte$", "recompte", "persones", "vacunats", "immunitz", "nombre"), label = "immunised count")
  c_pop_clean <- find_col(clean, c("poblacio", "population"), required = FALSE, label = "target population")
  c_year_clean <- find_col(clean, c("^any$", "year", "temporada"), required = FALSE, label = "campaign year")
  c_age_clean <- find_col(clean, c("grup.*edat", "edat.*grup", "^edat$"), required = FALSE, label = "age group")
  c_ses_clean <- find_col(clean, c("socio", "index.*socio"), required = FALSE, label = "socioeconomic index")
  c_sex_clean <- find_col(clean, c("^sexe$", "^sex$"), required = FALSE, label = "sex")

  c_campaign <- unname(lookup[c_campaign_clean])
  c_abs <- unname(lookup[c_abs_clean])
  c_count <- unname(lookup[c_count_clean])
  optional <- c(c_pop_clean, c_year_clean, c_age_clean, c_ses_clean, c_sex_clean)
  optional <- unname(lookup[optional[!is.na(optional)]])

  campaigns <- fetch_distinct_values(
    dataset_id = immun_id,
    field = c_campaign
  )
  rsv_campaigns <- campaigns[
    grepl("VRS|sincitial|nirsevimab", campaigns, ignore.case = TRUE, perl = TRUE)
  ]

  if (!length(rsv_campaigns)) {
    stop(
      "No RSV campaign was identified. Available campaigns: ",
      paste(campaigns, collapse = ", ")
    )
  }
  message("RSV campaigns selected: ", paste(rsv_campaigns, collapse = ", "))

  where_parts <- soql_in(c_campaign, rsv_campaigns)
  selected_immun_age_levels <- character()
  if (!is.na(c_age_clean)) {
    c_age <- unname(lookup[c_age_clean])
    immun_age_levels <- fetch_distinct_values(
      dataset_id = immun_id,
      field = c_age,
      where = where_parts
    )
    selected_immun_age_levels <- choose_infant_levels(
      immun_age_levels,
      label = "RSV immunisation age groups",
      fallback_pattern = age_pattern
    )
    message(
      "RSV infant age groups selected: ",
      paste(dQuote(selected_immun_age_levels), collapse = ", ")
    )
    where_parts <- paste(where_parts, soql_in(c_age, selected_immun_age_levels), sep = " and ")
    optional <- setdiff(optional, c_age)
  }

  # The age field is omitted after server-side filtering for the same reason as
  # in SIVIC: interval labels can contain commas that are not safely represented
  # by the portal's CSV exporter.
  fields <- unique(c(c_campaign, c_abs, c_count, optional))
  where <- where_parts

  download_filtered_query(
    dataset_id = immun_id,
    dest = dest,
    select = paste(fields, collapse = ","),
    where = where,
    order = paste(c_campaign, c_abs, sep = ","),
    chunk_size = chunk_size
  )
}

if (force_download || !file.exists(sivic_file) || file.info(sivic_file)$size < 20L) {
  message("Downloading only infant bronchiolitis records from SIVIC...")
  download_sivic_target(sivic_file)
} else {
  message("Using cached targeted SIVIC file: ", sivic_file)
}

if (force_download || !file.exists(immun_file) || file.info(immun_file)$size < 20L) {
  message("Downloading only RSV immunisation campaign records...")
  download_immun_target(immun_file)
} else {
  message("Using cached targeted immunisation file: ", immun_file)
}

# -----------------------------------------------------------------------------
# Local preparation.
# -----------------------------------------------------------------------------
message("Reading targeted SIVIC primary-care data...")
siv <- data.table::fread(sivic_file, encoding = "UTF-8", showProgress = TRUE)
orig <- names(siv)
setnames(siv, orig, clean_names_base(orig))

c_date <- find_col(names(siv), c("^data$", "^date$", "data_diagnostic", "dia"), label = "date")
c_abs <- find_col(names(siv), c("^abs$", "area_basica", "codi_abs"), label = "ABS")
c_diag <- find_col(names(siv), c("diagnostic", "diagnosi", "malaltia"), label = "diagnosis")
c_age <- find_col(names(siv), c("grup.*edat", "edat.*grup", "^edat$"), required = FALSE, label = "age group")
c_cases <- find_col(names(siv), c("^recompte$", "^casos$", "nombre.*cas", "num.*cas"), label = "cases")
c_pop <- find_col(names(siv), c("poblacio", "population"), label = "population")
c_ses <- find_col(names(siv), c("socio", "index.*socio"), required = FALSE, label = "socioeconomic index")
c_sex <- find_col(names(siv), c("^sexe$", "^sex$"), required = FALSE, label = "sex")

siv[, date := parse_date_flex(get(c_date))]
siv[, diagnosis := as.character(get(c_diag))]
siv[, age_group := if (!is.na(c_age)) as.character(get(c_age)) else "under 1 year"]
siv[, abs_code := as.character(get(c_abs))]
siv[, cases := as.numeric(get(c_cases))]
siv[, population := as.numeric(get(c_pop))]
siv[, ses := if (!is.na(c_ses)) suppressWarnings(as.numeric(get(c_ses))) else NA_real_]
siv[, sex := if (!is.na(c_sex)) as.character(get(c_sex)) else "All"]

bron <- siv[
  grepl("bronquiol", diagnosis, ignore.case = TRUE) &
    !is.na(date) & !is.na(abs_code)
]
if (!nrow(bron)) {
  stop("No infant bronchiolitis records were selected after the targeted download.")
}

bron[, week_date := floor_date(date, unit = "week", week_start = 1)]
bron[, season := epidemic_season(week_date)]
weekly <- bron[, .(
  y = sum(cases, na.rm = TRUE),
  population = max(population, na.rm = TRUE),
  ses = weighted_mean_safe(ses, pmax(population, 1))
), by = .(abs_code, week_date, season)]
weekly[!is.finite(population) | population <= 0, population := NA_real_]
weekly <- weekly[is.finite(population)]

message("Reading targeted immunisation data...")
imm <- data.table::fread(immun_file, encoding = "UTF-8", showProgress = TRUE)
orig <- names(imm)
setnames(imm, orig, clean_names_base(orig))

c_campaign <- find_col(names(imm), c("campanya", "campaign", "vacuna", "immun"), label = "campaign")
c_abs_i <- find_col(names(imm), c("^abs$", "area_basica", "codi_abs"), label = "ABS")
c_count_i <- find_col(
  names(imm),
  c("^recompte$", "recompte", "persones", "vacunats", "immunitz", "nombre"),
  label = "immunised count"
)
c_pop_i <- find_col(names(imm), c("poblacio", "population"), required = FALSE, label = "target population")
c_year_i <- find_col(names(imm), c("^any$", "year", "temporada"), required = FALSE, label = "campaign year")
c_age_i <- find_col(names(imm), c("grup.*edat", "edat.*grup", "^edat$"), required = FALSE, label = "age group")

imm[, campaign := as.character(get(c_campaign))]
imm[, abs_code := as.character(get(c_abs_i))]
imm[, immunised := as.numeric(get(c_count_i))]
imm[, target_population := if (!is.na(c_pop_i)) as.numeric(get(c_pop_i)) else NA_real_]
imm[, campaign_year := if (!is.na(c_year_i)) first_year_in_text(get(c_year_i)) else NA_integer_]
imm[, age_group := if (!is.na(c_age_i)) as.character(get(c_age_i)) else ""]

rsv <- imm[grepl("VRS|sincitial|nirsevimab", campaign, ignore.case = TRUE)]
if (!nrow(rsv)) {
  warning("No RSV campaign rows were identified. Immunisation intensity will be missing.")
  intensity <- data.table(
    abs_code = character(), season = character(), immunised = numeric(),
    reference_population = numeric(), immunisation_intensity_raw = numeric(),
    immunisation_intensity = numeric()
  )
  intensity_scale <- NA_real_
} else {
  # Age restriction was already applied in the Socrata query when an age
  # field was available. No second regex filter is needed here.
  rsv[is.na(campaign_year), campaign_year := first_year_in_text(campaign)]
  if (anyNA(rsv$campaign_year)) {
    warning("Some RSV campaign years could not be inferred and were removed.")
    rsv <- rsv[!is.na(campaign_year)]
  }
  # Exclude campaigns belonging to incomplete epidemic seasons.
  rsv <- rsv[campaign_year <= last_complete_season_start]
  rsv[, season := paste0(campaign_year, "-", campaign_year + 1L)]

  intensity <- rsv[, .(
    immunised = if (all(is.na(immunised))) NA_real_ else sum(immunised, na.rm = TRUE),
    reference_population = if (all(is.na(target_population))) {
      NA_real_
    } else {
      sum(target_population, na.rm = TRUE)
    }
  ), by = .(abs_code, season)]

  # The campaign numerator is accumulated over an open infant cohort, whereas
  # the published population is a reference population rather than a complete
  # eligible-campaign denominator. Their ratio is therefore an immunisation
  # intensity and is not constrained to [0, 1].
  intensity[, immunisation_intensity_raw := fifelse(
    is.finite(reference_population) & reference_population > 0 &
      is.finite(immunised) & immunised >= 0,
    immunised / reference_population,
    NA_real_
  )]

  intensity_scale_env <- suppressWarnings(as.numeric(
    Sys.getenv("BEACON_IMMUNISATION_INTENSITY_SCALE", NA_character_)
  ))
  if (is.finite(intensity_scale_env) && intensity_scale_env > 0) {
    intensity_scale <- intensity_scale_env
  } else {
    intensity_scale <- median(
      intensity$immunisation_intensity_raw[
        is.finite(intensity$immunisation_intensity_raw) &
          intensity$immunisation_intensity_raw > 0
      ],
      na.rm = TRUE
    )
  }
  if (!is.finite(intensity_scale) || intensity_scale <= 0) {
    stop(
      "Could not define a positive immunisation-intensity scale. Check the ",
      "immunised counts and reference-population field."
    )
  }

  # A value of 1 corresponds to the median observed ABS-season intensity in
  # programme seasons. Zero represents no recorded immunisation activity.
  intensity[, immunisation_intensity :=
    immunisation_intensity_raw / intensity_scale]
  intensity[, intensity_issue := fifelse(
    !is.finite(immunisation_intensity_raw),
    "missing_or_invalid_reference_population",
    fifelse(immunisation_intensity_raw < 0, "negative", "valid")
  )]

  invalid_intensity <- intensity[intensity_issue != "valid"]
  if (nrow(invalid_intensity)) {
    warning(
      nrow(invalid_intensity),
      " ABS-season immunisation-intensity estimates were invalid and were ",
      "set to NA. Inspect results/catalonia_immunisation_intensity_diagnostics.csv."
    )
  }
}

# Exclude incomplete seasons from the surveillance series before merging.
weekly <- weekly[as.integer(substr(season, 1, 4)) <= last_complete_season_start]
intensity <- intensity[as.integer(substr(season, 1, 4)) <= last_complete_season_start]

analysis <- merge(
  weekly,
  intensity[, .(
    abs_code, season, immunised, reference_population,
    immunisation_intensity_raw, immunisation_intensity
  )],
  by = c("abs_code", "season"),
  all.x = TRUE
)
analysis[, treated := as.integer(as.integer(substr(season, 1, 4)) >= 2023)]
analysis[treated == 0 & is.na(immunisation_intensity), immunisation_intensity := 0]
analysis[treated == 0 & is.na(immunisation_intensity_raw), immunisation_intensity_raw := 0]
analysis[, week_in_season := as.integer(difftime(
  week_date,
  as.Date(paste0(substr(season, 1, 4), "-07-01")),
  units = "weeks"
)) + 1L]
setorder(analysis, abs_code, season, week_in_season)
# Lag WITHIN an ABS-season series, not within an ABS. Lagging by ABS alone
# carries the last week of June into the first week of the following July, and
# the dynamic g-formula in beacon_affine_v2.stan restarts the recursion at the
# first week of each series, so the two constructions must agree.
analysis[, lag_log_y := shift(log1p(y)), by = .(abs_code, season)]
# The first week of each series has no predecessor. It is set to zero here and
# to lag_init in Stan; both must use the same convention.
analysis[is.na(lag_log_y), lag_log_y := 0]
analysis[, ses_z := as.numeric(scale(ses))]
analysis[!is.finite(ses_z), ses_z := 0]

# ---------------------------------------------------------------------------
# Overlap / positivity diagnostic.
#
# beta_intensity is identified by between-area variation in intensity WITHIN
# programme seasons, because the season random effects absorb anything that is
# programme-wide. If nearly every ABS sits at the season median, that contrast
# is close to degenerate and the posterior will revert towards its prior. This
# has to be established before any effect estimate is interpreted, so it is
# written out here rather than after fitting.
# ---------------------------------------------------------------------------
overlap_desc <- overlap_diagnostic(analysis, area_col = "abs_code",
                                   season_col = "season")
fwrite(overlap_desc, "results/catalonia_overlap_diagnostic.csv")
if (nrow(overlap_desc)) {
  message("Overlap in scaled immunisation intensity by programme season:")
  print(overlap_desc)
  if (any(overlap_desc$prop_within_5pct > 0.5, na.rm = TRUE)) {
    warning(
      "In at least one programme season more than half of ABSs lie within 5% ",
      "of the season median intensity. The dose-response contrast is close to ",
      "degenerate in that season and the intervention posterior should be ",
      "read as weakly identified rather than as evidence of no effect."
    )
  }
}

fwrite(weekly, "data/derived/catalonia_bronchiolitis_weekly.csv")
fwrite(intensity, "data/derived/catalonia_rsv_immunisation_intensity.csv")
fwrite(analysis, "data/derived/catalonia_beacon_analysis.csv")
if (exists("invalid_intensity")) {
  fwrite(
    invalid_intensity,
    "results/catalonia_immunisation_intensity_diagnostics.csv"
  )
} else {
  fwrite(
    data.table(
      abs_code = character(), season = character(), immunised = numeric(),
      reference_population = numeric(), immunisation_intensity_raw = numeric(),
      immunisation_intensity = numeric(), intensity_issue = character()
    ),
    "results/catalonia_immunisation_intensity_diagnostics.csv"
  )
}

intensity_by_abs <- analysis[, .(
  immunisation_intensity = if (all(is.na(immunisation_intensity))) {
    NA_real_
  } else {
    mean(immunisation_intensity, na.rm = TRUE)
  },
  immunisation_intensity_raw = if (all(is.na(immunisation_intensity_raw))) {
    NA_real_
  } else {
    mean(immunisation_intensity_raw, na.rm = TRUE)
  },
  area_population = median(population, na.rm = TRUE),
  treated = max(treated, na.rm = TRUE)
), by = .(abs_code, season)]
intensity_by_abs[!is.finite(area_population), area_population := NA_real_]
fwrite(
  intensity_by_abs,
  "results/catalonia_immunisation_intensity_by_abs.csv"
)

weekly_desc <- analysis[
  week_in_season >= 1 & week_in_season <= 52,
  .(
    total_cases = sum(y, na.rm = TRUE),
    total_population = sum(population, na.rm = TRUE),
    n_abs = uniqueN(abs_code),
    treated = max(treated, na.rm = TRUE)
  ),
  by = .(season, week_date, week_in_season)
]
weekly_desc[, incidence_per_100k := fifelse(
  total_population > 0,
  1e5 * total_cases / total_population,
  NA_real_
)]
setorder(weekly_desc, week_date)
fwrite(weekly_desc, "results/catalonia_descriptive_weekly.csv")

season_summary <- weekly_desc[, .(
  total_cases = sum(total_cases, na.rm = TRUE),
  person_weeks = sum(total_population, na.rm = TRUE),
  incidence_per_100k_pw = fifelse(
    sum(total_population, na.rm = TRUE) > 0,
    1e5 * sum(total_cases, na.rm = TRUE) / sum(total_population, na.rm = TRUE),
    NA_real_
  ),
  mean_population = mean(total_population, na.rm = TRUE),
  n_weeks = .N,
  n_abs = max(n_abs, na.rm = TRUE),
  first_week = min(week_date),
  last_week = max(week_date),
  treated = max(treated, na.rm = TRUE)
), by = season]

peak_summary <- weekly_desc[order(season, -incidence_per_100k), .SD[1], by = season][,
  .(
    season,
    peak_week_date = week_date,
    peak_week_in_season = week_in_season,
    peak_cases = total_cases,
    peak_incidence_per_100k = incidence_per_100k
  )
]

intensity_season <- intensity_by_abs[, .(
  mean_immunisation_intensity = if (all(is.na(immunisation_intensity))) {
    NA_real_
  } else {
    weighted_mean_safe(
      immunisation_intensity,
      fifelse(is.finite(area_population) & area_population > 0, area_population, 1)
    )
  },
  median_immunisation_intensity = if (all(is.na(immunisation_intensity))) {
    NA_real_
  } else {
    median(immunisation_intensity, na.rm = TRUE)
  },
  intensity_q25 = if (all(is.na(immunisation_intensity))) {
    NA_real_
  } else {
    quantile(immunisation_intensity, 0.25, na.rm = TRUE)
  },
  intensity_q75 = if (all(is.na(immunisation_intensity))) {
    NA_real_
  } else {
    quantile(immunisation_intensity, 0.75, na.rm = TRUE)
  },
  mean_raw_immunisation_intensity = if (all(is.na(immunisation_intensity_raw))) {
    NA_real_
  } else {
    weighted_mean_safe(
      immunisation_intensity_raw,
      fifelse(is.finite(area_population) & area_population > 0, area_population, 1)
    )
  }
), by = season]

season_summary <- merge(season_summary, peak_summary, by = "season", all.x = TRUE)
season_summary <- merge(season_summary, intensity_season, by = "season", all.x = TRUE)
setorder(season_summary, season)
fwrite(season_summary, "results/catalonia_season_summary.csv")

overall_programme_intensity <- intensity_by_abs[
  treated == 1 & is.finite(immunisation_intensity),
  weighted_mean_safe(
    immunisation_intensity,
    fifelse(is.finite(area_population) & area_population > 0, area_population, 1)
  )
]
if (!length(overall_programme_intensity) || !is.finite(overall_programme_intensity)) {
  overall_programme_intensity <- NA_real_
}

meta <- data.table(
  item = c(
    "download_date", "n_rows_analysis", "n_abs", "n_seasons", "first_week",
    "last_week", "total_cases", "n_programme_seasons",
    "mean_programme_immunisation_intensity", "immunisation_intensity_scale",
    "immunisation_intensity_definition", "infant_age_pattern",
    "last_complete_season"
  ),
  value = c(
    as.character(Sys.Date()), nrow(analysis), uniqueN(analysis$abs_code),
    uniqueN(analysis$season), as.character(min(analysis$week_date)),
    as.character(max(analysis$week_date)), sum(analysis$y, na.rm = TRUE),
    uniqueN(analysis[treated == 1, season]), overall_programme_intensity,
    intensity_scale,
    "(immunised/reference_population)/median_treated_ABS_season_ratio",
    age_pattern, last_complete_season
  )
)
fwrite(meta, "results/catalonia_data_metadata.csv")

message("Prepared targeted analysis and descriptive outputs.")

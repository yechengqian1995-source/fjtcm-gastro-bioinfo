#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(curl)
  library(xml2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

inventory_path <- file.path(project_dir, "data", "data_inventory.tsv")
out_dir <- file.path(project_dir, "data", "geo_metadata")
raw_dir <- file.path(out_dir, "raw_geo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

inventory <- read_tsv(inventory_path, show_col_types = FALSE)
geo_accessions <- inventory %>%
  filter(str_detect(accession, "^GSE[0-9]+$")) %>%
  pull(accession) %>%
  unique()

collapse_value <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(unique(x), collapse = " || ")
}

sanitize_key <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

strip_quotes <- function(x) {
  x <- str_trim(as.character(x))
  x <- str_replace_all(x, '^"', "")
  x <- str_replace_all(x, '"$', "")
  x
}

split_geo_line <- function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  strip_quotes(parts)
}

series_prefix <- function(accession) {
  num <- as.integer(str_remove(accession, "^GSE"))
  paste0("GSE", floor(num / 1000), "nnn")
}

matrix_listing_url <- function(accession) {
  paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    series_prefix(accession), "/",
    accession, "/matrix/"
  )
}

find_series_matrix_files <- function(accession) {
  listing_url <- matrix_listing_url(accession)
  html <- read_html(listing_url)
  hrefs <- xml_attr(xml_find_all(html, ".//a"), "href")
  hrefs <- hrefs[!is.na(hrefs)]
  files <- hrefs[str_detect(hrefs, "series_matrix.*\\.txt\\.gz$")]
  files <- unique(basename(files))
  if (!length(files)) {
    stop("No series_matrix.txt.gz file found at ", listing_url)
  }
  paste0(listing_url, files)
}

download_matrix_file <- function(url) {
  local_path <- file.path(raw_dir, basename(url))
  curl_download(url, local_path, mode = "wb", quiet = TRUE)
  local_path
}

empty_parse_result <- function(accession, matrix_url = NA_character_, matrix_file = NA_character_) {
  list(
    accession = accession,
    matrix_url = matrix_url,
    matrix_file = matrix_file,
    n_samples = 0L,
    series = tibble(accession = character(), matrix_file = character(), field = character(), value = character()),
    samples = tibble(accession = character(), matrix_file = character(), sample_index = integer(), gsm = character()),
    characteristics = tibble(
      accession = character(), matrix_file = character(), sample_index = integer(), gsm = character(),
      characteristic_index = integer(), raw_characteristic = character(), key = character(), value = character()
    )
  )
}

parse_series_metadata <- function(lines, accession, matrix_file) {
  series_lines <- lines[str_detect(lines, "^!Series_")]
  if (!length(series_lines)) {
    return(tibble(accession = character(), matrix_file = character(), field = character(), value = character()))
  }

  bind_rows(lapply(series_lines, function(line) {
    parts <- split_geo_line(line)
    field <- sanitize_key(str_remove(parts[[1]], "^!Series_"))
    values <- parts[-1]
    tibble(
      accession = accession,
      matrix_file = basename(matrix_file),
      field = field,
      value = values
    )
  })) %>%
    filter(!is.na(value), nzchar(value)) %>%
    group_by(accession, matrix_file, field) %>%
    summarise(value = collapse_value(value), .groups = "drop")
}

parse_sample_metadata <- function(lines, accession, matrix_file) {
  sample_lines <- lines[str_detect(lines, "^!Sample_")]
  if (!length(sample_lines)) {
    return(empty_parse_result(accession, matrix_file = basename(matrix_file))$samples)
  }

  sample_long <- bind_rows(lapply(sample_lines, function(line) {
    parts <- split_geo_line(line)
    field <- sanitize_key(str_remove(parts[[1]], "^!Sample_"))
    values <- parts[-1]
    tibble(
      accession = accession,
      matrix_file = basename(matrix_file),
      sample_index = seq_along(values),
      field = field,
      value = values
    )
  })) %>%
    filter(!is.na(value), nzchar(value))

  if (!nrow(sample_long)) {
    return(empty_parse_result(accession, matrix_file = basename(matrix_file))$samples)
  }

  sample_wide <- sample_long %>%
    group_by(accession, matrix_file, sample_index, field) %>%
    summarise(value = collapse_value(value), .groups = "drop") %>%
    pivot_wider(names_from = field, values_from = value)

  if (!"geo_accession" %in% names(sample_wide)) {
    sample_wide$geo_accession <- NA_character_
  }

  sample_wide %>%
    mutate(
      gsm = if_else(
        !is.na(geo_accession) & nzchar(geo_accession),
        geo_accession,
        paste0(accession, "_matrix_", dense_rank(matrix_file), "_sample_", sample_index)
      ),
      .before = 4
    )
}

parse_characteristics <- function(lines, accession, matrix_file, sample_meta) {
  char_lines <- lines[str_detect(lines, "^!Sample_characteristics_ch1")]
  if (!length(char_lines) || !nrow(sample_meta)) {
    return(empty_parse_result(accession, matrix_file = basename(matrix_file))$characteristics)
  }

  gsm_map <- sample_meta %>% select(sample_index, gsm)

  bind_rows(lapply(seq_along(char_lines), function(i) {
    parts <- split_geo_line(char_lines[[i]])
    values <- parts[-1]
    tibble(
      accession = accession,
      matrix_file = basename(matrix_file),
      sample_index = seq_along(values),
      characteristic_index = i,
      raw_characteristic = values
    )
  })) %>%
    left_join(gsm_map, by = "sample_index") %>%
    mutate(
      key = if_else(
        str_detect(raw_characteristic, ":"),
        str_trim(str_replace(raw_characteristic, ":.*$", "")),
        "characteristics_ch1"
      ),
      value = if_else(
        str_detect(raw_characteristic, ":"),
        str_trim(str_replace(raw_characteristic, "^[^:]+:", "")),
        raw_characteristic
      ),
      key = sanitize_key(key)
    ) %>%
    filter(!is.na(value), nzchar(value))
}

parse_matrix_file <- function(accession, matrix_url, matrix_file) {
  con <- gzfile(matrix_file, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  if (!length(lines)) {
    return(empty_parse_result(accession, matrix_url, basename(matrix_file)))
  }

  series_meta <- parse_series_metadata(lines, accession, matrix_file)
  sample_meta <- parse_sample_metadata(lines, accession, matrix_file)
  characteristics <- parse_characteristics(lines, accession, matrix_file, sample_meta)

  list(
    accession = accession,
    matrix_url = matrix_url,
    matrix_file = basename(matrix_file),
    n_samples = n_distinct(sample_meta$gsm),
    series = series_meta,
    samples = sample_meta,
    characteristics = characteristics
  )
}

fetch_one <- function(accession) {
  message("Freezing GEO metadata from series matrix: ", accession)
  tryCatch({
    matrix_urls <- find_series_matrix_files(accession)
    parsed <- lapply(matrix_urls, function(matrix_url) {
      matrix_file <- download_matrix_file(matrix_url)
      parse_matrix_file(accession, matrix_url, matrix_file)
    })

    list(
      accession = accession,
      ok = TRUE,
      error = NA_character_,
      matrix_urls = matrix_urls,
      n_matrix_files = length(matrix_urls),
      n_samples = sum(vapply(parsed, `[[`, integer(1), "n_samples")),
      parsed = parsed
    )
  }, error = function(e) {
    warning("Failed ", accession, ": ", conditionMessage(e))
    list(
      accession = accession,
      ok = FALSE,
      error = conditionMessage(e),
      matrix_urls = character(),
      n_matrix_files = 0L,
      n_samples = 0L,
      parsed = list(empty_parse_result(accession))
    )
  })
}

results <- lapply(geo_accessions, fetch_one)
parsed_results <- unlist(lapply(results, `[[`, "parsed"), recursive = FALSE)

status <- tibble(
  accession = vapply(results, `[[`, character(1), "accession"),
  fetch_ok = vapply(results, `[[`, logical(1), "ok"),
  n_matrix_files = vapply(results, `[[`, integer(1), "n_matrix_files"),
  n_geo_samples_parsed = vapply(results, `[[`, integer(1), "n_samples"),
  matrix_urls = vapply(lapply(results, `[[`, "matrix_urls"), collapse_value, character(1)),
  error = vapply(results, `[[`, character(1), "error")
)

series_meta <- bind_rows(lapply(parsed_results, `[[`, "series"))
sample_meta <- bind_rows(lapply(parsed_results, `[[`, "samples"))
characteristics_long <- bind_rows(lapply(parsed_results, `[[`, "characteristics")) %>%
  filter(!is.na(key), !is.na(value), nzchar(key), nzchar(value))

characteristics_wide <- characteristics_long %>%
  group_by(accession, gsm, key) %>%
  summarise(value = collapse_value(value), .groups = "drop") %>%
  pivot_wider(names_from = key, values_from = value)

metadata_master <- sample_meta %>%
  left_join(characteristics_wide, by = c("accession", "gsm"), suffix = c("", "_characteristic")) %>%
  mutate(
    gse202379_metadata_swap_watch = accession == "GSE202379" & gsm %in% c("GSM6112262", "GSM6112263")
  )

keywords <- c(
  disease = "disease|diagnosis|condition|group|phenotype|status|histology|nash|nafld|masld|mash",
  fibrosis = "fibrosis|fibrotic|fibro|stage|f[0-4]",
  nas = "\\bnas\\b|nafld_activity|activity_score|activity",
  steatosis = "steatosis|fat|lipid",
  inflammation = "inflammation|inflammatory|lobular",
  ballooning = "balloon",
  bmi = "\\bbmi\\b|body_mass|obes",
  diabetes = "diabetes|t2d|diabetic",
  age = "\\bage\\b",
  sex = "\\bsex\\b|gender",
  donor = "donor|patient|subject|individual|case|participant",
  region = "region|section|roi|zone|area|slide|spot",
  platform = "platform|instrument|array|sequenc|library",
  tissue = "tissue|source|biopsy|liver"
)

candidate_fields <- if (nrow(metadata_master)) {
  bind_rows(lapply(unique(metadata_master$accession), function(acc) {
    dat <- metadata_master %>% filter(accession == acc)
    cols <- setdiff(names(dat), c("accession", "matrix_file", "sample_index", "gsm"))
    bind_rows(lapply(names(keywords), function(keyword_name) {
      pattern <- keywords[[keyword_name]]
      matched_cols <- cols[str_detect(cols, regex(pattern, ignore_case = TRUE))]
      if (!length(matched_cols)) {
        return(tibble(
          accession = acc,
          keyword = keyword_name,
          candidate_column = NA_character_,
          n_nonmissing = 0L,
          examples = NA_character_
        ))
      }
      bind_rows(lapply(matched_cols, function(col) {
        vals <- dat[[col]]
        vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
        tibble(
          accession = acc,
          keyword = keyword_name,
          candidate_column = col,
          n_nonmissing = length(vals),
          examples = collapse_value(head(unique(as.character(vals)), 5))
        )
      }))
    }))
  }))
} else {
  tibble(
    accession = character(),
    keyword = character(),
    candidate_column = character(),
    n_nonmissing = integer(),
    examples = character()
  )
}

sample_counts <- sample_meta %>%
  summarise(n_geo_samples = n_distinct(gsm), .by = accession) %>%
  right_join(tibble(accession = geo_accessions), by = "accession") %>%
  mutate(n_geo_samples = replace_na(n_geo_samples, 0L)) %>%
  left_join(inventory %>% select(accession, layer, species, planned_role, priority, verification_status), by = "accession") %>%
  left_join(status %>% select(accession, fetch_ok, n_matrix_files), by = "accession")

metadata_alerts <- bind_rows(
  tibble(
    accession = "GSE202379",
    alert_level = "hard_gate",
    alert = "GEO notes a 2025-06-06 metadata swap for GSM6112262 and GSM6112263; freeze and manually confirm updated sample labels before disease-group or donor-level analyses."
  ),
  tibble(
    accession = "GSE212837",
    alert_level = "hard_gate",
    alert = "Split human and mouse data and model donor/region nesting; regions, nuclei, and samples are not independent patients."
  ),
  tibble(
    accession = "GSE292268",
    alert_level = "scope_limit",
    alert = "Human MASLD CosMx spatial dataset is F2-only and small; use for localization, not fibrosis severity gradients."
  ),
  tibble(
    accession = "GSE248077",
    alert_level = "scope_limit",
    alert = "Mouse Visium dataset lacks a GEO-linked PMID at this phase; use only as supplementary mouse context until citation is verified."
  )
)

write_tsv(status, file.path(out_dir, "geo_fetch_status.tsv"))
write_tsv(series_meta, file.path(out_dir, "series_metadata_long.tsv"))
write_tsv(sample_meta, file.path(out_dir, "sample_metadata_raw.tsv"))
write_tsv(characteristics_long, file.path(out_dir, "sample_characteristics_long.tsv"))
write_tsv(characteristics_wide, file.path(out_dir, "sample_characteristics_wide.tsv"))
write_tsv(metadata_master, file.path(out_dir, "metadata_master_initial.tsv"))
write_tsv(candidate_fields, file.path(out_dir, "metadata_candidate_fields.tsv"))
write_tsv(sample_counts, file.path(out_dir, "geo_sample_counts.tsv"))
write_tsv(metadata_alerts, file.path(out_dir, "metadata_hard_gates.tsv"))

log_path <- file.path(project_dir, "docs", "phase1_geo_metadata_freeze_log.md")
cat(
  "# Phase 1 GEO Metadata Freeze Log\n\n",
  "Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "Metadata source: GEO FTP `series_matrix.txt.gz` files parsed directly, rather than GEOquery GSMList parsing.\n\n",
  "Accessions attempted: ", paste(geo_accessions, collapse = ", "), "\n\n",
  "Output directory: `data/geo_metadata/`\n\n",
  "Generated files:\n\n",
  "- `geo_fetch_status.tsv`\n",
  "- `series_metadata_long.tsv`\n",
  "- `sample_metadata_raw.tsv`\n",
  "- `sample_characteristics_long.tsv`\n",
  "- `sample_characteristics_wide.tsv`\n",
  "- `metadata_master_initial.tsv`\n",
  "- `metadata_candidate_fields.tsv`\n",
  "- `geo_sample_counts.tsv`\n",
  "- `metadata_hard_gates.tsv`\n\n",
  "Hard gates before biological modeling:\n\n",
  "- Confirm sample-level fibrosis/NAS/disease fields from the frozen metadata and any supplementary phenotype files.\n",
  "- Resolve the `GSE202379` GSM6112262/GSM6112263 metadata-swap warning before group comparisons.\n",
  "- Treat donors/patients as statistical units; do not treat nuclei, spots, regions, slides, or sections as independent patients.\n",
  "- Use these outputs for metadata readiness only. Expression matrices are not yet approved for discovery or validation analysis.\n",
  file = log_path,
  sep = ""
)

message("Done. Metadata written to ", out_dir)

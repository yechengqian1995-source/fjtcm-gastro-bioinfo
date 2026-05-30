#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(curl)
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

raw_dir <- file.path(project_dir, "data", "expression_raw")
proc_dir <- file.path(project_dir, "data", "expression_processed")
result_dir <- file.path(project_dir, "results", "bulk")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_tsv(file.path(project_dir, "data", "processed", "metadata_freeze_v1.tsv"), show_col_types = FALSE)

download_if_missing <- function(url, dest) {
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading ", basename(dest))
    curl_download(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

gene_symbol_from_entrez <- function(ids) {
  ids <- as.character(ids)
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ids,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  out <- unname(mapped[ids])
  out[is.na(out) | !nzchar(out)] <- ids[is.na(out) | !nzchar(out)]
  out
}

gene_symbol_from_ensembl <- function(ids) {
  ids <- str_replace(as.character(ids), "\\..*$", "")
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  out <- unname(mapped[ids])
  out[is.na(out) | !nzchar(out)] <- ids[is.na(out) | !nzchar(out)]
  out
}

read_gzip_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

fread_maybe_gz <- function(path, header = "auto") {
  if (str_detect(path, "\\.gz$")) {
    lines <- read_gzip_lines(path)
    data.table::fread(text = paste(lines, collapse = "\n"), header = header, check.names = FALSE)
  } else {
    data.table::fread(path, header = header, check.names = FALSE)
  }
}

collapse_duplicate_symbols <- function(mat) {
  symbol <- rownames(mat)
  keep <- !is.na(symbol) & nzchar(symbol)
  mat <- mat[keep, , drop = FALSE]
  symbol <- symbol[keep]
  split_idx <- split(seq_along(symbol), symbol)
  collapsed <- vapply(split_idx, function(idx) {
    if (length(idx) == 1) return(mat[idx, ])
    colSums(mat[idx, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(mat)))
  collapsed <- t(collapsed)
  colnames(collapsed) <- colnames(mat)
  collapsed
}

make_logcpm <- function(counts) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  edgeR::cpm(edgeR::DGEList(counts = counts), log = TRUE, prior.count = 1)
}

sample_key_source_for_accession <- function(acc) {
  case_when(
    acc == "GSE130970" ~ "geo_title",
    acc == "GSE162694" ~ "geo_title_548nash_code",
    acc == "GSE167523" ~ "geo_description",
    acc == "GSE126848" ~ "geo_description_numeric_padded_4",
    TRUE ~ "geo_gsm"
  )
}

make_metadata_key <- function(acc, dat) {
  if (acc == "GSE130970") {
    return(as.character(dat$title))
  }
  if (acc == "GSE162694") {
    return(str_extract(as.character(dat$title), "548nash[0-9]+"))
  }
  if (acc == "GSE167523") {
    return(as.character(dat$description))
  }
  if (acc == "GSE126848") {
    raw_key <- as.character(dat$description)
    return(if_else(
      !is.na(raw_key) & str_detect(raw_key, "^[0-9]{1,3}$"),
      str_pad(raw_key, width = 4, pad = "0"),
      raw_key
    ))
  }
  as.character(dat$gsm)
}

make_expression_key <- function(acc, sample_ids) {
  sample_ids <- as.character(sample_ids)
  if (acc == "GSE126848") {
    return(if_else(
      !is.na(sample_ids) & str_detect(sample_ids, "^[0-9]{1,3}$"),
      str_pad(sample_ids, width = 4, pad = "0"),
      sample_ids
    ))
  }
  sample_ids
}

parse_fibrosis_stage <- function(x) {
  raw <- str_trim(str_to_lower(as.character(x)))
  case_when(
    is.na(x) | raw == "" ~ NA_real_,
    raw %in% c("normal liver histology", "normal", "n", "control") ~ 0,
    str_detect(raw, "^f[0-4]$") ~ as.numeric(str_remove(raw, "^f")),
    str_detect(raw, "^[0-4](\\.0)?$") ~ suppressWarnings(as.numeric(raw)),
    TRUE ~ suppressWarnings(as.numeric(x))
  )
}

load_gse167523_xlsx_pheno <- function() {
  xlsx_path <- file.path(raw_dir, "GSE167523_Sample_phenotype_correspondence.xlsx")
  if (!file.exists(xlsx_path)) {
    return(tibble(gsm = character(), disease_subtype_xlsx = character(), sample_id_xlsx = character()))
  }
  openxlsx::read.xlsx(xlsx_path, sheet = 1) %>%
    transmute(
      gsm = as.character(.data[["GEO.accession.ID"]]),
      disease_subtype_xlsx = as.character(.data[["NAFL/NASH"]]),
      sample_id_xlsx = as.character(.data[["Sample.ID"]])
    )
}

standardize_pheno <- function(acc, sample_ids) {
  sample_ids <- as.character(sample_ids)
  sample_key_source <- sample_key_source_for_accession(acc)

  mdat <- metadata %>%
    filter(accession == acc)
  mdat$metadata_sample_key <- make_metadata_key(acc, mdat)
  mdat$metadata_key_source <- sample_key_source
  mdat <- mdat %>%
    add_count(metadata_sample_key, name = "metadata_key_count") %>%
    distinct(metadata_sample_key, .keep_all = TRUE)

  dat <- tibble(
    accession = acc,
    sample_id = sample_ids,
    sample_key = make_expression_key(acc, sample_ids),
    sample_key_source = sample_key_source
  ) %>%
    left_join(mdat, by = c("sample_key" = "metadata_sample_key")) %>%
    mutate(
      metadata_match_status = if_else(!is.na(gsm), "matched", "unmatched"),
      metadata_key_count = replace_na(metadata_key_count, 0L)
    )

  if (acc == "GSE167523") {
    dat <- dat %>%
      left_join(load_gse167523_xlsx_pheno(), by = "gsm")
  } else {
    dat <- dat %>%
      mutate(
        disease_subtype_xlsx = NA_character_,
        sample_id_xlsx = NA_character_
      )
  }

  dat <- dat %>%
    mutate(
      title = as.character(title),
      description = as.character(description),
      disease_label_raw = as.character(disease_label_raw),
      disease_subtype_raw = coalesce(disease_subtype_xlsx, as.character(disease_subtype), as.character(disease_state)),
      fibrosis_stage_raw = as.character(fibrosis_stage_raw),
      nas_raw = as.character(nas_raw),
      age_raw = as.character(age_raw),
      sex_raw = as.character(sex_raw),
      bmi_raw = as.character(bmi_raw),
      diabetes_raw = as.character(diabetes_raw),
      donor_or_patient_raw = as.character(donor_or_patient_raw),
      platform_raw = as.character(platform_raw),
      fibrosis_stage = parse_fibrosis_stage(fibrosis_stage_raw),
      fibrosis_recoding_note = case_when(
        str_to_lower(fibrosis_stage_raw) == "normal liver histology" ~ "recoded_normal_liver_histology_to_F0",
        !is.na(fibrosis_stage_raw) & is.na(fibrosis_stage) ~ "unparsed_fibrosis_stage",
        TRUE ~ "none"
      ),
      nas_score = suppressWarnings(as.numeric(nas_raw)),
      advanced_fibrosis = case_when(
        !is.na(fibrosis_stage) & fibrosis_stage >= 3 ~ "advanced",
        !is.na(fibrosis_stage) & fibrosis_stage < 3 ~ "non_advanced",
        TRUE ~ NA_character_
      ),
      include_bulk_primary = metadata_match_status == "matched" & !is.na(fibrosis_stage),
      include_disease_state_secondary = metadata_match_status == "matched" &
        (!is.na(disease_label_raw) | !is.na(disease_subtype_raw))
    ) %>%
    arrange(match(sample_id, sample_ids))

  mapping_qc <- dat %>%
    transmute(
      accession = acc,
      sample_id,
      sample_key,
      sample_key_source,
      metadata_match_status,
      metadata_key_count,
      gsm,
      title,
      description,
      disease_label_raw,
      disease_subtype_raw,
      fibrosis_stage_raw,
      fibrosis_stage,
      nas_raw,
      nas_score
    )
  write_tsv(mapping_qc, file.path(result_dir, paste0(acc, "_sample_mapping_qc.tsv")))

  dat %>%
    transmute(
      accession = acc,
      sample_id,
      sample_key,
      sample_key_source,
      metadata_match_status,
      metadata_key_count,
      gsm,
      title,
      description,
      disease_label_raw,
      disease_subtype_raw,
      fibrosis_stage_raw,
      fibrosis_recoding_note,
      nas_raw,
      age_raw,
      sex_raw,
      bmi_raw,
      diabetes_raw,
      donor_or_patient_raw,
      platform_raw,
      fibrosis_stage,
      nas_score,
      advanced_fibrosis,
      include_bulk_primary,
      include_disease_state_secondary
    )
}

write_dataset <- function(acc, counts, pheno, matrix_type, source_url) {
  common <- intersect(colnames(counts), pheno$sample_id)
  counts <- counts[, common, drop = FALSE]
  pheno <- pheno %>% filter(sample_id %in% common) %>% arrange(match(sample_id, common))
  logcpm <- make_logcpm(counts)

  saveRDS(
    list(
      accession = acc,
      counts = counts,
      logcpm = logcpm,
      pheno = pheno,
      matrix_type = matrix_type,
      source_url = source_url
    ),
    file.path(proc_dir, paste0(acc, "_bulk_expression.rds"))
  )
  write_tsv(pheno, file.path(proc_dir, paste0(acc, "_phenotype.tsv")))
  write_tsv(
    tibble(
      accession = acc,
      n_genes = nrow(counts),
      n_samples = ncol(counts),
      n_with_fibrosis = sum(!is.na(pheno$fibrosis_stage)),
      n_advanced = sum(pheno$advanced_fibrosis == "advanced", na.rm = TRUE),
      n_non_advanced = sum(pheno$advanced_fibrosis == "non_advanced", na.rm = TRUE),
      matrix_type = matrix_type,
      source_url = source_url
    ),
    file.path(result_dir, paste0(acc, "_preprocess_summary.tsv"))
  )
}

read_gse135251 <- function() {
  acc <- "GSE135251"
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135251/suppl/GSE135251_RAW.tar"
  tar_path <- download_if_missing(url, file.path(raw_dir, "GSE135251_RAW.tar"))
  extract_dir <- file.path(raw_dir, "GSE135251_RAW")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(list.files(extract_dir, pattern = "\\.counts\\.txt\\.gz$", full.names = TRUE))) {
    message("Extracting GSE135251 RAW tar")
    untar(tar_path, exdir = extract_dir)
  }
  files <- list.files(extract_dir, pattern = "^GSM.*\\.counts\\.txt\\.gz$", full.names = TRUE)
  if (!length(files)) stop("No GSE135251 count files found after extraction.")

  read_one <- function(path) {
    gsm <- str_extract(basename(path), "^GSM[0-9]+")
    dat <- fread_maybe_gz(path, header = FALSE)
    if (ncol(dat) < 2) stop("Unexpected count file format: ", path)
    setnames(dat, c(names(dat)[1], names(dat)[ncol(dat)]), c("gene_id", gsm))
    dat <- dat[, .(gene_id = as.character(gene_id), value = as.numeric(get(gsm)))]
    setnames(dat, "value", gsm)
    dat
  }
  lst <- lapply(files, read_one)
  merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), lst)
  mat <- as.matrix(merged[, -1, with = FALSE])
  rownames(mat) <- merged$gene_id
  mat[is.na(mat)] <- 0
  htseq_summary_rows <- str_detect(rownames(mat), "^__")
  if (any(htseq_summary_rows)) {
    mat <- mat[!htseq_summary_rows, , drop = FALSE]
  }
  rownames(mat) <- if (all(str_detect(rownames(mat), "^ENSG"))) {
    gene_symbol_from_ensembl(rownames(mat))
  } else if (all(str_detect(rownames(mat), "^[0-9]+$"))) {
    gene_symbol_from_entrez(rownames(mat))
  } else {
    toupper(rownames(mat))
  }
  mat <- collapse_duplicate_symbols(mat)
  pheno <- standardize_pheno(acc, colnames(mat))
  write_dataset(acc, mat, pheno, "raw_counts_individual_gsm_files", url)
}

read_matrix_with_gene_col <- function(path, header = TRUE) {
  dat <- fread_maybe_gz(path, header = header)
  gene_col <- names(dat)[1]
  genes <- as.character(dat[[gene_col]])
  mat <- as.matrix(dat[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  mat[is.na(mat)] <- 0
  mat
}

read_gse130970 <- function() {
  acc <- "GSE130970"
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130970/suppl/GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"
  path <- download_if_missing(url, file.path(raw_dir, basename(url)))
  mat <- read_matrix_with_gene_col(path)
  rownames(mat) <- gene_symbol_from_entrez(rownames(mat))
  mat <- collapse_duplicate_symbols(mat)
  pheno <- standardize_pheno(acc, colnames(mat))
  write_dataset(acc, mat, pheno, "salmon_tximport_counts_entrez_gene_id", url)
}

read_gse162694 <- function() {
  acc <- "GSE162694"
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162694/suppl/GSE162694_raw_counts.csv.gz"
  path <- download_if_missing(url, file.path(raw_dir, basename(url)))
  mat <- read_matrix_with_gene_col(path)
  rownames(mat) <- if (all(str_detect(rownames(mat), "^ENSG"))) {
    gene_symbol_from_ensembl(rownames(mat))
  } else if (all(str_detect(rownames(mat), "^[0-9]+$"))) {
    gene_symbol_from_entrez(rownames(mat))
  } else {
    toupper(rownames(mat))
  }
  mat <- collapse_duplicate_symbols(mat)
  pheno <- standardize_pheno(acc, colnames(mat))
  write_dataset(acc, mat, pheno, "raw_counts", url)
}

read_gse167523 <- function() {
  acc <- "GSE167523"
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167523/suppl/GSE167523_Raw_gene_counts_matrix.txt.gz"
  path <- download_if_missing(url, file.path(raw_dir, basename(url)))
  mat <- read_matrix_with_gene_col(path)
  rownames(mat) <- if (all(str_detect(rownames(mat), "^ENSG"))) {
    gene_symbol_from_ensembl(rownames(mat))
  } else if (all(str_detect(rownames(mat), "^[0-9]+$"))) {
    gene_symbol_from_entrez(rownames(mat))
  } else {
    toupper(rownames(mat))
  }
  mat <- collapse_duplicate_symbols(mat)
  pheno <- standardize_pheno(acc, colnames(mat))
  write_dataset(acc, mat, pheno, "raw_gene_counts_matrix", url)
}

read_gse126848 <- function() {
  acc <- "GSE126848"
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126848/suppl/GSE126848_Gene_counts_raw.txt.gz"
  path <- download_if_missing(url, file.path(raw_dir, basename(url)))
  mat <- read_matrix_with_gene_col(path)
  rownames(mat) <- if (all(str_detect(rownames(mat), "^ENSG"))) {
    gene_symbol_from_ensembl(rownames(mat))
  } else if (all(str_detect(rownames(mat), "^[0-9]+$"))) {
    gene_symbol_from_entrez(rownames(mat))
  } else {
    toupper(rownames(mat))
  }
  mat <- collapse_duplicate_symbols(mat)
  pheno <- standardize_pheno(acc, colnames(mat))
  write_dataset(acc, mat, pheno, "raw_gene_counts", url)
}

read_gse135251()
read_gse130970()
read_gse162694()
read_gse167523()
read_gse126848()

summaries <- list.files(result_dir, pattern = "_preprocess_summary\\.tsv$", full.names = TRUE) %>%
  discard(~ basename(.x) == "bulk_preprocess_summary.tsv") %>%
  lapply(read_tsv, show_col_types = FALSE) %>%
  bind_rows()
write_tsv(summaries, file.path(result_dir, "bulk_preprocess_summary.tsv"))

cat(
  "# Bulk Preprocess Summary\n\n",
  "Generated by `scripts/04_bulk_preprocess.R`.\n\n",
  "## Datasets\n\n",
  paste0(
    "- `", summaries$accession, "`: ", summaries$n_samples, " samples, ",
    summaries$n_genes, " genes, ", summaries$n_with_fibrosis,
    " samples with candidate fibrosis stage.\n",
    collapse = ""
  ),
  "\n## Boundary\n\n",
  "These are expression preprocessing outputs. They do not yet constitute biological results. Sample overlap, endpoint lock, and signature scoring are handled in later scripts.\n",
  file = file.path(project_dir, "docs", "bulk_preprocess_summary.md"),
  sep = ""
)

message("Done. Bulk expression objects written to ", proc_dir)

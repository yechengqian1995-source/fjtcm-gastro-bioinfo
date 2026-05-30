#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

sig_path <- file.path(project_dir, "data", "signatures", "signature_v1_locked.tsv")
coverage_path <- file.path(project_dir, "results", "bulk", "bulk_signature_gene_coverage.tsv")
random_summary_path <- file.path(project_dir, "results", "bulk", "bulk_random_gene_set_summary.tsv")
out_dir <- file.path(project_dir, "manuscript_outputs", "supplementary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

signature <- read_tsv(sig_path, show_col_types = FALSE) %>%
  mutate(
    gene_symbol = toupper(gene_symbol),
    include_in_primary_score = as.logical(include_in_primary_score),
    is_ecm_structural_gene = as.logical(is_ecm_structural_gene),
    include_in_ecm_excluded_score = as.logical(include_in_ecm_excluded_score)
  )

coverage <- read_tsv(coverage_path, show_col_types = FALSE) %>%
  filter(score_name == "primary") %>%
  select(accession, present_genes, missing_genes) %>%
  pivot_longer(c(present_genes, missing_genes), names_to = "status_field", values_to = "genes") %>%
  mutate(genes = str_split(coalesce(genes, ""), ";")) %>%
  unnest(genes) %>%
  mutate(
    genes = toupper(str_trim(genes)),
    detected = status_field == "present_genes"
  ) %>%
  filter(genes != "") %>%
  select(gene_symbol = genes, accession, detected) %>%
  distinct()

coverage_wide <- coverage %>%
  mutate(accession = paste0("detected_", accession)) %>%
  pivot_wider(names_from = accession, values_from = detected, values_fill = FALSE)

random_n <- if (file.exists(random_summary_path)) {
  random_summary <- read_tsv(random_summary_path, show_col_types = FALSE)
  paste(sort(unique(random_summary$n_random_sets)), collapse = "; ")
} else {
  "not available"
}

supp <- signature %>%
  left_join(coverage_wide, by = "gene_symbol") %>%
  mutate(
    across(starts_with("detected_"), ~ if_else(is.na(.x), FALSE, .x)),
    primary_score_member = if_else(include_in_primary_score, "yes", "no"),
    ecm_excluded_score_member = if_else(include_in_ecm_excluded_score, "yes", "no"),
    structural_ecm = if_else(is_ecm_structural_gene, "yes", "no")
  ) %>%
  transmute(
    gene_symbol,
    source_category,
    evidence_tier,
    evidence_note,
    primary_score_member,
    ecm_excluded_score_member,
    structural_ecm,
    lock_status,
    lock_date,
    detected_GSE135251,
    detected_GSE130970,
    detected_GSE162694,
    detected_GSE126848,
    detected_GSE167523,
    hgnc_status
  ) %>%
  arrange(source_category, gene_symbol)

summary <- tibble(
  item = c(
    "locked_gene_count",
    "primary_score_genes",
    "ecm_excluded_score_genes",
    "structural_ecm_genes",
    "lock_date",
    "random_gene_sets_per_fibrosis_mapped_cohort"
  ),
  value = c(
    nrow(supp),
    sum(supp$primary_score_member == "yes"),
    sum(supp$ecm_excluded_score_member == "yes"),
    sum(supp$structural_ecm == "yes"),
    paste(unique(supp$lock_date), collapse = "; "),
    random_n
  )
)

tsv_out <- file.path(out_dir, "Supplementary_Table_1_locked_97_gene_program.tsv")
xlsx_out <- file.path(out_dir, "Supplementary_Table_1_locked_97_gene_program.xlsx")
write_tsv(supp, tsv_out)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "locked_97_gene_program")
  openxlsx::writeData(wb, "locked_97_gene_program", supp)
  openxlsx::freezePane(wb, "locked_97_gene_program", firstRow = TRUE)
  openxlsx::setColWidths(wb, "locked_97_gene_program", cols = 1:ncol(supp), widths = "auto")

  openxlsx::addWorksheet(wb, "table_summary")
  openxlsx::writeData(wb, "table_summary", summary)
  openxlsx::freezePane(wb, "table_summary", firstRow = TRUE)
  openxlsx::setColWidths(wb, "table_summary", cols = 1:2, widths = "auto")

  openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
} else {
  write_tsv(summary, file.path(out_dir, "Supplementary_Table_1_locked_97_gene_program_summary.tsv"))
}

message("Supplementary Table 1 written:")
message(tsv_out)
if (file.exists(xlsx_out)) message(xlsx_out)

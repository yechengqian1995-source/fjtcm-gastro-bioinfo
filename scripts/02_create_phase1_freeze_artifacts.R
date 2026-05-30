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

metadata_dir <- file.path(project_dir, "data", "geo_metadata")
processed_dir <- file.path(project_dir, "data", "processed")
docs_dir <- file.path(project_dir, "docs")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_tsv(file.path(metadata_dir, "metadata_master_initial.tsv"), show_col_types = FALSE)
sample_counts <- read_tsv(file.path(metadata_dir, "geo_sample_counts.tsv"), show_col_types = FALSE)
candidate_fields <- read_tsv(file.path(metadata_dir, "metadata_candidate_fields.tsv"), show_col_types = FALSE)
hard_gates <- read_tsv(file.path(metadata_dir, "metadata_hard_gates.tsv"), show_col_types = FALSE)

choose_first <- function(dat, cols) {
  available <- cols[cols %in% names(dat)]
  if (!length(available)) return(rep(NA_character_, nrow(dat)))
  out <- rep(NA_character_, nrow(dat))
  for (col in available) {
    vals <- as.character(dat[[col]])
    use <- (is.na(out) | !nzchar(out)) & !is.na(vals) & nzchar(vals)
    out[use] <- vals[use]
  }
  out
}

metadata_freeze <- metadata
metadata_freeze$disease_label_raw <- choose_first(metadata_freeze, c("disease", "disease_status", "disease_state", "disease_subtype", "diagnosis", "group", "group_in_paper", "source_name_ch1", "title"))
metadata_freeze$fibrosis_stage_raw <- choose_first(metadata_freeze, c("fibrosis_stage", "fibrosis", "stage"))
metadata_freeze$nas_raw <- choose_first(metadata_freeze, c("nas_score", "nafld_activity_score", "nas"))
metadata_freeze$steatosis_raw <- choose_first(metadata_freeze, c("steatosis_grade", "steatosis", "fat"))
metadata_freeze$inflammation_raw <- choose_first(metadata_freeze, c("lobular_inflammation_grade", "lobular_inflammation_severity", "inflammation"))
metadata_freeze$ballooning_raw <- choose_first(metadata_freeze, c("cytological_ballooning_grade", "ballooning_intensity"))
metadata_freeze$donor_or_patient_raw <- choose_first(metadata_freeze, c("patient_id", "other_id", "title"))
metadata_freeze$age_raw <- choose_first(metadata_freeze, c("age", "age_y", "age_at_biopsy"))
metadata_freeze$sex_raw <- choose_first(metadata_freeze, c("sex", "gender"))
metadata_freeze$bmi_raw <- choose_first(metadata_freeze, c("bmi", "body_mass_index_kg_m2"))
metadata_freeze$diabetes_raw <- choose_first(metadata_freeze, c("diabetes"))
metadata_freeze$tissue_raw <- choose_first(metadata_freeze, c("tissue", "tissue_type", "source_name_ch1"))
metadata_freeze$platform_raw <- choose_first(metadata_freeze, c("platform_id", "instrument_model"))

metadata_freeze <- metadata_freeze %>%
  mutate(
    statistical_unit_rule = case_when(
      str_detect(accession, "GSE202379|GSE212837|GSE189600|GSE289173") ~ "donor/patient-level pseudobulk or mixed model; cells/nuclei are visualization units only",
      str_detect(accession, "GSE292268|GSE248077") ~ "patient/section-level summaries or mixed model; spots/ROIs/sections are nested observations",
      TRUE ~ "sample or participant-level analysis after duplicate/overlap audit"
    ),
    manual_correction_flag = case_when(
      accession == "GSE202379" & gsm %in% c("GSM6112262", "GSM6112263") ~ "GEO 2025-06-06 metadata-swap watch: manually confirm updated labels before analysis",
      accession == "GSE212837" ~ "split human/mouse and resolve donor-region nesting before analysis",
      accession == "GSE292268" ~ "F2-only human spatial dataset: localization only",
      accession == "GSE248077" ~ "mouse-only spatial context; citation still requires manual verification",
      TRUE ~ "none"
    ),
    metadata_freeze_version = "phase1_v1_2026-05-23",
    freeze_source = "GEO FTP series_matrix.txt.gz"
  ) %>%
  select(
    metadata_freeze_version, accession, matrix_file, sample_index, gsm, title,
    disease_label_raw, fibrosis_stage_raw, nas_raw, steatosis_raw,
    inflammation_raw, ballooning_raw, donor_or_patient_raw, age_raw, sex_raw,
    bmi_raw, diabetes_raw, tissue_raw, platform_raw, statistical_unit_rule,
    manual_correction_flag, freeze_source, everything()
  )

write_tsv(metadata_freeze, file.path(processed_dir, "metadata_freeze_v1.tsv"))

field_qc <- candidate_fields %>%
  mutate(has_candidate = !is.na(candidate_column) & n_nonmissing > 0) %>%
  summarise(
    candidate_columns = paste(unique(candidate_column[has_candidate]), collapse = "; "),
    n_candidate_columns = n_distinct(candidate_column[has_candidate]),
    n_nonmissing_max = if (any(has_candidate)) max(n_nonmissing[has_candidate], na.rm = TRUE) else 0L,
    example_values = paste(unique(examples[has_candidate & !is.na(examples)]), collapse = " | "),
    .by = c(accession, keyword)
  ) %>%
  arrange(accession, keyword)
write_tsv(field_qc, file.path(metadata_dir, "metadata_field_qc.tsv"))

accession_qc <- metadata_freeze %>%
  summarise(
    n_geo_records = n_distinct(gsm),
    n_records_with_disease_label = sum(!is.na(disease_label_raw) & nzchar(disease_label_raw)),
    n_records_with_fibrosis = sum(!is.na(fibrosis_stage_raw) & nzchar(fibrosis_stage_raw)),
    n_records_with_nas = sum(!is.na(nas_raw) & nzchar(nas_raw)),
    n_records_with_donor_or_patient = sum(!is.na(donor_or_patient_raw) & nzchar(donor_or_patient_raw)),
    n_manual_correction_flags = sum(manual_correction_flag != "none"),
    .by = accession
  ) %>%
  left_join(sample_counts %>% select(accession, layer, species, planned_role, priority, verification_status), by = "accession") %>%
  arrange(desc(priority == "high"), accession)
write_tsv(accession_qc, file.path(metadata_dir, "metadata_accession_qc.tsv"))

md_path <- file.path(docs_dir, "metadata_freeze_qc.md")
cat(
  "# Metadata Freeze QC\n\n",
  "Freeze version: `phase1_v1_2026-05-23`\n\n",
  "Source: GEO FTP `series_matrix.txt.gz` files parsed directly by `scripts/00_freeze_geo_metadata.R`, then standardized by `scripts/02_create_phase1_freeze_artifacts.R`.\n\n",
  "## Output Tables\n\n",
  "- `data/processed/metadata_freeze_v1.tsv`: one row per GEO sample record or matrix sample record.\n",
  "- `data/geo_metadata/metadata_accession_qc.tsv`: accession-level metadata completeness summary.\n",
  "- `data/geo_metadata/metadata_field_qc.tsv`: keyword-based candidate field report.\n",
  "- `data/geo_metadata/metadata_hard_gates.tsv`: non-negotiable metadata and scope gates.\n\n",
  "## Accessions Parsed\n\n",
  paste0(
    "- `", accession_qc$accession, "`: ",
    accession_qc$n_geo_records, " GEO sample records; ",
    accession_qc$n_records_with_fibrosis, " records with a candidate fibrosis field; ",
    accession_qc$n_records_with_nas, " records with a candidate NAS field; ",
    accession_qc$n_records_with_donor_or_patient, " records with a candidate donor or patient field.\n",
    collapse = ""
  ),
  "\n## Hard Gates\n\n",
  paste0("- `", hard_gates$accession, "` [", hard_gates$alert_level, "]: ", hard_gates$alert, "\n", collapse = ""),
  "\n## Interpretation Boundary\n\n",
  "This freeze verifies source-level sample metadata availability. It does not yet verify expression matrices, cell annotations, supplementary phenotype files, sample overlap, independent donor counts, or statistical model readiness.\n",
  file = md_path,
  sep = ""
)

message("Done. Freeze artifacts written to data/processed and data/geo_metadata.")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

meta_path <- file.path(project_dir, "data", "processed", "metadata_freeze_v1.tsv")
result_dir <- file.path(project_dir, "results", "single_spatial")
fig_dir <- file.path(project_dir, "figures", "single_spatial_readiness")
doc_dir <- file.path(project_dir, "docs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_tsv(meta_path, show_col_types = FALSE)

target_acc <- c("GSE202379", "GSE212837", "GSE289173", "GSE292268", "GSE248077", "GSE189600")
supp_cols <- names(metadata)[str_detect(names(metadata), "^supplementary_file")]

links <- metadata %>%
  filter(accession %in% target_acc) %>%
  select(accession, gsm, title, organism_ch1, platform_id, disease_label_raw, fibrosis_stage_raw,
         tissue_raw, gse202379_metadata_swap_watch, all_of(supp_cols)) %>%
  pivot_longer(all_of(supp_cols), names_to = "supplementary_field", values_to = "url") %>%
  filter(!is.na(url), url != "NA", url != "") %>%
  mutate(
    file_name = basename(url),
    file_type = case_when(
      str_detect(file_name, "\\.h5ad(\\.gz)?$") ~ "h5ad",
      str_detect(file_name, "\\.h5(\\.gz)?$") ~ "h5",
      str_detect(file_name, "barcodes\\.tsv\\.gz$") ~ "10x_barcodes",
      str_detect(file_name, "features\\.tsv\\.gz$|genes\\.tsv\\.gz$") ~ "10x_features",
      str_detect(file_name, "matrix\\.mtx\\.gz$") ~ "10x_matrix",
      str_detect(file_name, "counts\\.csv\\.gz$|norm_counts\\.csv\\.gz$|raw_counts\\.csv\\.gz$") ~ "counts_csv",
      str_detect(file_name, "tx_file\\.csv\\.gz$") ~ "cosmx_tx_file",
      str_detect(file_name, "jpg\\.gz$|png\\.gz$|tif(f)?\\.gz$") ~ "image",
      TRUE ~ "other"
    )
  )
write_tsv(links, file.path(result_dir, "single_spatial_supplementary_links.tsv"))

file_summary <- links %>%
  count(accession, file_type, name = "n_files") %>%
  arrange(accession, file_type)
write_tsv(file_summary, file.path(result_dir, "single_spatial_file_type_summary.tsv"))

sample_summary <- metadata %>%
  filter(accession %in% target_acc) %>%
  group_by(accession) %>%
  summarise(
    n_geo_records = n(),
    organisms = paste(sort(unique(na.omit(organism_ch1))), collapse = ";"),
    platforms = paste(sort(unique(na.omit(platform_id))), collapse = ";"),
    n_human_records = sum(organism_ch1 == "Homo sapiens", na.rm = TRUE),
    n_mouse_records = sum(organism_ch1 == "Mus musculus", na.rm = TRUE),
    n_with_fibrosis_stage = sum(!is.na(fibrosis_stage_raw)),
    n_gse202379_swap_watch = sum(gse202379_metadata_swap_watch, na.rm = TRUE),
    disease_labels = paste(sort(unique(na.omit(disease_label_raw))), collapse = ";"),
    .groups = "drop"
  ) %>%
  left_join(
    links %>% count(accession, name = "n_supplementary_links"),
    by = "accession"
  ) %>%
  mutate(n_supplementary_links = replace_na(n_supplementary_links, 0L))
write_tsv(sample_summary, file.path(result_dir, "single_spatial_sample_summary.tsv"))

readiness <- sample_summary %>%
  mutate(
    expression_links_available = n_supplementary_links > 0,
    human_records_available = n_human_records > 0,
    mouse_or_species_split_needed = n_mouse_records > 0 & n_human_records > 0,
    fibrosis_stage_available = n_with_fibrosis_stage > 0,
    gse202379_swap_clear = accession != "GSE202379" | n_gse202379_swap_watch == 0,
    localization_only = accession %in% c("GSE292268", "GSE248077"),
    donor_or_section_nesting_required = accession %in% c("GSE212837", "GSE292268", "GSE248077", "GSE202379", "GSE289173"),
    ready_for_biological_claims = case_when(
      accession == "GSE202379" ~ FALSE,
      accession == "GSE212837" ~ FALSE,
      accession == "GSE292268" ~ FALSE,
      TRUE ~ FALSE
    )
  ) %>%
  select(
    accession,
    expression_links_available,
    human_records_available,
    mouse_or_species_split_needed,
    fibrosis_stage_available,
    gse202379_swap_clear,
    donor_or_section_nesting_required,
    localization_only,
    ready_for_biological_claims
  ) %>%
  pivot_longer(-accession, names_to = "gate", values_to = "status") %>%
  mutate(
    gate_label = recode(
      gate,
      expression_links_available = "Expression links",
      human_records_available = "Human records",
      mouse_or_species_split_needed = "Species split needed",
      fibrosis_stage_available = "Fibrosis stage field",
      gse202379_swap_clear = "GSE202379 swap clear",
      donor_or_section_nesting_required = "Nesting required",
      localization_only = "Localization only",
      ready_for_biological_claims = "Ready for claims"
    ),
    gate_label = factor(gate_label, levels = c(
      "Expression links", "Human records", "Species split needed",
      "Fibrosis stage field", "GSE202379 swap clear", "Nesting required",
      "Localization only", "Ready for claims"
    )),
    display = case_when(
      gate %in% c("mouse_or_species_split_needed", "donor_or_section_nesting_required", "localization_only") & status ~ "requires handling",
      status ~ "yes",
      TRUE ~ "no"
    ),
    fill_group = case_when(
      gate == "ready_for_biological_claims" & !status ~ "blocked",
      gate %in% c("mouse_or_species_split_needed", "donor_or_section_nesting_required", "localization_only") & status ~ "caution",
      status ~ "available",
      TRUE ~ "not_available"
    )
  )
write_tsv(readiness, file.path(result_dir, "single_spatial_readiness_gates.tsv"))

theme_readiness <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "#333333"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold")
  )

p <- ggplot(readiness, aes(x = accession, y = gate_label, fill = fill_group)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = display), size = 2.4, color = "#111111") +
  scale_fill_manual(
    values = c(
      available = "#007A87",
      caution = "#E69F00",
      not_available = "#B8B8B8",
      blocked = "#D55E00"
    ),
    breaks = c("available", "caution", "not_available", "blocked"),
    labels = c("Available", "Requires handling", "Not available", "Blocked")
  ) +
  labs(
    title = "Single-nucleus and spatial readiness gates",
    subtitle = "Expression links exist, but claim-level use is blocked until nesting, species, and metadata gates are resolved",
    fill = "Gate status"
  ) +
  theme_readiness

ggsave(file.path(fig_dir, "single_spatial_readiness_gates.pdf"), p, width = 7.2, height = 3.8, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "single_spatial_readiness_gates.png"), p, width = 7.2, height = 3.8, units = "in", dpi = 600, bg = "white")

cat(
  "# Single-Nucleus and Spatial Readiness Gate\n\n",
  "Generated by `scripts/07_single_spatial_readiness.R`.\n\n",
  "## Outputs\n\n",
  "- `results/single_spatial/single_spatial_supplementary_links.tsv`\n",
  "- `results/single_spatial/single_spatial_file_type_summary.tsv`\n",
  "- `results/single_spatial/single_spatial_sample_summary.tsv`\n",
  "- `results/single_spatial/single_spatial_readiness_gates.tsv`\n",
  "- `figures/single_spatial_readiness/single_spatial_readiness_gates.png` and `.pdf`\n\n",
  "## Chief-Agent Decision\n\n",
  "Single-nucleus and spatial expression links are present in GEO metadata, but these datasets are not yet released for biological claim figures in this local project state.\n\n",
  "Key gates:\n\n",
  "- `GSE202379` remains blocked for disease-group statistics until the `GSM6112262` and `GSM6112263` metadata-swap issue is manually resolved.\n",
  "- `GSE212837` requires human/mouse split and donor-region nesting before HSC-lineage inference.\n",
  "- `GSE292268` is suitable only for localization because the local metadata indicate limited human MASLD spatial records and F2-only context from the project inventory.\n",
  "- `GSE248077` is mouse supplementary spatial context only.\n",
  "- `GSE289173` is a promising external single-nucleus candidate but fibrosis and access details need manual verification.\n\n",
  "## Interpretation Boundary\n\n",
  "No single-cell, single-nucleus, or spatial result should be described as completed validation until expression matrices are downloaded, parsed, QC'd, annotated, and analyzed with donor/section-aware statistical units.\n",
  file = file.path(doc_dir, "single_spatial_readiness_gate.md"),
  sep = ""
)

message("Done. Single/spatial readiness outputs written to ", result_dir)

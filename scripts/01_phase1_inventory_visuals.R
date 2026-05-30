#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

inventory_path <- file.path(project_dir, "data", "data_inventory.tsv")
metadata_dir <- file.path(project_dir, "data", "geo_metadata")
fig_dir <- file.path(project_dir, "figures", "phase1")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

inventory <- read_tsv(inventory_path, show_col_types = FALSE)

sample_counts_path <- file.path(metadata_dir, "geo_sample_counts.tsv")
field_path <- file.path(metadata_dir, "metadata_candidate_fields.tsv")

read_tsv_if_nonempty <- function(path, fallback) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(fallback)
  }
  dat <- tryCatch(read_tsv(path, show_col_types = FALSE), error = function(e) fallback)
  dat
}

sample_counts <- if (file.exists(sample_counts_path) && file.info(sample_counts_path)$size > 0) {
  read_tsv(sample_counts_path, show_col_types = FALSE)
} else {
  tibble(accession = inventory$accession, n_geo_samples = NA_integer_) %>%
    left_join(inventory %>% select(accession, layer, species, planned_role, priority, verification_status), by = "accession")
}

field_report <- read_tsv_if_nonempty(
  field_path,
  tibble(accession = character(), keyword = character(), candidate_column = character(), n_nonmissing = integer(), examples = character())
)

required_field_cols <- c("accession", "keyword", "candidate_column", "n_nonmissing", "examples")
missing_field_cols <- setdiff(required_field_cols, names(field_report))
if (length(missing_field_cols)) {
  for (col in missing_field_cols) field_report[[col]] <- NA
}
field_report <- field_report %>% select(all_of(required_field_cols))

hard_gate_path <- file.path(metadata_dir, "metadata_hard_gates.tsv")
hard_gates <- read_tsv_if_nonempty(
  hard_gate_path,
  tibble(accession = character(), alert_level = character(), alert = character())
)

status_path <- file.path(metadata_dir, "geo_fetch_status.tsv")
fetch_status <- read_tsv_if_nonempty(
  status_path,
  tibble(accession = sample_counts$accession, fetch_ok = NA, n_matrix_files = NA_integer(), n_geo_samples_parsed = sample_counts$n_geo_samples)
)

theme_bio <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

priority_cols <- c(high = "#D55E00", medium = "#0072B2", low = "#999999")
layer_cols <- c(
  "bulk RNA-seq" = "#0072B2",
  "microarray" = "#56B4E9",
  "snRNA-seq" = "#D55E00",
  "snRNA/snATAC/bulk" = "#CC79A7",
  "spatial transcriptomics" = "#009E73",
  "spatial multi-omics" = "#7B3294",
  "reference/experiment" = "#E69F00",
  "template" = "#4D4D4D"
)

keyword_levels <- c("disease", "fibrosis", "nas", "steatosis", "inflammation", "ballooning", "bmi", "diabetes", "age", "sex", "donor", "region", "platform", "tissue")

p1 <- inventory %>%
  mutate(layer = factor(layer, levels = names(layer_cols))) %>%
  count(layer, priority, name = "n") %>%
  ggplot(aes(x = layer, y = n, fill = priority)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = priority_cols, na.value = "grey70") +
  labs(
    title = "Phase 1 Candidate Dataset Inventory",
    subtitle = "Counts by assay layer and current priority; not biological results",
    x = NULL,
    y = "Number of entries",
    fill = "Priority"
  ) +
  theme_bio()
ggsave(file.path(fig_dir, "phase1_dataset_inventory.png"), p1, width = 8, height = 4.8, dpi = 300)
ggsave(file.path(fig_dir, "phase1_dataset_inventory.pdf"), p1, width = 8, height = 4.8)

p2 <- sample_counts %>%
  mutate(
    accession = factor(accession, levels = accession[order(n_geo_samples, decreasing = FALSE)]),
    layer = factor(layer, levels = names(layer_cols))
  ) %>%
  ggplot(aes(x = accession, y = n_geo_samples, fill = layer)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = layer_cols, na.value = "grey70") +
  labs(
    title = "GEO Sample Counts From Frozen Metadata",
    subtitle = "Counts reflect GEO sample records, not independent patients or donors",
    x = NULL,
    y = "GEO sample records",
    fill = "Layer"
  ) +
  theme_bio()
ggsave(file.path(fig_dir, "phase1_geo_sample_counts.png"), p2, width = 8, height = 5.6, dpi = 300)
ggsave(file.path(fig_dir, "phase1_geo_sample_counts.pdf"), p2, width = 8, height = 5.6)

field_presence <- field_report %>%
  filter(!is.na(candidate_column), n_nonmissing > 0) %>%
  distinct(accession, keyword) %>%
  mutate(present = TRUE) %>%
  right_join(
    expand_grid(
      accession = unique(sample_counts$accession),
      keyword = keyword_levels
    ),
    by = c("accession", "keyword")
  ) %>%
  mutate(
    present = replace_na(present, FALSE),
    accession = factor(accession, levels = rev(unique(sample_counts$accession))),
    keyword = factor(keyword, levels = keyword_levels)
  )

p3 <- ggplot(field_presence, aes(x = keyword, y = accession, fill = present)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#E6E6E6"), labels = c(`TRUE` = "Candidate field found", `FALSE` = "Not detected")) +
  labs(
    title = "Candidate Metadata Field Availability",
    subtitle = "Automated keyword screen; each field still requires manual interpretation",
    x = "Metadata keyword group",
    y = NULL,
    fill = NULL
  ) +
  theme_bio() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "phase1_metadata_field_availability.png"), p3, width = 9.5, height = 5.8, dpi = 300)
ggsave(file.path(fig_dir, "phase1_metadata_field_availability.pdf"), p3, width = 9.5, height = 5.8)

verification_levels <- c("metadata_partly_verified", "partially_verified", "needs_manual_verification")

p4 <- inventory %>%
  mutate(
    verification_status = factor(verification_status, levels = verification_levels),
    priority = factor(priority, levels = c("high", "medium", "low"))
  ) %>%
  count(priority, verification_status, name = "n") %>%
  ggplot(aes(x = verification_status, y = priority, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n), size = 3.4, fontface = "bold") +
  scale_fill_gradient(low = "#F2F2F2", high = "#0072B2") +
  labs(
    title = "Dataset Priority by Verification Status",
    subtitle = "Readiness matrix for Phase 1 curation; unresolved entries remain gated",
    x = "Verification status",
    y = "Priority",
    fill = "Entries"
  ) +
  theme_bio() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "phase1_priority_verification_matrix.png"), p4, width = 7.4, height = 4.6, dpi = 300)
ggsave(file.path(fig_dir, "phase1_priority_verification_matrix.pdf"), p4, width = 7.4, height = 4.6)

alert_counts <- hard_gates %>%
  count(alert_level, name = "n") %>%
  mutate(alert_level = if_else(is.na(alert_level), "none", alert_level))

p5 <- ggplot(alert_counts, aes(x = reorder(alert_level, n), y = n, fill = alert_level)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.2) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.4, fontface = "bold") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(hard_gate = "#D55E00", scope_limit = "#E69F00", none = "#999999"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Phase 1 Metadata Gates",
    subtitle = "Hard gates and scope limits to resolve before biological result figures",
    x = NULL,
    y = "Number of gate notes"
  ) +
  theme_bio()
ggsave(file.path(fig_dir, "phase1_metadata_gate_counts.png"), p5, width = 6.4, height = 3.8, dpi = 300)
ggsave(file.path(fig_dir, "phase1_metadata_gate_counts.pdf"), p5, width = 6.4, height = 3.8)

fetch_plot <- fetch_status %>%
  select(accession, fetch_ok, n_matrix_files, n_geo_samples_parsed) %>%
  right_join(sample_counts %>% select(accession, layer, priority), by = "accession") %>%
  mutate(
    fetch_label = case_when(
      isTRUE(fetch_ok) ~ "parsed",
      isFALSE(fetch_ok) ~ "failed",
      TRUE ~ "not checked"
    ),
    accession = factor(accession, levels = accession[order(replace_na(n_geo_samples_parsed, 0), decreasing = FALSE)])
  )

p6 <- ggplot(fetch_plot, aes(x = accession, y = replace_na(n_geo_samples_parsed, 0), fill = fetch_label)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c(parsed = "#009E73", failed = "#D55E00", `not checked` = "#999999")) +
  labs(
    title = "GEO Series Matrix Parsing Status",
    subtitle = "Parsed records are GEO sample records, not independent biological units",
    x = NULL,
    y = "Parsed GEO sample records",
    fill = "Status"
  ) +
  theme_bio()
ggsave(file.path(fig_dir, "phase1_geo_parse_status.png"), p6, width = 8, height = 5.6, dpi = 300)
ggsave(file.path(fig_dir, "phase1_geo_parse_status.pdf"), p6, width = 8, height = 5.6)

summary_path <- file.path(project_dir, "docs", "phase1_visual_outputs.md")
cat(
  "# Phase 1 Visual Outputs\n\n",
  "Generated with publication-style biomedical visualization principles: data-dense, explicit labels, and conservative evidence boundaries.\n\n",
  "## Files\n\n",
  "- `figures/phase1/phase1_dataset_inventory.png`\n",
  "- `figures/phase1/phase1_dataset_inventory.pdf`\n",
  "- `figures/phase1/phase1_geo_sample_counts.png`\n",
  "- `figures/phase1/phase1_geo_sample_counts.pdf`\n",
  "- `figures/phase1/phase1_metadata_field_availability.png`\n",
  "- `figures/phase1/phase1_metadata_field_availability.pdf`\n",
  "- `figures/phase1/phase1_priority_verification_matrix.png`\n",
  "- `figures/phase1/phase1_priority_verification_matrix.pdf`\n",
  "- `figures/phase1/phase1_metadata_gate_counts.png`\n",
  "- `figures/phase1/phase1_metadata_gate_counts.pdf`\n",
  "- `figures/phase1/phase1_geo_parse_status.png`\n",
  "- `figures/phase1/phase1_geo_parse_status.pdf`\n\n",
  "## Interpretation Boundary\n\n",
  "These figures summarize project readiness and metadata availability only. They do not provide biological findings, differential expression results, or validation evidence.\n",
  file = summary_path,
  sep = ""
)

message("Done. Figures written to ", fig_dir)

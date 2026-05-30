#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

result_dir <- file.path(project_dir, "results", "bulk")
figure_dir <- file.path(project_dir, "figures", "phase2_bulk")
source_dir <- file.path(result_dir, "figure_source_data")
doc_dir <- file.path(project_dir, "docs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

theme_publication <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "#333333", margin = margin(b = 6)),
      plot.caption = element_text(size = base_size - 1, color = "#555555", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#111111"),
      axis.line = element_line(linewidth = 0.35, color = "#111111"),
      axis.ticks = element_line(linewidth = 0.25, color = "#111111"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

score_label <- c(
  primary = "Primary locked score",
  ecm_excluded = "ECM-excluded score",
  structural_ecm = "Structural ECM-only",
  hsc_context = "HSC context"
)

cohort_role <- c(
  GSE135251 = "Discovery",
  GSE130970 = "Primary validation",
  GSE162694 = "Primary validation",
  GSE167523 = "Disease-state context",
  GSE126848 = "Disease-state context"
)

cohort_order <- c("GSE135251", "GSE130970", "GSE162694", "GSE167523", "GSE126848")

palette_status <- c(
  "Available" = "#007A87",
  "Partial" = "#E69F00",
  "Not available" = "#B8B8B8"
)
palette_score <- c(
  "Primary locked score" = "#0072B2",
  "ECM-excluded score" = "#009E73",
  "Structural ECM-only" = "#D55E00",
  "HSC context" = "#7B3294"
)
palette_advanced <- c(
  "non_advanced" = "#4D4D4D",
  "advanced" = "#D55E00"
)

endpoint <- read_tsv(file.path(result_dir, "bulk_endpoint_summary.tsv"), show_col_types = FALSE) %>%
  mutate(
    accession = factor(accession, levels = cohort_order),
    cohort_role = recode(as.character(accession), !!!cohort_role)
  )
coverage <- read_tsv(file.path(result_dir, "bulk_signature_gene_coverage.tsv"), show_col_types = FALSE) %>%
  mutate(
    accession = factor(accession, levels = cohort_order),
    score_label = recode(score_name, !!!score_label)
  )
scores <- read_tsv(file.path(result_dir, "bulk_program_scores.tsv"), show_col_types = FALSE) %>%
  mutate(
    accession = factor(accession, levels = cohort_order),
    cohort_role = recode(as.character(accession), !!!cohort_role),
    advanced_fibrosis = factor(advanced_fibrosis, levels = c("non_advanced", "advanced"))
  )
models <- read_tsv(file.path(result_dir, "bulk_model_results.tsv"), show_col_types = FALSE) %>%
  mutate(
    accession = factor(accession, levels = cohort_order),
    score_label = recode(score_name, !!!score_label)
  )
random_controls <- read_tsv(file.path(result_dir, "bulk_random_gene_set_controls.tsv"), show_col_types = FALSE) %>%
  mutate(accession = factor(accession, levels = cohort_order))
random_summary <- read_tsv(file.path(result_dir, "bulk_random_gene_set_summary.tsv"), show_col_types = FALSE) %>%
  mutate(accession = factor(accession, levels = cohort_order))

mapping_files <- list.files(result_dir, pattern = "_sample_mapping_qc\\.tsv$", full.names = TRUE)
mapping <- map_dfr(mapping_files, ~ read_tsv(.x, col_types = cols(.default = col_character()), show_col_types = FALSE)) %>%
  mutate(accession = factor(accession, levels = cohort_order))

mapping_summary <- mapping %>%
  group_by(accession) %>%
  summarise(
    n_mapping_samples = n(),
    n_matched = sum(metadata_match_status == "matched", na.rm = TRUE),
    mapping_fraction = n_matched / n_mapping_samples,
    .groups = "drop"
  )

primary_cov <- coverage %>%
  filter(score_name == "primary") %>%
  select(accession, primary_coverage = coverage_fraction)

readiness <- endpoint %>%
  left_join(mapping_summary, by = "accession") %>%
  left_join(primary_cov, by = "accession") %>%
  mutate(
    expression_loaded = n_samples > 0,
    sample_id_matched = mapping_fraction >= 0.95,
    fibrosis_endpoint_mapped = n_with_fibrosis > 0,
    advanced_endpoint_available = n_advanced > 0 & n_non_advanced > 0,
    nas_mapped = n_with_nas > 0,
    age_sex_available = n_with_age > 0 & n_with_sex > 0,
    signature_coverage_available = primary_coverage >= 0.90,
    disease_state_context = n_with_fibrosis == 0 & str_detect(disease_labels, "NAFL|NASH|healthy|obese")
  ) %>%
  select(
    accession, expression_loaded, sample_id_matched, signature_coverage_available,
    fibrosis_endpoint_mapped, advanced_endpoint_available, nas_mapped,
    age_sex_available, disease_state_context
  ) %>%
  pivot_longer(-accession, names_to = "metric", values_to = "available") %>%
  mutate(
    metric = recode(
      metric,
      expression_loaded = "Expression matrix loaded",
      sample_id_matched = "Sample IDs matched",
      signature_coverage_available = "Locked-gene coverage",
      fibrosis_endpoint_mapped = "Fibrosis endpoint",
      advanced_endpoint_available = "Advanced fibrosis split",
      nas_mapped = "NAS available",
      age_sex_available = "Age/sex covariates",
      disease_state_context = "Disease-state context"
    ),
    status = if_else(available, "Available", "Not available"),
    metric = factor(metric, levels = c(
      "Expression matrix loaded", "Sample IDs matched", "Locked-gene coverage",
      "Fibrosis endpoint", "Advanced fibrosis split", "NAS available",
      "Age/sex covariates", "Disease-state context"
    ))
  )

write_tsv(readiness, file.path(source_dir, "phase2_bulk_readiness_heatmap_source.tsv"))

p_readiness <- ggplot(readiness, aes(x = accession, y = metric, fill = status)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = if_else(status == "Available", "yes", "no")), size = 2.7, color = "#111111") +
  scale_fill_manual(values = palette_status) +
  labs(
    title = "Bulk cohort readiness for Phase 2 scoring",
    subtitle = "Readiness/QC view before biological interpretation",
    x = NULL,
    y = NULL,
    fill = "Status",
    caption = "Fibrosis-unmapped cohorts are used only as disease-state context unless additional phenotype fields are verified."
  ) +
  theme_publication(8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_readiness, "phase2_bulk_readiness_heatmap", 6.8, 3.8)

coverage_plot_data <- coverage %>%
  filter(score_name %in% c("primary", "ecm_excluded", "structural_ecm")) %>%
  mutate(
    score_label = factor(score_label, levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only")),
    coverage_label = paste0(n_present, "/", n_signature_genes)
  )
write_tsv(coverage_plot_data, file.path(source_dir, "phase2_signature_coverage_source.tsv"))

p_coverage <- ggplot(coverage_plot_data, aes(x = accession, y = coverage_fraction, fill = score_label)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68, color = "#222222", linewidth = 0.15) +
  geom_text(aes(label = coverage_label), position = position_dodge(width = 0.78), vjust = -0.35, size = 2.3) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = palette_score) +
  labs(
    title = "Coverage of locked prior genes",
    subtitle = "Coverage is a technical prerequisite, not validation of the program",
    x = NULL,
    y = "Gene coverage",
    fill = "Score set"
  ) +
  theme_publication(8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_coverage, "phase2_signature_gene_coverage", 6.8, 3.6)

fibrosis_scores <- scores %>%
  filter(!is.na(fibrosis_stage), !is.na(primary_score_z)) %>%
  mutate(
    fibrosis_stage_factor = factor(paste0("F", fibrosis_stage), levels = paste0("F", 0:4)),
    cohort_panel = paste0(accession, " (n=", ave(sample_id, accession, FUN = length), ")")
  )
write_tsv(fibrosis_scores, file.path(source_dir, "phase2_primary_score_by_fibrosis_source.tsv"))

p_stage <- ggplot(fibrosis_scores, aes(x = fibrosis_stage_factor, y = primary_score_z, color = advanced_fibrosis)) +
  geom_violin(aes(fill = advanced_fibrosis), alpha = 0.12, color = NA, trim = FALSE, linewidth = 0.2) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "#222222", linewidth = 0.3) +
  geom_jitter(width = 0.13, height = 0, alpha = 0.68, size = 0.85) +
  facet_wrap(~ accession, nrow = 1, scales = "free_x") +
  scale_color_manual(values = palette_advanced, na.value = "#777777") +
  scale_fill_manual(values = palette_advanced, na.value = "#777777") +
  labs(
    title = "Locked program score across fibrosis stages",
    subtitle = "Three cohorts with mapped fibrosis endpoints; score is standardized within cohort",
    x = "Fibrosis stage",
    y = "Primary score, within-cohort z",
    color = "Fibrosis group",
    fill = "Fibrosis group"
  ) +
  theme_publication(8)
save_plot(p_stage, "phase2_primary_score_by_fibrosis_stage", 7.2, 3.7)

advanced_scores <- scores %>%
  filter(!is.na(advanced_fibrosis)) %>%
  select(accession, sample_id, advanced_fibrosis, primary_score_z, ecm_excluded_score_z, structural_ecm_score_z) %>%
  pivot_longer(
    c(primary_score_z, ecm_excluded_score_z, structural_ecm_score_z),
    names_to = "score_name",
    values_to = "score_z"
  ) %>%
  mutate(
    score_name = str_remove(score_name, "_score_z"),
    score_label = factor(recode(score_name, !!!score_label), levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only"))
  )
write_tsv(advanced_scores, file.path(source_dir, "phase2_score_by_advanced_fibrosis_source.tsv"))

p_advanced <- ggplot(advanced_scores, aes(x = advanced_fibrosis, y = score_z, color = advanced_fibrosis)) +
  geom_boxplot(width = 0.34, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(width = 0.12, alpha = 0.62, size = 0.75) +
  facet_grid(score_label ~ accession, scales = "free_y") +
  scale_color_manual(values = palette_advanced) +
  labs(
    title = "Score distributions by advanced fibrosis endpoint",
    subtitle = "Primary and sensitivity scores use the same locked prior gene list",
    x = NULL,
    y = "Score, within-cohort z",
    color = "Fibrosis group"
  ) +
  theme_publication(7.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_advanced, "phase2_score_by_advanced_fibrosis", 7.2, 5.8)

forest_data <- models %>%
  filter(
    endpoint == "advanced_fibrosis_F3_F4_vs_F0_F2",
    model_type == "logistic_base",
    score_name %in% c("primary", "ecm_excluded", "structural_ecm")
  ) %>%
  mutate(
    score_label = factor(score_label, levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only")),
    accession = factor(accession, levels = rev(c("GSE135251", "GSE130970", "GSE162694"))),
    estimate_label = paste0("OR ", signif(estimate, 2), " [", signif(conf_low, 2), ", ", signif(conf_high, 2), "]")
  )
write_tsv(forest_data, file.path(source_dir, "phase2_advanced_fibrosis_forest_source.tsv"))

p_forest <- ggplot(forest_data, aes(x = estimate, y = accession, color = score_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#666666", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, linewidth = 0.45, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  scale_x_log10(breaks = c(1, 2, 4, 8, 16, 32), limits = c(0.8, 75)) +
  scale_color_manual(values = palette_score) +
  labs(
    title = "Within-cohort association with advanced fibrosis",
    subtitle = "Logistic models use continuous within-cohort z scores; no cutoff optimization was used",
    x = "Odds ratio per 1 SD score, log scale",
    y = NULL,
    color = "Score set",
    caption = "These estimates are association-level evidence; sample-overlap audit remains required before strong external-validation language."
  ) +
  theme_publication(8)
save_plot(p_forest, "phase2_advanced_fibrosis_forest", 7.5, 3.6)

endpoint_dot_data <- models %>%
  filter(
    endpoint %in% c("ordinal_fibrosis_stage", "NAS"),
    model_type == "spearman",
    score_name %in% c("primary", "ecm_excluded", "structural_ecm", "hsc_context")
  ) %>%
  mutate(
    endpoint_label = recode(
      endpoint,
      ordinal_fibrosis_stage = "Ordinal fibrosis stage",
      NAS = "NAS"
    ),
    score_label = factor(
      score_label,
      levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only", "HSC context")
    ),
    accession = factor(accession, levels = c("GSE135251", "GSE130970", "GSE162694")),
    sig_label = case_when(
      p_value < 0.001 ~ "P < 0.001",
      p_value < 0.01 ~ "P < 0.01",
      p_value < 0.05 ~ "P < 0.05",
      TRUE ~ "P >= 0.05"
    ),
    sig_label = factor(sig_label, levels = c("P < 0.001", "P < 0.01", "P < 0.05", "P >= 0.05"))
  )
write_tsv(endpoint_dot_data, file.path(source_dir, "phase2_endpoint_spearman_dotplot_source.tsv"))

p_endpoint <- ggplot(endpoint_dot_data, aes(x = accession, y = score_label)) +
  geom_point(aes(size = abs(estimate), fill = estimate, shape = sig_label), color = "#111111", stroke = 0.22) +
  facet_wrap(~ endpoint_label, nrow = 1) +
  scale_fill_gradient2(
    low = "#2D2A7B",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-0.65, 0.65),
    oob = squish,
    name = "Spearman rho"
  ) +
  scale_size_continuous(range = c(1.2, 5.2), limits = c(0, 0.65), name = "|rho|") +
  scale_shape_manual(values = c("P < 0.001" = 21, "P < 0.01" = 22, "P < 0.05" = 24, "P >= 0.05" = 1), drop = FALSE) +
  labs(
    title = "Ordinal endpoint correlations across fibrosis-mapped cohorts",
    subtitle = "Spearman correlations for fibrosis stage and NAS; scores are standardized within cohort",
    x = NULL,
    y = NULL,
    shape = "Nominal P"
  ) +
  theme_publication(8) +
  theme(panel.grid.major.x = element_line(color = "#F1F1F1", linewidth = 0.25))
save_plot(p_endpoint, "phase2_endpoint_spearman_dotplot", 7.4, 3.6)

random_plot_data <- random_controls %>%
  mutate(accession = factor(accession, levels = c("GSE135251", "GSE130970", "GSE162694")))
write_tsv(random_plot_data, file.path(source_dir, "phase2_random_gene_set_controls_source.tsv"))

p_random <- ggplot(random_plot_data, aes(x = beta)) +
  geom_density(fill = "#B8B8B8", color = "#4D4D4D", linewidth = 0.35, alpha = 0.75) +
  geom_vline(aes(xintercept = observed_beta), color = "#D55E00", linewidth = 0.65) +
  facet_wrap(~ accession, nrow = 1, scales = "free_y") +
  labs(
    title = "Expression-matched random gene-set controls",
    subtitle = "Orange line marks the observed primary-score beta in the advanced-fibrosis model",
    x = "Log-odds beta per 1 SD random score",
    y = "Density"
  ) +
  theme_publication(8)
save_plot(p_random, "phase2_random_gene_set_controls", 7.2, 3.3)

secondary_context <- scores %>%
  filter(accession %in% c("GSE126848", "GSE167523")) %>%
  mutate(
    context_group = case_when(
      accession == "GSE167523" & !is.na(disease_subtype_raw) ~ disease_subtype_raw,
      accession == "GSE126848" & !is.na(disease_label_raw) ~ disease_label_raw,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(context_group), !is.na(primary_score_z))
write_tsv(secondary_context, file.path(source_dir, "phase2_secondary_disease_state_context_source.tsv"))

p_secondary <- ggplot(secondary_context, aes(x = context_group, y = primary_score_z, color = context_group)) +
  geom_boxplot(width = 0.38, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(width = 0.12, alpha = 0.65, size = 0.8) +
  facet_wrap(~ accession, scales = "free_x") +
  scale_color_manual(values = c(
    healthy = "#0072B2", obese = "#E69F00", NAFLD = "#009E73",
    NASH = "#D55E00", NAFL = "#56B4E9"
  ), na.value = "#777777") +
  labs(
    title = "Secondary disease-state context cohorts",
    subtitle = "These cohorts lack mapped fibrosis stage in the current local files",
    x = NULL,
    y = "Primary score, within-cohort z",
    color = "Group"
  ) +
  theme_publication(8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_secondary, "phase2_secondary_disease_state_context", 6.8, 3.5)

figure_manifest <- tibble(
  file_stem = c(
    "phase2_bulk_readiness_heatmap",
    "phase2_signature_gene_coverage",
    "phase2_primary_score_by_fibrosis_stage",
    "phase2_score_by_advanced_fibrosis",
    "phase2_advanced_fibrosis_forest",
    "phase2_endpoint_spearman_dotplot",
    "phase2_random_gene_set_controls",
    "phase2_secondary_disease_state_context"
  ),
  role = c(
    "QC/readiness",
    "QC/signature coverage",
    "Primary bulk association visualization",
    "Sensitivity distribution visualization",
    "Within-cohort effect estimates",
    "Ordinal endpoint correlations",
    "Negative-control benchmark",
    "Secondary context only"
  ),
  boundary = c(
    "Shows whether cohorts are ready for scoring and endpoint modeling.",
    "Shows gene availability, not biological validation.",
    "Shows association with fibrosis stage in mapped cohorts.",
    "Compares primary, ECM-excluded, and structural ECM-only scores.",
    "Shows association-level ORs without cutoff optimization.",
    "Shows Spearman correlations with fibrosis stage and NAS.",
    "Benchmarks observed signal against expression-matched random gene sets.",
    "Does not support fibrosis interpretation."
  )
)
write_tsv(figure_manifest, file.path(source_dir, "phase2_figure_manifest.tsv"))

cat(
  "# Phase 2 Bulk Figure Manifest\n\n",
  "Generated by `scripts/06_phase2_bulk_visuals.R`.\n\n",
  "## Figure Files\n\n",
  paste0(
    "- `figures/phase2_bulk/", figure_manifest$file_stem, ".png` and `.pdf`: ",
    figure_manifest$role, ". ",
    figure_manifest$boundary,
    "\n",
    collapse = ""
  ),
  "\n## Claim Boundary\n\n",
  "These figures support cohort-level association and technical readiness claims only. They should not be described as proof of mechanical memory, causal HSC activation, diagnostic utility, prognostic utility, or therapeutic target validity.\n",
  file = file.path(doc_dir, "phase2_bulk_figure_manifest.md"),
  sep = ""
)

message("Done. Phase 2 bulk figures written to ", figure_dir)

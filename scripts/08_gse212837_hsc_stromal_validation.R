#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(curl)
  library(Matrix)
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

raw_dir <- file.path(project_dir, "data", "single_cell_raw", "GSE212837")
result_dir <- file.path(project_dir, "results", "single_cell")
fig_dir <- file.path(project_dir, "figures", "single_cell")
doc_dir <- file.path(project_dir, "docs")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

links_path <- file.path(project_dir, "results", "single_spatial", "single_spatial_supplementary_links.tsv")
if (!file.exists(links_path)) {
  links_path <- file.path(project_dir, "data", "geo_metadata", "single_spatial_supplementary_links.tsv")
}
if (!file.exists(links_path)) stop("Missing single/spatial supplementary links table.")

signature <- read_tsv(file.path(project_dir, "data", "signatures", "signature_v1_locked.tsv"), show_col_types = FALSE) %>%
  mutate(
    gene_symbol = toupper(gene_symbol),
    include_in_primary_score = as.logical(include_in_primary_score),
    include_in_ecm_excluded_score = as.logical(include_in_ecm_excluded_score),
    is_ecm_structural_gene = as.logical(is_ecm_structural_gene)
  )

gene_sets <- list(
  primary = signature %>% filter(include_in_primary_score) %>% pull(gene_symbol) %>% unique(),
  ecm_excluded = signature %>% filter(include_in_ecm_excluded_score) %>% pull(gene_symbol) %>% unique(),
  structural_ecm = signature %>% filter(is_ecm_structural_gene) %>% pull(gene_symbol) %>% unique(),
  hsc_context = signature %>% filter(str_detect(source_category, "HSC|portal_fibroblast")) %>% pull(gene_symbol) %>% unique(),
  stromal_marker = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "COL6A3", "PDGFRB", "RGS5", "ACTA2", "TAGLN", "DES", "LRAT", "RBP1"),
  quiescent_hsc = c("LRAT", "RBP1", "DES", "REELIN", "DCN"),
  activated_hsc = c("COL1A1", "COL1A2", "ACTA2", "TAGLN", "TIMP1", "MMP2", "POSTN", "LUM"),
  immune_marker = c("PTPRC", "LST1", "C1QA", "C1QB", "C1QC", "CD68", "TYROBP", "CD3D", "NKG7", "MS4A1"),
  endothelial_marker = c("PECAM1", "VWF", "KDR", "FLT1", "RAMP2", "CLEC4G", "STAB2"),
  hepatocyte_marker = c("ALB", "APOA1", "APOB", "TTR", "TF", "FGB", "CYP3A4"),
  cholangiocyte_marker = c("KRT19", "KRT7", "EPCAM", "SOX9", "MUC1")
)

download_if_missing <- function(url, dest) {
  url <- str_replace(url, "^ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading ", basename(dest))
    curl_download(url, dest, mode = "wb", quiet = FALSE)
  }
  dest
}

parse_sample_info <- function(title) {
  title <- as.character(title)
  disease_group <- case_when(
    str_detect(title, regex("Control", ignore_case = TRUE)) ~ "Control",
    str_detect(title, regex("NASH", ignore_case = TRUE)) ~ "NASH",
    TRUE ~ "Other"
  )
  donor_num <- str_match(title, regex("Liver\\s+([0-9]+)", ignore_case = TRUE))[, 2]
  donor_id <- case_when(
    disease_group == "Control" & !is.na(donor_num) ~ paste0("Control_", donor_num),
    disease_group == "NASH" & !is.na(donor_num) ~ paste0("NASH_", donor_num),
    TRUE ~ str_replace_all(title, "[^A-Za-z0-9]+", "_")
  )
  region <- str_match(title, "(R[0-9]+C[0-9]+)") [, 2]
  region[is.na(region)] <- "region_unspecified"
  tibble(disease_group = disease_group, donor_id = donor_id, region = region)
}

read_10x_sample <- function(sample_links) {
  gsm <- unique(sample_links$gsm)
  if (length(gsm) != 1) stop("Expected one GSM per sample group.")
  files <- sample_links %>%
    mutate(local_path = file.path(raw_dir, file_name))
  walk2(files$url, files$local_path, download_if_missing)

  matrix_path <- files$local_path[files$file_type == "10x_matrix"][[1]]
  features_path <- files$local_path[files$file_type == "10x_features"][[1]]
  barcodes_path <- files$local_path[files$file_type == "10x_barcodes"][[1]]

  mat <- Matrix::readMM(gzfile(matrix_path))
  features <- read_tsv(features_path, col_names = FALSE, show_col_types = FALSE)
  barcodes <- read_tsv(barcodes_path, col_names = FALSE, show_col_types = FALSE)[[1]]

  gene_symbols <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  gene_symbols <- toupper(make.unique(as.character(gene_symbols)))
  rownames(mat) <- gene_symbols
  colnames(mat) <- paste(gsm, barcodes, sep = "_")
  as(mat, "dgCMatrix")
}

score_gene_set <- function(mat, genes, totals = NULL) {
  if (is.null(totals)) totals <- Matrix::colSums(mat)
  totals <- as.numeric(totals)
  totals[!is.finite(totals) | totals <= 0] <- NA_real_
  present <- intersect(toupper(genes), rownames(mat))
  if (!length(present)) return(rep(NA_real_, ncol(mat)))
  gene_counts <- as.numeric(Matrix::colSums(mat[present, , drop = FALSE]))
  score <- log1p(gene_counts / totals * 1e4) / length(unique(toupper(genes)))
  score[!is.finite(score)] <- NA_real_
  score
}

count_detected <- function(mat, genes) {
  present <- intersect(toupper(genes), rownames(mat))
  if (!length(present)) return(rep(0L, ncol(mat)))
  Matrix::colSums(mat[present, , drop = FALSE] > 0)
}

process_sample <- function(sample_links) {
  gsm <- unique(sample_links$gsm)
  title <- unique(sample_links$title)
  info <- parse_sample_info(title)
  message("Processing ", gsm, " - ", title)

  mat <- read_10x_sample(sample_links)
  n_count <- Matrix::colSums(mat)
  n_feature <- Matrix::colSums(mat > 0)
  keep <- n_count >= 500 & n_feature >= 200
  mat <- mat[, keep, drop = FALSE]
  n_count <- n_count[keep]
  n_feature <- n_feature[keep]

  score_df <- tibble(
    cell_id = colnames(mat),
    n_count = as.numeric(n_count),
    n_feature = as.numeric(n_feature)
  )
  for (nm in names(gene_sets)) {
    score_df[[paste0(nm, "_score")]] <- score_gene_set(mat, gene_sets[[nm]], n_count)
  }

  stromal_tail_threshold <- suppressWarnings(
    quantile(score_df$stromal_marker_score, probs = 0.99, na.rm = TRUE, names = FALSE)
  )
  if (!is.finite(stromal_tail_threshold)) stromal_tail_threshold <- Inf

  score_df <- score_df %>%
    mutate(
      stromal_genes_detected = count_detected(mat, gene_sets$stromal_marker),
      activated_hsc_genes_detected = count_detected(mat, gene_sets$activated_hsc),
      quiescent_hsc_genes_detected = count_detected(mat, gene_sets$quiescent_hsc),
      max_nonstromal_score = pmax(immune_marker_score, endothelial_marker_score, hepatocyte_marker_score, cholangiocyte_marker_score, na.rm = TRUE),
      stromal_candidate = stromal_genes_detected >= 2 &
        stromal_marker_score > max_nonstromal_score,
      stromal_enriched_candidate = stromal_genes_detected >= 1 &
        !is.na(stromal_marker_score) &
        stromal_marker_score > 0 &
        stromal_marker_score >= stromal_tail_threshold,
      hsc_state_candidate = case_when(
        stromal_candidate & activated_hsc_score >= quiescent_hsc_score ~ "activated_stromal_like",
        stromal_candidate & activated_hsc_score < quiescent_hsc_score ~ "quiescent_stromal_like",
        TRUE ~ "non_stromal_or_unassigned"
      ),
      gsm = gsm,
      title = title,
      disease_group = info$disease_group[[1]],
      donor_id = info$donor_id[[1]],
      region = info$region[[1]]
    )

  sample_qc <- score_df %>%
    summarise(
      gsm = first(gsm),
      title = first(title),
      disease_group = first(disease_group),
      donor_id = first(donor_id),
      region = first(region),
      n_cells_qc = n(),
      median_n_count = median(n_count),
      median_n_feature = median(n_feature),
      n_stromal_candidate = sum(stromal_candidate),
      stromal_fraction = n_stromal_candidate / n_cells_qc,
      n_stromal_enriched_candidate = sum(stromal_enriched_candidate),
      stromal_enriched_fraction = n_stromal_enriched_candidate / n_cells_qc,
      stromal_enriched_threshold = stromal_tail_threshold,
      n_activated_stromal_like = sum(hsc_state_candidate == "activated_stromal_like"),
      n_quiescent_stromal_like = sum(hsc_state_candidate == "quiescent_stromal_like")
    )

  stromal_df <- score_df %>% filter(stromal_candidate)
  if (nrow(stromal_df)) {
    region_scores <- stromal_df %>%
      summarise(
        gsm = first(gsm),
        title = first(title),
        disease_group = first(disease_group),
        donor_id = first(donor_id),
        region = first(region),
        n_stromal_candidate = n(),
        primary_score = mean(primary_score, na.rm = TRUE),
        ecm_excluded_score = mean(ecm_excluded_score, na.rm = TRUE),
        structural_ecm_score = mean(structural_ecm_score, na.rm = TRUE),
        hsc_context_score = mean(hsc_context_score, na.rm = TRUE),
        stromal_marker_score = mean(stromal_marker_score, na.rm = TRUE),
        activated_hsc_score = mean(activated_hsc_score, na.rm = TRUE),
        quiescent_hsc_score = mean(quiescent_hsc_score, na.rm = TRUE)
      )
  } else {
    region_scores <- tibble(
      gsm = gsm,
      title = title,
      disease_group = info$disease_group[[1]],
      donor_id = info$donor_id[[1]],
      region = info$region[[1]],
      n_stromal_candidate = 0L,
      primary_score = NA_real_,
      ecm_excluded_score = NA_real_,
      structural_ecm_score = NA_real_,
      hsc_context_score = NA_real_,
      stromal_marker_score = NA_real_,
      activated_hsc_score = NA_real_,
      quiescent_hsc_score = NA_real_
    )
  }

  enriched_df <- score_df %>% filter(stromal_enriched_candidate)
  if (nrow(enriched_df)) {
    enriched_region_scores <- enriched_df %>%
      summarise(
        gsm = first(gsm),
        title = first(title),
        disease_group = first(disease_group),
        donor_id = first(donor_id),
        region = first(region),
        n_stromal_enriched_candidate = n(),
        primary_score = mean(primary_score, na.rm = TRUE),
        ecm_excluded_score = mean(ecm_excluded_score, na.rm = TRUE),
        structural_ecm_score = mean(structural_ecm_score, na.rm = TRUE),
        hsc_context_score = mean(hsc_context_score, na.rm = TRUE),
        stromal_marker_score = mean(stromal_marker_score, na.rm = TRUE),
        activated_hsc_score = mean(activated_hsc_score, na.rm = TRUE),
        quiescent_hsc_score = mean(quiescent_hsc_score, na.rm = TRUE)
      )
  } else {
    enriched_region_scores <- tibble(
      gsm = gsm,
      title = title,
      disease_group = info$disease_group[[1]],
      donor_id = info$donor_id[[1]],
      region = info$region[[1]],
      n_stromal_enriched_candidate = 0L,
      primary_score = NA_real_,
      ecm_excluded_score = NA_real_,
      structural_ecm_score = NA_real_,
      hsc_context_score = NA_real_,
      stromal_marker_score = NA_real_,
      activated_hsc_score = NA_real_,
      quiescent_hsc_score = NA_real_
    )
  }

  list(
    sample_qc = sample_qc,
    region_scores = region_scores,
    enriched_region_scores = enriched_region_scores
  )
}

links <- read_tsv(links_path, show_col_types = FALSE) %>%
  filter(accession == "GSE212837", organism_ch1 == "Homo sapiens", file_type %in% c("10x_matrix", "10x_features", "10x_barcodes"))

human_samples <- links %>%
  count(gsm, title, organism_ch1, name = "n_files") %>%
  filter(n_files == 3) %>%
  arrange(gsm)

if (!nrow(human_samples)) stop("No complete human 10x samples found for GSE212837.")

write_tsv(human_samples, file.path(result_dir, "gse212837_human_sample_manifest.tsv"))

outputs <- vector("list", nrow(human_samples))
for (i in seq_len(nrow(human_samples))) {
  sample_links <- links %>% filter(gsm == human_samples$gsm[[i]])
  outputs[[i]] <- process_sample(sample_links)
  gc()
}

sample_qc <- map_dfr(outputs, "sample_qc")
region_scores <- map_dfr(outputs, "region_scores")
enriched_region_scores <- map_dfr(outputs, "enriched_region_scores")

donor_scores <- region_scores %>%
  group_by(donor_id, disease_group) %>%
  summarise(
    n_regions = n(),
    n_stromal_candidate = sum(n_stromal_candidate),
    primary_score = mean(primary_score, na.rm = TRUE),
    ecm_excluded_score = mean(ecm_excluded_score, na.rm = TRUE),
    structural_ecm_score = mean(structural_ecm_score, na.rm = TRUE),
    hsc_context_score = mean(hsc_context_score, na.rm = TRUE),
    activated_hsc_score = mean(activated_hsc_score, na.rm = TRUE),
    quiescent_hsc_score = mean(quiescent_hsc_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    primary_score_z = as.numeric(scale(primary_score)),
    ecm_excluded_score_z = as.numeric(scale(ecm_excluded_score)),
    structural_ecm_score_z = as.numeric(scale(structural_ecm_score)),
    hsc_context_score_z = as.numeric(scale(hsc_context_score))
  )

enriched_donor_scores <- enriched_region_scores %>%
  group_by(donor_id, disease_group) %>%
  summarise(
    n_regions = n(),
    n_stromal_enriched_candidate = sum(n_stromal_enriched_candidate),
    primary_score = mean(primary_score, na.rm = TRUE),
    ecm_excluded_score = mean(ecm_excluded_score, na.rm = TRUE),
    structural_ecm_score = mean(structural_ecm_score, na.rm = TRUE),
    hsc_context_score = mean(hsc_context_score, na.rm = TRUE),
    activated_hsc_score = mean(activated_hsc_score, na.rm = TRUE),
    quiescent_hsc_score = mean(quiescent_hsc_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    primary_score_z = as.numeric(scale(primary_score)),
    ecm_excluded_score_z = as.numeric(scale(ecm_excluded_score)),
    structural_ecm_score_z = as.numeric(scale(structural_ecm_score)),
    hsc_context_score_z = as.numeric(scale(hsc_context_score))
  )

safe_wilcox <- function(dat, value_col) {
  tmp <- dat %>% filter(!is.na(.data[[value_col]]), disease_group %in% c("Control", "NASH"))
  if (n_distinct(tmp$disease_group) < 2 || min(table(tmp$disease_group)) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(as.formula(paste(value_col, "~ disease_group")), data = tmp, exact = FALSE)$p.value)
}

stats <- tibble(
  score_name = c("primary", "ecm_excluded", "structural_ecm", "hsc_context"),
  donor_level_wilcox_p = c(
    safe_wilcox(donor_scores, "primary_score"),
    safe_wilcox(donor_scores, "ecm_excluded_score"),
    safe_wilcox(donor_scores, "structural_ecm_score"),
    safe_wilcox(donor_scores, "hsc_context_score")
  )
) %>%
  mutate(
    interpretation_boundary = "donor-level marker-gated stromal validation; not proof of HSC mechanical memory"
  )

enriched_stats <- tibble(
  score_name = c("primary", "ecm_excluded", "structural_ecm", "hsc_context"),
  donor_level_wilcox_p = c(
    safe_wilcox(enriched_donor_scores, "primary_score"),
    safe_wilcox(enriched_donor_scores, "ecm_excluded_score"),
    safe_wilcox(enriched_donor_scores, "structural_ecm_score"),
    safe_wilcox(enriched_donor_scores, "hsc_context_score")
  )
) %>%
  mutate(
    interpretation_boundary = "donor-level stromal-marker-enriched sensitivity analysis; not definitive HSC annotation or proof of mechanical memory"
  )

write_tsv(sample_qc, file.path(result_dir, "gse212837_human_sample_qc.tsv"))
write_tsv(region_scores, file.path(result_dir, "gse212837_human_stromal_region_scores.tsv"))
write_tsv(donor_scores, file.path(result_dir, "gse212837_human_stromal_donor_scores.tsv"))
write_tsv(stats, file.path(result_dir, "gse212837_human_stromal_donor_stats.tsv"))
write_tsv(enriched_region_scores, file.path(result_dir, "gse212837_human_stromal_enriched_region_scores.tsv"))
write_tsv(enriched_donor_scores, file.path(result_dir, "gse212837_human_stromal_enriched_donor_scores.tsv"))
write_tsv(enriched_stats, file.path(result_dir, "gse212837_human_stromal_enriched_donor_stats.tsv"))

theme_publication <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "#333333"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#111111"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.25)
  )

score_plot_data <- donor_scores %>%
  filter(n_stromal_candidate > 0) %>%
  select(donor_id, disease_group, n_regions, n_stromal_candidate,
         primary_score_z, ecm_excluded_score_z, structural_ecm_score_z, hsc_context_score_z) %>%
  pivot_longer(ends_with("_score_z"), names_to = "score_name", values_to = "score_z") %>%
  mutate(
    score_name = recode(
      score_name,
      primary_score_z = "Primary locked score",
      ecm_excluded_score_z = "ECM-excluded score",
      structural_ecm_score_z = "Structural ECM-only",
      hsc_context_score_z = "HSC context"
    ),
    score_name = factor(score_name, levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only", "HSC context")),
    disease_group = factor(disease_group, levels = c("Control", "NASH"))
  )

p_score <- ggplot(score_plot_data, aes(x = disease_group, y = score_z, color = disease_group)) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(aes(size = n_stromal_candidate), width = 0.08, alpha = 0.75) +
  facet_wrap(~ score_name, nrow = 1) +
  scale_color_manual(values = c(Control = "#0072B2", NASH = "#D55E00")) +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    title = "Donor-level stromal program scores in GSE212837",
    subtitle = "Marker-gated human stromal nuclei; donors are the statistical units",
    x = NULL,
    y = "Donor-level score, z",
    color = "Group",
    size = "Stromal nuclei"
  ) +
  theme_publication
ggsave(file.path(fig_dir, "gse212837_donor_stromal_scores.pdf"), p_score, width = 7.2, height = 3.4, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse212837_donor_stromal_scores.png"), p_score, width = 7.2, height = 3.4, units = "in", dpi = 600, bg = "white")

enriched_score_plot_data <- enriched_donor_scores %>%
  select(donor_id, disease_group, n_regions, n_stromal_enriched_candidate,
         primary_score_z, ecm_excluded_score_z, structural_ecm_score_z, hsc_context_score_z) %>%
  pivot_longer(ends_with("_score_z"), names_to = "score_name", values_to = "score_z") %>%
  mutate(
    score_name = recode(
      score_name,
      primary_score_z = "Primary locked score",
      ecm_excluded_score_z = "ECM-excluded score",
      structural_ecm_score_z = "Structural ECM-only",
      hsc_context_score_z = "HSC context"
    ),
    score_name = factor(score_name, levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only", "HSC context")),
    disease_group = factor(disease_group, levels = c("Control", "NASH"))
  )

p_enriched_score <- ggplot(enriched_score_plot_data, aes(x = disease_group, y = score_z, color = disease_group)) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(aes(size = n_stromal_enriched_candidate), width = 0.08, alpha = 0.78) +
  facet_wrap(~ score_name, nrow = 1) +
  scale_color_manual(values = c(Control = "#0072B2", NASH = "#D55E00")) +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    title = "Donor-level stromal-enriched sensitivity scores in GSE212837",
    subtitle = "Top 1% stromal-marker-enriched nuclei per sample; donors are the statistical units",
    x = NULL,
    y = "Donor-level score, z",
    color = "Group",
    size = "Enriched nuclei"
  ) +
  theme_publication
ggsave(file.path(fig_dir, "gse212837_donor_stromal_enriched_sensitivity_scores.pdf"), p_enriched_score, width = 7.2, height = 3.4, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse212837_donor_stromal_enriched_sensitivity_scores.png"), p_enriched_score, width = 7.2, height = 3.4, units = "in", dpi = 600, bg = "white")

p_qc <- sample_qc %>%
  mutate(
    sample_label = paste0(donor_id, " ", region),
    disease_group = factor(disease_group, levels = c("Control", "NASH"))
  ) %>%
  ggplot(aes(x = reorder(sample_label, stromal_fraction), y = stromal_fraction, fill = disease_group)) +
  geom_col(color = "#222222", linewidth = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.5), breaks = seq(0, 0.03, 0.01), limits = c(0, 0.03)) +
  scale_fill_manual(values = c(Control = "#0072B2", NASH = "#D55E00")) +
  labs(
    title = "Marker-gated stromal nuclei fraction",
    subtitle = "QC view for GSE212837 human samples and regions",
    x = NULL,
    y = "Stromal candidate fraction",
    fill = "Group"
  ) +
  theme_publication
ggsave(file.path(fig_dir, "gse212837_stromal_fraction_qc.pdf"), p_qc, width = 6.4, height = 5.2, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse212837_stromal_fraction_qc.png"), p_qc, width = 6.4, height = 5.2, units = "in", dpi = 600, bg = "white")

p_enriched_qc <- sample_qc %>%
  mutate(
    sample_label = paste0(donor_id, " ", region),
    disease_group = factor(disease_group, levels = c("Control", "NASH"))
  ) %>%
  ggplot(aes(x = reorder(sample_label, stromal_enriched_fraction), y = stromal_enriched_fraction, fill = disease_group)) +
  geom_col(color = "#222222", linewidth = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_fill_manual(values = c(Control = "#0072B2", NASH = "#D55E00")) +
  labs(
    title = "Stromal-marker-enriched nuclei fraction",
    subtitle = "Top 1% per-sample marker-tail sensitivity gate",
    x = NULL,
    y = "Enriched candidate fraction",
    fill = "Group"
  ) +
  theme_publication
ggsave(file.path(fig_dir, "gse212837_stromal_enriched_fraction_qc.pdf"), p_enriched_qc, width = 6.4, height = 5.2, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse212837_stromal_enriched_fraction_qc.png"), p_enriched_qc, width = 6.4, height = 5.2, units = "in", dpi = 600, bg = "white")

cat(
  "# GSE212837 Donor-Aware HSC/Stromal Validation\n\n",
  "Generated by `scripts/08_gse212837_hsc_stromal_validation.R`.\n\n",
  "## Outputs\n\n",
  "- `results/single_cell/gse212837_human_sample_manifest.tsv`\n",
  "- `results/single_cell/gse212837_human_sample_qc.tsv`\n",
  "- `results/single_cell/gse212837_human_stromal_region_scores.tsv`\n",
  "- `results/single_cell/gse212837_human_stromal_donor_scores.tsv`\n",
  "- `results/single_cell/gse212837_human_stromal_donor_stats.tsv`\n",
  "- `results/single_cell/gse212837_human_stromal_enriched_region_scores.tsv`\n",
  "- `results/single_cell/gse212837_human_stromal_enriched_donor_scores.tsv`\n",
  "- `results/single_cell/gse212837_human_stromal_enriched_donor_stats.tsv`\n",
  "- `figures/single_cell/gse212837_donor_stromal_scores.png` and `.pdf`\n",
  "- `figures/single_cell/gse212837_stromal_fraction_qc.png` and `.pdf`\n\n",
  "- `figures/single_cell/gse212837_donor_stromal_enriched_sensitivity_scores.png` and `.pdf`\n",
  "- `figures/single_cell/gse212837_stromal_enriched_fraction_qc.png` and `.pdf`\n\n",
  "## Analysis Boundary\n\n",
  "This is a donor-aware marker-gated stromal validation. It uses nuclei-level marker scores only to define a candidate stromal compartment, then aggregates scores to donor-level summaries. It does not constitute a full single-nucleus atlas, definitive HSC annotation, or proof of mechanical memory. Because the high-confidence marker gate recovered no control stromal candidates, a separate top 1% stromal-marker-enriched per-sample sensitivity analysis is reported for donor-level score comparison.\n\n",
  "## High-Confidence Marker-Gated Donor-Level Statistics\n\n",
  paste0(
    "- `", stats$score_name, "` Wilcoxon P = ",
    ifelse(is.na(stats$donor_level_wilcox_p), "NA", signif(stats$donor_level_wilcox_p, 3)),
    ".\n",
    collapse = ""
  ),
  "\n## Stromal-Marker-Enriched Sensitivity Donor-Level Statistics\n\n",
  paste0(
    "- `", enriched_stats$score_name, "` Wilcoxon P = ",
    ifelse(is.na(enriched_stats$donor_level_wilcox_p), "NA", signif(enriched_stats$donor_level_wilcox_p, 3)),
    ".\n",
    collapse = ""
  ),
  file = file.path(doc_dir, "gse212837_hsc_stromal_validation_summary.md"),
  sep = ""
)

message("Done. GSE212837 donor-aware stromal outputs written to ", result_dir)

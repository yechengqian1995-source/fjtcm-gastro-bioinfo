#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(Seurat)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

raw_dir <- file.path(project_dir, "data", "single_cell_raw", "GSE212837")
result_dir <- file.path(project_dir, "results", "single_cell")
fig_dir <- file.path(project_dir, "figures", "single_cell")
doc_dir <- file.path(project_dir, "docs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(212837)

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

sample_manifest <- read_tsv(file.path(result_dir, "gse212837_human_sample_manifest.tsv"), show_col_types = FALSE)

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

score_gene_set <- function(mat, base_symbols, genes, totals = NULL) {
  if (is.null(totals)) totals <- Matrix::colSums(mat)
  totals <- as.numeric(totals)
  totals[!is.finite(totals) | totals <= 0] <- NA_real_
  rows <- which(base_symbols %in% toupper(genes))
  if (!length(rows)) return(rep(NA_real_, ncol(mat)))
  gene_counts <- as.numeric(Matrix::colSums(mat[rows, , drop = FALSE]))
  score <- log1p(gene_counts / totals * 1e4) / length(unique(toupper(genes)))
  score[!is.finite(score)] <- NA_real_
  score
}

count_detected <- function(mat, base_symbols, genes) {
  rows <- which(base_symbols %in% toupper(genes))
  if (!length(rows)) return(rep(0L, ncol(mat)))
  Matrix::colSums(mat[rows, , drop = FALSE] > 0)
}

read_10x_paths <- function(gsm) {
  files <- list.files(raw_dir, pattern = paste0("^", gsm, ".*_(barcodes|features|matrix)\\.tsv\\.gz$|^", gsm, ".*_matrix\\.mtx\\.gz$"), full.names = TRUE)
  barcode <- files[str_detect(basename(files), "barcodes\\.tsv\\.gz$")]
  feature <- files[str_detect(basename(files), "features\\.tsv\\.gz$")]
  matrix <- files[str_detect(basename(files), "matrix\\.mtx\\.gz$")]
  if (length(barcode) != 1 || length(feature) != 1 || length(matrix) != 1) {
    stop("Could not resolve 10x files for ", gsm)
  }
  list(barcode = barcode, feature = feature, matrix = matrix)
}

read_sketch_sample <- function(gsm, title, max_random_per_sample = 700) {
  message("Embedding sketch: ", gsm, " - ", title)
  paths <- read_10x_paths(gsm)
  mat <- Matrix::readMM(gzfile(paths$matrix))
  features <- read_tsv(paths$feature, col_names = FALSE, show_col_types = FALSE)
  barcodes <- read_tsv(paths$barcode, col_names = FALSE, show_col_types = FALSE)[[1]]
  base_symbols <- toupper(as.character(if (ncol(features) >= 2) features[[2]] else features[[1]]))
  rownames(mat) <- make.unique(base_symbols)
  colnames(mat) <- paste(gsm, barcodes, sep = "_")

  n_count <- Matrix::colSums(mat)
  n_feature <- Matrix::colSums(mat > 0)
  keep <- n_count >= 500 & n_feature >= 200
  mat <- mat[, keep, drop = FALSE]
  n_count <- n_count[keep]
  n_feature <- n_feature[keep]

  info <- parse_sample_info(title)
  meta <- tibble(
    cell_id = colnames(mat),
    gsm = gsm,
    title = title,
    disease_group = info$disease_group[[1]],
    donor_id = info$donor_id[[1]],
    region = info$region[[1]],
    n_count = as.numeric(n_count),
    n_feature = as.numeric(n_feature)
  )

  for (nm in names(gene_sets)) {
    meta[[paste0(nm, "_score")]] <- score_gene_set(mat, base_symbols, gene_sets[[nm]], n_count)
  }
  meta <- meta %>%
    mutate(
      stromal_genes_detected = count_detected(mat, base_symbols, gene_sets$stromal_marker),
      activated_hsc_genes_detected = count_detected(mat, base_symbols, gene_sets$activated_hsc),
      quiescent_hsc_genes_detected = count_detected(mat, base_symbols, gene_sets$quiescent_hsc),
      max_nonstromal_score = pmax(immune_marker_score, endothelial_marker_score, hepatocyte_marker_score, cholangiocyte_marker_score, na.rm = TRUE),
      stromal_candidate = stromal_genes_detected >= 2 & stromal_marker_score > max_nonstromal_score
    )
  tail_threshold <- suppressWarnings(quantile(meta$stromal_marker_score, probs = 0.99, na.rm = TRUE, names = FALSE))
  if (!is.finite(tail_threshold)) tail_threshold <- Inf
  meta <- meta %>%
    mutate(
      stromal_enriched_candidate = stromal_genes_detected >= 1 &
        !is.na(stromal_marker_score) &
        stromal_marker_score > 0 &
        stromal_marker_score >= tail_threshold,
      sketch_priority = stromal_candidate | stromal_enriched_candidate
    )

  priority_idx <- which(meta$sketch_priority)
  background_idx <- setdiff(seq_len(nrow(meta)), priority_idx)
  if (length(background_idx) > max_random_per_sample) {
    background_idx <- sample(background_idx, max_random_per_sample)
  }
  selected_idx <- sort(unique(c(priority_idx, background_idx)))

  mat_sel <- mat[, selected_idx, drop = FALSE]
  meta_sel <- meta[selected_idx, , drop = FALSE] %>%
    mutate(sketch_selection = if_else(sketch_priority, "marker-priority", "balanced-background"))

  rm(mat)
  gc(verbose = FALSE)
  list(counts = mat_sel, meta = meta_sel, base_symbols = base_symbols)
}

sample_objects <- pmap(
  list(sample_manifest$gsm, sample_manifest$title),
  ~ read_sketch_sample(..1, ..2)
)

feature_union <- reduce(map(sample_objects, ~ rownames(.x$counts)), union)
align_counts <- function(obj, feature_union) {
  mat <- obj$counts
  if (identical(rownames(mat), feature_union)) return(mat)
  out <- Matrix(0, nrow = length(feature_union), ncol = ncol(mat), sparse = TRUE)
  rownames(out) <- feature_union
  colnames(out) <- colnames(mat)
  idx <- match(rownames(mat), feature_union)
  out[idx, ] <- mat
  out
}

counts <- Reduce(Matrix::cbind2, map(sample_objects, ~ align_counts(.x, feature_union)))
meta <- bind_rows(map(sample_objects, "meta")) %>%
  as.data.frame()
rownames(meta) <- meta$cell_id

stopifnot(identical(colnames(counts), rownames(meta)))

obj <- CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 0, min.features = 0, project = "GSE212837_sketch")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.45, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

umap <- Embeddings(obj, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id")
names(umap)[2:3] <- c("UMAP_1", "UMAP_2")

meta_df <- obj@meta.data %>%
  as.data.frame()
if (!"cell_id" %in% names(meta_df)) {
  meta_df <- rownames_to_column(meta_df, "cell_id")
} else {
  meta_df$cell_id <- as.character(meta_df$cell_id)
}

cell_scores <- meta_df %>%
  left_join(umap, by = "cell_id") %>%
  mutate(
    broad_compartment = case_when(
      stromal_candidate ~ "High-confidence stromal",
      stromal_enriched_candidate ~ "Stromal-marker enriched",
      hepatocyte_marker_score >= pmax(immune_marker_score, endothelial_marker_score, cholangiocyte_marker_score, stromal_marker_score, na.rm = TRUE) ~ "Hepatocyte-like",
      immune_marker_score >= pmax(hepatocyte_marker_score, endothelial_marker_score, cholangiocyte_marker_score, stromal_marker_score, na.rm = TRUE) ~ "Immune-like",
      endothelial_marker_score >= pmax(hepatocyte_marker_score, immune_marker_score, cholangiocyte_marker_score, stromal_marker_score, na.rm = TRUE) ~ "Endothelial-like",
      cholangiocyte_marker_score >= pmax(hepatocyte_marker_score, immune_marker_score, endothelial_marker_score, stromal_marker_score, na.rm = TRUE) ~ "Cholangiocyte-like",
      TRUE ~ "Unassigned"
    ),
    disease_group = factor(disease_group, levels = c("Control", "NASH")),
    broad_compartment = factor(
      broad_compartment,
      levels = c("Hepatocyte-like", "Immune-like", "Endothelial-like", "Cholangiocyte-like", "High-confidence stromal", "Stromal-marker enriched", "Unassigned")
    )
  )

obj$broad_compartment <- cell_scores$broad_compartment[match(colnames(obj), cell_scores$cell_id)]
obj$disease_group <- cell_scores$disease_group[match(colnames(obj), cell_scores$cell_id)]

write_tsv(
  cell_scores %>%
    select(cell_id, gsm, disease_group, donor_id, region, n_count, n_feature,
           seurat_clusters, UMAP_1, UMAP_2, broad_compartment, sketch_selection,
           primary_score, ecm_excluded_score, structural_ecm_score, hsc_context_score,
           stromal_marker_score, activated_hsc_score, quiescent_hsc_score,
           stromal_candidate, stromal_enriched_candidate),
  file.path(result_dir, "gse212837_sketch_umap_cell_scores.tsv")
)

score_qc_columns <- c(
  "primary_score", "ecm_excluded_score", "structural_ecm_score", "hsc_context_score",
  "stromal_marker_score", "activated_hsc_score", "quiescent_hsc_score"
)
score_finite_qc <- map_dfr(score_qc_columns, function(score_col) {
  x <- cell_scores[[score_col]]
  tibble(
    score_column = score_col,
    n_values = length(x),
    n_finite = sum(is.finite(x)),
    n_na = sum(is.na(x)),
    n_nan = sum(is.nan(x)),
    n_inf = sum(is.infinite(x)),
    n_zero = sum(is.finite(x) & x == 0),
    n_nonzero = sum(is.finite(x) & x != 0),
    min = if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_,
    q25 = if (any(is.finite(x))) unname(quantile(x[is.finite(x)], 0.25, names = FALSE)) else NA_real_,
    median = if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_,
    q75 = if (any(is.finite(x))) unname(quantile(x[is.finite(x)], 0.75, names = FALSE)) else NA_real_,
    max = if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  )
})
write_tsv(score_finite_qc, file.path(result_dir, "gse212837_sketch_score_finite_qc.tsv"))

embedding_summary <- cell_scores %>%
  count(disease_group, donor_id, region, broad_compartment, name = "n_sketch_nuclei") %>%
  arrange(disease_group, donor_id, region, broad_compartment)
write_tsv(embedding_summary, file.path(result_dir, "gse212837_sketch_umap_summary.tsv"))

saveRDS(
  obj,
  file.path(result_dir, "gse212837_sketch_seurat_object.rds")
)

theme_umap <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "#333333"),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.14, "in"),
    plot.margin = margin(8, 8, 8, 8)
  )

compartment_colors <- c(
  "Hepatocyte-like" = "#0072B2",
  "Immune-like" = "#D55E00",
  "Endothelial-like" = "#009E73",
  "Cholangiocyte-like" = "#CC79A7",
  "High-confidence stromal" = "#E69F00",
  "Stromal-marker enriched" = "#7B3294",
  "Unassigned" = "#9E9E9E"
)

plot_umap_category <- function(data, color_col, title, subtitle, colors, out_stem) {
  p <- ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
    geom_point(size = 0.2, alpha = 0.78, stroke = 0) +
    scale_color_manual(values = colors, na.value = "#BDBDBD") +
    guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    theme_umap
  ggsave(file.path(fig_dir, paste0(out_stem, ".pdf")), p, width = 4.2, height = 3.6, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(out_stem, ".png")), p, width = 4.2, height = 3.6, units = "in", dpi = 600, bg = "white")
  p
}

plot_umap_score <- function(data, score_col, title, subtitle, out_stem) {
  p <- ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = .data[[score_col]])) +
    geom_point(size = 0.2, alpha = 0.84, stroke = 0) +
    scale_color_gradientn(colors = c("#CFCFCF", "#56B4E9", "#E69F00", "#C44E29"), na.value = "#EFEFEF") +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2", color = "Score") +
    theme_umap
  ggsave(file.path(fig_dir, paste0(out_stem, ".pdf")), p, width = 4.2, height = 3.6, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(out_stem, ".png")), p, width = 4.2, height = 3.6, units = "in", dpi = 600, bg = "white")
  p
}

p_compartment <- plot_umap_category(
  cell_scores,
  "broad_compartment",
  "GSE212837 sketch UMAP by marker context",
  "Balanced nuclei sketch with stromal-marker-priority retention",
  compartment_colors,
  "gse212837_sketch_umap_marker_context"
)

p_group <- plot_umap_category(
  cell_scores,
  "disease_group",
  "GSE212837 sketch UMAP by donor group",
  "Control and NASH nuclei shown for localization only",
  c(Control = "#0072B2", NASH = "#D55E00"),
  "gse212837_sketch_umap_disease_group"
)

p_cluster <- plot_umap_category(
  cell_scores %>% mutate(seurat_clusters = factor(seurat_clusters)),
  "seurat_clusters",
  "GSE212837 sketch UMAP by Seurat cluster",
  "Unsupervised clusters from the balanced sketch object",
  setNames(scales::hue_pal()(n_distinct(cell_scores$seurat_clusters)), sort(unique(as.character(cell_scores$seurat_clusters)))),
  "gse212837_sketch_umap_seurat_cluster"
)

p_primary <- plot_umap_score(
  cell_scores,
  "primary_score",
  "Primary program score on sketch UMAP",
  "Cell-level score map; donor summaries remain the inferential layer",
  "gse212837_sketch_umap_primary_score"
)

p_ecm_excluded <- plot_umap_score(
  cell_scores,
  "ecm_excluded_score",
  "ECM-excluded score on sketch UMAP",
  "Non-structural mechanostress-associated score context",
  "gse212837_sketch_umap_ecm_excluded_score"
)

p_stromal <- plot_umap_score(
  cell_scores,
  "stromal_marker_score",
  "Stromal-marker score on sketch UMAP",
  "Marker-tail retention supports stromal-context visualization",
  "gse212837_sketch_umap_stromal_marker_score"
)

p_hsc <- plot_umap_score(
  cell_scores,
  "hsc_context_score",
  "HSC-context score on sketch UMAP",
  "Context score projected onto the same sketch embedding",
  "gse212837_sketch_umap_hsc_context_score"
)

marker_panel <- c("ALB", "PTPRC", "PECAM1", "KRT19", "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "RGS5", "ACTA2", "TAGLN", "LRAT", "RBP1")
marker_panel <- marker_panel[marker_panel %in% rownames(obj)]
if (length(marker_panel) >= 4) {
  p_dot <- DotPlot(
    obj,
    features = marker_panel,
    group.by = "broad_compartment",
    cols = c("#D9D9D9", "#D55E00"),
    dot.scale = 5.2
  ) +
    coord_flip() +
    labs(
      title = "Canonical marker expression by UMAP marker context",
      subtitle = "Dot size shows detection fraction; color shows average scaled expression",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 9, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8.5, color = "#333333"),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(face = "bold")
    )
  ggsave(file.path(fig_dir, "gse212837_sketch_marker_dotplot.pdf"), p_dot, width = 7.6, height = 5.0, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, "gse212837_sketch_marker_dotplot.png"), p_dot, width = 7.6, height = 5.0, units = "in", dpi = 600, bg = "white")
} else {
  p_dot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Marker genes not available for dot plot") +
    theme_void(base_family = "Arial")
}

violin_data <- cell_scores %>%
  filter(!is.na(broad_compartment)) %>%
  mutate(
    broad_compartment = factor(broad_compartment),
    disease_group = factor(disease_group, levels = c("Control", "NASH"))
  )
p_violin <- ggplot(violin_data, aes(x = broad_compartment, y = stromal_marker_score, fill = broad_compartment)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.78) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "#222222", linewidth = 0.25) +
  facet_wrap(~ disease_group, nrow = 1) +
  scale_fill_manual(values = compartment_colors, na.value = "#BDBDBD") +
  labs(
    title = "Stromal-marker score across marker-context compartments",
    subtitle = "Nucleus-level distributions are descriptive; donor-level summaries are reported separately",
    x = NULL,
    y = "Stromal-marker score",
    fill = NULL
  ) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "#333333"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(fig_dir, "gse212837_sketch_stromal_score_violin.pdf"), p_violin, width = 7.0, height = 4.1, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse212837_sketch_stromal_score_violin.png"), p_violin, width = 7.0, height = 4.1, units = "in", dpi = 600, bg = "white")

composite <- (p_cluster + p_compartment + p_group) / (p_primary + p_ecm_excluded + p_hsc) / (p_dot + p_violin + plot_layout(widths = c(1.15, 1))) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 13, color = "#0B2545"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

ggsave(file.path(fig_dir, "gse212837_sketch_umap_composite.pdf"), composite, width = 11.2, height = 11.4, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse212837_sketch_umap_composite.png"), composite, width = 11.2, height = 11.4, units = "in", dpi = 600, bg = "white")

figure_manifest <- tribble(
  ~figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "gse212837_sketch_umap_seurat_cluster", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "cluster UMAP", "descriptive sketch embedding",
  "gse212837_sketch_umap_marker_context", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "marker-context UMAP", "descriptive sketch embedding",
  "gse212837_sketch_umap_disease_group", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "donor-group UMAP", "localization only, not cell-level group inference",
  "gse212837_sketch_umap_primary_score", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "primary-score feature map", "cell-level score map with donor-level inference retained separately",
  "gse212837_sketch_umap_ecm_excluded_score", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "ECM-excluded feature map", "non-structural score context",
  "gse212837_sketch_umap_stromal_marker_score", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "stromal-marker feature map", "marker-context visualization",
  "gse212837_sketch_umap_hsc_context_score", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "HSC-context feature map", "context visualization",
  "gse212837_sketch_marker_dotplot", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv; results/single_cell/gse212837_sketch_seurat_object.rds", "scripts/16_gse212837_single_nucleus_embedding.R", "canonical marker dot plot", "marker support for broad UMAP context labels",
  "gse212837_sketch_stromal_score_violin", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "stromal score violin", "descriptive nucleus-level distribution"
)
write_tsv(figure_manifest, file.path(result_dir, "gse212837_sketch_umap_figure_manifest.tsv"))

cat(
  "# GSE212837 Single-Nucleus Sketch UMAP\n\n",
  "Generated by `scripts/16_gse212837_single_nucleus_embedding.R`.\n\n",
  "The workflow reads public 10x matrices, applies the same nucleus-level QC thresholds used in the donor-aware stromal analysis, retains a balanced random sketch from each sample, and prioritizes high-confidence or top-tail stromal-marker nuclei so that rare stromal contexts are visible on the embedding.\n\n",
  "## Outputs\n\n",
  "- `results/single_cell/gse212837_sketch_umap_cell_scores.tsv`\n",
  "- `results/single_cell/gse212837_sketch_umap_summary.tsv`\n",
  "- `results/single_cell/gse212837_sketch_umap_figure_manifest.tsv`\n",
  "- `results/single_cell/gse212837_sketch_seurat_object.rds`\n",
  "- `figures/single_cell/gse212837_sketch_umap_composite.png` and `.pdf`\n",
  "- Individual UMAP panels in `figures/single_cell/`.\n\n",
  "## Sketch Summary\n\n",
  "- Sketch nuclei: ", nrow(cell_scores), ".\n",
  "- Donors represented: ", n_distinct(cell_scores$donor_id), ".\n",
  "- Samples or regions represented: ", n_distinct(cell_scores$gsm), ".\n",
  "- High-confidence stromal nuclei retained: ", sum(cell_scores$stromal_candidate), ".\n",
  "- Stromal-marker-enriched nuclei retained: ", sum(cell_scores$stromal_enriched_candidate), ".\n\n",
  "## Interpretation\n\n",
  "The UMAP panels are descriptive cell-state localization figures. Donor-level summaries remain the analysis layer for disease-group comparisons.\n",
  file = file.path(doc_dir, "gse212837_sketch_umap_summary.md"),
  sep = ""
)

message("Done. GSE212837 sketch UMAP outputs written to ", fig_dir)

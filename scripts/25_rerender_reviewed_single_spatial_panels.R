#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

single_result_dir <- file.path(project_dir, "results", "single_cell")
single_fig_dir <- file.path(project_dir, "figures", "single_cell")
spatial_result_dir <- file.path(project_dir, "results", "spatial")
spatial_fig_dir <- file.path(project_dir, "figures", "spatial")
dir.create(single_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_umap_reviewed <- theme_classic(base_size = 8.8, base_family = "") +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8.4, color = "#333333"),
    axis.title = element_text(face = "bold", size = 8.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7.4),
    legend.key.height = unit(0.15, "in"),
    plot.margin = margin(8, 9, 8, 8)
  )

theme_spatial_reviewed <- theme_classic(base_size = 8.8, base_family = "") +
  theme(
    plot.title = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 8.4, color = "#333333"),
    axis.title = element_text(face = "bold", size = 8.5),
    axis.text = element_text(color = "#111111", size = 7.8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7.4),
    panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "#F2F2F2", linewidth = 0.2),
    plot.margin = margin(8, 9, 8, 8)
  )

score_palette <- c("#253494", "#2C7FB8", "#41B6C4", "#FED976", "#F03B20", "#99000D")
neighbor_palette <- c("#5B5B5B", "#74ADD1", "#1A9850", "#FEE08B", "#D73027")

save_plot_pair <- function(plot, dir, stem, width, height) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::pdf)
  ggsave(file.path(dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

plot_umap_category <- function(data, color_col, title, subtitle, colors, out_stem) {
  p <- ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
    geom_point(size = 0.26, alpha = 0.88, stroke = 0) +
    scale_color_manual(values = colors, na.value = "#8A8A8A") +
    guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    theme_umap_reviewed
  save_plot_pair(p, single_fig_dir, out_stem, 4.35, 3.75)
  p
}

plot_umap_score <- function(data, score_col, title, subtitle, out_stem) {
  p <- ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = .data[[score_col]])) +
    geom_point(size = 0.27, alpha = 0.92, stroke = 0) +
    scale_color_gradientn(colors = score_palette, na.value = "#D9D9D9", name = "Score") +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2", color = "Score") +
    theme_umap_reviewed
  save_plot_pair(p, single_fig_dir, out_stem, 4.35, 3.75)
  p
}

cell_scores <- read_tsv(file.path(single_result_dir, "gse212837_sketch_umap_cell_scores.tsv"), show_col_types = FALSE) %>%
  mutate(
    disease_group = factor(disease_group, levels = c("Control", "NASH")),
    broad_compartment = factor(
      broad_compartment,
      levels = c("Hepatocyte-like", "Immune-like", "Endothelial-like", "Cholangiocyte-like", "High-confidence stromal", "Stromal-marker enriched", "Unassigned")
    )
  )

compartment_colors <- c(
  "Hepatocyte-like" = "#0072B2",
  "Immune-like" = "#D55E00",
  "Endothelial-like" = "#009E73",
  "Cholangiocyte-like" = "#CC79A7",
  "High-confidence stromal" = "#E69F00",
  "Stromal-marker enriched" = "#7B3294",
  "Unassigned" = "#777777"
)

plot_umap_category(
  cell_scores,
  "broad_compartment",
  "GSE212837 sketch UMAP by marker context",
  "Balanced nuclei sketch with stromal-marker-priority retention",
  compartment_colors,
  "gse212837_sketch_umap_marker_context"
)

plot_umap_category(
  cell_scores,
  "disease_group",
  "GSE212837 sketch UMAP by donor group",
  "Control and NASH nuclei shown for localization only",
  c(Control = "#0072B2", NASH = "#D55E00"),
  "gse212837_sketch_umap_disease_group"
)

cluster_levels <- sort(unique(as.character(cell_scores$seurat_clusters)))
plot_umap_category(
  cell_scores %>% mutate(seurat_clusters = factor(seurat_clusters, levels = cluster_levels)),
  "seurat_clusters",
  "GSE212837 sketch UMAP by Seurat cluster",
  "Unsupervised clusters from the balanced sketch object",
  setNames(scales::hue_pal()(length(cluster_levels)), cluster_levels),
  "gse212837_sketch_umap_seurat_cluster"
)

plot_umap_score(
  cell_scores,
  "primary_score",
  "Primary program score on sketch UMAP",
  "Cell-level score map; donor summaries remain the inferential layer",
  "gse212837_sketch_umap_primary_score"
)

plot_umap_score(
  cell_scores,
  "ecm_excluded_score",
  "ECM-excluded score on sketch UMAP",
  "Non-structural locked stromal program context",
  "gse212837_sketch_umap_ecm_excluded_score"
)

plot_umap_score(
  cell_scores,
  "stromal_marker_score",
  "Stromal-marker score on sketch UMAP",
  "Marker-tail retention supports stromal-context visualization",
  "gse212837_sketch_umap_stromal_marker_score"
)

rds_path <- file.path(single_result_dir, "gse212837_sketch_seurat_object.rds")
marker_panel <- c("ALB", "PTPRC", "PECAM1", "KRT19", "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "RGS5", "ACTA2", "TAGLN", "LRAT", "RBP1")
if (file.exists(rds_path) && requireNamespace("Seurat", quietly = TRUE)) {
  suppressPackageStartupMessages(library(Seurat))
  obj <- readRDS(rds_path)
  obj$broad_compartment <- factor(
    obj$broad_compartment,
    levels = c("Hepatocyte-like", "Immune-like", "Endothelial-like", "Cholangiocyte-like", "High-confidence stromal", "Stromal-marker enriched", "Unassigned")
  )
  marker_panel <- marker_panel[marker_panel %in% rownames(obj)]
  p_dot <- DotPlot(
    obj,
    features = marker_panel,
    group.by = "broad_compartment",
    cols = c("#E6E6E6", "#B40426"),
    dot.scale = 5.4
  ) +
    coord_flip() +
    labs(
      title = "Canonical marker expression by UMAP marker context",
      subtitle = "Dot size shows detection fraction; color shows average scaled expression",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 9.2, base_family = "") +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8.5, color = "#333333"),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8.2),
      axis.text.y = element_text(size = 8.2),
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 7.5)
    )
  save_plot_pair(p_dot, single_fig_dir, "gse212837_sketch_marker_dotplot", 7.8, 5.1)
}

spatial_scores <- read_tsv(file.path(spatial_result_dir, "gse292268_spatial_cell_scores.tsv"), show_col_types = FALSE)
coverage <- read_tsv(file.path(spatial_result_dir, "gse292268_panel_coverage.tsv"), show_col_types = FALSE)
niche_summary <- read_tsv(file.path(spatial_result_dir, "gse292268_spatial_niche_summary.tsv"), show_col_types = FALSE)
bin_summary <- read_tsv(file.path(spatial_result_dir, "gse292268_spatial_stellate_neighbor_bin_summary.tsv"), show_col_types = FALSE)

coverage_plot_data <- coverage %>%
  mutate(score_name = recode(
    score_name,
    primary = "Primary locked score",
    ecm_excluded = "ECM-excluded score",
    structural_ecm = "Structural ECM-only",
    hsc_context = "HSC context",
    stromal_marker = "Stromal markers"
  )) %>%
  mutate(score_name = factor(score_name, levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only", "HSC context", "Stromal markers")))

p_cov <- ggplot(coverage_plot_data, aes(x = score_name, y = coverage_fraction)) +
  geom_col(width = 0.62, fill = "#2A9D8F", color = "#1A1A1A", linewidth = 0.22) +
  geom_text(aes(label = paste0(n_present_genes, "/", n_total_genes)), vjust = -0.35, size = 2.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "GSE292268 CosMx panel coverage",
    subtitle = "Exact feature overlap with the locked stromal program",
    x = NULL,
    y = "Coverage"
  ) +
  theme_spatial_reviewed +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot_pair(p_cov, spatial_fig_dir, "gse292268_panel_coverage", 5.8, 3.55)

map_data <- spatial_scores %>%
  filter(qc_pass, is.finite(CenterX_global_px), is.finite(CenterY_global_px)) %>%
  group_by(gsm, slide_id) %>%
  mutate(
    primary_score_z_clip = pmax(pmin(primary_score_z, 2.5), -2.5),
    x_norm = rescale(CenterX_global_px),
    y_norm = rescale(CenterY_global_px),
    slide_label = paste0(gsm, " ", slide_id)
  ) %>%
  ungroup()

p_map <- ggplot(map_data, aes(x = x_norm, y = y_norm, color = primary_score_z_clip)) +
  geom_point(size = 0.105, alpha = 0.96, stroke = 0) +
  facet_wrap(~ slide_label, nrow = 1) +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_gradientn(colors = score_palette, limits = c(-2.5, 2.5), oob = squish, name = "Score z") +
  labs(
    title = "Spatial localization of the locked stromal program in GSE292268",
    subtitle = "Cell-level CosMx map, z-scored within slide; all QC-passing cells are shown",
    x = "Normalized slide X coordinate",
    y = "Normalized slide Y coordinate"
  ) +
  theme_spatial_reviewed +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "#F1F1F1", color = NA),
    panel.border = element_rect(fill = NA, color = "#9A9A9A", linewidth = 0.3)
  )
save_plot_pair(p_map, spatial_fig_dir, "gse292268_program_localization_map", 8.4, 3.95)

neighbor_map_data <- map_data %>%
  mutate(
    stellate_neighbors_clip = pmin(stellate_neighbors, quantile(stellate_neighbors, 0.995, na.rm = TRUE)),
    stellate_neighbors_clip = if_else(is.finite(stellate_neighbors_clip), stellate_neighbors_clip, NA_real_)
  )

p_neighbor_map <- ggplot(neighbor_map_data, aes(x = x_norm, y = y_norm, color = stellate_neighbors_clip)) +
  geom_point(size = 0.105, alpha = 0.96, stroke = 0) +
  facet_wrap(~ slide_label, nrow = 1) +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_gradientn(colors = neighbor_palette, name = "Stellate\nneighbors", na.value = "#D9D9D9") +
  labs(
    title = "Spatial distribution of Stellate-cell-neighbor abundance",
    subtitle = "Cell-level neighbor metadata projected on matched slide coordinates",
    x = "Normalized slide X coordinate",
    y = "Normalized slide Y coordinate"
  ) +
  theme_spatial_reviewed +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "#F1F1F1", color = NA),
    panel.border = element_rect(fill = NA, color = "#9A9A9A", linewidth = 0.3)
  )
save_plot_pair(p_neighbor_map, spatial_fig_dir, "gse292268_stellate_neighbor_spatial_map", 8.4, 3.95)

bubble_data <- niche_summary %>%
  mutate(slide_label = paste0(gsm, " ", slide_id))

p_bubble <- ggplot(bubble_data, aes(x = mean_stellate_neighbors, y = mean_primary_score_z)) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.25, linetype = "dashed") +
  geom_point(aes(size = n_cells, color = mean_hsc_context_score_z), alpha = 0.9) +
  geom_text(aes(label = niche), size = 2.25, vjust = -0.75, check_overlap = TRUE) +
  facet_wrap(~ slide_label, nrow = 1, scales = "free_x") +
  scale_size_continuous(range = c(2.4, 6.4), labels = comma) +
  scale_color_gradientn(colors = score_palette, name = "HSC context z") +
  labs(
    title = "Niche-level localization relative to stellate-cell-neighbor-rich regions",
    subtitle = "Each point is one spatial niche; descriptive localization, not donor-level inference",
    x = "Mean Stellate-cell neighbors",
    y = "Mean primary score z",
    size = "Cells"
  ) +
  theme_spatial_reviewed
save_plot_pair(p_bubble, spatial_fig_dir, "gse292268_niche_stellate_localization_bubble", 7.6, 3.95)

bin_plot_data <- bin_summary %>%
  select(gsm, slide_id, stellate_neighbor_bin, n_cells, starts_with("median_")) %>%
  pivot_longer(starts_with("median_"), names_to = "score_name", values_to = "median_score_z") %>%
  mutate(
    score_name = recode(
      score_name,
      median_primary_score_z = "Primary locked score",
      median_ecm_excluded_score_z = "ECM-excluded score",
      median_structural_ecm_score_z = "Structural ECM-only",
      median_hsc_context_score_z = "HSC context"
    ),
    score_name = factor(score_name, levels = c("Primary locked score", "ECM-excluded score", "Structural ECM-only", "HSC context")),
    slide_label = paste0(gsm, " ", slide_id)
  )

p_bin <- ggplot(bin_plot_data, aes(x = stellate_neighbor_bin, y = median_score_z, group = score_name, color = score_name)) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.25, linetype = "dashed") +
  geom_line(linewidth = 0.46) +
  geom_point(size = 2.0) +
  facet_wrap(~ slide_label, nrow = 1) +
  scale_color_manual(values = c(
    "Primary locked score" = "#0072B2",
    "ECM-excluded score" = "#2A9D8F",
    "Structural ECM-only" = "#D55E00",
    "HSC context" = "#7B3294"
  )) +
  labs(
    title = "Score gradient across Stellate-neighbor quartiles",
    subtitle = "Within-slide descriptive bins; cells are not treated as independent donors",
    x = "Stellate-neighbor bin",
    y = "Median score z",
    color = "Score"
  ) +
  theme_spatial_reviewed +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot_pair(p_bin, spatial_fig_dir, "gse292268_stellate_neighbor_score_gradient", 8.2, 3.95)

rerender_manifest <- tibble(
  panel_group = c("single_nucleus", "spatial"),
  source_inputs = c(
    "results/single_cell/gse212837_sketch_umap_cell_scores.tsv; results/single_cell/gse212837_sketch_seurat_object.rds",
    "results/spatial/gse292268_spatial_cell_scores.tsv; results/spatial/gse292268_panel_coverage.tsv; results/spatial/gse292268_spatial_niche_summary.tsv; results/spatial/gse292268_spatial_stellate_neighbor_bin_summary.tsv"
  ),
  outputs = c(
    "figures/single_cell/gse212837_sketch_umap_*.png/pdf; figures/single_cell/gse212837_sketch_marker_dotplot.png/pdf",
    "figures/spatial/gse292268_*.png/pdf"
  ),
  revision_note = c(
    "Improved UMAP score-map contrast and replaced stronger mechanostress wording with locked stromal-program wording.",
    "Improved spatial point visibility and replaced stronger mechanical-stress wording in coverage/localization labels."
  )
)
write_tsv(rerender_manifest, file.path(project_dir, "docs", "figure_team_reviewed_single_spatial_rerender_manifest.tsv"))

message("Reviewed single-nucleus and spatial source panels rerendered.")

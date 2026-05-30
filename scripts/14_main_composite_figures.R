#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(magick)
  library(tibble)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

out_dir <- file.path(project_dir, "figures", "main_text_composites")
doc_dir <- file.path(project_dir, "docs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)

panel_label <- function(img, label) {
  image_annotate(
    img,
    text = label,
    size = 76,
    weight = 700,
    color = "#0B2545",
    gravity = "northwest",
    location = "+38+26",
    font = "Arial"
  )
}

read_panel <- function(path, label, width = 2500, height = 1700, border = FALSE) {
  full_path <- file.path(project_dir, path)
  if (!file.exists(full_path)) stop("Missing panel: ", full_path)
  img <- image_read(full_path)
  img <- image_trim(img)
  img <- image_resize(img, paste0(width - 110, "x", height - 110, ">"))
  canvas <- image_blank(width, height, color = "white")
  canvas <- image_composite(canvas, img, gravity = "center")
  if (border) {
    canvas <- image_border(canvas, color = "#D8DEE6", geometry = "2x2")
  }
  panel_label(canvas, label)
}

row_join <- function(images) {
  image_append(image_join(images), stack = FALSE)
}

col_join <- function(images) {
  image_append(image_join(images), stack = TRUE)
}

extend_to_width <- function(img, target_width) {
  info <- image_info(img)
  if (info$width >= target_width) return(img)
  image_extent(img, geometry = paste0(target_width, "x", info$height), gravity = "center", color = "white")
}

col_join_center <- function(images) {
  widths <- vapply(images, function(x) image_info(x)$width, numeric(1))
  target_width <- max(widths)
  col_join(lapply(images, extend_to_width, target_width = target_width))
}

save_composite <- function(img, stem, caption, source_rows) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  image_write(img, png_path, format = "png", density = "600x600")
  image_write(img, pdf_path, format = "pdf")
  tibble(
    figure = stem,
    output_png = png_path,
    output_pdf = pdf_path,
    caption = caption
  ) %>%
    write_tsv(file.path(out_dir, paste0(stem, "_caption.tsv")))
  source_rows %>%
    mutate(figure = stem, output_png = png_path, output_pdf = pdf_path) %>%
    select(figure, panel, source_figure, source_data, generating_script, role, wording_note, output_png, output_pdf) %>%
    write_tsv(file.path(out_dir, paste0(stem, "_source_manifest.tsv")))
  invisible(list(png = png_path, pdf = pdf_path))
}

bulk_panels <- tribble(
  ~panel, ~source_figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "A", "figures/phase2_bulk/phase2_primary_score_by_fibrosis_stage.png", "results/bulk/figure_source_data/phase2_primary_score_by_fibrosis_stage_source.tsv", "scripts/06_phase2_bulk_visuals.R", "fibrosis-stage distribution", "within-cohort score distribution by fibrosis stage",
  "B", "figures/phase2_bulk/phase2_score_by_advanced_fibrosis.png", "results/bulk/figure_source_data/phase2_score_by_advanced_fibrosis_source.tsv", "scripts/06_phase2_bulk_visuals.R", "advanced-fibrosis score distributions", "primary and sensitivity score distributions",
  "C", "figures/phase2_bulk/phase2_advanced_fibrosis_forest.png", "results/bulk/figure_source_data/phase2_advanced_fibrosis_forest_source.tsv", "scripts/06_phase2_bulk_visuals.R", "advanced-fibrosis logistic forest plot", "OR per 1-SD score increase",
  "D", "figures/phase2_bulk/phase2_endpoint_spearman_dotplot.png", "results/bulk/figure_source_data/phase2_endpoint_spearman_dotplot_source.tsv", "scripts/06_phase2_bulk_visuals.R", "fibrosis and NAS Spearman dot plot", "ordinal endpoint correlations",
  "E", "figures/phase2_bulk/phase2_random_gene_set_controls.png", "results/bulk/figure_source_data/phase2_random_gene_set_controls_source.tsv", "scripts/06_phase2_bulk_visuals.R", "expression-matched random gene-set benchmark", "negative-control benchmark",
  "F", "figures/phase2_bulk/phase2_signature_gene_coverage.png", "results/bulk/figure_source_data/phase2_signature_gene_coverage_source.tsv", "scripts/06_phase2_bulk_visuals.R", "locked-gene coverage", "technical feature coverage"
)
bulk_img <- col_join(list(
  row_join(list(read_panel(bulk_panels$source_figure[1], "A"), read_panel(bulk_panels$source_figure[2], "B"))),
  row_join(list(read_panel(bulk_panels$source_figure[3], "C"), read_panel(bulk_panels$source_figure[4], "D"))),
  row_join(list(read_panel(bulk_panels$source_figure[5], "E"), read_panel(bulk_panels$source_figure[6], "F")))
))
save_composite(
  bulk_img,
  "Figure_1_bulk_standard_association",
  "Figure 1. Standard bulk transcriptomic visualization of the locked mechanical-stress-associated program across fibrosis-mapped MASLD/MASH cohorts.",
  bulk_panels
)

single_umap_panels <- tribble(
  ~panel, ~source_figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "A", "figures/single_cell/gse212837_sketch_umap_seurat_cluster.png", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "Seurat cluster UMAP", "balanced single-nucleus sketch",
  "B", "figures/single_cell/gse212837_sketch_umap_marker_context.png", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "marker-context UMAP", "broad marker context labels",
  "C", "figures/single_cell/gse212837_sketch_umap_disease_group.png", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "disease-group UMAP", "visualization only; donor-level summaries used for group comparisons",
  "D", "figures/single_cell/gse212837_sketch_umap_primary_score.png", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "primary-score feature map", "program-score localization on UMAP",
  "E", "figures/single_cell/gse212837_sketch_umap_ecm_excluded_score.png", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "ECM-excluded feature map", "non-structural score context",
  "F", "figures/single_cell/gse212837_sketch_marker_dotplot.png", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv; results/single_cell/gse212837_sketch_seurat_object.rds", "scripts/16_gse212837_single_nucleus_embedding.R", "canonical marker dot plot", "marker support for broad context labels"
)
single_umap_img <- col_join_center(list(
  row_join(list(
    read_panel(single_umap_panels$source_figure[1], "A", width = 1800, height = 1550),
    read_panel(single_umap_panels$source_figure[2], "B", width = 1800, height = 1550),
    read_panel(single_umap_panels$source_figure[3], "C", width = 1800, height = 1550)
  )),
  row_join(list(
    read_panel(single_umap_panels$source_figure[4], "D", width = 1700, height = 1650),
    read_panel(single_umap_panels$source_figure[5], "E", width = 1700, height = 1650),
    read_panel(single_umap_panels$source_figure[6], "F", width = 2500, height = 1650)
  ))
))
save_composite(
  single_umap_img,
  "Figure_2_single_nucleus_standard_umap",
  "Figure 2. Standard single-nucleus sketch UMAP and marker visualization in GSE212837.",
  single_umap_panels
)

single_donor_panels <- tribble(
  ~panel, ~source_figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "A", "figures/single_cell/gse212837_donor_stromal_scores.png", "results/single_cell/gse212837_human_stromal_donor_scores.tsv", "scripts/08_gse212837_hsc_stromal_validation.R", "high-confidence stromal donor scores", "marker-gated stromal context",
  "B", "figures/single_cell/gse212837_stromal_fraction_qc.png", "results/single_cell/gse212837_human_sample_qc.tsv", "scripts/08_gse212837_hsc_stromal_validation.R", "high-confidence stromal fraction QC", "gate behavior by donor-region",
  "C", "figures/single_cell/gse212837_donor_stromal_enriched_sensitivity_scores.png", "results/single_cell/gse212837_human_stromal_enriched_donor_scores.tsv", "scripts/08_gse212837_hsc_stromal_validation.R", "stromal-marker-enriched donor scores", "donor-level sensitivity summary",
  "D", "figures/single_cell/gse212837_stromal_enriched_fraction_qc.png", "results/single_cell/gse212837_human_sample_qc.tsv", "scripts/08_gse212837_hsc_stromal_validation.R", "stromal-marker-enriched fraction QC", "top-tail marker enrichment QC",
  "E", "figures/single_cell/gse212837_sketch_stromal_score_violin.png", "results/single_cell/gse212837_sketch_umap_cell_scores.tsv", "scripts/16_gse212837_single_nucleus_embedding.R", "stromal-marker score violin", "nucleus-level descriptive distribution"
)
single_donor_img <- col_join(list(
  row_join(list(
    read_panel(single_donor_panels$source_figure[1], "A", width = 3400, height = 1750),
    read_panel(single_donor_panels$source_figure[2], "B", width = 1800, height = 1750)
  )),
  row_join(list(
    read_panel(single_donor_panels$source_figure[3], "C", width = 3400, height = 1750),
    read_panel(single_donor_panels$source_figure[4], "D", width = 1800, height = 1750)
  )),
  row_join(list(read_panel(single_donor_panels$source_figure[5], "E", width = 5200, height = 1550)))
))
save_composite(
  single_donor_img,
  "Figure_3_single_nucleus_donor_aware_stromal",
  "Figure 3. Donor-aware stromal summaries and stromal-marker score distributions in GSE212837.",
  single_donor_panels
)

spatial_panels <- tribble(
  ~panel, ~source_figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "A", "figures/spatial/gse292268_program_localization_map.png", "results/spatial/gse292268_spatial_cell_scores.tsv", "scripts/09_gse292268_spatial_localization.R", "cell-level program feature map", "spatial localization on true coordinates",
  "B", "figures/spatial/gse292268_stellate_neighbor_spatial_map.png", "results/spatial/gse292268_spatial_cell_scores.tsv", "scripts/09_gse292268_spatial_localization.R", "Stellate-neighbor spatial map", "neighbor metadata on true coordinates",
  "C", "figures/spatial/gse292268_niche_stellate_localization_bubble.png", "results/spatial/gse292268_spatial_niche_summary.tsv", "scripts/09_gse292268_spatial_localization.R", "niche-level localization bubble plot", "Stellate-neighbor context",
  "D", "figures/spatial/gse292268_stellate_neighbor_score_gradient.png", "results/spatial/gse292268_spatial_stellate_neighbor_bin_summary.tsv", "scripts/09_gse292268_spatial_localization.R", "Stellate-neighbor quartile gradient", "within-slide descriptive bins",
  "E", "figures/spatial/gse292268_panel_coverage.png", "results/spatial/gse292268_panel_coverage.tsv", "scripts/09_gse292268_spatial_localization.R", "CosMx panel coverage", "feature availability"
)
spatial_img <- col_join(list(
  row_join(list(
    read_panel(spatial_panels$source_figure[1], "A", width = 3000, height = 1900),
    read_panel(spatial_panels$source_figure[2], "B", width = 3000, height = 1900)
  )),
  row_join(list(
    read_panel(spatial_panels$source_figure[3], "C", width = 3000, height = 1750),
    read_panel(spatial_panels$source_figure[4], "D", width = 3000, height = 1750)
  )),
  row_join(list(read_panel(spatial_panels$source_figure[5], "E", width = 6000, height = 1250)))
))
save_composite(
  spatial_img,
  "Figure_4_spatial_standard_localization",
  "Figure 4. Standard spatial transcriptomic localization of program activity and Stellate-cell-neighbor context in GSE292268.",
  spatial_panels
)

network_panels <- tribble(
  ~panel, ~source_figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "A", "figures/network/network_prioritization_string_network.png", "results/network/network_prioritization_string_filtered_edges.tsv", "scripts/10_network_prioritization.R", "STRING network", "high-confidence functional association context",
  "B", "figures/network/network_prioritization_lollipop.png", "results/network/network_prioritization_centrality_table.tsv", "scripts/10_network_prioritization.R", "candidate-node lollipop ranking", "centrality-based prioritization",
  "C", "figures/network/network_prioritization_algorithm_upset.png", "results/network/network_prioritization_algorithm_overlap.tsv", "scripts/10_network_prioritization.R", "centrality overlap", "pre-specified top-quartile overlap"
)
network_img <- col_join(list(
  row_join(list(
    read_panel(network_panels$source_figure[1], "A", width = 3400, height = 2350),
    read_panel(network_panels$source_figure[2], "B", width = 2600, height = 2350)
  )),
  row_join(list(
    read_panel(network_panels$source_figure[3], "C", width = 6000, height = 1700)
  ))
))
save_composite(
  network_img,
  "Figure_5_network_standard_prioritization",
  "Figure 5. Standard network visualization and centrality-based candidate prioritization within the locked program.",
  network_panels
)

coexp_panels <- tribble(
  ~panel, ~source_figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "A", "figures/network/bulk_locked_gene_fibrosis_consistency_heatmap.png", "results/network/bulk_locked_gene_consistency_summary.tsv", "scripts/11_bulk_coexpression_network_context.R", "gene-level fibrosis association heatmap", "cross-cohort gene-level context",
  "B", "figures/network/bulk_locked_gene_coexpression_heatmap.png", "results/network/bulk_locked_gene_coexpression_edge_summary.tsv", "scripts/11_bulk_coexpression_network_context.R", "coexpression heatmap", "sample-level gene-gene coordination",
  "C", "figures/network/bulk_locked_gene_coexpression_network.png", "results/network/bulk_locked_gene_coexpression_edge_summary.tsv; results/network/network_prioritization_gene_evidence_matrix.tsv", "scripts/11_bulk_coexpression_network_context.R", "reproducible coexpression network", "positive replicated coexpression context",
  "D", "figures/network/network_prioritization_integrated_candidate_evidence_matrix.png", "results/network/network_prioritization_gene_evidence_matrix.tsv", "scripts/11_bulk_coexpression_network_context.R", "integrated evidence matrix", "candidate prioritization summary"
)
coexp_img <- col_join(list(
  row_join(list(
    read_panel(coexp_panels$source_figure[1], "A", width = 3100, height = 2350),
    read_panel(coexp_panels$source_figure[2], "B", width = 3100, height = 2350)
  )),
  row_join(list(
    read_panel(coexp_panels$source_figure[3], "C", width = 3100, height = 2200),
    read_panel(coexp_panels$source_figure[4], "D", width = 3100, height = 2200)
  ))
))
save_composite(
  coexp_img,
  "Figure_6_gene_level_coexpression_standard_context",
  "Figure 6. Standard gene-level fibrosis association and reproducible coexpression context across fibrosis-mapped cohorts.",
  coexp_panels
)

figure7_panels <- tribble(
  ~panel, ~source_figure, ~source_data, ~generating_script, ~role, ~wording_note,
  "A", "figures/benchmark_robustness_sensitivity/Figure_7_benchmark_robustness_sensitivity.png", "results/benchmark_robustness_sensitivity/fibrosis_stromal_benchmark_results.tsv", "scripts/22_benchmark_robustness_sensitivity_analyses.R", "fibrosis and proxy benchmark associations", "benchmark comparator analysis, not superiority testing",
  "B", "figures/benchmark_robustness_sensitivity/Figure_7_benchmark_robustness_sensitivity.png", "results/benchmark_robustness_sensitivity/bulk_confounding_sensitivity_models.tsv", "scripts/22_benchmark_robustness_sensitivity_analyses.R", "proxy-adjusted sensitivity models", "bulk proxy adjustment with collinearity limitations",
  "C", "figures/benchmark_robustness_sensitivity/Figure_7_benchmark_robustness_sensitivity.png", "results/benchmark_robustness_sensitivity/small_sample_bootstrap_iterations.tsv; results/benchmark_robustness_sensitivity/small_sample_bootstrap_stability_summary.tsv", "scripts/22_benchmark_robustness_sensitivity_analyses.R", "small-cohort stratified bootstrap", "direction stability check, not definitive external validation",
  "D", "figures/benchmark_robustness_sensitivity/Figure_7_benchmark_robustness_sensitivity.png", "results/benchmark_robustness_sensitivity/small_sample_leave_one_out_models.tsv; results/benchmark_robustness_sensitivity/small_sample_leave_one_out_summary.tsv", "scripts/22_benchmark_robustness_sensitivity_analyses.R", "leave-one-out small-cohort stability", "single-sample influence check, not definitive external validation",
  "E", "figures/benchmark_robustness_sensitivity/Figure_7_benchmark_robustness_sensitivity.png", "results/benchmark_robustness_sensitivity/bulk_score_proxy_correlation_matrix.tsv", "scripts/22_benchmark_robustness_sensitivity_analyses.R", "score and proxy correlation context", "bulk score-proxy overlap context"
)
if (file.exists(file.path(project_dir, figure7_panels$source_figure[[1]]))) {
  file.copy(
    file.path(project_dir, "figures", "benchmark_robustness_sensitivity", "Figure_7_benchmark_robustness_sensitivity.png"),
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity.png"),
    overwrite = TRUE
  )
  file.copy(
    file.path(project_dir, "figures", "benchmark_robustness_sensitivity", "Figure_7_benchmark_robustness_sensitivity.pdf"),
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity.pdf"),
    overwrite = TRUE
  )
  file.copy(
    file.path(project_dir, "figures", "benchmark_robustness_sensitivity", "Figure_7_benchmark_robustness_sensitivity_caption.tsv"),
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity_caption.tsv"),
    overwrite = TRUE
  )
  write_tsv(figure7_panels, file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity_source_manifest.tsv"))
  file.copy(
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity.png"),
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity_benchmark_sensitivity.png"),
    overwrite = TRUE
  )
  file.copy(
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity.pdf"),
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity_benchmark_sensitivity.pdf"),
    overwrite = TRUE
  )
  file.copy(
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity_caption.tsv"),
    file.path(out_dir, "Figure_7_benchmark_robustness_sensitivity_caption.tsv"),
    overwrite = TRUE
  )
}

combined_manifest <- bind_rows(
  bulk_panels %>% mutate(figure = "Figure_1_bulk_standard_association"),
  single_umap_panels %>% mutate(figure = "Figure_2_single_nucleus_standard_umap"),
  single_donor_panels %>% mutate(figure = "Figure_3_single_nucleus_donor_aware_stromal"),
  spatial_panels %>% mutate(figure = "Figure_4_spatial_standard_localization"),
  network_panels %>% mutate(figure = "Figure_5_network_standard_prioritization"),
  coexp_panels %>% mutate(figure = "Figure_6_gene_level_coexpression_standard_context"),
  figure7_panels %>% mutate(figure = "Figure_7_benchmark_robustness_sensitivity")
) %>%
  select(figure, panel, source_figure, source_data, generating_script, role, wording_note)

write_tsv(combined_manifest, file.path(out_dir, "main_text_composite_source_manifest.tsv"))

manifest_md <- c(
  "# Main Text Composite Figure Manifest",
  "",
  "Generated by reproducible R composite assembly from R-generated publication-style source plots.",
  "",
  "The composite step arranges R-generated panels into manuscript-ready multi-panel figures. It does not recompute statistics and does not create new analytical results.",
  "",
  "The composite PNG files are manuscript-facing 600 ppi raster exports. The companion composite PDFs are raster-based PDF wrappers unless source-level vector recomposition is performed; they should not be described as fully vector figures.",
  "",
  "## Composite Figures",
  "",
  "- `Figure_1_bulk_standard_association`: standard bulk score distributions, forest plot, endpoint dot plot, random controls, and coverage.",
  "- `Figure_2_single_nucleus_standard_umap`: Seurat sketch UMAP, marker-context UMAP, score feature maps, and marker dot plot.",
  "- `Figure_3_single_nucleus_donor_aware_stromal`: donor-aware stromal summaries and stromal-marker score distributions.",
  "- `Figure_4_spatial_standard_localization`: CosMx coordinate maps, Stellate-neighbor context, niche plot, gradient plot, and panel coverage.",
  "- `Figure_5_network_standard_prioritization`: STRING network and centrality prioritization.",
  "- `Figure_6_gene_level_coexpression_standard_context`: gene-level fibrosis association and coexpression context.",
  "- `Figure_7_benchmark_robustness_sensitivity`: benchmark, proxy-adjustment, small-cohort stability, and score-proxy correlation analyses.",
  "",
  "Source-data traceability is stored in `figures/main_text_composites/main_text_composite_source_manifest.tsv`.",
  "",
  "Python scripts are document-packaging utilities only and are not part of the statistical or visualization analysis workflow."
)
writeLines(manifest_md, file.path(doc_dir, "main_text_composite_figure_manifest.md"))

message("Main-text composite figures written to: ", out_dir)

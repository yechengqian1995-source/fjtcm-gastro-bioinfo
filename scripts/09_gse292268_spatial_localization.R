#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(curl)
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

raw_dir <- file.path(project_dir, "data", "spatial_raw", "GSE292268")
result_dir <- file.path(project_dir, "results", "spatial")
fig_dir <- file.path(project_dir, "figures", "spatial")
doc_dir <- file.path(project_dir, "docs")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

links_path <- file.path(project_dir, "results", "single_spatial", "single_spatial_supplementary_links.tsv")
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
  stromal_marker = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "COL6A3", "PDGFRB", "RGS5", "ACTA2", "TAGLN", "DES", "LRAT", "RBP1")
)

download_if_missing <- function(url, dest) {
  url <- str_replace(url, "^ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading ", basename(dest))
    curl_download(url, dest, mode = "wb", quiet = FALSE)
  }
  dest
}

expand_supplementary_links <- function(links_path) {
  read_tsv(links_path, show_col_types = FALSE) %>%
    filter(accession == "GSE292268", organism_ch1 == "Homo sapiens") %>%
    mutate(url_split = str_split(url, "\\s+\\|\\|\\s+")) %>%
    select(accession, gsm, title, organism_ch1, platform_id, disease_label_raw, fibrosis_stage_raw, tissue_raw, url_split) %>%
    unnest(url_split) %>%
    transmute(
      accession,
      gsm,
      title,
      organism_ch1,
      platform_id,
      disease_label_raw,
      fibrosis_stage_raw,
      tissue_raw,
      url = str_trim(url_split),
      url_https = str_replace(str_trim(url_split), "^ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov"),
      file_name = basename(url_https),
      slide_id = str_match(file_name, "Slide[0-9]+")[, 1],
      file_kind = case_when(
        str_detect(file_name, "exprMat_file") ~ "exprMat",
        str_detect(file_name, "metadata_file") ~ "metadata",
        str_detect(file_name, "fov_positions_file") ~ "fov_positions",
        str_detect(file_name, "SeuratObject") ~ "seurat_object",
        str_detect(file_name, "tx_file") ~ "transcript_file",
        str_detect(file_name, "polygons") ~ "polygons",
        TRUE ~ "other"
      ),
      local_path = file.path(raw_dir, file_name)
    )
}

score_gene_set <- function(count_mat, gene_cols, genes, total_count) {
  present <- intersect(toupper(genes), gene_cols)
  if (!length(present)) return(rep(NA_real_, nrow(count_mat)))
  idx <- match(present, gene_cols)
  totals <- total_count
  totals[totals <= 0] <- NA_real_
  norm <- sweep(count_mat[, idx, drop = FALSE], 1, totals, "/") * 1e4
  rowMeans(log1p(norm), na.rm = TRUE)
}

score_z <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

theme_publication <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "#333333"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#111111"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "#F2F2F2", linewidth = 0.2)
  )

process_slide <- function(slide_files) {
  gsm <- unique(slide_files$gsm)
  slide_id <- unique(slide_files$slide_id)
  title <- unique(slide_files$title)
  if (length(gsm) != 1 || length(slide_id) != 1) stop("Expected one GSM and one slide per slide group.")

  expr_path <- slide_files$local_path[slide_files$file_kind == "exprMat"][[1]]
  meta_path <- slide_files$local_path[slide_files$file_kind == "metadata"][[1]]

  message("Processing ", gsm, " ", slide_id)
  expr <- read_csv(expr_path, show_col_types = FALSE, progress = FALSE)
  meta <- read_csv(meta_path, show_col_types = FALSE, progress = FALSE)

  key_cols <- c("fov", "cell_ID")
  if (!all(key_cols %in% names(expr))) stop("Expression matrix missing fov/cell_ID columns.")
  if (!all(key_cols %in% names(meta))) stop("Metadata missing fov/cell_ID columns.")

  gene_cols <- setdiff(names(expr), key_cols)
  gene_cols <- toupper(gene_cols)
  names(expr) <- c(key_cols, gene_cols)
  count_mat <- as.matrix(expr[, gene_cols, drop = FALSE])
  storage.mode(count_mat) <- "numeric"

  total_count <- rowSums(count_mat, na.rm = TRUE)
  n_feature <- rowSums(count_mat > 0, na.rm = TRUE)

  score_df <- expr[, key_cols] %>%
    mutate(
      cell_key = paste(fov, cell_ID, sep = "_"),
      total_count = total_count,
      n_feature = n_feature,
      primary_score = score_gene_set(count_mat, gene_cols, gene_sets$primary, total_count),
      ecm_excluded_score = score_gene_set(count_mat, gene_cols, gene_sets$ecm_excluded, total_count),
      structural_ecm_score = score_gene_set(count_mat, gene_cols, gene_sets$structural_ecm, total_count),
      hsc_context_score = score_gene_set(count_mat, gene_cols, gene_sets$hsc_context, total_count),
      stromal_marker_score = score_gene_set(count_mat, gene_cols, gene_sets$stromal_marker, total_count)
    )

  niche_col <- names(meta)[grepl("^spatialclust_.*assignments$", names(meta))][1]
  stellate_neighbor_col <- names(meta)[grepl("neighbours_Stellate[.]cells", names(meta))][1]
  portal_endo_col <- names(meta)[grepl("neighbours_Portal[.]endothelial[.]cells", names(meta))][1]
  periportal_lsec_col <- names(meta)[grepl("neighbours_Periportal[.]LSECs", names(meta))][1]
  inflammatory_mac_col <- names(meta)[grepl("neighbours_Inflammatory[.]macrophages", names(meta))][1]
  noninflammatory_mac_col <- names(meta)[grepl("neighbours_Non[.]inflammatory[.]macrophages", names(meta))][1]
  hep_cols <- names(meta)[grepl("neighbours_Hep[.][0-9]+$", names(meta))]

  required_meta <- c("fov", "cell_ID", "cell", "cell_id", "CenterX_global_px", "CenterY_global_px",
                     "qcFlagsCellCounts", "qcFlagsCellPropNeg", "qcFlagsCellComplex", "qcFlagsCellArea", "qcFlagsFOV")
  optional_meta <- c(niche_col, stellate_neighbor_col, portal_endo_col, periportal_lsec_col,
                     inflammatory_mac_col, noninflammatory_mac_col, hep_cols)
  meta_small <- meta %>%
    select(any_of(c(required_meta, optional_meta))) %>%
    mutate(cell_key = paste(fov, cell_ID, sep = "_"))

  qc_cols <- intersect(c("qcFlagsCellCounts", "qcFlagsCellPropNeg", "qcFlagsCellComplex", "qcFlagsCellArea", "qcFlagsFOV"), names(meta_small))
  if (length(qc_cols)) {
    qc_mat <- as.data.frame(meta_small[, qc_cols, drop = FALSE])
    meta_small$qc_pass <- apply(qc_mat, 1, function(x) all(is.na(x) | x == "Pass"))
  } else {
    meta_small$qc_pass <- TRUE
  }

  if (!is.na(stellate_neighbor_col)) {
    meta_small$stellate_neighbors <- meta_small[[stellate_neighbor_col]]
  } else {
    meta_small$stellate_neighbors <- NA_real_
  }
  if (!is.na(portal_endo_col)) {
    meta_small$portal_endothelial_neighbors <- meta_small[[portal_endo_col]]
  } else {
    meta_small$portal_endothelial_neighbors <- NA_real_
  }
  if (!is.na(periportal_lsec_col)) {
    meta_small$periportal_lsec_neighbors <- meta_small[[periportal_lsec_col]]
  } else {
    meta_small$periportal_lsec_neighbors <- NA_real_
  }
  if (!is.na(inflammatory_mac_col)) {
    meta_small$inflammatory_macrophage_neighbors <- meta_small[[inflammatory_mac_col]]
  } else {
    meta_small$inflammatory_macrophage_neighbors <- NA_real_
  }
  if (!is.na(noninflammatory_mac_col)) {
    meta_small$noninflammatory_macrophage_neighbors <- meta_small[[noninflammatory_mac_col]]
  } else {
    meta_small$noninflammatory_macrophage_neighbors <- NA_real_
  }
  if (length(hep_cols)) {
    meta_small$hepatocyte_neighbors <- rowSums(as.matrix(meta_small[, hep_cols, drop = FALSE]), na.rm = TRUE)
  } else {
    meta_small$hepatocyte_neighbors <- NA_real_
  }
  if (!is.na(niche_col)) {
    meta_small$niche <- meta_small[[niche_col]]
  } else {
    meta_small$niche <- NA_character_
  }

  cell_scores <- meta_small %>%
    select(
      cell_key, any_of(c("cell", "cell_id", "fov", "cell_ID", "CenterX_global_px", "CenterY_global_px")),
      qc_pass, niche, stellate_neighbors, portal_endothelial_neighbors, periportal_lsec_neighbors,
      inflammatory_macrophage_neighbors, noninflammatory_macrophage_neighbors, hepatocyte_neighbors
    ) %>%
    left_join(score_df %>% select(-fov, -cell_ID), by = "cell_key") %>%
    mutate(
      gsm = gsm,
      slide_id = slide_id,
      title = title,
      primary_score_z = score_z(primary_score),
      ecm_excluded_score_z = score_z(ecm_excluded_score),
      structural_ecm_score_z = score_z(structural_ecm_score),
      hsc_context_score_z = score_z(hsc_context_score),
      stromal_marker_score_z = score_z(stromal_marker_score)
    )

  slide_qc <- cell_scores %>%
    summarise(
      gsm = first(gsm),
      slide_id = first(slide_id),
      title = first(title),
      n_cells = n(),
      n_qc_pass = sum(qc_pass, na.rm = TRUE),
      qc_pass_fraction = n_qc_pass / n_cells,
      median_total_count = median(total_count, na.rm = TRUE),
      median_n_feature = median(n_feature, na.rm = TRUE),
      median_stellate_neighbors = median(stellate_neighbors, na.rm = TRUE),
      mean_primary_score = mean(primary_score, na.rm = TRUE),
      mean_ecm_excluded_score = mean(ecm_excluded_score, na.rm = TRUE),
      mean_hsc_context_score = mean(hsc_context_score, na.rm = TRUE)
    )

  niche_summary <- cell_scores %>%
    filter(qc_pass, !is.na(niche)) %>%
    group_by(gsm, slide_id, niche) %>%
    summarise(
      n_cells = n(),
      mean_primary_score_z = mean(primary_score_z, na.rm = TRUE),
      mean_ecm_excluded_score_z = mean(ecm_excluded_score_z, na.rm = TRUE),
      mean_structural_ecm_score_z = mean(structural_ecm_score_z, na.rm = TRUE),
      mean_hsc_context_score_z = mean(hsc_context_score_z, na.rm = TRUE),
      mean_stromal_marker_score_z = mean(stromal_marker_score_z, na.rm = TRUE),
      mean_stellate_neighbors = mean(stellate_neighbors, na.rm = TRUE),
      median_stellate_neighbors = median(stellate_neighbors, na.rm = TRUE),
      mean_portal_endothelial_neighbors = mean(portal_endothelial_neighbors, na.rm = TRUE),
      mean_periportal_lsec_neighbors = mean(periportal_lsec_neighbors, na.rm = TRUE),
      mean_hepatocyte_neighbors = mean(hepatocyte_neighbors, na.rm = TRUE),
      .groups = "drop"
    )

  bin_summary <- cell_scores %>%
    filter(qc_pass, !is.na(stellate_neighbors)) %>%
    group_by(gsm, slide_id) %>%
    mutate(stellate_neighbor_bin = ntile(stellate_neighbors, 4)) %>%
    ungroup() %>%
    mutate(
      stellate_neighbor_bin = factor(
        stellate_neighbor_bin,
        levels = 1:4,
        labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")
      )
    ) %>%
    group_by(gsm, slide_id, stellate_neighbor_bin) %>%
    summarise(
      n_cells = n(),
      mean_stellate_neighbors = mean(stellate_neighbors, na.rm = TRUE),
      median_primary_score_z = median(primary_score_z, na.rm = TRUE),
      median_ecm_excluded_score_z = median(ecm_excluded_score_z, na.rm = TRUE),
      median_structural_ecm_score_z = median(structural_ecm_score_z, na.rm = TRUE),
      median_hsc_context_score_z = median(hsc_context_score_z, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    cell_scores = cell_scores,
    slide_qc = slide_qc,
    niche_summary = niche_summary,
    bin_summary = bin_summary,
    panel_genes = gene_cols
  )
}

manifest <- expand_supplementary_links(links_path)
target_manifest <- manifest %>%
  filter(file_kind %in% c("exprMat", "metadata", "fov_positions")) %>%
  arrange(gsm, slide_id, file_kind)

walk2(target_manifest$url_https, target_manifest$local_path, download_if_missing)

target_manifest <- target_manifest %>%
  mutate(
    local_path_relative = file.path("data", "spatial_raw", "GSE292268", basename(local_path)),
    exists_local = file.exists(local_path),
    file_size_mb = round(file.info(local_path)$size / 1024^2, 3)
  )
write_tsv(target_manifest, file.path(result_dir, "gse292268_spatial_file_manifest.tsv"))

slide_groups <- target_manifest %>%
  filter(file_kind %in% c("exprMat", "metadata")) %>%
  count(gsm, title, slide_id, file_kind) %>%
  pivot_wider(names_from = file_kind, values_from = n, values_fill = 0) %>%
  filter(exprMat == 1, metadata == 1)

if (!nrow(slide_groups)) stop("No complete GSE292268 exprMat plus metadata slide pairs found.")

outputs <- vector("list", nrow(slide_groups))
for (i in seq_len(nrow(slide_groups))) {
  slide_files <- target_manifest %>%
    filter(gsm == slide_groups$gsm[[i]], slide_id == slide_groups$slide_id[[i]])
  outputs[[i]] <- process_slide(slide_files)
  gc()
}

cell_scores <- map_dfr(outputs, "cell_scores")
slide_qc <- map_dfr(outputs, "slide_qc")
niche_summary <- map_dfr(outputs, "niche_summary")
bin_summary <- map_dfr(outputs, "bin_summary")
panel_genes <- unique(unlist(map(outputs, "panel_genes")))

coverage <- imap_dfr(gene_sets, function(genes, set_name) {
  present <- intersect(toupper(genes), panel_genes)
  missing <- setdiff(toupper(genes), panel_genes)
  tibble(
    score_name = set_name,
    n_total_genes = length(unique(toupper(genes))),
    n_present_genes = length(present),
    n_missing_genes = length(missing),
    coverage_fraction = n_present_genes / n_total_genes,
    present_genes = paste(sort(present), collapse = ";"),
    missing_genes = paste(sort(missing), collapse = ";"),
    coverage_boundary = "Exact panel-feature match only; ambiguous slash-combined probes were not decomposed."
  )
})

niche_cor <- niche_summary %>%
  group_by(gsm, slide_id) %>%
  summarise(
    n_niches = n(),
    spearman_primary_vs_stellate_neighbors = suppressWarnings(cor(mean_primary_score_z, mean_stellate_neighbors, method = "spearman", use = "pairwise.complete.obs")),
    spearman_ecm_excluded_vs_stellate_neighbors = suppressWarnings(cor(mean_ecm_excluded_score_z, mean_stellate_neighbors, method = "spearman", use = "pairwise.complete.obs")),
    spearman_hsc_context_vs_stellate_neighbors = suppressWarnings(cor(mean_hsc_context_score_z, mean_stellate_neighbors, method = "spearman", use = "pairwise.complete.obs")),
    interpretation_boundary = "Niche-level descriptive localization only; not a donor-level inferential test.",
    .groups = "drop"
  )

write_tsv(cell_scores, file.path(result_dir, "gse292268_spatial_cell_scores.tsv"))
write_tsv(slide_qc, file.path(result_dir, "gse292268_spatial_slide_qc.tsv"))
write_tsv(niche_summary, file.path(result_dir, "gse292268_spatial_niche_summary.tsv"))
write_tsv(bin_summary, file.path(result_dir, "gse292268_spatial_stellate_neighbor_bin_summary.tsv"))
write_tsv(coverage, file.path(result_dir, "gse292268_panel_coverage.tsv"))
write_tsv(niche_cor, file.path(result_dir, "gse292268_spatial_niche_stellate_correlations.tsv"))

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
  geom_col(width = 0.62, fill = "#2A9D8F", color = "#1A1A1A", linewidth = 0.2) +
  geom_text(aes(label = paste0(n_present_genes, "/", n_total_genes)), vjust = -0.35, size = 2.6) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "GSE292268 CosMx panel coverage",
    subtitle = "Exact feature overlap with the locked mechanical-stress-associated program",
    x = NULL,
    y = "Coverage"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "gse292268_panel_coverage.pdf"), p_cov, width = 5.6, height = 3.4, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse292268_panel_coverage.png"), p_cov, width = 5.6, height = 3.4, units = "in", dpi = 600, bg = "white")

map_data <- cell_scores %>%
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
  geom_point(size = 0.07, alpha = 0.92, stroke = 0) +
  facet_wrap(~ slide_label, nrow = 1) +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_gradientn(
    colors = c("#233B7A", "#2A9D8F", "#F4D35E", "#C44E29"),
    limits = c(-2.5, 2.5),
    oob = squish,
    name = "Score z"
  ) +
  labs(
    title = "Spatial localization of the locked program in GSE292268",
    subtitle = "Cell-level CosMx map, z-scored within slide; all QC-passing cells are shown",
    x = "Normalized slide X coordinate",
    y = "Normalized slide Y coordinate"
  ) +
  theme_publication +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "#F7F7F7", color = NA),
    panel.border = element_rect(fill = NA, color = "#B8B8B8", linewidth = 0.25)
  )
ggsave(file.path(fig_dir, "gse292268_program_localization_map.pdf"), p_map, width = 8.2, height = 3.8, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse292268_program_localization_map.png"), p_map, width = 8.2, height = 3.8, units = "in", dpi = 600, bg = "white")

neighbor_map_data <- map_data %>%
  mutate(
    stellate_neighbors_clip = pmin(stellate_neighbors, quantile(stellate_neighbors, 0.995, na.rm = TRUE)),
    stellate_neighbors_clip = if_else(is.finite(stellate_neighbors_clip), stellate_neighbors_clip, NA_real_)
  )

p_neighbor_map <- ggplot(neighbor_map_data, aes(x = x_norm, y = y_norm, color = stellate_neighbors_clip)) +
  geom_point(size = 0.07, alpha = 0.92, stroke = 0) +
  facet_wrap(~ slide_label, nrow = 1) +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_gradientn(
    colors = c("#D9D9D9", "#56B4E9", "#009E73", "#C44E29"),
    name = "Stellate\nneighbors",
    na.value = "#EFEFEF"
  ) +
  labs(
    title = "Spatial distribution of Stellate-cell-neighbor abundance",
    subtitle = "Cell-level neighbor metadata projected on matched slide coordinates",
    x = "Normalized slide X coordinate",
    y = "Normalized slide Y coordinate"
  ) +
  theme_publication +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "#F7F7F7", color = NA),
    panel.border = element_rect(fill = NA, color = "#B8B8B8", linewidth = 0.25)
  )
ggsave(file.path(fig_dir, "gse292268_stellate_neighbor_spatial_map.pdf"), p_neighbor_map, width = 8.2, height = 3.8, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse292268_stellate_neighbor_spatial_map.png"), p_neighbor_map, width = 8.2, height = 3.8, units = "in", dpi = 600, bg = "white")

bubble_data <- niche_summary %>%
  mutate(slide_label = paste0(gsm, " ", slide_id))

p_bubble <- ggplot(bubble_data, aes(x = mean_stellate_neighbors, y = mean_primary_score_z)) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.25, linetype = "dashed") +
  geom_point(aes(size = n_cells, color = mean_hsc_context_score_z), alpha = 0.88) +
  geom_text(aes(label = niche), size = 2.2, vjust = -0.75, check_overlap = TRUE) +
  facet_wrap(~ slide_label, nrow = 1, scales = "free_x") +
  scale_size_continuous(range = c(2.2, 6.2), labels = comma) +
  scale_color_gradientn(
    colors = c("#2D2A7B", "#2AA7A1", "#F6E05E", "#D55E00"),
    name = "HSC context z"
  ) +
  labs(
    title = "Niche-level localization relative to stellate-cell-neighbor-rich regions",
    subtitle = "Each point is one spatial niche; descriptive localization, not donor-level inference",
    x = "Mean Stellate-cell neighbors",
    y = "Mean primary score z",
    size = "Cells"
  ) +
  theme_publication
ggsave(file.path(fig_dir, "gse292268_niche_stellate_localization_bubble.pdf"), p_bubble, width = 7.4, height = 3.8, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse292268_niche_stellate_localization_bubble.png"), p_bubble, width = 7.4, height = 3.8, units = "in", dpi = 600, bg = "white")

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
  geom_line(linewidth = 0.42) +
  geom_point(size = 1.8) +
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
  theme_publication +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(fig_dir, "gse292268_stellate_neighbor_score_gradient.pdf"), p_bin, width = 8.0, height = 3.8, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "gse292268_stellate_neighbor_score_gradient.png"), p_bin, width = 8.0, height = 3.8, units = "in", dpi = 600, bg = "white")

coverage_lines <- paste0(
  "- `", coverage$score_name, "`: ",
  coverage$n_present_genes, "/", coverage$n_total_genes,
  " genes represented by exact panel-feature matches.\n",
  collapse = ""
)

slide_lines <- paste0(
  "- `", slide_qc$gsm, " ", slide_qc$slide_id, "`: ",
  slide_qc$n_cells, " cells; QC-pass fraction ",
  percent(slide_qc$qc_pass_fraction, accuracy = 0.1),
  "; median Stellate-cell neighbors ",
  signif(slide_qc$median_stellate_neighbors, 3), ".\n",
  collapse = ""
)

cor_lines <- paste0(
  "- `", niche_cor$gsm, " ", niche_cor$slide_id, "`: niche-level Spearman rho primary score vs Stellate-cell neighbors = ",
  signif(niche_cor$spearman_primary_vs_stellate_neighbors, 3),
  "; ECM-excluded score rho = ",
  signif(niche_cor$spearman_ecm_excluded_vs_stellate_neighbors, 3),
  "; HSC-context score rho = ",
  signif(niche_cor$spearman_hsc_context_vs_stellate_neighbors, 3),
  ".\n",
  collapse = ""
)

cat(
  "# GSE292268 Spatial Localization Summary\n\n",
  "Generated by `scripts/09_gse292268_spatial_localization.R`.\n\n",
  "## Outputs\n\n",
  "- `results/spatial/gse292268_spatial_file_manifest.tsv`\n",
  "- `results/spatial/gse292268_panel_coverage.tsv`\n",
  "- `results/spatial/gse292268_spatial_cell_scores.tsv`\n",
  "- `results/spatial/gse292268_spatial_slide_qc.tsv`\n",
  "- `results/spatial/gse292268_spatial_niche_summary.tsv`\n",
  "- `results/spatial/gse292268_spatial_stellate_neighbor_bin_summary.tsv`\n",
  "- `results/spatial/gse292268_spatial_niche_stellate_correlations.tsv`\n",
  "- `figures/spatial/gse292268_panel_coverage.png` and `.pdf`\n",
  "- `figures/spatial/gse292268_program_localization_map.png` and `.pdf`\n",
  "- `figures/spatial/gse292268_stellate_neighbor_spatial_map.png` and `.pdf`\n",
  "- `figures/spatial/gse292268_niche_stellate_localization_bubble.png` and `.pdf`\n",
  "- `figures/spatial/gse292268_stellate_neighbor_score_gradient.png` and `.pdf`\n\n",
  "## Panel Coverage\n\n",
  coverage_lines,
  "\n## Slide QC\n\n",
  slide_lines,
  "\n## Niche-Level Localization Readout\n\n",
  cor_lines,
  "\n## Analysis Boundary\n\n",
  "This spatial analysis uses GSE292268 human MASLD CosMx expression and metadata files for localization only. It does not test fibrosis severity, donor-level association, cell-cell communication mechanism, or mechanical memory. Spatial cells and niches are descriptive localization units, not independent patients. The appropriate wording is that the locked program shows spatial localization patterns relative to Stellate-cell-neighbor-rich niches, consistent with a fibrogenic stromal mechanostress hypothesis.\n",
  file = file.path(doc_dir, "gse292268_spatial_localization_summary.md"),
  sep = ""
)

message("Done. GSE292268 spatial localization outputs written to ", result_dir)

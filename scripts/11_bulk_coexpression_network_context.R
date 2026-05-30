#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(igraph)
  library(ggraph)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

expr_dir <- file.path(project_dir, "data", "expression_processed")
result_dir <- file.path(project_dir, "results", "network")
figure_dir <- file.path(project_dir, "figures", "network")
doc_dir <- file.path(project_dir, "docs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

fibrosis_cohorts <- c("GSE135251", "GSE130970", "GSE162694")

theme_publication <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size, color = "#333333"),
      plot.caption = element_text(size = base_size - 1, color = "#555555", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#111111"),
      axis.line = element_line(linewidth = 0.3, color = "#111111"),
      axis.ticks = element_line(linewidth = 0.25, color = "#111111"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

zscore <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

aggregate_duplicate_rows <- function(mat) {
  rownames(mat) <- toupper(rownames(mat))
  if (!anyDuplicated(rownames(mat))) return(mat)
  group <- rownames(mat)
  summed <- rowsum(mat, group = group, reorder = FALSE)
  counts <- as.numeric(table(factor(group, levels = rownames(summed))))
  summed / counts
}

safe_logistic_gene <- function(expr_z, advanced_binary) {
  dat <- tibble(expr_z = expr_z, advanced_binary = advanced_binary) %>% drop_na()
  if (nrow(dat) < 10 || length(unique(dat$advanced_binary)) < 2 || any(table(dat$advanced_binary) < 3)) {
    return(tibble(n_model = nrow(dat), beta = NA_real_, estimate_or = NA_real_, p_value = NA_real_))
  }
  fit <- tryCatch(
    suppressWarnings(glm(advanced_binary ~ expr_z, data = dat, family = binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(n_model = nrow(dat), beta = NA_real_, estimate_or = NA_real_, p_value = NA_real_))
  }
  coefs <- summary(fit)$coefficients
  if (!"expr_z" %in% rownames(coefs)) {
    return(tibble(n_model = nrow(dat), beta = NA_real_, estimate_or = NA_real_, p_value = NA_real_))
  }
  beta <- unname(coefs["expr_z", "Estimate"])
  tibble(
    n_model = nrow(dat),
    beta = beta,
    estimate_or = exp(beta),
    p_value = unname(coefs["expr_z", "Pr(>|z|)"])
  )
}

safe_spearman_gene <- function(expr_z, trait) {
  dat <- tibble(expr_z = expr_z, trait = trait) %>% drop_na()
  if (nrow(dat) < 10 || length(unique(dat$trait)) < 3) {
    return(tibble(n_model = nrow(dat), rho = NA_real_, p_value = NA_real_))
  }
  ct <- tryCatch(
    suppressWarnings(cor.test(dat$expr_z, dat$trait, method = "spearman", exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(ct)) return(tibble(n_model = nrow(dat), rho = NA_real_, p_value = NA_real_))
  tibble(n_model = nrow(dat), rho = unname(ct$estimate), p_value = ct$p.value)
}

safe_spearman_pair <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(tibble(n_pair = sum(ok), rho = NA_real_, p_value = NA_real_))
  ct <- tryCatch(
    suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(ct)) return(tibble(n_pair = sum(ok), rho = NA_real_, p_value = NA_real_))
  tibble(n_pair = sum(ok), rho = unname(ct$estimate), p_value = ct$p.value)
}

signature <- read_tsv(file.path(project_dir, "data", "signatures", "signature_v1_locked.tsv"), show_col_types = FALSE) %>%
  mutate(
    gene_symbol = toupper(gene_symbol),
    include_in_primary_score = as.logical(include_in_primary_score),
    include_in_ecm_excluded_score = as.logical(include_in_ecm_excluded_score),
    is_ecm_structural_gene = as.logical(is_ecm_structural_gene)
  ) %>%
  distinct(gene_symbol, .keep_all = TRUE)

primary_genes <- signature %>%
  filter(include_in_primary_score) %>%
  pull(gene_symbol) %>%
  unique()

spatial_present <- character()
coverage_path <- file.path(project_dir, "results", "spatial", "gse292268_panel_coverage.tsv")
if (file.exists(coverage_path)) {
  spatial_present <- read_tsv(coverage_path, show_col_types = FALSE) %>%
    filter(score_name == "primary") %>%
    pull(present_genes) %>%
    str_split(";") %>%
    unlist(use.names = FALSE) %>%
    toupper() %>%
    unique()
}

ppi_table_path <- file.path(result_dir, "network_prioritization_centrality_table.tsv")
if (!file.exists(ppi_table_path)) {
  ppi_table_path <- file.path(result_dir, "zhao_method_hub_centrality_table.tsv")
}
ppi_prior <- if (file.exists(ppi_table_path)) {
  read_tsv(ppi_table_path, show_col_types = FALSE) %>%
    select(gene_symbol, degree, n_algorithms_top_quartile, network_priority_group)
} else {
  tibble(gene_symbol = primary_genes, degree = NA_real_, n_algorithms_top_quartile = NA_real_, network_priority_group = "not_run")
}

expr_files <- file.path(expr_dir, paste0(fibrosis_cohorts, "_bulk_expression.rds"))
names(expr_files) <- fibrosis_cohorts
expr_files <- expr_files[file.exists(expr_files)]
if (!length(expr_files)) stop("No fibrosis-mapped expression objects found.")

gene_assoc_tables <- list()
edge_tables <- list()
cohort_gene_presence <- list()

for (acc in names(expr_files)) {
  obj <- readRDS(expr_files[[acc]])
  log_expr <- as.matrix(obj$logcpm)
  log_expr <- aggregate_duplicate_rows(log_expr)
  pheno <- obj$pheno %>%
    mutate(
      fibrosis_stage = suppressWarnings(as.numeric(fibrosis_stage)),
      advanced_binary = case_when(
        advanced_fibrosis == "advanced" ~ 1,
        advanced_fibrosis == "non_advanced" ~ 0,
        TRUE ~ NA_real_
      )
    )

  if (ncol(log_expr) != nrow(pheno)) {
    stop("Expression and phenotype dimensions do not match for ", acc)
  }

  present <- intersect(primary_genes, rownames(log_expr))
  cohort_gene_presence[[acc]] <- tibble(
    accession = acc,
    gene_symbol = primary_genes,
    present = gene_symbol %in% present
  )

  gene_assoc_tables[[acc]] <- map_dfr(present, function(gene) {
    expr_z <- zscore(as.numeric(log_expr[gene, ]))
    adv <- safe_logistic_gene(expr_z, pheno$advanced_binary)
    fib <- safe_spearman_gene(expr_z, pheno$fibrosis_stage)
    tibble(
      accession = acc,
      gene_symbol = gene,
      advanced_n = adv$n_model,
      advanced_beta = adv$beta,
      advanced_or = adv$estimate_or,
      advanced_p = adv$p_value,
      fibrosis_n = fib$n_model,
      fibrosis_rho = fib$rho,
      fibrosis_p = fib$p_value
    )
  }) %>%
    mutate(
      advanced_fdr = p.adjust(advanced_p, method = "BH"),
      fibrosis_fdr = p.adjust(fibrosis_p, method = "BH"),
      advanced_direction = case_when(
        advanced_beta > 0 ~ "positive",
        advanced_beta < 0 ~ "negative",
        TRUE ~ "not_estimable"
      ),
      fibrosis_direction = case_when(
        fibrosis_rho > 0 ~ "positive",
        fibrosis_rho < 0 ~ "negative",
        TRUE ~ "not_estimable"
      )
    )

  expr_sub <- t(log_expr[present, , drop = FALSE])
  pair_grid <- combn(present, 2, simplify = FALSE)
  edge_tables[[acc]] <- map_dfr(pair_grid, function(pair) {
    res <- safe_spearman_pair(expr_sub[, pair[[1]]], expr_sub[, pair[[2]]])
    tibble(
      accession = acc,
      gene_a = min(pair),
      gene_b = max(pair),
      n_pair = res$n_pair,
      rho = res$rho,
      p_value = res$p_value
    )
  }) %>%
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      edge_significant = !is.na(fdr) & fdr < 0.05 & abs(rho) >= 0.30,
      edge_direction = case_when(
        rho > 0 ~ "positive",
        rho < 0 ~ "negative",
        TRUE ~ "not_estimable"
      ),
      boundary = "Bulk coexpression edge among locked genes; association only, not direct interaction or causal regulation."
    )
}

gene_assoc <- bind_rows(gene_assoc_tables) %>%
  right_join(bind_rows(cohort_gene_presence), by = c("accession", "gene_symbol")) %>%
  left_join(signature %>% select(gene_symbol, source_category, evidence_tier, include_in_ecm_excluded_score, is_ecm_structural_gene), by = "gene_symbol") %>%
  mutate(
    boundary = "Gene-level bulk association within locked program; not a mechanism validation or feature-selection step."
  )

coexpression_edges <- bind_rows(edge_tables)

write_tsv(gene_assoc, file.path(result_dir, "bulk_locked_gene_fibrosis_associations.tsv"))
write_tsv(coexpression_edges, file.path(result_dir, "bulk_locked_gene_coexpression_edges.tsv"))

gene_consistency <- gene_assoc %>%
  group_by(gene_symbol) %>%
  summarise(
    n_cohorts_present = sum(present, na.rm = TRUE),
    n_advanced_positive = sum(advanced_beta > 0, na.rm = TRUE),
    n_advanced_positive_fdr05 = sum(advanced_beta > 0 & advanced_fdr < 0.05, na.rm = TRUE),
    n_fibrosis_positive = sum(fibrosis_rho > 0, na.rm = TRUE),
    n_fibrosis_positive_fdr05 = sum(fibrosis_rho > 0 & fibrosis_fdr < 0.05, na.rm = TRUE),
    median_advanced_beta = median(advanced_beta, na.rm = TRUE),
    median_fibrosis_rho = median(fibrosis_rho, na.rm = TRUE),
    min_advanced_fdr = suppressWarnings(min(advanced_fdr, na.rm = TRUE)),
    min_fibrosis_fdr = suppressWarnings(min(fibrosis_fdr, na.rm = TRUE)),
    source_category = first(na.omit(source_category)),
    evidence_tier = first(na.omit(evidence_tier)),
    include_in_ecm_excluded_score = first(na.omit(include_in_ecm_excluded_score)),
    is_ecm_structural_gene = first(na.omit(is_ecm_structural_gene)),
    .groups = "drop"
  ) %>%
  mutate(
    min_advanced_fdr = if_else(is.infinite(min_advanced_fdr), NA_real_, min_advanced_fdr),
    min_fibrosis_fdr = if_else(is.infinite(min_fibrosis_fdr), NA_real_, min_fibrosis_fdr),
    bulk_consistency_group = case_when(
      n_advanced_positive_fdr05 >= 2 | n_fibrosis_positive_fdr05 >= 2 ~ "replicated_positive_fdr05",
      n_advanced_positive >= 3 | n_fibrosis_positive >= 3 ~ "directionally_consistent_positive",
      TRUE ~ "limited_or_mixed"
    )
  ) %>%
  left_join(ppi_prior, by = "gene_symbol") %>%
  mutate(
    spatial_panel_present = gene_symbol %in% spatial_present,
    boundary = "Cross-cohort consistency summary for prioritization only; not used to alter the locked score."
  ) %>%
  arrange(
    desc(n_advanced_positive_fdr05 + n_fibrosis_positive_fdr05),
    desc(n_advanced_positive + n_fibrosis_positive),
    desc(n_algorithms_top_quartile),
    desc(abs(median_fibrosis_rho)),
    gene_symbol
  )

write_tsv(gene_consistency, file.path(result_dir, "bulk_locked_gene_consistency_summary.tsv"))

edge_summary <- coexpression_edges %>%
  group_by(gene_a, gene_b) %>%
  summarise(
    n_cohorts_tested = sum(!is.na(rho)),
    n_significant_edges = sum(edge_significant, na.rm = TRUE),
    n_positive_edges = sum(edge_significant & rho > 0, na.rm = TRUE),
    n_negative_edges = sum(edge_significant & rho < 0, na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    median_abs_rho = median(abs(rho), na.rm = TRUE),
    min_fdr = suppressWarnings(min(fdr, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    min_fdr = if_else(is.infinite(min_fdr), NA_real_, min_fdr),
    reproducible_edge = n_significant_edges >= 2 & n_positive_edges >= 2,
    boundary = "Reproducible bulk coexpression edge; not a direct interaction or causal edge."
  ) %>%
  arrange(desc(reproducible_edge), desc(n_significant_edges), desc(median_abs_rho), gene_a, gene_b)

write_tsv(edge_summary, file.path(result_dir, "bulk_locked_gene_coexpression_edge_summary.tsv"))

coexpression_degree <- edge_summary %>%
  filter(reproducible_edge) %>%
  select(gene_a, gene_b) %>%
  pivot_longer(everything(), values_to = "gene_symbol") %>%
  count(gene_symbol, name = "reproducible_coexpression_degree")

evidence_matrix <- gene_consistency %>%
  left_join(coexpression_degree, by = "gene_symbol") %>%
  mutate(
    reproducible_coexpression_degree = replace_na(reproducible_coexpression_degree, 0L),
    ppi_high_priority = n_algorithms_top_quartile >= 3,
    bulk_replicated_positive = bulk_consistency_group == "replicated_positive_fdr05",
    coexpression_reproducible = reproducible_coexpression_degree > 0,
    spatial_panel_present = spatial_panel_present,
    structural_ecm = is_ecm_structural_gene
  ) %>%
  select(
    gene_symbol, source_category, evidence_tier, structural_ecm,
    ppi_high_priority, bulk_replicated_positive, coexpression_reproducible,
    spatial_panel_present, reproducible_coexpression_degree,
    n_algorithms_top_quartile, n_advanced_positive_fdr05, n_fibrosis_positive_fdr05,
    median_advanced_beta, median_fibrosis_rho, boundary
  )
write_tsv(evidence_matrix, file.path(result_dir, "network_prioritization_gene_evidence_matrix.tsv"))

top_heatmap_genes <- gene_consistency %>%
  slice_head(n = min(32L, nrow(gene_consistency))) %>%
  pull(gene_symbol)

heatmap_data <- gene_assoc %>%
  filter(gene_symbol %in% top_heatmap_genes) %>%
  mutate(
    accession = factor(accession, levels = fibrosis_cohorts),
    gene_symbol = factor(gene_symbol, levels = rev(top_heatmap_genes)),
    fdr_for_size = pmin(-log10(fibrosis_fdr), 6),
    fdr_for_size = if_else(is.finite(fdr_for_size), fdr_for_size, 0)
  )

p_gene_heatmap <- ggplot(heatmap_data, aes(x = accession, y = gene_symbol)) +
  geom_tile(aes(fill = fibrosis_rho), color = "white", linewidth = 0.25) +
  geom_point(aes(size = fdr_for_size), shape = 21, fill = "white", color = "#111111", stroke = 0.15, alpha = 0.85) +
  scale_fill_gradient2(
    low = "#2D2A7B",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-0.65, 0.65),
    oob = squish,
    name = "Spearman rho"
  ) +
  scale_size_continuous(range = c(0.2, 2.8), breaks = c(0, 2, 4, 6), name = "-log10 FDR") +
  labs(
    title = "Gene-level fibrosis association across bulk cohorts",
    subtitle = "Top locked genes prioritized by cross-cohort direction, FDR support, and network context",
    x = NULL,
    y = NULL,
    caption = "Gene-level associations are exploratory and do not redefine the locked score."
  ) +
  theme_publication(8.2) +
  theme(axis.text.y = element_text(size = 7.2), axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_gene_heatmap, "bulk_locked_gene_fibrosis_consistency_heatmap", 6.2, 8.4)

top_corr_genes <- evidence_matrix %>%
  arrange(desc(bulk_replicated_positive), desc(coexpression_reproducible), desc(ppi_high_priority), desc(abs(median_fibrosis_rho))) %>%
  slice_head(n = min(24L, nrow(.))) %>%
  pull(gene_symbol)

corr_heatmap <- edge_summary %>%
  filter(gene_a %in% top_corr_genes, gene_b %in% top_corr_genes) %>%
  select(gene_a, gene_b, median_rho) %>%
  bind_rows(tibble(gene_a = top_corr_genes, gene_b = top_corr_genes, median_rho = 1)) %>%
  bind_rows(
    edge_summary %>%
      filter(gene_a %in% top_corr_genes, gene_b %in% top_corr_genes) %>%
      transmute(gene_a = gene_b, gene_b = gene_a, median_rho)
  ) %>%
  mutate(
    gene_a = factor(gene_a, levels = top_corr_genes),
    gene_b = factor(gene_b, levels = rev(top_corr_genes))
  )

p_corr <- ggplot(corr_heatmap, aes(x = gene_a, y = gene_b, fill = median_rho)) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_gradient2(
    low = "#2D2A7B",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-1, 1),
    oob = squish,
    name = "Median rho"
  ) +
  labs(
    title = "Bulk coexpression among prioritized locked genes",
    subtitle = "Median Spearman correlation across fibrosis-mapped cohorts",
    x = NULL,
    y = NULL,
    caption = "Coexpression is sample-level coordination, not direct interaction or causality."
  ) +
  theme_publication(8.0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.0), axis.text.y = element_text(size = 7.0))
save_plot(p_corr, "bulk_locked_gene_coexpression_heatmap", 7.4, 7.2)

network_edges <- edge_summary %>%
  filter(reproducible_edge, gene_a %in% top_corr_genes, gene_b %in% top_corr_genes) %>%
  slice_head(n = min(120L, nrow(.))) %>%
  transmute(from = gene_a, to = gene_b, median_rho, n_significant_edges, min_fdr)

if (nrow(network_edges) >= 1) {
  network_nodes <- evidence_matrix %>%
    filter(gene_symbol %in% unique(c(network_edges$from, network_edges$to))) %>%
    mutate(
      node_class = case_when(
        structural_ecm ~ "Structural ECM",
        ppi_high_priority & bulk_replicated_positive ~ "PPI + bulk",
        ppi_high_priority ~ "PPI-prioritized",
        bulk_replicated_positive ~ "Bulk replicated",
        TRUE ~ "Context"
      )
    )

  g_net <- graph_from_data_frame(network_edges, directed = FALSE, vertices = network_nodes %>% rename(name = gene_symbol))

  p_net <- ggraph(g_net, layout = "fr") +
    geom_edge_link(aes(width = abs(median_rho), alpha = n_significant_edges), color = "#666666") +
    geom_node_point(aes(size = reproducible_coexpression_degree + 1, fill = node_class), shape = 21, color = "#111111", stroke = 0.25) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.4, max.overlaps = Inf) +
    scale_edge_width(range = c(0.2, 1.1), guide = "none") +
    scale_edge_alpha(range = c(0.3, 0.8), guide = "none") +
    scale_size_continuous(range = c(2.4, 7.0), name = "Replicated edge degree") +
    scale_fill_manual(values = c(
      "Structural ECM" = "#D55E00",
      "PPI + bulk" = "#7B3294",
      "PPI-prioritized" = "#0072B2",
      "Bulk replicated" = "#2A9D8F",
      "Context" = "#999999"
    )) +
    labs(
      title = "Reproducible bulk coexpression context",
      subtitle = "Edges require positive coexpression with FDR < 0.05 and rho >= 0.30 in at least two cohorts",
      fill = "Node evidence",
      caption = "Network is a candidate context map only; it does not validate regulatory or physical interactions."
    ) +
    theme_void(base_family = "Arial", base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, color = "#333333"),
      plot.caption = element_text(size = 7, color = "#555555", hjust = 0),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 8)
    )
  save_plot(p_net, "bulk_locked_gene_coexpression_network", 7.4, 5.4)
}

matrix_top_genes <- evidence_matrix %>%
  arrange(desc(bulk_replicated_positive), desc(coexpression_reproducible), desc(ppi_high_priority), desc(spatial_panel_present), desc(abs(median_fibrosis_rho))) %>%
  slice_head(n = min(30L, nrow(.))) %>%
  pull(gene_symbol)

evidence_plot_data <- evidence_matrix %>%
  filter(gene_symbol %in% matrix_top_genes) %>%
  mutate(
    structural_ecm = if_else(structural_ecm, "Structural ECM", "Non-structural or remodeling"),
    ppi_high_priority = if_else(ppi_high_priority, "Present", "Absent"),
    bulk_replicated_positive = if_else(bulk_replicated_positive, "Present", "Absent"),
    coexpression_reproducible = if_else(coexpression_reproducible, "Present", "Absent"),
    spatial_panel_present = if_else(spatial_panel_present, "Present", "Absent")
  ) %>%
  select(gene_symbol, structural_ecm, ppi_high_priority, bulk_replicated_positive, coexpression_reproducible, spatial_panel_present) %>%
  pivot_longer(-gene_symbol, names_to = "evidence_layer", values_to = "status") %>%
  mutate(
    evidence_layer = recode(
      evidence_layer,
      structural_ecm = "Structural ECM flag",
      ppi_high_priority = "PPI high priority",
      bulk_replicated_positive = "Bulk FDR support",
      coexpression_reproducible = "Replicated coexpression",
      spatial_panel_present = "CosMx panel present"
    ),
    evidence_layer = factor(evidence_layer, levels = c("Bulk FDR support", "Replicated coexpression", "PPI high priority", "CosMx panel present", "Structural ECM flag")),
    gene_symbol = factor(gene_symbol, levels = rev(matrix_top_genes)),
    status_plot = case_when(
      evidence_layer == "Structural ECM flag" & status == "Structural ECM" ~ "Structural ECM",
      evidence_layer == "Structural ECM flag" ~ "Non-structural",
      status == "Present" ~ "Present",
      TRUE ~ "Absent"
    )
  )

p_evidence <- ggplot(evidence_plot_data, aes(x = evidence_layer, y = gene_symbol, fill = status_plot)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Present" = "#2A9D8F",
    "Absent" = "#E8E8E8",
    "Structural ECM" = "#D55E00",
    "Non-structural" = "#B8B8B8"
  )) +
  labs(
    title = "Integrated candidate evidence matrix",
    subtitle = "Bulk, coexpression, PPI, and CosMx panel support for prioritized locked genes",
    x = NULL,
    y = NULL,
    fill = "Status",
    caption = "This matrix prioritizes candidates for follow-up; it is not multi-omics proof of mechanical memory."
  ) +
  theme_publication(8.0) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.text.y = element_text(size = 7.0))
save_plot(p_evidence, "network_prioritization_integrated_candidate_evidence_matrix", 6.6, 7.0)

fig_manifest <- tibble(
  file_stem = c(
    "bulk_locked_gene_fibrosis_consistency_heatmap",
    "bulk_locked_gene_coexpression_heatmap",
    "bulk_locked_gene_coexpression_network",
    "network_prioritization_integrated_candidate_evidence_matrix"
  ),
  source_data = c(
    "results/network/bulk_locked_gene_fibrosis_associations.tsv; results/network/bulk_locked_gene_consistency_summary.tsv",
    "results/network/bulk_locked_gene_coexpression_edge_summary.tsv",
    "results/network/bulk_locked_gene_coexpression_edge_summary.tsv; results/network/network_prioritization_gene_evidence_matrix.tsv",
    "results/network/network_prioritization_gene_evidence_matrix.tsv"
  ),
  role = c(
    "gene-level fibrosis association context",
    "inter-gene coexpression context",
    "reproducible coexpression network context",
    "integrated prioritization matrix"
  ),
  boundary = c(
    "Exploratory gene-level association; not signature redefinition.",
    "Coexpression among samples; not direct interaction.",
    "Candidate context network; not validated regulatory mechanism.",
    "Candidate prioritization only; not proof of mechanical memory."
  )
)
write_tsv(fig_manifest, file.path(result_dir, "bulk_coexpression_figure_manifest.tsv"))

summary_stats <- tibble(
  n_primary_genes = length(primary_genes),
  n_genes_present_all_three = sum(gene_consistency$n_cohorts_present == length(fibrosis_cohorts)),
  n_bulk_replicated_positive = sum(gene_consistency$bulk_consistency_group == "replicated_positive_fdr05"),
  n_directionally_consistent_positive = sum(gene_consistency$bulk_consistency_group == "directionally_consistent_positive"),
  n_reproducible_coexpression_edges = sum(edge_summary$reproducible_edge, na.rm = TRUE),
  n_evidence_matrix_genes = nrow(evidence_matrix),
  boundary = "Summary for bulk gene-level/coexpression context only."
)
write_tsv(summary_stats, file.path(result_dir, "bulk_coexpression_network_summary.tsv"))

top_gene_lines <- gene_consistency %>%
  slice_head(n = min(15L, nrow(.))) %>%
  mutate(line = paste0(
    "- `", gene_symbol, "`: bulk group ", bulk_consistency_group,
    "; advanced FDR-positive cohorts ", n_advanced_positive_fdr05,
    "; fibrosis FDR-positive cohorts ", n_fibrosis_positive_fdr05,
    "; median fibrosis rho ", signif(median_fibrosis_rho, 3),
    "; PPI priority ", ifelse(is.na(n_algorithms_top_quartile), "not run", paste0(n_algorithms_top_quartile, "/5")),
    "."
  )) %>%
  pull(line) %>%
  paste(collapse = "\n")

cat(
  "# Bulk Gene-Level and Coexpression Network Context\n\n",
  "Generated by `scripts/11_bulk_coexpression_network_context.R` on 2026-05-23.\n\n",
  "## Analysis Boundary\n\n",
  "This module evaluates gene-level fibrosis associations and inter-gene coexpression among the already locked primary program genes. It does not discover a new signature or change the locked score.\n\n",
  "## Scope\n\n",
  "- Fibrosis-mapped cohorts: `", paste(fibrosis_cohorts, collapse = "`, `"), "`.\n",
  "- Input universe: locked primary score genes only.\n",
  "- Gene-level readout: univariate advanced-fibrosis logistic association and Spearman association with ordinal fibrosis stage within each cohort.\n",
  "- Coexpression readout: Spearman gene-pair correlation within each cohort; reproducible edges require positive rho >= 0.30 and FDR < 0.05 in at least two cohorts.\n\n",
  "## Summary\n\n",
  "- Primary genes evaluated: ", summary_stats$n_primary_genes, ".\n",
  "- Genes present in all three fibrosis-mapped cohorts: ", summary_stats$n_genes_present_all_three, ".\n",
  "- Genes with replicated positive FDR support: ", summary_stats$n_bulk_replicated_positive, ".\n",
  "- Directionally consistent positive genes without replicated FDR support: ", summary_stats$n_directionally_consistent_positive, ".\n",
  "- Reproducible positive coexpression edges: ", summary_stats$n_reproducible_coexpression_edges, ".\n\n",
  "## Top Candidate Context Genes\n\n",
  top_gene_lines,
  "\n\n## Outputs\n\n",
  "- `results/network/bulk_locked_gene_fibrosis_associations.tsv`\n",
  "- `results/network/bulk_locked_gene_consistency_summary.tsv`\n",
  "- `results/network/bulk_locked_gene_coexpression_edges.tsv`\n",
  "- `results/network/bulk_locked_gene_coexpression_edge_summary.tsv`\n",
  "- `results/network/network_prioritization_gene_evidence_matrix.tsv`\n",
  "- `results/network/bulk_coexpression_figure_manifest.tsv`\n",
  "- `figures/network/bulk_locked_gene_fibrosis_consistency_heatmap.png` and `.pdf`\n",
  "- `figures/network/bulk_locked_gene_coexpression_heatmap.png` and `.pdf`\n",
  "- `figures/network/bulk_locked_gene_coexpression_network.png` and `.pdf` if reproducible edges are available.\n",
  "- `figures/network/network_prioritization_integrated_candidate_evidence_matrix.png` and `.pdf`\n",
  file = file.path(doc_dir, "bulk_coexpression_network_context_summary.md"),
  sep = ""
)

message("Done. Bulk gene-level and coexpression context outputs written to ", result_dir)

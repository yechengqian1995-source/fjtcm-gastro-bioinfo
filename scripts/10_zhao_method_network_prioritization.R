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
  library(UpSetR)
  library(httr2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

result_dir <- file.path(project_dir, "results", "network")
figure_dir <- file.path(project_dir, "figures", "network")
doc_dir <- file.path(project_dir, "docs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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
      panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.25),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
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

fetch_string_edges <- function(genes, required_score = 700L) {
  api_url <- "https://string-db.org/api/tsv/network"
  identifiers <- paste(unique(genes), collapse = "\r")
  req <- request(api_url) %>%
    req_user_agent("masld_mechanical_memory/1.0") %>%
    req_url_query(
      identifiers = identifiers,
      species = 9606,
      required_score = required_score,
      caller_identity = "masld_mechanical_memory"
    ) %>%
    req_timeout(120)

  resp <- req_perform(req)
  txt <- resp_body_string(resp)
  if (!nzchar(txt)) stop("STRING API returned empty response.")

  edges <- read_tsv(I(txt), show_col_types = FALSE, progress = FALSE)
  if (!nrow(edges)) stop("STRING API returned no edges at the selected threshold.")
  edges
}

string_status <- tibble(
  source = "STRING REST API",
  species = "Homo sapiens 9606",
  required_score = 700L,
  evidence_channel_note = "STRING functional association network; edges may include multiple evidence channels and are not interpreted as direct physical binding.",
  query_status = "not_run",
  error_message = NA_character_
)

string_raw <- tryCatch(
  {
    string_status$query_status <- "success"
    fetch_string_edges(primary_genes, required_score = 700L)
  },
  error = function(e) {
    string_status$query_status <<- "failed"
    string_status$error_message <<- conditionMessage(e)
    tibble()
  }
)

if (!nrow(string_raw)) {
  write_tsv(string_status, file.path(result_dir, "network_prioritization_string_query_status.tsv"))
  cat(
    "# STRING Network Candidate Prioritization\n\n",
    "The STRING high-confidence network query failed or returned no edges. No PPI/network-prioritized candidate nodes were generated.\n\n",
    "This missing network result does not affect the locked bulk score, single-nucleus stromal context, or spatial localization analyses.\n",
    file = file.path(doc_dir, "network_prioritization_summary.md"),
    sep = ""
  )
  stop("STRING network unavailable; status file and summary written.")
}

write_tsv(string_status, file.path(result_dir, "network_prioritization_string_query_status.tsv"))
write_tsv(string_raw, file.path(result_dir, "network_prioritization_string_raw_edges.tsv"))

edge_cols_needed <- c("preferredName_A", "preferredName_B", "score")
if (!all(edge_cols_needed %in% names(string_raw))) {
  stop("STRING response did not include required columns: preferredName_A, preferredName_B, score.")
}

edges <- string_raw %>%
  transmute(
    from = toupper(preferredName_A),
    to = toupper(preferredName_B),
    string_score = as.numeric(score),
    nscore = suppressWarnings(as.numeric(.data[["nscore"]])),
    fscore = suppressWarnings(as.numeric(.data[["fscore"]])),
    pscore = suppressWarnings(as.numeric(.data[["pscore"]])),
    ascore = suppressWarnings(as.numeric(.data[["ascore"]])),
    escore = suppressWarnings(as.numeric(.data[["escore"]])),
    dscore = suppressWarnings(as.numeric(.data[["dscore"]])),
    tscore = suppressWarnings(as.numeric(.data[["tscore"]]))
  ) %>%
  filter(from %in% primary_genes, to %in% primary_genes, from != to) %>%
  mutate(
    gene_a = pmin(from, to),
    gene_b = pmax(from, to)
  ) %>%
  distinct(gene_a, gene_b, .keep_all = TRUE) %>%
  select(from = gene_a, to = gene_b, string_score, nscore, fscore, pscore, ascore, escore, dscore, tscore)

if (!nrow(edges)) stop("No STRING edges remained after filtering to locked primary genes.")

write_tsv(
  edges %>%
    mutate(boundary = "STRING network edge used for candidate prioritization only; not direct mechanistic validation."),
  file.path(result_dir, "network_prioritization_string_filtered_edges.tsv")
)

g <- graph_from_data_frame(edges %>% select(from, to, string_score), directed = FALSE, vertices = tibble(name = primary_genes))
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(string_score = "max"))

calc_mcc_like <- function(graph) {
  out <- setNames(rep(0, vcount(graph)), V(graph)$name)
  clq <- max_cliques(graph, min = 2)
  if (!length(clq)) return(out)
  for (clique in clq) {
    genes <- as_ids(clique)
    score <- factorial(max(length(genes) - 1, 1))
    out[genes] <- out[genes] + score
  }
  out
}

calc_mnc_like <- function(graph) {
  setNames(map_dbl(V(graph), function(v) {
    nbr <- neighbors(graph, v)
    if (length(nbr) == 0) return(0)
    sub <- induced_subgraph(graph, nbr)
    if (vcount(sub) == 0) return(0)
    max(components(sub)$csize)
  }), V(graph)$name)
}

safe_metric <- function(x) {
  x[!is.finite(x)] <- 0
  x
}

degree_score <- safe_metric(degree(g, mode = "all", normalized = FALSE))
mcc_score <- safe_metric(calc_mcc_like(g))
mnc_score <- safe_metric(calc_mnc_like(g))
closeness_score <- safe_metric(closeness(g, mode = "all", normalized = TRUE))
eigen_score <- safe_metric(eigen_centrality(g, directed = FALSE, scale = TRUE)$vector)

node_metrics <- tibble(
  gene_symbol = V(g)$name,
  degree = as.numeric(degree_score[gene_symbol]),
  mcc_like = as.numeric(mcc_score[gene_symbol]),
  mnc_like = as.numeric(mnc_score[gene_symbol]),
  closeness = as.numeric(closeness_score[gene_symbol]),
  eigenvector_epc_surrogate = as.numeric(eigen_score[gene_symbol])
) %>%
  left_join(
    signature %>%
      select(gene_symbol, source_category, evidence_tier, include_in_ecm_excluded_score, is_ecm_structural_gene),
    by = "gene_symbol"
  )

rank_desc <- function(x) {
  rank(-x, ties.method = "min", na.last = "keep")
}

metric_names <- c("degree", "mcc_like", "mnc_like", "closeness", "eigenvector_epc_surrogate")
top_n <- max(1L, ceiling(0.25 * nrow(node_metrics)))

node_metrics <- node_metrics %>%
  mutate(
    degree_rank = rank_desc(degree),
    mcc_like_rank = rank_desc(mcc_like),
    mnc_like_rank = rank_desc(mnc_like),
    closeness_rank = rank_desc(closeness),
    eigenvector_epc_surrogate_rank = rank_desc(eigenvector_epc_surrogate),
    degree_top_quartile = degree_rank <= top_n,
    mcc_like_top_quartile = mcc_like_rank <= top_n,
    mnc_like_top_quartile = mnc_like_rank <= top_n,
    closeness_top_quartile = closeness_rank <= top_n,
    eigenvector_epc_surrogate_top_quartile = eigenvector_epc_surrogate_rank <= top_n,
    n_algorithms_top_quartile = rowSums(across(ends_with("_top_quartile")), na.rm = TRUE),
    composite_rank_score = rowMeans(across(ends_with("_rank")), na.rm = TRUE),
    network_priority_group = case_when(
      n_algorithms_top_quartile >= 3 ~ "high_priority_3plus_algorithms",
      n_algorithms_top_quartile == 2 ~ "moderate_priority_2_algorithms",
      n_algorithms_top_quartile == 1 ~ "single_algorithm_priority",
      TRUE ~ "not_prioritized"
    ),
    boundary = "Network-prioritized candidate node; not used to redefine the locked program score."
  ) %>%
  arrange(desc(n_algorithms_top_quartile), composite_rank_score, desc(degree), gene_symbol)

write_tsv(node_metrics, file.path(result_dir, "network_prioritization_centrality_table.tsv"))

algorithm_overlap <- node_metrics %>%
  transmute(
    gene_symbol,
    Degree = as.integer(degree_top_quartile),
    MCC_like = as.integer(mcc_like_top_quartile),
    MNC_like = as.integer(mnc_like_top_quartile),
    Closeness = as.integer(closeness_top_quartile),
    EPC_surrogate = as.integer(eigenvector_epc_surrogate_top_quartile),
    n_algorithms_top_quartile,
    network_priority_group
  )

write_tsv(algorithm_overlap, file.path(result_dir, "network_prioritization_algorithm_overlap.tsv"))

network_summary <- tibble(
  n_locked_primary_genes = length(primary_genes),
  n_nodes_in_graph = vcount(g),
  n_edges_in_graph = ecount(g),
  string_required_score = 700L,
  top_quartile_n_per_algorithm = top_n,
  n_high_priority_3plus = sum(node_metrics$n_algorithms_top_quartile >= 3),
  n_moderate_priority_2 = sum(node_metrics$n_algorithms_top_quartile == 2),
  n_structural_ecm_high_priority = sum(node_metrics$n_algorithms_top_quartile >= 3 & node_metrics$is_ecm_structural_gene, na.rm = TRUE),
  n_ecm_excluded_high_priority = sum(node_metrics$n_algorithms_top_quartile >= 3 & node_metrics$include_in_ecm_excluded_score, na.rm = TRUE),
  boundary = "PPI/network centrality used for candidate prioritization only; it does not alter the locked score."
)
write_tsv(network_summary, file.path(result_dir, "network_prioritization_summary.tsv"))

overlap_plot_data <- algorithm_overlap %>%
  mutate(
    combination = pmap_chr(
      select(., Degree, MCC_like, MNC_like, Closeness, EPC_surrogate),
      function(Degree, MCC_like, MNC_like, Closeness, EPC_surrogate) {
        active <- c("Degree", "MCC-like", "MNC-like", "Closeness", "EPC surrogate")[
          c(Degree, MCC_like, MNC_like, Closeness, EPC_surrogate) == 1
        ]
        if (!length(active)) "No top-quartile metric" else paste(active, collapse = " + ")
      }
    )
  ) %>%
  filter(combination != "No top-quartile metric") %>%
  count(combination, n_algorithms_top_quartile, name = "n_genes", sort = TRUE) %>%
  mutate(
    combination_label = str_wrap(combination, width = 36),
    combination_label = factor(combination_label, levels = rev(unique(combination_label))),
    n_algorithms_top_quartile = factor(n_algorithms_top_quartile, levels = sort(unique(n_algorithms_top_quartile)))
  )

p_overlap <- ggplot(overlap_plot_data, aes(x = n_genes, y = combination_label, fill = n_algorithms_top_quartile)) +
  geom_col(width = 0.68, color = "#111111", linewidth = 0.18) +
  geom_text(aes(label = n_genes), hjust = -0.25, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  scale_fill_manual(values = c("1" = "#8DA0CB", "2" = "#66C2A5", "3" = "#FC8D62", "4" = "#E78AC3", "5" = "#A6D854")) +
  labs(
    title = "Overlap among top-quartile centrality metrics",
    subtitle = "Each bar shows the number of locked genes selected by the indicated metric combination",
    x = "Gene count",
    y = NULL,
    fill = "Metrics",
    caption = "Prespecified centrality metrics were used for candidate prioritization only."
  ) +
  theme_publication(8) +
  theme(panel.grid.major.y = element_blank())
save_plot(p_overlap, "network_prioritization_algorithm_upset", 6.6, 4.2)

lollipop_top_n <- min(30L, nrow(node_metrics))
lollipop_data <- node_metrics %>%
  slice_head(n = lollipop_top_n) %>%
  mutate(
    gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)),
    ecm_flag = if_else(is_ecm_structural_gene, "Structural ECM", "ECM-excluded/non-structural"),
    label_metric = paste0(n_algorithms_top_quartile, "/5")
  )

p_lollipop <- ggplot(lollipop_data, aes(x = n_algorithms_top_quartile, y = gene_symbol)) +
  geom_segment(aes(x = 0, xend = n_algorithms_top_quartile, yend = gene_symbol), color = "#777777", linewidth = 0.35) +
  geom_point(aes(size = degree, fill = ecm_flag), shape = 21, color = "#111111", stroke = 0.25) +
  geom_text(aes(label = label_metric), hjust = -0.35, size = 2.4) +
  scale_x_continuous(breaks = 0:5, limits = c(0, 5.65), expand = expansion(mult = c(0, 0.02))) +
  scale_size_continuous(range = c(1.8, 5.6)) +
  scale_fill_manual(values = c("ECM-excluded/non-structural" = "#0072B2", "Structural ECM" = "#D55E00")) +
  labs(
    title = "STRING network candidate prioritization",
    subtitle = "Genes ranked by top-quartile membership across five pre-specified centrality metrics",
    x = "Number of centrality metrics in top quartile",
    y = NULL,
    size = "STRING degree",
    fill = "Gene class",
    caption = "Candidate prioritization only; the locked score is not redefined by this network analysis."
  ) +
  theme_publication(8)
save_plot(p_lollipop, "network_prioritization_lollipop", 6.9, 5.4)

plot_node_candidates <- node_metrics %>%
  filter(n_algorithms_top_quartile >= 2 | degree_rank <= 25) %>%
  arrange(desc(n_algorithms_top_quartile), composite_rank_score, desc(degree), gene_symbol)
plot_nodes <- plot_node_candidates %>%
  slice_head(n = min(45L, nrow(plot_node_candidates))) %>%
  pull(gene_symbol)

g_plot <- induced_subgraph(g, vids = V(g)[name %in% plot_nodes])
V(g_plot)$n_algorithms_top_quartile <- node_metrics$n_algorithms_top_quartile[match(V(g_plot)$name, node_metrics$gene_symbol)]
V(g_plot)$is_ecm_structural_gene <- node_metrics$is_ecm_structural_gene[match(V(g_plot)$name, node_metrics$gene_symbol)]
V(g_plot)$network_priority_group <- node_metrics$network_priority_group[match(V(g_plot)$name, node_metrics$gene_symbol)]

p_network <- ggraph(g_plot, layout = "fr") +
  geom_edge_link(aes(width = string_score), alpha = 0.28, color = "#777777") +
  geom_node_point(
    aes(size = n_algorithms_top_quartile, fill = is_ecm_structural_gene),
    shape = 21,
    color = "#111111",
    stroke = 0.25
  ) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.4, max.overlaps = Inf) +
  scale_edge_width(range = c(0.15, 0.8), guide = "none") +
  scale_size_continuous(range = c(2.2, 6.4), breaks = 0:5) +
  scale_fill_manual(
    values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
    labels = c(`FALSE` = "ECM-excluded/non-structural", `TRUE` = "Structural ECM"),
    na.value = "#999999"
  ) +
  labs(
    title = "High-confidence STRING subnetwork of prioritized candidates",
    subtitle = "Displayed nodes are high/moderate-priority candidates plus top-degree context nodes",
    size = "Top-quartile metrics",
    fill = "Gene class",
    caption = "Edges represent STRING functional associations and are not interpreted as direct physical or causal interactions."
  ) +
  theme_void(base_family = "Arial", base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "#333333"),
    plot.caption = element_text(size = 7, color = "#555555", hjust = 0),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )
save_plot(p_network, "network_prioritization_string_network", 7.4, 5.4)

top_candidates <- node_metrics %>%
  filter(n_algorithms_top_quartile >= 3) %>%
  arrange(composite_rank_score, desc(degree), gene_symbol)

top_lines <- if (nrow(top_candidates)) {
  paste0(
    "- `", top_candidates$gene_symbol, "`: ",
    top_candidates$n_algorithms_top_quartile, "/5 centrality metrics; degree ",
    top_candidates$degree, "; class ",
    ifelse(top_candidates$is_ecm_structural_gene, "structural ECM", "ECM-excluded/non-structural"),
    ".\n",
    collapse = ""
  )
} else {
  "- No gene met the >=3 of 5 top-quartile prioritization rule.\n"
}

cat(
  "# STRING Network Candidate Prioritization Summary\n\n",
  "Generated by `scripts/10_network_prioritization.R` on 2026-05-23.\n\n",
  "## Analysis Scope\n\n",
  "This analysis uses STRING functional associations and prespecified centrality metrics as a candidate-prioritization layer. The locked primary score was not changed. STRING edges are functional association edges and are not interpreted as direct physical binding or causal regulation.\n\n",
  "## Network Parameters\n\n",
  "- Species: Homo sapiens, NCBI taxonomy 9606.\n",
  "- STRING required score: 700.\n",
  "- Input universe: locked primary program genes only.\n",
  "- Prioritization rule: top quartile within each of five pre-specified centrality metrics; high-priority candidates require top-quartile status in at least 3 of 5 metrics.\n",
  "- Metrics: degree, MCC-like maximal-clique centrality, MNC-like maximum-neighborhood-component size, closeness, and eigenvector centrality as an EPC surrogate.\n\n",
  "## Summary\n\n",
  "- Locked primary genes queried: ", network_summary$n_locked_primary_genes, ".\n",
  "- Nodes in STRING graph: ", network_summary$n_nodes_in_graph, ".\n",
  "- Edges in filtered graph: ", network_summary$n_edges_in_graph, ".\n",
  "- High-priority candidates by >=3 algorithms: ", network_summary$n_high_priority_3plus, ".\n",
  "- Moderate-priority candidates by 2 algorithms: ", network_summary$n_moderate_priority_2, ".\n\n",
  "## High-Priority Candidate Nodes\n\n",
  top_lines,
  "\n## Outputs\n\n",
  "- `results/network/network_prioritization_string_query_status.tsv`\n",
  "- `results/network/network_prioritization_string_raw_edges.tsv`\n",
  "- `results/network/network_prioritization_string_filtered_edges.tsv`\n",
  "- `results/network/network_prioritization_centrality_table.tsv`\n",
  "- `results/network/network_prioritization_algorithm_overlap.tsv`\n",
  "- `results/network/network_prioritization_summary.tsv`\n",
  "- `figures/network/network_prioritization_algorithm_upset.png` and `.pdf`\n",
  "- `figures/network/network_prioritization_lollipop.png` and `.pdf`\n",
  "- `figures/network/network_prioritization_string_network.png` and `.pdf`\n",
  file = file.path(doc_dir, "network_prioritization_summary.md"),
  sep = ""
)

message("Done. STRING network candidate-prioritization outputs written to ", result_dir)

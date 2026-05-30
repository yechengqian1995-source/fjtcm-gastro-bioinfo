#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1]]
}
positional_args <- args[!str_starts(args, "--")]
project_dir <- if (length(positional_args) >= 1) positional_args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

expr_dir <- file.path(project_dir, "data", "expression_processed")
sig_dir <- file.path(project_dir, "data", "signatures")
result_dir <- file.path(project_dir, "results", "bulk")
doc_dir <- file.path(project_dir, "docs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

random_seed <- as.integer(get_opt("--seed", Sys.getenv("RANDOM_SET_SEED", "20260523")))
n_random_sets <- as.integer(get_opt("--n-random", Sys.getenv("N_RANDOM_SETS", "10000")))
if (is.na(n_random_sets) || n_random_sets < 1L) stop("n_random_sets must be a positive integer.")
set.seed(random_seed)

signature <- read_tsv(file.path(sig_dir, "signature_v1_locked.tsv"), show_col_types = FALSE) %>%
  mutate(
    gene_symbol = toupper(gene_symbol),
    include_in_primary_score = as.logical(include_in_primary_score),
    include_in_ecm_excluded_score = as.logical(include_in_ecm_excluded_score),
    is_ecm_structural_gene = as.logical(is_ecm_structural_gene)
  )

gene_sets <- list(
  primary = signature %>%
    filter(include_in_primary_score) %>%
    pull(gene_symbol),
  ecm_excluded = signature %>%
    filter(include_in_ecm_excluded_score) %>%
    pull(gene_symbol),
  structural_ecm = signature %>%
    filter(is_ecm_structural_gene) %>%
    pull(gene_symbol),
  hsc_context = signature %>%
    filter(str_detect(source_category, "HSC|portal_fibroblast")) %>%
    pull(gene_symbol)
)
gene_sets <- lapply(gene_sets, unique)

rank_percentile_matrix <- function(log_expr) {
  ranks <- apply(log_expr, 2, rank, ties.method = "average", na.last = "keep")
  ranks <- as.matrix(ranks)
  rownames(ranks) <- rownames(log_expr)
  colnames(ranks) <- colnames(log_expr)
  (ranks - 1) / pmax(nrow(log_expr) - 1, 1)
}

score_from_pct <- function(rank_pct, genes) {
  present <- intersect(unique(toupper(genes)), rownames(rank_pct))
  if (!length(present)) {
    return(rep(NA_real_, ncol(rank_pct)))
  }
  colMeans(rank_pct[present, , drop = FALSE], na.rm = TRUE)
}

score_from_ranks <- function(log_expr, genes) {
  rank_pct <- rank_percentile_matrix(log_expr)
  score_from_pct(rank_pct, genes)
}

zscore <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

standardize_sex <- function(x) {
  raw <- str_to_lower(str_trim(as.character(x)))
  case_when(
    raw %in% c("m", "male") ~ "male",
    raw %in% c("f", "female") ~ "female",
    TRUE ~ NA_character_
  )
}

numeric_or_na <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_logistic <- function(dat, score_col, endpoint_col = "advanced_binary", covariates = character()) {
  use_cols <- c(endpoint_col, score_col, covariates)
  model_dat <- dat %>%
    select(all_of(use_cols)) %>%
    drop_na()

  if (nrow(model_dat) < 10 || length(unique(model_dat[[endpoint_col]])) < 2) {
    return(NULL)
  }
  if (any(table(model_dat[[endpoint_col]]) < 3)) {
    return(NULL)
  }

  rhs <- paste(c(score_col, covariates), collapse = " + ")
  frm <- as.formula(paste(endpoint_col, "~", rhs))
  fit <- tryCatch(
    suppressWarnings(glm(frm, data = model_dat, family = binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  coefs <- summary(fit)$coefficients
  if (!score_col %in% rownames(coefs)) return(NULL)
  beta <- unname(coefs[score_col, "Estimate"])
  se <- unname(coefs[score_col, "Std. Error"])
  p <- unname(coefs[score_col, "Pr(>|z|)"])
  tibble(
    n_model = nrow(model_dat),
    n_advanced = sum(model_dat[[endpoint_col]] == 1),
    n_non_advanced = sum(model_dat[[endpoint_col]] == 0),
    effect_scale = "OR_per_1SD_score",
    estimate = exp(beta),
    conf_low = exp(beta - 1.96 * se),
    conf_high = exp(beta + 1.96 * se),
    beta = beta,
    beta_se = se,
    p_value = p,
    convergence_note = if_else(isTRUE(fit$converged), "glm_converged", "glm_not_converged")
  )
}

safe_spearman <- function(dat, score_col, endpoint_col) {
  model_dat <- dat %>%
    select(all_of(c(endpoint_col, score_col))) %>%
    drop_na()
  if (nrow(model_dat) < 10 || length(unique(model_dat[[endpoint_col]])) < 3) {
    return(NULL)
  }
  ct <- tryCatch(
    suppressWarnings(cor.test(model_dat[[score_col]], model_dat[[endpoint_col]], method = "spearman", exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(ct)) return(NULL)
  tibble(
    n_model = nrow(model_dat),
    n_advanced = NA_integer_,
    n_non_advanced = NA_integer_,
    effect_scale = paste0("spearman_rho_with_", endpoint_col),
    estimate = unname(ct$estimate),
    conf_low = NA_real_,
    conf_high = NA_real_,
    beta = NA_real_,
    beta_se = NA_real_,
    p_value = ct$p.value,
    convergence_note = "nonparametric_spearman"
  )
}

make_coverage <- function(acc, log_expr) {
  imap_dfr(gene_sets, function(genes, score_name) {
    present <- intersect(genes, rownames(log_expr))
    tibble(
      accession = acc,
      score_name = score_name,
      n_signature_genes = length(unique(genes)),
      n_present = length(present),
      coverage_fraction = if_else(length(unique(genes)) > 0, length(present) / length(unique(genes)), NA_real_),
      present_genes = paste(sort(present), collapse = ";"),
      missing_genes = paste(sort(setdiff(unique(genes), present)), collapse = ";")
    )
  })
}

make_expression_bins <- function(log_expr, available_genes) {
  gene_mean <- rowMeans(log_expr, na.rm = TRUE)
  tibble(
    gene = rownames(log_expr),
    mean_expr = gene_mean
  ) %>%
    filter(gene %in% available_genes, is.finite(mean_expr)) %>%
    mutate(expr_bin = ntile(mean_expr, 10))
}

sample_matched_genes <- function(sig_present, bins, universe, n_sets) {
  sig_bins <- bins %>%
    filter(gene %in% sig_present) %>%
    count(expr_bin, name = "n")
  candidates <- bins %>%
    filter(gene %in% universe, !gene %in% sig_present)

  map(seq_len(n_sets), function(iter) {
    sampled <- map2(sig_bins$expr_bin, sig_bins$n, function(bin, n) {
      pool <- candidates %>% filter(expr_bin == bin) %>% pull(gene)
      if (length(pool) < n) {
        pool <- candidates$gene
      }
      sample(pool, size = n, replace = length(pool) < n)
    }) %>%
      unlist(use.names = FALSE)
    unique(sampled)
  })
}

run_random_controls <- function(acc, log_expr, rank_pct, pheno, observed_beta) {
  pheno_model <- pheno %>%
    mutate(advanced_binary = case_when(
      advanced_fibrosis == "advanced" ~ 1,
      advanced_fibrosis == "non_advanced" ~ 0,
      TRUE ~ NA_real_
    ))
  if (is.na(observed_beta) ||
      sum(!is.na(pheno_model$advanced_binary)) < 10 ||
      length(unique(na.omit(pheno_model$advanced_binary))) < 2) {
    return(tibble())
  }

  sig_present <- intersect(gene_sets$primary, rownames(log_expr))
  if (length(sig_present) < 10) return(tibble())
  bins <- make_expression_bins(log_expr, rownames(log_expr))
  random_sets <- sample_matched_genes(sig_present, bins, rownames(log_expr), n_random_sets)

  map_dfr(seq_along(random_sets), function(i) {
    random_score <- score_from_pct(rank_pct, random_sets[[i]])
    tmp <- pheno_model %>%
      mutate(random_score_z = zscore(random_score))
    res <- safe_logistic(tmp, "random_score_z")
    if (is.null(res)) return(tibble())
    res %>%
      transmute(
        accession = acc,
        random_set_id = i,
        n_genes = length(random_sets[[i]]),
        beta = beta,
        estimate_or = estimate,
        p_value = p_value,
        observed_beta = observed_beta,
        abs_beta_ge_observed = abs(beta) >= abs(observed_beta)
      )
  })
}

expr_files <- list.files(expr_dir, pattern = "_bulk_expression\\.rds$", full.names = TRUE)
if (!length(expr_files)) stop("No processed bulk expression RDS files found.")

score_tables <- list()
coverage_tables <- list()
model_tables <- list()
random_tables <- list()
endpoint_tables <- list()

for (file in expr_files) {
  obj <- readRDS(file)
  acc <- obj$accession
  log_expr <- as.matrix(obj$logcpm)
  rownames(log_expr) <- toupper(rownames(log_expr))
  rank_pct <- rank_percentile_matrix(log_expr)
  pheno <- obj$pheno %>%
    mutate(
      sample_id = as.character(sample_id),
      age_numeric = numeric_or_na(age_raw),
      sex_binary = standardize_sex(sex_raw),
      advanced_binary = case_when(
        advanced_fibrosis == "advanced" ~ 1,
        advanced_fibrosis == "non_advanced" ~ 0,
        TRUE ~ NA_real_
      )
    )

  scores <- imap_dfc(gene_sets, function(genes, nm) {
    tibble(!!paste0(nm, "_score") := score_from_pct(rank_pct, genes))
  })
  scores <- scores %>%
    mutate(across(everything(), zscore, .names = "{.col}_z"))

  score_dat <- bind_cols(pheno, scores) %>%
    mutate(accession = acc)
  score_tables[[acc]] <- score_dat
  coverage_tables[[acc]] <- make_coverage(acc, log_expr)

  endpoint_tables[[acc]] <- score_dat %>%
    summarise(
      accession = acc,
      n_samples = n(),
      n_with_fibrosis = sum(!is.na(fibrosis_stage)),
      n_advanced = sum(advanced_fibrosis == "advanced", na.rm = TRUE),
      n_non_advanced = sum(advanced_fibrosis == "non_advanced", na.rm = TRUE),
      n_with_nas = sum(!is.na(nas_score)),
      n_with_age = sum(!is.na(age_numeric)),
      n_with_sex = sum(!is.na(sex_binary)),
      disease_labels = paste(sort(unique(na.omit(c(disease_label_raw, disease_subtype_raw)))), collapse = ";")
    )

  score_cols <- names(scores)[str_detect(names(scores), "_score_z$")]
  for (score_col in score_cols) {
    score_name <- str_remove(score_col, "_score_z$")

    base <- safe_logistic(score_dat, score_col)
    if (!is.null(base)) {
      model_tables[[paste(acc, score_name, "advanced_base", sep = "_")]] <- base %>%
        mutate(
          accession = acc,
          score_name = score_name,
          endpoint = "advanced_fibrosis_F3_F4_vs_F0_F2",
          model_type = "logistic_base"
        )
    }

    covs <- character()
    if (sum(!is.na(score_dat$age_numeric)) >= 20) covs <- c(covs, "age_numeric")
    if (n_distinct(na.omit(score_dat$sex_binary)) == 2 && sum(!is.na(score_dat$sex_binary)) >= 20) {
      score_dat <- score_dat %>% mutate(sex_binary = factor(sex_binary))
      covs <- c(covs, "sex_binary")
    }
    if (length(covs)) {
      adj <- safe_logistic(score_dat, score_col, covariates = covs)
      if (!is.null(adj)) {
        model_tables[[paste(acc, score_name, "advanced_age_sex", sep = "_")]] <- adj %>%
          mutate(
            accession = acc,
            score_name = score_name,
            endpoint = "advanced_fibrosis_F3_F4_vs_F0_F2",
            model_type = paste0("logistic_adjusted_", paste(covs, collapse = "_"))
          )
      }
    }

    stage_res <- safe_spearman(score_dat, score_col, "fibrosis_stage")
    if (!is.null(stage_res)) {
      model_tables[[paste(acc, score_name, "fibrosis_spearman", sep = "_")]] <- stage_res %>%
        mutate(
          accession = acc,
          score_name = score_name,
          endpoint = "ordinal_fibrosis_stage",
          model_type = "spearman"
        )
    }

    nas_res <- safe_spearman(score_dat, score_col, "nas_score")
    if (!is.null(nas_res)) {
      model_tables[[paste(acc, score_name, "nas_spearman", sep = "_")]] <- nas_res %>%
        mutate(
          accession = acc,
          score_name = score_name,
          endpoint = "NAS",
          model_type = "spearman"
        )
    }
  }

  primary_base_beta <- model_tables[[paste(acc, "primary", "advanced_base", sep = "_")]]
  observed_beta <- if (!is.null(primary_base_beta)) primary_base_beta$beta[[1]] else NA_real_
  random_tables[[acc]] <- run_random_controls(acc, log_expr, rank_pct, score_dat, observed_beta)
}

program_scores <- bind_rows(score_tables)
coverage <- bind_rows(coverage_tables)
endpoint_summary <- bind_rows(endpoint_tables)
model_results <- bind_rows(model_tables) %>%
  select(accession, score_name, endpoint, model_type, n_model, n_advanced, n_non_advanced,
         effect_scale, estimate, conf_low, conf_high, beta, beta_se, p_value, convergence_note) %>%
  group_by(endpoint, model_type) %>%
  mutate(p_fdr_within_endpoint_model = p.adjust(p_value, method = "BH")) %>%
  ungroup()
random_controls <- bind_rows(random_tables)

random_summary <- random_controls %>%
  group_by(accession) %>%
  summarise(
    n_random_sets = n(),
    observed_beta = first(observed_beta),
    n_abs_beta_ge_observed = sum(abs_beta_ge_observed, na.rm = TRUE),
    empirical_two_sided_p = mean(abs_beta_ge_observed, na.rm = TRUE),
    empirical_two_sided_p_plus1 = (sum(abs_beta_ge_observed, na.rm = TRUE) + 1) / (n() + 1),
    empirical_p_resolution = 1 / n(),
    random_beta_median = median(beta, na.rm = TRUE),
    random_beta_q025 = quantile(beta, 0.025, na.rm = TRUE),
    random_beta_q975 = quantile(beta, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(program_scores, file.path(result_dir, "bulk_program_scores.tsv"))
write_tsv(coverage, file.path(result_dir, "bulk_signature_gene_coverage.tsv"))
write_tsv(endpoint_summary, file.path(result_dir, "bulk_endpoint_summary.tsv"))
write_tsv(model_results, file.path(result_dir, "bulk_model_results.tsv"))
write_tsv(random_controls, file.path(result_dir, "bulk_random_gene_set_controls.tsv"))
write_tsv(random_summary, file.path(result_dir, "bulk_random_gene_set_summary.tsv"))

advanced_model_text <- model_results %>%
  filter(endpoint == "advanced_fibrosis_F3_F4_vs_F0_F2", model_type == "logistic_base", score_name %in% c("primary", "ecm_excluded", "structural_ecm")) %>%
  arrange(accession, score_name) %>%
  mutate(
    line = paste0(
      "- `", accession, "` ", score_name, ": n=", n_model,
      ", OR per 1 SD=", signif(estimate, 3),
      " (95% CI ", signif(conf_low, 3), "-", signif(conf_high, 3),
      "), P=", signif(p_value, 3), "."
    )
  ) %>%
  pull(line) %>%
  paste(collapse = "\n")

if (!nzchar(advanced_model_text)) {
  advanced_model_text <- "- No advanced-fibrosis model passed minimum endpoint and event-count checks."
}

coverage_text <- coverage %>%
  filter(score_name %in% c("primary", "ecm_excluded", "structural_ecm")) %>%
  mutate(line = paste0(
    "- `", accession, "` ", score_name, ": ", n_present, "/", n_signature_genes,
    " genes present (", scales::percent(coverage_fraction, accuracy = 0.1), ")."
  )) %>%
  pull(line) %>%
  paste(collapse = "\n")

cat(
  "# Bulk Phase 2 Results Summary\n\n",
  "Generated by `scripts/05_bulk_score_validate.R` on ", as.character(Sys.Date()), ".\n\n",
  "Random-control seed: ", random_seed, ". Random gene sets per fibrosis-mapped cohort: ", n_random_sets, ".\n\n",
  "## Analysis Boundary\n\n",
  "Scores quantify a pre-specified mechanical-stress-associated transcriptional memory-like program. These retrospective public-data results support association-level interpretation only. They do not prove mechanical memory, HSC causality, diagnosis, prognosis, or therapeutic target validity.\n\n",
  "## Output Files\n\n",
  "- `results/bulk/bulk_program_scores.tsv`\n",
  "- `results/bulk/bulk_signature_gene_coverage.tsv`\n",
  "- `results/bulk/bulk_endpoint_summary.tsv`\n",
  "- `results/bulk/bulk_model_results.tsv`\n",
  "- `results/bulk/bulk_random_gene_set_controls.tsv`\n",
  "- `results/bulk/bulk_random_gene_set_summary.tsv`\n\n",
  "## Signature Coverage\n\n",
  coverage_text,
  "\n\n## Primary Advanced-Fibrosis Association Models\n\n",
  advanced_model_text,
  "\n\n## Interpretation Rules\n\n",
  "- Full-score signals must be compared with ECM-excluded and structural ECM-only scores.\n",
  "- If ECM-excluded associations weaken or disappear, interpretation should be downgraded to matrix or fibrotic-stromal abundance-associated signal.\n",
  "- Cohorts without mapped fibrosis stage are secondary disease-state context only.\n",
  "- Cross-cohort claims remain conditional on sample-overlap auditing and endpoint harmonization.\n",
  file = file.path(doc_dir, "bulk_phase2_results_summary.md"),
  sep = ""
)

message("Done. Bulk score and validation outputs written to ", result_dir)

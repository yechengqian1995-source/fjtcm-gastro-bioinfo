#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1 && !grepl("^--", args[[1]])) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

expr_dir <- file.path(project_dir, "data", "expression_processed")
sig_dir <- file.path(project_dir, "data", "signatures")
result_dir <- file.path(project_dir, "results", "benchmark_robustness_sensitivity")
figure_dir <- file.path(project_dir, "figures", "benchmark_robustness_sensitivity")
supp_dir <- file.path(project_dir, "manuscript_outputs", "supplementary_tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260524)
n_boot <- 2000L
fibrosis_cohorts <- c("GSE135251", "GSE130970", "GSE162694")
small_cohorts <- c("GSE130970", "GSE162694")

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "")
}

write_tsv_base <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

to_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

standardize_sex <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("m", "male")] <- "male"
  out[raw %in% c("f", "female")] <- "female"
  out
}

rank_percentile_matrix <- function(log_expr) {
  ranks <- apply(log_expr, 2, rank, ties.method = "average", na.last = "keep")
  ranks <- as.matrix(ranks)
  rownames(ranks) <- rownames(log_expr)
  colnames(ranks) <- colnames(log_expr)
  (ranks - 1) / pmax(nrow(log_expr) - 1, 1)
}

score_from_pct <- function(rank_pct, genes) {
  present <- intersect(unique(toupper(genes)), rownames(rank_pct))
  if (!length(present)) return(rep(NA_real_, ncol(rank_pct)))
  colMeans(rank_pct[present, , drop = FALSE], na.rm = TRUE)
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 10 || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    return(c(estimate = NA_real_, p_value = NA_real_, n = sum(keep)))
  }
  ct <- tryCatch(
    suppressWarnings(stats::cor.test(x[keep], y[keep], method = method, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(ct)) return(c(estimate = NA_real_, p_value = NA_real_, n = sum(keep)))
  c(estimate = unname(ct$estimate), p_value = ct$p.value, n = sum(keep))
}

fit_logistic_target <- function(dat, target_col, covariates = character(), endpoint_col = "advanced_binary") {
  cols <- unique(c(endpoint_col, target_col, covariates))
  model_dat <- dat[, cols, drop = FALSE]
  model_dat <- model_dat[stats::complete.cases(model_dat), , drop = FALSE]
  if (nrow(model_dat) < 10) return(NULL)
  if (length(unique(model_dat[[endpoint_col]])) < 2) return(NULL)
  if (any(table(model_dat[[endpoint_col]]) < 3)) return(NULL)
  for (nm in covariates) {
    if (is.character(model_dat[[nm]])) model_dat[[nm]] <- factor(model_dat[[nm]])
    if (is.factor(model_dat[[nm]]) && length(unique(model_dat[[nm]])) < 2) return(NULL)
    if (is.numeric(model_dat[[nm]]) && stats::sd(model_dat[[nm]], na.rm = TRUE) == 0) return(NULL)
  }
  rhs <- paste(c(target_col, covariates), collapse = " + ")
  frm <- stats::as.formula(paste(endpoint_col, "~", rhs))
  fit <- tryCatch(
    suppressWarnings(stats::glm(frm, data = model_dat, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  coefs <- summary(fit)$coefficients
  if (!target_col %in% rownames(coefs)) return(NULL)
  beta <- unname(coefs[target_col, "Estimate"])
  se <- unname(coefs[target_col, "Std. Error"])
  p <- unname(coefs[target_col, "Pr(>|z|)"])
  if (!is.finite(beta) || !is.finite(se)) return(NULL)

  vif <- NA_real_
  pred_cols <- c(target_col, covariates)
  numeric_pred <- pred_cols[vapply(model_dat[, pred_cols, drop = FALSE], is.numeric, logical(1))]
  if (length(pred_cols) >= 2 && target_col %in% numeric_pred) {
    mm <- stats::model.matrix(stats::as.formula(paste("~", paste(covariates, collapse = " + "))), data = model_dat)
    if (ncol(mm) > 1) {
      aux <- tryCatch(summary(stats::lm(model_dat[[target_col]] ~ mm[, -1, drop = FALSE]))$r.squared, error = function(e) NA_real_)
      if (is.finite(aux) && aux < 1) vif <- 1 / (1 - aux)
    }
  }

  data.frame(
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
    target_vif = vif,
    convergence_note = ifelse(isTRUE(fit$converged), "glm_converged", "glm_not_converged"),
    stringsAsFactors = FALSE
  )
}

bootstrap_score <- function(dat, score_col, accession, score_name, n_iter = n_boot) {
  model_dat <- dat[, c("advanced_binary", score_col), drop = FALSE]
  model_dat <- model_dat[stats::complete.cases(model_dat), , drop = FALSE]
  cases <- which(model_dat$advanced_binary == 1)
  controls <- which(model_dat$advanced_binary == 0)
  if (length(cases) < 3 || length(controls) < 3) return(data.frame())
  rows <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    idx <- c(sample(cases, length(cases), replace = TRUE), sample(controls, length(controls), replace = TRUE))
    tmp <- model_dat[idx, , drop = FALSE]
    fit <- tryCatch(
      suppressWarnings(stats::glm(stats::as.formula(paste("advanced_binary ~", score_col)), data = tmp, family = stats::binomial())),
      error = function(e) NULL
    )
    if (!is.null(fit) && score_col %in% names(stats::coef(fit)) && is.finite(stats::coef(fit)[[score_col]])) {
      beta <- unname(stats::coef(fit)[[score_col]])
      rows[[i]] <- data.frame(
        accession = accession,
        score_name = score_name,
        bootstrap_id = i,
        beta = beta,
        estimate_or = exp(beta),
        converged = isTRUE(fit$converged),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

leave_one_out_score <- function(dat, score_col, accession, score_name) {
  model_dat <- dat[, c("sample_id", "advanced_binary", score_col), drop = FALSE]
  model_dat <- model_dat[stats::complete.cases(model_dat), , drop = FALSE]
  if (nrow(model_dat) < 10 || any(table(model_dat$advanced_binary) < 3)) return(data.frame())
  rows <- vector("list", nrow(model_dat))
  for (i in seq_len(nrow(model_dat))) {
    tmp <- model_dat[-i, , drop = FALSE]
    res <- fit_logistic_target(tmp, score_col)
    if (!is.null(res)) {
      rows[[i]] <- data.frame(
        accession = accession,
        score_name = score_name,
        removed_sample_id = model_dat$sample_id[[i]],
        removed_endpoint = model_dat$advanced_binary[[i]],
        res,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

signature <- read_tsv_base(file.path(sig_dir, "signature_v1_locked.tsv"))
signature$gene_symbol <- toupper(signature$gene_symbol)

locked_sets <- list(
  primary = unique(signature$gene_symbol[signature$include_in_primary_score == TRUE]),
  ecm_excluded = unique(signature$gene_symbol[signature$include_in_ecm_excluded_score == TRUE]),
  structural_ecm = unique(signature$gene_symbol[signature$is_ecm_structural_gene == TRUE]),
  hsc_context = unique(signature$gene_symbol[grepl("HSC|portal_fibroblast", signature$source_category)])
)

benchmark_sets <- list(
  core_fibrosis_stromal_benchmark = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "DCN", "LUM", "FN1", "POSTN", "SPARC", "THBS1", "THBS2", "MMP2", "MMP14",
    "TIMP1", "TIMP2", "LOX", "LOXL2", "ACTA2", "TGFB1", "SERPINE1", "CTGF",
    "PDGFRB", "VIM", "VCAN", "ITGB1", "ITGAV"
  ),
  hsc_myofibroblast_benchmark = c(
    "ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "COL3A1", "PDGFRB", "RGS5",
    "CSPG4", "DES", "DCN", "LUM", "THY1", "COL6A3", "SPARC", "MMP2", "TIMP1",
    "PDGFRA", "FAP"
  ),
  inflammatory_activity_benchmark = c(
    "PTPRC", "LST1", "AIF1", "CD68", "CD14", "ITGAM", "FCGR3A", "CSF1R", "TYROBP",
    "NFKB1", "RELA", "IL1B", "IL6", "TNF", "CCL2", "CXCL8", "CXCL10", "CCL5", "NLRP3"
  ),
  myeloid_immune_proxy = c(
    "PTPRC", "LST1", "AIF1", "CD68", "CD14", "ITGAM", "FCGR3A", "CSF1R", "TYROBP",
    "MS4A7", "LILRB4", "CTSS", "SPI1"
  ),
  hepatocyte_identity_proxy = c(
    "ALB", "APOA1", "APOB", "TTR", "TF", "FGA", "FGB", "APOC3", "CYP2E1", "CYP3A4",
    "PCK1", "HNF4A", "KRT18", "SERPINA1", "HPD", "ARG1"
  )
)
benchmark_sets <- lapply(benchmark_sets, function(x) unique(toupper(x)))
all_sets <- c(locked_sets, benchmark_sets)

gene_set_rows <- do.call(rbind, lapply(names(all_sets), function(nm) {
  role <- ifelse(nm %in% names(locked_sets), "locked_program_or_internal_sensitivity", "transparent_benchmark_or_proxy")
  note <- switch(
    nm,
    core_fibrosis_stromal_benchmark = "Core fibrosis and stromal-remodeling comparator compiled from canonical fibrogenic ECM, remodeling, TGF-beta, integrin and myofibroblast markers; not claimed as a single published signature.",
    hsc_myofibroblast_benchmark = "HSC and myofibroblast marker-context comparator; used as stromal abundance sensitivity, not as direct deconvolution.",
    inflammatory_activity_benchmark = "Inflammatory and myeloid-activation marker proxy used to evaluate inflammation-related bulk confounding.",
    myeloid_immune_proxy = "Myeloid and pan-immune marker proxy used to evaluate immune-cell composition sensitivity.",
    hepatocyte_identity_proxy = "Hepatocyte identity proxy used as an alternative tissue-composition control.",
    "Locked program component from the predefined 97-gene table."
  )
  data.frame(
    gene_set = nm,
    gene_symbol = all_sets[[nm]],
    set_size = length(all_sets[[nm]]),
    role = role,
    source_note = note,
    stringsAsFactors = FALSE
  )
}))
write_tsv_base(gene_set_rows, file.path(result_dir, "fibrosis_stromal_benchmark_gene_sets.tsv"))
write_tsv_base(gene_set_rows, file.path(supp_dir, "Supplementary_Table_2_benchmark_and_proxy_gene_sets.tsv"))

score_tables <- list()
coverage_tables <- list()
benchmark_models <- list()
confound_models <- list()
correlation_tables <- list()
bootstrap_tables <- list()
bootstrap_summary <- list()
loo_tables <- list()
loo_summary <- list()

expr_files <- list.files(expr_dir, pattern = "_bulk_expression\\.rds$", full.names = TRUE)
for (file in expr_files) {
  obj <- readRDS(file)
  acc <- as.character(obj$accession)
  if (!acc %in% fibrosis_cohorts) next

  log_expr <- as.matrix(obj$logcpm)
  rownames(log_expr) <- toupper(rownames(log_expr))
  rank_pct <- rank_percentile_matrix(log_expr)

  pheno <- obj$pheno
  pheno$sample_id <- as.character(pheno$sample_id)
  if (all(colnames(log_expr) %in% pheno$sample_id)) {
    pheno <- pheno[match(colnames(log_expr), pheno$sample_id), , drop = FALSE]
  }
  pheno$advanced_binary <- ifelse(pheno$advanced_fibrosis == "advanced", 1,
                                  ifelse(pheno$advanced_fibrosis == "non_advanced", 0, NA))
  pheno$fibrosis_stage <- to_numeric(pheno$fibrosis_stage)
  pheno$nas_score <- to_numeric(pheno$nas_score)
  pheno$age_numeric <- to_numeric(pheno$age_raw)
  pheno$sex_binary <- standardize_sex(pheno$sex_raw)

  score_df <- pheno
  for (nm in names(all_sets)) {
    raw_score <- score_from_pct(rank_pct, all_sets[[nm]])
    score_df[[paste0(nm, "_score")]] <- raw_score
    score_df[[paste0(nm, "_score_z")]] <- zscore(raw_score)
    present <- intersect(all_sets[[nm]], rownames(log_expr))
    coverage_tables[[paste(acc, nm, sep = "__")]] <- data.frame(
      accession = acc,
      gene_set = nm,
      n_genes = length(all_sets[[nm]]),
      n_present = length(present),
      coverage_fraction = ifelse(length(all_sets[[nm]]) > 0, length(present) / length(all_sets[[nm]]), NA_real_),
      present_genes = paste(sort(present), collapse = ";"),
      missing_genes = paste(sort(setdiff(all_sets[[nm]], present)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  score_tables[[acc]] <- score_df

  benchmark_score_names <- c(
    "primary", "ecm_excluded", "structural_ecm", "hsc_context",
    "core_fibrosis_stromal_benchmark", "hsc_myofibroblast_benchmark",
    "inflammatory_activity_benchmark", "myeloid_immune_proxy", "hepatocyte_identity_proxy"
  )

  for (score_name in benchmark_score_names) {
    score_col <- paste0(score_name, "_score_z")
    if (!score_col %in% names(score_df)) next
    base <- fit_logistic_target(score_df, score_col)
    if (!is.null(base)) {
      benchmark_models[[paste(acc, score_name, "adv", sep = "__")]] <- data.frame(
        accession = acc,
        score_name = score_name,
        endpoint = "advanced_fibrosis_F3_F4_vs_F0_F2",
        model_type = "logistic_base",
        base,
        stringsAsFactors = FALSE
      )
    }
    st <- safe_cor(score_df[[score_col]], score_df$fibrosis_stage, method = "spearman")
    if (is.finite(st[["estimate"]])) {
      benchmark_models[[paste(acc, score_name, "stage", sep = "__")]] <- data.frame(
        accession = acc,
        score_name = score_name,
        endpoint = "ordinal_fibrosis_stage",
        model_type = "spearman",
        n_model = as.integer(st[["n"]]),
        n_advanced = NA_integer_,
        n_non_advanced = NA_integer_,
        effect_scale = "spearman_rho_with_fibrosis_stage",
        estimate = st[["estimate"]],
        conf_low = NA_real_,
        conf_high = NA_real_,
        beta = NA_real_,
        beta_se = NA_real_,
        p_value = st[["p_value"]],
        target_vif = NA_real_,
        convergence_note = "nonparametric_spearman",
        stringsAsFactors = FALSE
      )
    }
    ns <- safe_cor(score_df[[score_col]], score_df$nas_score, method = "spearman")
    if (is.finite(ns[["estimate"]])) {
      benchmark_models[[paste(acc, score_name, "nas", sep = "__")]] <- data.frame(
        accession = acc,
        score_name = score_name,
        endpoint = "NAS",
        model_type = "spearman",
        n_model = as.integer(ns[["n"]]),
        n_advanced = NA_integer_,
        n_non_advanced = NA_integer_,
        effect_scale = "spearman_rho_with_NAS",
        estimate = ns[["estimate"]],
        conf_low = NA_real_,
        conf_high = NA_real_,
        beta = NA_real_,
        beta_se = NA_real_,
        p_value = ns[["p_value"]],
        target_vif = NA_real_,
        convergence_note = "nonparametric_spearman",
        stringsAsFactors = FALSE
      )
    }
  }

  confounder_map <- list(
    primary_alone = character(),
    primary_plus_structural_ecm = "structural_ecm_score_z",
    primary_plus_core_fibrosis_benchmark = "core_fibrosis_stromal_benchmark_score_z",
    primary_plus_hsc_myofibroblast_benchmark = "hsc_myofibroblast_benchmark_score_z",
    primary_plus_inflammatory_activity = "inflammatory_activity_benchmark_score_z",
    primary_plus_myeloid_immune_proxy = "myeloid_immune_proxy_score_z",
    primary_plus_hepatocyte_identity_proxy = "hepatocyte_identity_proxy_score_z",
    ecm_excluded_alone = character(),
    ecm_excluded_plus_core_fibrosis_benchmark = "core_fibrosis_stromal_benchmark_score_z",
    ecm_excluded_plus_hsc_myofibroblast_benchmark = "hsc_myofibroblast_benchmark_score_z",
    ecm_excluded_plus_inflammatory_activity = "inflammatory_activity_benchmark_score_z"
  )
  if (sum(!is.na(score_df$nas_score)) >= 20) {
    confounder_map$primary_plus_NAS_activity <- "nas_score"
    confounder_map$ecm_excluded_plus_NAS_activity <- "nas_score"
  }
  if (sum(!is.na(score_df$age_numeric)) >= 20) {
    confounder_map$primary_plus_age <- "age_numeric"
  }
  if (length(unique(na.omit(score_df$sex_binary))) == 2 && sum(!is.na(score_df$sex_binary)) >= 20) {
    confounder_map$primary_plus_age_sex <- c("age_numeric", "sex_binary")
  }

  for (model_name in names(confounder_map)) {
    target <- ifelse(grepl("^ecm_excluded", model_name), "ecm_excluded_score_z", "primary_score_z")
    covs <- confounder_map[[model_name]]
    if (any(!covs %in% names(score_df))) next
    res <- fit_logistic_target(score_df, target, covariates = covs)
    if (!is.null(res)) {
      confound_models[[paste(acc, model_name, sep = "__")]] <- data.frame(
        accession = acc,
        target_score = sub("_score_z$", "", target),
        endpoint = "advanced_fibrosis_F3_F4_vs_F0_F2",
        model_type = model_name,
        covariates = ifelse(length(covs), paste(covs, collapse = ";"), "none"),
        res,
        stringsAsFactors = FALSE
      )
    }
  }

  corr_cols <- paste0(benchmark_score_names, "_score_z")
  corr_cols <- corr_cols[corr_cols %in% names(score_df)]
  pair_rows <- list()
  k <- 1L
  for (i in seq_along(corr_cols)) {
    for (j in seq_along(corr_cols)) {
      cr <- safe_cor(score_df[[corr_cols[[i]]]], score_df[[corr_cols[[j]]]], method = "spearman")
      pair_rows[[k]] <- data.frame(
        accession = acc,
        score_x = sub("_score_z$", "", corr_cols[[i]]),
        score_y = sub("_score_z$", "", corr_cols[[j]]),
        spearman_rho = cr[["estimate"]],
        p_value = cr[["p_value"]],
        n_pair = as.integer(cr[["n"]]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  correlation_tables[[acc]] <- do.call(rbind, pair_rows)

  if (acc %in% small_cohorts) {
    for (score_name in c("primary", "ecm_excluded")) {
      score_col <- paste0(score_name, "_score_z")
      bt <- bootstrap_score(score_df, score_col, acc, score_name, n_iter = n_boot)
      bootstrap_tables[[paste(acc, score_name, sep = "__")]] <- bt
      if (nrow(bt)) {
        bootstrap_summary[[paste(acc, score_name, sep = "__")]] <- data.frame(
          accession = acc,
          score_name = score_name,
          n_bootstrap_requested = n_boot,
          n_bootstrap_success = nrow(bt),
          beta_median = stats::median(bt$beta, na.rm = TRUE),
          beta_q025 = unname(stats::quantile(bt$beta, 0.025, na.rm = TRUE)),
          beta_q975 = unname(stats::quantile(bt$beta, 0.975, na.rm = TRUE)),
          or_median = stats::median(bt$estimate_or, na.rm = TRUE),
          or_q025 = unname(stats::quantile(bt$estimate_or, 0.025, na.rm = TRUE)),
          or_q975 = unname(stats::quantile(bt$estimate_or, 0.975, na.rm = TRUE)),
          fraction_beta_positive = mean(bt$beta > 0, na.rm = TRUE),
          fraction_glm_converged = mean(bt$converged, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }

      loo <- leave_one_out_score(score_df, score_col, acc, score_name)
      loo_tables[[paste(acc, score_name, sep = "__")]] <- loo
      if (nrow(loo)) {
        loo_summary[[paste(acc, score_name, sep = "__")]] <- data.frame(
          accession = acc,
          score_name = score_name,
          n_leave_one_out_models = nrow(loo),
          n_expected = sum(stats::complete.cases(score_df[, c("advanced_binary", score_col)])),
          beta_min = min(loo$beta, na.rm = TRUE),
          beta_median = stats::median(loo$beta, na.rm = TRUE),
          beta_max = max(loo$beta, na.rm = TRUE),
          or_min = min(loo$estimate, na.rm = TRUE),
          or_median = stats::median(loo$estimate, na.rm = TRUE),
          or_max = max(loo$estimate, na.rm = TRUE),
          fraction_beta_positive = mean(loo$beta > 0, na.rm = TRUE),
          fraction_p_lt_0_05 = mean(loo$p_value < 0.05, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

all_scores <- do.call(rbind, score_tables)
all_coverage <- do.call(rbind, coverage_tables)
all_benchmark_models <- do.call(rbind, benchmark_models)
all_confound_models <- do.call(rbind, confound_models)
all_correlations <- do.call(rbind, correlation_tables)
all_bootstrap <- do.call(rbind, bootstrap_tables)
all_bootstrap_summary <- do.call(rbind, bootstrap_summary)
all_loo <- do.call(rbind, loo_tables)
all_loo_summary <- do.call(rbind, loo_summary)

if (!is.null(all_benchmark_models) && nrow(all_benchmark_models)) {
  all_benchmark_models$p_fdr_within_endpoint_model <- ave(
    all_benchmark_models$p_value,
    paste(all_benchmark_models$endpoint, all_benchmark_models$model_type),
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
}

write_tsv_base(all_scores, file.path(result_dir, "benchmark_robustness_sensitivity_scores.tsv"))
write_tsv_base(all_coverage, file.path(result_dir, "benchmark_and_proxy_gene_set_coverage.tsv"))
write_tsv_base(all_benchmark_models, file.path(result_dir, "fibrosis_stromal_benchmark_results.tsv"))
write_tsv_base(all_confound_models, file.path(result_dir, "bulk_confounding_sensitivity_models.tsv"))
write_tsv_base(all_correlations, file.path(result_dir, "bulk_score_proxy_correlation_matrix.tsv"))
write_tsv_base(all_bootstrap, file.path(result_dir, "small_sample_bootstrap_iterations.tsv"))
write_tsv_base(all_bootstrap_summary, file.path(result_dir, "small_sample_bootstrap_stability_summary.tsv"))
write_tsv_base(all_loo, file.path(result_dir, "small_sample_leave_one_out_models.tsv"))
write_tsv_base(all_loo_summary, file.path(result_dir, "small_sample_leave_one_out_summary.tsv"))

write_tsv_base(all_benchmark_models, file.path(supp_dir, "Supplementary_Table_3_benchmark_model_results.tsv"))
write_tsv_base(all_confound_models, file.path(supp_dir, "Supplementary_Table_4_bulk_confounding_sensitivity.tsv"))
write_tsv_base(all_bootstrap_summary, file.path(supp_dir, "Supplementary_Table_5_small_sample_stability_summary.tsv"))
write_tsv_base(all_correlations, file.path(supp_dir, "Supplementary_Table_6_score_proxy_correlations.tsv"))

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

advanced_bench <- all_benchmark_models[
  all_benchmark_models$endpoint == "advanced_fibrosis_F3_F4_vs_F0_F2" &
    all_benchmark_models$model_type == "logistic_base",
]
advanced_bench <- advanced_bench[advanced_bench$score_name %in% c(
  "primary", "ecm_excluded", "structural_ecm", "core_fibrosis_stromal_benchmark",
  "hsc_myofibroblast_benchmark", "inflammatory_activity_benchmark"
), ]

summary_lines <- c(
  "# Benchmark and Robustness Sensitivity Analysis Summary",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Scope",
  "",
  "This analysis added transparent fibrosis/stromal comparator scores, small-cohort stability checks, and bulk proxy-confounding sensitivity analyses. Comparator gene sets were used as transparent biological benchmarks and proxies, not as claimed externally validated signatures.",
  "",
  "## Benchmark Advanced-Fibrosis Models",
  ""
)
for (i in seq_len(nrow(advanced_bench))) {
  r <- advanced_bench[i, ]
  summary_lines <- c(summary_lines, paste0(
    "- ", r$accession, " ", r$score_name, ": OR ", format_num(r$estimate),
    " (95% CI ", format_num(r$conf_low), "-", format_num(r$conf_high),
    "), P=", format_num(r$p_value), "."
  ))
}
summary_lines <- c(
  summary_lines,
  "",
  "## Small-Cohort Stability",
  ""
)
for (i in seq_len(nrow(all_bootstrap_summary))) {
  r <- all_bootstrap_summary[i, ]
  summary_lines <- c(summary_lines, paste0(
    "- ", r$accession, " ", r$score_name, ": bootstrap median OR ",
    format_num(r$or_median), " (2.5%-97.5% ", format_num(r$or_q025),
    "-", format_num(r$or_q975), "), positive beta fraction ",
    format_num(r$fraction_beta_positive), "."
  ))
}
for (i in seq_len(nrow(all_loo_summary))) {
  r <- all_loo_summary[i, ]
  summary_lines <- c(summary_lines, paste0(
    "- ", r$accession, " ", r$score_name, ": leave-one-out OR range ",
    format_num(r$or_min), "-", format_num(r$or_max), ", positive beta fraction ",
    format_num(r$fraction_beta_positive), "."
  ))
}
summary_lines <- c(
  summary_lines,
  "",
  "## Interpretation Boundary",
  "",
  "- These analyses test robustness to selected transparent comparators and bulk marker proxies.",
  "- They do not prove cell-intrinsic mechanotransduction, stiffness exposure, or mechanical memory.",
  "- Estimated marker-proxy adjustment cannot fully separate confounding from biological mediation in bulk tissue.",
  "- Small-cohort bootstrap and leave-one-out analyses evaluate directional stability, not definitive external validation."
)
writeLines(summary_lines, file.path(result_dir, "benchmark_robustness_sensitivity_analysis_summary.md"), useBytes = TRUE)

panel_label <- function(label) {
  graphics::mtext(label, side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
}

short_score_label <- function(x) {
  y <- as.character(x)
  y[y == "primary"] <- "Primary"
  y[y == "ecm_excluded"] <- "ECM-excl."
  y[y == "structural_ecm"] <- "Structural ECM"
  y[y == "core_fibrosis_stromal_benchmark"] <- "Fibrosis/stromal"
  y[y == "hsc_myofibroblast_benchmark"] <- "HSC/myofib."
  y[y == "inflammatory_activity_benchmark"] <- "Inflammatory"
  y[y == "myeloid_immune_proxy"] <- "Myeloid"
  y[y == "hepatocyte_identity_proxy"] <- "Hepatocyte"
  y
}

short_model_label <- function(x) {
  y <- as.character(x)
  y[y == "primary_alone"] <- "Primary alone"
  y[y == "primary_plus_core_fibrosis_benchmark"] <- "+ fibrosis/stromal"
  y[y == "primary_plus_hsc_myofibroblast_benchmark"] <- "+ HSC/myofib."
  y[y == "primary_plus_inflammatory_activity"] <- "+ inflammatory"
  y[y == "primary_plus_myeloid_immune_proxy"] <- "+ myeloid"
  y[y == "primary_plus_hepatocyte_identity_proxy"] <- "+ hepatocyte"
  y[y == "primary_plus_NAS_activity"] <- "+ NAS"
  y
}

draw_forest <- function(dat, title, score_order = NULL) {
  dat <- dat[is.finite(dat$estimate) & is.finite(dat$conf_low) & is.finite(dat$conf_high), ]
  if (!nrow(dat)) {
    plot.new(); title(title); text(0.5, 0.5, "No estimable models"); return(invisible(NULL))
  }
  if (!is.null(score_order)) {
    dat$score_name <- factor(dat$score_name, levels = score_order)
    dat <- dat[order(dat$accession, dat$score_name), ]
  }
  labels <- paste(dat$accession, short_score_label(dat$score_name), sep = " | ")
  y <- seq_len(nrow(dat))
  xlim <- range(log(c(dat$conf_low, dat$conf_high, 1)), finite = TRUE)
  xlim <- xlim + c(-0.15, 0.15)
  graphics::plot(log(dat$estimate), y, xlim = xlim, ylim = c(0.5, length(y) + 0.8),
                 yaxt = "n", ylab = "", xlab = "Log odds ratio per 1 SD score",
                 pch = 19, col = "#1f4e79", main = title)
  graphics::segments(log(dat$conf_low), y, log(dat$conf_high), y, col = "#4f6d7a", lwd = 1.5)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  graphics::axis(2, at = y, labels = labels, las = 2, cex.axis = 0.68)
}

draw_adjusted <- function(dat) {
  dat <- dat[dat$target_score == "primary" & is.finite(dat$estimate), ]
  keep <- c(
    "primary_alone", "primary_plus_core_fibrosis_benchmark",
    "primary_plus_hsc_myofibroblast_benchmark", "primary_plus_inflammatory_activity",
    "primary_plus_myeloid_immune_proxy", "primary_plus_hepatocyte_identity_proxy",
    "primary_plus_NAS_activity"
  )
  dat <- dat[dat$model_type %in% keep, ]
  dat$model_type <- factor(dat$model_type, levels = keep)
  dat <- dat[order(dat$accession, dat$model_type), ]
  labels <- paste(dat$accession, short_model_label(dat$model_type), sep = " | ")
  y <- seq_len(nrow(dat))
  xlim <- range(log(c(dat$conf_low, dat$conf_high, 1)), finite = TRUE)
  xlim <- xlim + c(-0.15, 0.15)
  graphics::plot(log(dat$estimate), y, xlim = xlim, ylim = c(0.5, length(y) + 0.8),
                 yaxt = "n", ylab = "", xlab = "Adjusted log odds ratio",
                 pch = 19, col = "#7c2d12", main = "Primary score sensitivity models")
  graphics::segments(log(dat$conf_low), y, log(dat$conf_high), y, col = "#9a3412", lwd = 1.5)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  graphics::axis(2, at = y, labels = labels, las = 2, cex.axis = 0.66)
}

draw_bootstrap <- function(dat) {
  dat <- dat[is.finite(dat$estimate_or) & dat$estimate_or > 0, ]
  dat <- dat[dat$estimate_or < stats::quantile(dat$estimate_or, 0.995, na.rm = TRUE), ]
  groups <- paste(dat$accession, dat$score_name, sep = "\n")
  graphics::boxplot(log(dat$estimate_or) ~ groups, col = "#d8e2dc", border = "#264653",
                    ylab = "Bootstrap log odds ratio", las = 2,
                    main = "Bootstrap stability in small cohorts", cex.axis = 0.65)
  graphics::abline(h = 0, lty = 2, col = "grey50")
}

draw_loo <- function(dat) {
  dat <- dat[dat$score_name == "primary" & is.finite(dat$estimate), ]
  if (!nrow(dat)) {
    plot.new(); title("Leave-one-out primary score"); text(0.5, 0.5, "No estimable models"); return(invisible(NULL))
  }
  dat$group <- paste(dat$accession, dat$score_name, sep = " ")
  groups <- unique(dat$group)
  graphics::plot(NA, xlim = c(0.5, length(groups) + 0.5),
                 ylim = range(log(dat$estimate), finite = TRUE) + c(-0.2, 0.2),
                 xaxt = "n", xlab = "", ylab = "Leave-one-out log odds ratio",
                 main = "Single-sample deletion influence")
  graphics::abline(h = 0, lty = 2, col = "grey50")
  for (i in seq_along(groups)) {
    vals <- log(dat$estimate[dat$group == groups[[i]]])
    graphics::stripchart(vals, at = i, vertical = TRUE, add = TRUE, method = "jitter",
                         pch = 16, col = "#38664166", cex = 0.7)
    graphics::segments(i - 0.18, stats::median(vals), i + 0.18, stats::median(vals), lwd = 2, col = "#386641")
  }
  graphics::axis(1, at = seq_along(groups), labels = gsub(" ", "\n", groups), cex.axis = 0.75)
}

draw_correlation_heatmap <- function(corr_df, accession = "GSE162694") {
  use_scores <- c(
    "primary", "ecm_excluded", "structural_ecm", "core_fibrosis_stromal_benchmark",
    "hsc_myofibroblast_benchmark", "inflammatory_activity_benchmark",
    "myeloid_immune_proxy", "hepatocyte_identity_proxy"
  )
  df <- corr_df[corr_df$accession == accession & corr_df$score_x %in% use_scores & corr_df$score_y %in% use_scores, ]
  mat <- matrix(NA_real_, length(use_scores), length(use_scores), dimnames = list(use_scores, use_scores))
  for (i in seq_len(nrow(df))) mat[df$score_x[[i]], df$score_y[[i]]] <- df$spearman_rho[[i]]
  graphics::par(mar = c(7.5, 7.8, 3.2, 1.5))
  pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
  graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, ]),
                  col = pal, zlim = c(-1, 1), xaxt = "n", yaxt = "n",
                  xlab = "", ylab = "", main = paste0(accession, " score/proxy correlations"))
  graphics::axis(1, at = seq_len(ncol(mat)), labels = short_score_label(colnames(mat)), las = 2, cex.axis = 0.62)
  graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(short_score_label(rownames(mat))), las = 2, cex.axis = 0.62)
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      graphics::text(j, nrow(mat) - i + 1, sprintf("%.2f", mat[i, j]), cex = 0.55)
    }
  }
}

figure7_stem <- "Figure_7_benchmark_robustness_sensitivity"
png_path <- file.path(figure_dir, paste0(figure7_stem, ".png"))
pdf_path <- file.path(figure_dir, paste0(figure7_stem, ".pdf"))
for (device in c("png", "pdf")) {
  if (device == "png") {
    grDevices::png(png_path, width = 8400, height = 5200, res = 600, type = "cairo")
  } else {
    grDevices::pdf(pdf_path, width = 14, height = 8.7, onefile = TRUE)
  }
  graphics::layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 2, byrow = TRUE), widths = c(1.18, 1.18, 0.90))
  graphics::par(mar = c(4.2, 9.8, 3.2, 1.1), cex = 0.85)
  draw_forest(advanced_bench, "Benchmark score associations",
              score_order = c("primary", "ecm_excluded", "structural_ecm",
                              "core_fibrosis_stromal_benchmark",
                              "hsc_myofibroblast_benchmark",
                              "inflammatory_activity_benchmark"))
  panel_label("A")
  graphics::par(mar = c(4.2, 9.8, 3.2, 1.1), cex = 0.85)
  draw_adjusted(all_confound_models)
  panel_label("B")
  graphics::par(mar = c(5.2, 4.5, 3.2, 1.1), cex = 0.85)
  draw_bootstrap(all_bootstrap)
  panel_label("C")
  graphics::par(mar = c(5.0, 4.5, 3.2, 1.1), cex = 0.85)
  draw_loo(all_loo)
  panel_label("D")
  draw_correlation_heatmap(all_correlations, accession = "GSE162694")
  panel_label("E")
  grDevices::dev.off()
}

caption <- data.frame(
  figure = "Figure 7",
  caption = paste(
    "Benchmark and robustness analyses for fibrosis and proxy sensitivity.",
    "(A) Advanced-fibrosis logistic associations for the locked primary program, ECM-excluded score, structural ECM-only score and transparent fibrosis, stromal and inflammatory benchmark or proxy scores.",
    "(B) Primary-score association estimates after selected one-proxy adjustments for core fibrosis, myofibroblast, inflammatory, myeloid, hepatocyte and NAS-related signals.",
    "(C) Stratified bootstrap distributions of logistic effect estimates in the two smaller fibrosis-mapped cohorts.",
    "(D) Leave-one-out estimates for the primary score in the two smaller cohorts.",
    "(E) Representative score and proxy correlation heatmap in GSE162694.",
    "Benchmark and proxy analyses are sensitivity analyses and do not establish cell-intrinsic mechanotransduction or mechanical memory."
  ),
  stringsAsFactors = FALSE
)
write_tsv_base(caption, file.path(figure_dir, paste0(figure7_stem, "_caption.tsv")))

figure7_manifest <- data.frame(
  figure = figure7_stem,
  panel = c("A", "B", "C", "D", "E"),
  source_figure = c(
    paste0("figures/benchmark_robustness_sensitivity/", figure7_stem, ".png"),
    paste0("figures/benchmark_robustness_sensitivity/", figure7_stem, ".png"),
    paste0("figures/benchmark_robustness_sensitivity/", figure7_stem, ".png"),
    paste0("figures/benchmark_robustness_sensitivity/", figure7_stem, ".png"),
    paste0("figures/benchmark_robustness_sensitivity/", figure7_stem, ".png")
  ),
  source_data = c(
    "results/benchmark_robustness_sensitivity/fibrosis_stromal_benchmark_results.tsv",
    "results/benchmark_robustness_sensitivity/bulk_confounding_sensitivity_models.tsv",
    "results/benchmark_robustness_sensitivity/small_sample_bootstrap_iterations.tsv; results/benchmark_robustness_sensitivity/small_sample_bootstrap_stability_summary.tsv",
    "results/benchmark_robustness_sensitivity/small_sample_leave_one_out_models.tsv; results/benchmark_robustness_sensitivity/small_sample_leave_one_out_summary.tsv",
    "results/benchmark_robustness_sensitivity/bulk_score_proxy_correlation_matrix.tsv"
  ),
  generating_script = rep("scripts/22_benchmark_robustness_sensitivity_analyses.R", 5),
  role = c(
    "fibrosis and proxy benchmark associations",
    "proxy-adjusted sensitivity models",
    "small-cohort stratified bootstrap",
    "leave-one-out small-cohort stability",
    "score and proxy correlation context"
  ),
  wording_note = c(
    "benchmark comparator analysis, not superiority testing",
    "bulk proxy adjustment; collinearity and mediation cannot be fully separated",
    "direction stability check, not definitive external validation",
    "single-sample influence check, not definitive external validation",
    "bulk score-proxy overlap context"
  ),
  stringsAsFactors = FALSE
)
write_tsv_base(figure7_manifest, file.path(figure_dir, paste0(figure7_stem, "_source_manifest.tsv")))

legacy_png <- file.path(figure_dir, "Figure_7_benchmark_robustness_sensitivity_benchmark_sensitivity.png")
legacy_pdf <- file.path(figure_dir, "Figure_7_benchmark_robustness_sensitivity_benchmark_sensitivity.pdf")
legacy_caption <- file.path(figure_dir, "Figure_7_benchmark_robustness_sensitivity_caption.tsv")
file.copy(png_path, legacy_png, overwrite = TRUE)
file.copy(pdf_path, legacy_pdf, overwrite = TRUE)
file.copy(file.path(figure_dir, paste0(figure7_stem, "_caption.tsv")), legacy_caption, overwrite = TRUE)

message("Done. Benchmark robustness outputs written to ", result_dir)

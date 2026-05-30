#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

sig_dir <- file.path(project_dir, "data", "signatures")
dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

add_genes <- function(genes, source_category, evidence_tier, evidence_note, include_primary = TRUE, is_ecm_structural = FALSE) {
  tibble(
    gene_symbol = genes,
    source_category = source_category,
    evidence_tier = evidence_tier,
    evidence_note = evidence_note,
    include_in_primary_score = include_primary,
    is_ecm_structural_gene = is_ecm_structural
  )
}

prior <- bind_rows(
  add_genes(
    c("YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "CTGF", "CYR61", "AMOTL2", "ANKRD1"),
    "YAP_TAZ_TEAD_mechanotransduction",
    "A/B",
    "Mechanotransduction and stiffness-responsive transcriptional context; liver/HSC relevance requires dataset-level confirmation."
  ),
  add_genes(
    c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "ATF3", "MAFF", "MAFG"),
    "AP1_chromatin_priming_context",
    "A/B",
    "AP-1 and immediate-early transcriptional programs linked to stiffness response and HSC activation context."
  ),
  add_genes(
    c("RHOA", "ROCK1", "ROCK2", "MYL9", "MYLK", "ACTA2", "TAGLN", "CNN1", "VCL", "ZYX", "PXN", "TLN1", "VASP"),
    "RhoA_ROCK_actomyosin_focal_adhesion",
    "A/B",
    "Actomyosin contractility and focal-adhesion mechanotransduction; may overlap with myofibroblast abundance."
  ),
  add_genes(
    c("ITGA1", "ITGA2", "ITGA5", "ITGAV", "ITGB1", "ITGB3", "PTK2", "SRC", "FERMT2", "PARVA", "PARVB"),
    "integrin_FAK_adhesion",
    "A/B",
    "Integrin and focal adhesion signaling in matrix sensing."
  ),
  add_genes(
    c("MRTFA", "MRTFB", "SRF", "LMNA", "LIMD1", "LATS1", "LATS2", "RAC1", "CDC42"),
    "nuclear_mechanotransduction_MRTF_SRF_Hippo_context",
    "B",
    "Nuclear and cytoskeletal mechanotransduction prior; not liver-specific alone."
  ),
  add_genes(
    c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SERPINE1", "THBS1"),
    "TGF_beta_SMAD_fibrogenic_context",
    "A/B",
    "Fibrogenic signaling and mechanotransduction-adjacent context; not direct memory evidence."
  ),
  add_genes(
    c("RELA", "NFKB1", "NFKBIA", "IL6", "CXCL8", "CCL2"),
    "NFkB_inflammatory_mechanostress_context",
    "B",
    "Inflammatory mechanostress context; must be separated from inflammation confounding."
  ),
  add_genes(
    c("LRAT", "RBP1", "DES", "PDGFRB", "PDGFRA", "COL15A1", "COLEC11"),
    "HSC_portal_fibroblast_context_markers",
    "A",
    "HSC or portal fibroblast context markers for lineage interpretation; not necessarily primary mechanostress genes.",
    include_primary = FALSE
  ),
  add_genes(
    c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "FN1", "LAMA2", "LAMB1", "LUM", "DCN", "POSTN"),
    "ECM_structural_fibrosis_context",
    "D",
    "Structural ECM and fibrosis abundance context; excluded from ECM-excluded sensitivity score.",
    include_primary = TRUE,
    is_ecm_structural = TRUE
  ),
  add_genes(
    c("TIMP1", "TIMP2", "MMP2", "MMP14", "LOX", "LOXL2", "PLOD2", "SPARC", "VCAN", "THY1"),
    "ECM_remodeling_HSC_activation_context",
    "A/B",
    "ECM remodeling, crosslinking, and HSC activation context; may overlap with fibrosis burden."
  ),
  add_genes(
    c("PFN1", "CFL1", "ACTB", "FLNA", "FLNB", "CAV1"),
    "Zhao_2026_template_overlap_and_cytoskeleton_context",
    "C",
    "External breast-cancer mechanical-memory template overlap or cytoskeletal context; not sufficient alone for MASLD mechanism."
  )
) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  mutate(
    gene_symbol = str_to_upper(gene_symbol),
    include_in_ecm_excluded_score = include_in_primary_score & !is_ecm_structural_gene,
    hgnc_status = "requires_symbol_audit_before_submission",
    lock_status = "phase2_v1_locked_prior",
    lock_date = "2026-05-23"
  ) %>%
  arrange(source_category, gene_symbol)

write_tsv(prior, file.path(sig_dir, "signature_v0_prior.tsv"))
write_tsv(prior %>% filter(include_in_primary_score), file.path(sig_dir, "signature_v1_locked.tsv"))
write_tsv(
  prior %>%
    count(source_category, evidence_tier, include_in_primary_score, is_ecm_structural_gene, name = "n_genes"),
  file.path(sig_dir, "signature_source_manifest.tsv")
)
write_tsv(
  prior %>% filter(!include_in_primary_score) %>% mutate(exclusion_reason = "context marker retained for annotation, excluded from primary score"),
  file.path(sig_dir, "signature_exclusion_log.tsv")
)

summary_path <- file.path(project_dir, "docs", "signature_freeze_summary.md")
cat(
  "# Signature Freeze Summary\n\n",
  "Freeze date: 2026-05-23\n\n",
  "The Phase 2 prior signature is locked before expression-result inspection. It uses external mechanotransduction, HSC fibrosis, ECM remodeling, and bounded Zhao 2026 template-overlap categories. Validation cohorts must not be used for gene selection, weighting, or cutoff tuning.\n\n",
  "## Output Files\n\n",
  "- `data/signatures/signature_v0_prior.tsv`\n",
  "- `data/signatures/signature_v1_locked.tsv`\n",
  "- `data/signatures/signature_source_manifest.tsv`\n",
  "- `data/signatures/signature_exclusion_log.tsv`\n\n",
  "## Counts\n\n",
  "- Prior genes: ", nrow(prior), "\n",
  "- Primary score genes: ", sum(prior$include_in_primary_score), "\n",
  "- ECM-excluded score genes: ", sum(prior$include_in_ecm_excluded_score), "\n",
  "- Structural ECM genes flagged: ", sum(prior$is_ecm_structural_gene), "\n\n",
  "## Boundary\n\n",
  "This is a mechanical-stress-associated memory-like program prior, not proof of HSC mechanical memory. Direct mechanical-memory language requires stiffness-history experiments.\n",
  file = summary_path,
  sep = ""
)

message("Done. Signature files written to ", sig_dir)

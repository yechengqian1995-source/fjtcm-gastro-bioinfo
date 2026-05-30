# Single-Nucleus Processed Source Tables

This directory contains processed source tables used for the single-nucleus UMAP and donor-aware stromal panels.
Single-cell raw matrices and the Seurat RDS object are intentionally excluded from this review package to keep the package small and to avoid redistributing raw single-cell data.
The sketch Seurat object can be regenerated from public GSE212837 inputs with `scripts/08_gse212837_hsc_stromal_validation.R` followed by `scripts/16_gse212837_single_nucleus_embedding.R`.

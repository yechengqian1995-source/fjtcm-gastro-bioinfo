# Script roadmap for the reproducibility archive

This folder contains the reproducible analysis and figure-generation scripts for the manuscript-supporting archive.

Important path note: these scripts were written for the original project-root structure with top-level `data`, `results`, `figures`, `docs`, `manuscript_outputs` and `scripts` folders. Full source-level regeneration may require downloading the original public raw data and restoring that project layout.

Key figure scripts:

| Script | Purpose |
|---|---|
| `06_phase2_bulk_visuals.R` | Bulk association source panels. |
| `08_gse212837_hsc_stromal_validation.R` | Donor-aware single-nucleus stromal summaries. |
| `09_gse292268_spatial_localization.R` | CosMx spatial localization results. |
| `10_network_prioritization.R` | Network-prioritization panels. |
| `11_bulk_coexpression_network_context.R` | Gene-level association and coexpression context. |
| `14_main_composite_figures.R` | Current composite assembly of manuscript figures. |
| `16_gse212837_single_nucleus_embedding.R` | Single-nucleus sketch UMAP and marker visualization. |
| `21_build_locked_gene_supplementary_table.R` | Supplementary Table 1 builder. |
| `22_editorial_strengthening_analyses.R` | Benchmark, proxy-adjustment, bootstrap and leave-one-out sensitivity analyses. |
| `25_rerender_reviewed_single_spatial_panels.R` | Final reviewed rerender for the current Figure 2 and Figure 4 source panels. |

The single-nucleus raw matrices and sketch Seurat object are not packaged. Regenerate them from public GSE212837 inputs or use the original project-root object when full single-nucleus reproduction is required.

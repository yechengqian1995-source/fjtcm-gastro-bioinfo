# Benchmark and Robustness Sensitivity Analysis Summary

Generated: 2026-05-27

## Scope

This analysis added transparent fibrosis/stromal comparator scores, small-cohort stability checks, and bulk proxy-confounding sensitivity analyses. Comparator gene sets were used as transparent biological benchmarks and proxies, not as claimed externally validated signatures.

## Benchmark Advanced-Fibrosis Models

- GSE130970 primary: OR 5.52 (95% CI 2.07-14.7), P=0.000627.
- GSE130970 ecm_excluded: OR 5.78 (95% CI 2.11-15.8), P=0.000645.
- GSE130970 structural_ecm: OR 3.91 (95% CI 1.75-8.73), P=0.000886.
- GSE130970 core_fibrosis_stromal_benchmark: OR 5.91 (95% CI 2.22-15.7), P=0.000373.
- GSE130970 hsc_myofibroblast_benchmark: OR 3.41 (95% CI 1.68-6.95), P=0.000712.
- GSE130970 inflammatory_activity_benchmark: OR 2.63 (95% CI 1.36-5.10), P=0.00409.
- GSE135251 primary: OR 2.96 (95% CI 2.03-4.32), P=0.0000000190.
- GSE135251 ecm_excluded: OR 2.95 (95% CI 2.01-4.31), P=0.0000000254.
- GSE135251 structural_ecm: OR 2.20 (95% CI 1.58-3.07), P=0.00000300.
- GSE135251 core_fibrosis_stromal_benchmark: OR 2.56 (95% CI 1.77-3.69), P=0.000000500.
- GSE135251 hsc_myofibroblast_benchmark: OR 2.85 (95% CI 1.92-4.21), P=0.000000161.
- GSE135251 inflammatory_activity_benchmark: OR 1.65 (95% CI 1.20-2.26), P=0.00198.
- GSE162694 primary: OR 8.87 (95% CI 3.52-22.3), P=0.00000361.
- GSE162694 ecm_excluded: OR 8.19 (95% CI 3.33-20.1), P=0.00000456.
- GSE162694 structural_ecm: OR 5.03 (95% CI 2.38-10.6), P=0.0000232.
- GSE162694 core_fibrosis_stromal_benchmark: OR 7.70 (95% CI 3.35-17.7), P=0.00000156.
- GSE162694 hsc_myofibroblast_benchmark: OR 6.05 (95% CI 2.96-12.4), P=0.000000765.
- GSE162694 inflammatory_activity_benchmark: OR 3.99 (95% CI 2.18-7.31), P=0.00000736.

## Small-Cohort Stability

- GSE130970 primary: bootstrap median OR 6.12 (2.5%-97.5% 2.46-24.6), positive beta fraction 1.00.
- GSE130970 ecm_excluded: bootstrap median OR 6.13 (2.5%-97.5% 2.60-22.5), positive beta fraction 1.00.
- GSE162694 primary: bootstrap median OR 9.28 (2.5%-97.5% 4.70-34.9), positive beta fraction 1.00.
- GSE162694 ecm_excluded: bootstrap median OR 8.57 (2.5%-97.5% 4.33-29.7), positive beta fraction 1.00.
- GSE130970 primary: leave-one-out OR range 5.02-8.13, positive beta fraction 1.00.
- GSE130970 ecm_excluded: leave-one-out OR range 5.22-8.29, positive beta fraction 1.00.
- GSE162694 primary: leave-one-out OR range 8.21-11.6, positive beta fraction 1.00.
- GSE162694 ecm_excluded: leave-one-out OR range 7.55-10.7, positive beta fraction 1.00.

## Interpretation Boundary

- These analyses test robustness to selected transparent comparators and bulk marker proxies.
- They do not prove cell-intrinsic mechanotransduction, stiffness exposure, or mechanical memory.
- Estimated marker-proxy adjustment cannot fully separate confounding from biological mediation in bulk tissue.
- Small-cohort bootstrap and leave-one-out analyses evaluate directional stability, not definitive external validation.

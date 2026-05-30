# Reproducibility guide

## Scope

This archive supports a public-dataset bioinformatics study of a mechanostress/stromal transcriptomic program in MASLD-associated liver fibrosis.

## Contents

- `data/gene_sets/`: locked 97-gene program and benchmark/proxy gene sets.
- `data/supplementary_tables/`: Supplementary Tables 1-6 used by the manuscript.
- `data/processed_results/results/`: processed bulk, benchmark, sensitivity, network, single-cell/single-nucleus, spatial, and figure-source tables.
- `scripts/`: analysis and figure-generation scripts.
- `renv.lock`, `requirements.txt`, and `sessionInfo_2026-05-27.txt`: environment records.

## Suggested script order

The scripts are numbered to reflect the intended workflow. Start with gene-set freezing and GEO metadata inventory, then bulk preprocessing/scoring, benchmark and sensitivity analyses, single-nucleus and spatial analyses, network context, and final figure assembly.

## Raw-data policy

Raw public datasets are not redistributed. Download raw inputs from the original repositories listed in the manuscript. Large intermediate single-cell objects are intentionally excluded from this public archive.

## Version to cite

Use release `v1.0.3` for the manuscript submission archive. After a GitHub release and Zenodo archive are created, replace the placeholders in the manuscript with the real GitHub URL, commit hash, and Zenodo DOI.


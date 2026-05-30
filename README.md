# MASLD mechanostress/stromal transcriptomic program

This repository contains the manuscript-supporting analysis scripts, locked gene sets, processed result tables, figure source data, and reproducibility notes for a public-dataset bioinformatics study of a mechanostress/stromal transcriptomic program in MASLD-associated liver fibrosis.

## Repository contents

- `data/gene_sets/`: locked 97-gene program and benchmark/proxy gene sets.
- `data/supplementary_tables/`: Supplementary Tables 1-6.
- `data/processed_results/results/`: processed tables supporting bulk, benchmark, sensitivity, network, single-cell/single-nucleus, spatial, and figure-source analyses.
- `scripts/`: analysis and figure-generation scripts.
- `docs/`: reproducibility guide and raw-data redistribution notice.
- `renv.lock`, `requirements.txt`, `sessionInfo_2026-05-27.txt`: environment records.

## Data sources

Raw datasets analysed in the manuscript are publicly available from the original repositories under the accession numbers reported in the manuscript. Raw sequencing files, raw single-cell/single-nucleus matrices, large Seurat/RDS objects, protected clinical records, and copyrighted article PDFs are not redistributed here.

## Reproducibility

The scripts are numbered to reflect the intended analysis order. Processed result tables and figure source data are provided to support inspection of manuscript-level analyses and figures. Full regeneration of all single-cell and spatial outputs may require downloading the original public raw data and rebuilding large intermediate objects.

## Version

Manuscript submission archive: `v1.0.1`

GitHub repository: `https://github.com/yechengqian1995-source/fjtcm-gastro-bioinfo`

Commit: see the GitHub release tag for the exact archived commit.

Zenodo DOI: to be assigned by Zenodo after release archiving.

## License

Code is released under the MIT License. Processed tables and documentation are provided for academic reuse with citation of the associated manuscript and archived repository record.



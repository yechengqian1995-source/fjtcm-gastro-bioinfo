# GitHub upload checklist

Repository URL currently shown by GitHub:
https://github.com/yechengqian1995-source/fjtcm-gastro-bioinfo

## Upload source folder

Upload the contents of this local folder, not a zip file:
E:\网毒\masld_mechanical_memory\repository_deposit_v1.0.0

## Required top-level files and folders

- README.md
- LICENSE
- CITATION.cff
- .zenodo.json
- renv.lock
- requirements.txt
- sessionInfo_2026-05-27.txt
- R_ENVIRONMENT_FREEZE_NOTE.md
- data/
- scripts/
- docs/

## After upload

1. Commit message: Initial reproducibility archive for MASLD mechanostress study
2. After commit, copy the short commit hash from the GitHub page.
3. Create a release named v1.0.0.
4. Connect the GitHub repository to Zenodo and archive release v1.0.0.
5. Replace manuscript placeholders with the final GitHub URL, release tag, commit hash, and Zenodo DOI.

## Do not upload separately

- Word manuscript drafts
- reviewer comments
- team logs
- raw single-cell matrices
- large Seurat/RDS objects
- downloaded article PDFs
- local work logs or private WeChat folders

# AGENTS.md

**GOAL:** Identify what clusters 6 and 7 actually are, so the owner can decide whether to keep or remove cluster 6 and whether to fund follow-up sequencing for cluster 7.

See `spec.md` for all details.

## MUST-NOTs

- Never remove a cluster based on low gene count alone; check `% mito` and markers first.
- Never call a cluster novel just because it sits apart on the UMAP.
- Always check lineage co-expression on the raw counts layer, not `.X`.
- Never invent thresholds or numbers not in the data.

## Operating environment

- Use the existing virtual environment in `.venv/`.
- Run Python with `.venv/bin/python`.
- The dataset is `data/pbmc3k.h5ad`.
- Load it with `scanpy.read_h5ad`.
- Clusters are in `adata.obs['leiden']`; these are numbered labels, not cell-type names.
- UMAP coordinates are in `adata.obsm['X_umap']`.
- Per-cell QC columns are in `adata.obs`: `n_genes`, `total_counts`, and `pct_mito`.
- Raw counts are in `adata.layers['counts']`; `.X` contains log-normalised expression.
- Outputs belong in `results/`.

## Version control

Commit the current state before any big change. Commit again whenever something starts working. Use short, clear commit messages.

## Interpretation rule

Never name or judge a cluster without first reporting its top marker genes and its QC numbers—`n_genes`, total counts, and `% mito`—compared to the other clusters.

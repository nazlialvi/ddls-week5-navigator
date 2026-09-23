# PBMC Single-Cell Data Navigator

An interactive single-cell navigator that identifies two disputed clusters in a PBMC dataset and checks the analysis.

**Live app: LINK**

## The question

The owner wanted to remove cluster 6 as junk because of its low gene count and was considering funding follow-up sequencing for cluster 7 because it was small, distinct, and had high counts.

## The answer

**Cluster 6 — 13 cells: platelets. Recommendation: keep.** Its top markers include `PF4`, `PPBP`, `GNG11`, `SDPR`, and `NRGN`; these are expressed in approximately 92–100% of cluster 6 cells versus approximately 0.1–2.6% of other cells. Its median mitochondrial percentage is **1.57%**, the lowest of all clusters, so these are not dying cells. The low median gene count (**350**) is expected for platelets, which have no nucleus. Most cluster 6 cells have 0–1 counts of `LYZ`/`CD74`, consistent with ambient RNA; two individual cells have higher `LYZ`/`CD74` and may be single doublets, but the cluster as a whole is a coherent platelet population.

**Cluster 7 — 10 cells: proliferating cytotoxic NK/T cells. Recommendation: do not fund follow-up sequencing.** This is a known cell state, not a novel type. Its top-ranked markers are mostly housekeeping and cell-cycle genes; `MKI67` is expressed in **70%** of cluster 7 cells versus **0.1%** of other cells, and `TOP2A` in **50%** versus **0.1%**. Cytotoxic genes including `NKG7`, `GZMA`, and `CST7` are enriched. Its median total counts are **8,508** versus **2,194** for all other cells, while its median mitochondrial percentage is normal at **2.03%**. **0 cells** co-express NK/T with B-cell or myeloid markers, so the tested doublet signature is absent. Follow-up sequencing is therefore not funded.

## Cluster summary

| Cluster | Cells | Suggested cell type |
|---:|---:|---|
| 0 | 1,197 | CD4 T cells |
| 1 | 489 | CD14+ monocytes |
| 2 | 445 | NK / CD8 T cells |
| 3 | 347 | B cells |
| 4 | 163 | FCGR3A+ monocytes |
| 5 | 36 | dendritic cells |
| 6 | 13 | platelets |
| 7 | 10 | proliferating cytotoxic NK/T cells |

## What the navigator does

- Displays the UMAP coloured by cluster, QC metric, or curated gene.
- Shows a summary table for all clusters.
- Provides a cluster picker and **Run** button. The selected cluster is highlighted on the UMAP, and the results show its QC values versus other cells plus a gene-enrichment table.

## How to run

Open `index.html` directly in a browser, or use the live link above.

To rebuild the project:

```bash
python analysis.py
python build_app.py
python verify.py
node test_run.js
```

With the project environment:

```bash
.venv/bin/python analysis.py
.venv/bin/python build_app.py
.venv/bin/python verify.py
node test_run.js
```

The dataset at `data/pbmc3k.h5ad` is not included in this repository. It comes from the DDLS course portal.

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

## Limitations

The clusters are very small (13 and 10 cells). The doublet check only tested NK/T versus B-cell/myeloid co-expression. The identities are based on canonical markers, not experimental validation.

## How it was checked

`verify.py` and `test_run.js` reproduce the app's numbers in plain Python and Node. By inspecting the app, I found and fixed shifted cell-type labels, an empty UMAP, and a broken Run button. `results/findings.json` is the single source of truth for suggested cell types.

## AI disclosure

Pi (gpt-5.6-luna via the DDLS portal) wrote the code. Claude helped interpret marker genes and guide the workflow. I reviewed the specification and checked the app and every number myself.

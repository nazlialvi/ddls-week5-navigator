# PBMC Single-Cell Data Navigator

A browser-based navigator for exploring a processed PBMC single-cell RNA-seq dataset.

## The question

The data owner wanted to remove cluster 6 as junk because of its low gene count, and was considering funding follow-up sequencing for cluster 7 because it was small, distinct, and had high counts.

## The answer

Cluster 6 is **platelets**: its top markers include **PF4, PPBP, GNG11, SDPR, and NRGN**; `PPBP` is expressed in **13/13 cells**, and its median mitochondrial percentage is **1.57%**, the lowest cluster value, so these are not dying cells. Their low RNA signal is expected because platelets have no nucleus. **Keep cluster 6.** Cluster 7 is **proliferating NK/T-like cells**: its top markers are mostly housekeeping and cell-cycle genes; `MKI67` is expressed in **70%** of cells and `NKG7` in **80%**. **0 cells** co-express a strong NK/T marker with a strong B-cell or myeloid marker, so the group does not show the tested doublet signature. This is a known cell state rather than a novel cell type, so **do not fund follow-up sequencing**.

## What the navigator does

- Displays a UMAP of all cells, coloured by Leiden cluster, QC metric, or curated gene.
- Provides a summary table for every cluster with cell counts, markers, suggested type, and QC medians.
- Provides a cluster picker and **Run** button to highlight a cluster and show its markers, QC, and curated-gene expression percentages.

## How to run

Open `index.html` directly in a browser. No local server is required; Plotly is loaded from a CDN.

Live app: LINK

The dataset is not stored in this repository. Obtain `data/pbmc3k.h5ad` from the DDLS course portal.

To rebuild the analysis outputs and app from the dataset:

```bash
python analysis.py
python build_app.py
python verify.py
```

Use the project environment when available:

```bash
.venv/bin/python analysis.py
.venv/bin/python build_app.py
.venv/bin/python verify.py
```

## Checking the numbers

`verify.py` reproduces the app's cluster counts, QC medians, marker summaries, and selected expression percentages in plain Python. A bug involving shifted cell-type labels was found by inspecting the app and fixed. `results/findings.json` is now the single source of truth for the suggested cell-type mapping, and the CSV and app use that mapping.

## AI disclosure

Pi (gpt-5.6-luna via the DDLS portal) wrote the code. Claude was used to help interpret marker genes and guide the workflow. I reviewed the specification and checked the app and numbers myself.

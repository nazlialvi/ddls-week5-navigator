# Specification

## Decision

The owner needs two practical decisions before the figure deadline and next sequencing budget are committed:

1. **Cluster 6:** retain the 13 cells or remove them as junk.
2. **Cluster 7:** fund or do not fund another sequencing run for its 10 cells.

This work must first settle the cluster 6 question. The owner believes cluster 6 is junk/dead or empty material because it has about 350 genes detected per cell. The owner is also considering whether cluster 7 is unusual enough to justify follow-up sequencing because it is a tiny, apparently distinctive group with about 2,363 genes per cell. The transcript does not provide a specific biological claim that cluster 6 or 7 is a novel cell type; it says the owner wants the evidence checked rather than assumed.

## Dataset

The file is `data/pbmc3k.h5ad`, a processed single-cell RNA-seq dataset of human PBMCs.

The inspected dataset contains:

- 2,700 cells
- 13,714 genes
- 8 Leiden clusters, labelled `0` through `7`
- one row per cell in the AnnData observation matrix
- genes represented in the variables dimension
- `.obs` columns: `n_genes`, `total_counts`, `pct_mito`, `leiden`
- `.layers['counts']`: raw UMI counts
- `.obsm['X_umap']`: UMAP coordinates
- `.X`: log-normalised expression

The cluster sizes are:

| Cluster | Cells | Genes per cell | Mitochondrial signal |
|---|---:|---:|---:|
| 0 | 1,197 | 809 | 1.8% |
| 1 | 489 | 850 | 2.3% |
| 2 | 445 | 828 | 2.3% |
| 3 | 347 | 673 | 2.1% |
| 4 | 163 | 1,263 | 2.4% |
| 5 | 36 | 1,570 | 2.0% |
| 6 | 13 | 350 | 1.6% |
| 7 | 10 | 2,363 | 2.0% |

The transcript and file do not establish cell-type identities for the numbered clusters.

## Evidence that would settle the decision

- Top marker genes for the relevant cluster, reported before naming or judging it.
- Comparison of `n_genes`, total counts, and `% mito` with the other clusters.
- Whether the gene pattern is coherent rather than blank or empty.
- Co-expression of two lineage programs, checked using the raw counts in `.layers['counts']`, to identify possible mixed-lineage/doublet evidence.
- For cluster 7, the UMAP location, counts, QC, and top genes are relevant to the funding decision.

## Traps

- Low gene count does not by itself mean junk.
- High counts or many detected genes may indicate a doublet rather than a novel population.
- Sitting apart on the UMAP does not by itself establish a novel cell type.
- The owner's cluster number may be slightly off, so verify the relevant cells/cluster rather than assuming the number is exact.
- Cluster numbers are labels, not identities.
- A small cluster should not be treated as reliable without checking coherent expression and QC.

## Done

Done means reporting the relevant cluster's identity only after the required evidence, its top marker genes, and its QC values compared with the other clusters, then stating the cell type those markers point to. For the owner's immediate call, the output must also make an explicit retain/remove decision for cluster 6 and, subsequently, a fund/do-not-fund decision for cluster 7, with concise evidence-based reasons.

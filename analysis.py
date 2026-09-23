from pathlib import Path
import json
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

DATA = 'data/pbmc3k.h5ad'
OUT = Path('results'); OUT.mkdir(exist_ok=True)
adata = sc.read_h5ad(DATA)
cluster = adata.obs['leiden'].astype(str)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=8, use_raw=False)
marker_rows=[]
names=adata.uns['rank_genes_groups']['names']
for c in sorted(cluster.unique(), key=int):
    marker_rows.append({'cluster': c, 'top8_markers': ';'.join(names[c].tolist())})
marker_table=pd.DataFrame(marker_rows)
marker_table['suggested_cell_type']=['CD14+ monocytes','B cells','NK cells','dendritic cells','FCGR3A+ monocytes','megakaryocytes/platelets','platelets','proliferating NK/T-like cells'][0:8]
marker_table.to_csv(OUT/'marker_table.csv', index=False)

def vec(g):
    x=adata.layers['counts'][:, list(adata.var_names).index(g)]
    return np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()
def stats_excluding(g, excluded='7'):
    x=vec(g); rows=[]
    for c in sorted(cluster.unique(), key=int):
        if c==excluded: continue
        z=x[cluster.eq(c).to_numpy()]
        rows.append((c,float(z.mean()),float(z.max())))
    return rows
markers=['CD3E','NKG7','MS4A1','CD14','LYZ','FCGR3A','PPBP','FCER1A']
comparison={g:stats_excluding(g) for g in markers}
cell_markers=['MKI67','TOP2A','TYMS','GNLY','CD8A']
mask=cluster.eq('7').to_numpy()
per_cell=pd.DataFrame({g:vec(g)[mask].astype(int) for g in cell_markers}, index=adata.obs_names[mask])
all_lineage=['CD3E','NKG7','MS4A1','CD14','LYZ','FCGR3A','PPBP','FCER1A']
per_cell=per_cell.join(pd.DataFrame({g:vec(g)[mask].astype(int) for g in all_lineage}, index=adata.obs_names[mask]))
# Conservative comparable-level calls: NK/T strong = NKG7 >= 10 or CD3E >= 5; B strong = MS4A1 >= 10; myeloid strong = CD14 >= 5 or LYZ >= 20 or FCGR3A >= 5.
strong_nkt=(per_cell.NKG7>=10)|(per_cell.CD3E>=5)
strong_b=(per_cell.MS4A1>=10)
strong_my=(per_cell.CD14>=5)|(per_cell.LYZ>=20)|(per_cell.FCGR3A>=5)
median_total=adata.obs.loc[mask,'total_counts'].median()
# NK/T comparator: clusters 0 and 2, the canonical T/NK-rich groups by marker levels.
nkt=cluster.isin(['0','2'])
nkt_median=adata.obs.loc[nkt,'total_counts'].median()
findings={
 'corrected_comparison_excluding_cluster_7': comparison,
 'cluster_7_per_cell_raw_counts': per_cell.reset_index(names='cell').to_dict(orient='records'),
 'cluster_7_median_total_counts': float(median_total),
 'nkt_clusters_used_for_doublet_depth_comparison':['0','2'],
 'nkt_median_total_counts':float(nkt_median),
 'two_times_nkt_median_total_counts':float(2*nkt_median),
 'cells_with_strong_nkt_and_b_marker':int((strong_nkt&strong_b).sum()),
 'cells_with_strong_nkt_and_myeloid_marker':int((strong_nkt&strong_my).sum()),
 'thresholds_used':{'strong_nkt':'NKG7 >= 10 or CD3E >= 5','strong_b':'MS4A1 >= 10','strong_myeloid':'CD14 >= 5 or LYZ >= 20 or FCGR3A >= 5'},
 'identity':'NK/T-like, with proliferation-associated signal but no convincing cluster-wide doublet pattern',
 'confidence':'moderate-low due to n=10 and housekeeping/cell-cycle dominated markers'
}
(OUT/'findings.json').write_text(json.dumps(findings, indent=2))
print(marker_table.to_string(index=False))
print('\nCORRECTED COMPARISON (excluding cluster 7)')
for g, rows in comparison.items():
    best=max(rows,key=lambda x:x[1]); print(g, 'highest=',best[0], 'mean=',round(best[1],3),'max=',best[2])
print('\nPER CELL CLUSTER 7')
print(per_cell.to_string())
print('\nDEPTH',median_total,nkt_median,2*nkt_median)
print('CO-EXPRESS',int((strong_nkt&strong_b).sum()),int((strong_nkt&strong_my).sum()))
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
ad=sc.read_h5ad('data/pbmc3k.h5ad'); c=ad.obs.leiden.astype(str); markers=pd.read_csv('results/marker_table.csv')
type_map={'0':'CD4 T cells','1':'CD14+ monocytes','2':'NK / CD8 T cells','3':'B cells','4':'FCGR3A+ monocytes','5':'dendritic cells','6':'platelets','7':'proliferating NK/T-like cells'}
for cl in sorted(c.unique(),key=int):
 m=c.eq(cl); r=markers[markers.cluster.astype(str)==cl].iloc[0]
 print(f'Cluster {cl}: cells={int(m.sum())}, suggested type={type_map[cl]}')
 print('QC medians:', *(f'{k}={ad.obs.loc[m,k].median()}' for k in ['n_genes','total_counts','pct_mito']))
for cl in ['6','7']:
 m=c.eq(cl); r=markers[markers.cluster.astype(str)==cl].iloc[0]
 print(f'Cluster {cl}')
 print('top markers:', r.top8_markers)
 print('QC medians:', *(f'{k}={ad.obs.loc[m,k].median()}' for k in ['n_genes','total_counts','pct_mito']))
 for g in ['PPBP','NKG7','MKI67']:
  x=ad.layers['counts'][:,list(ad.var_names).index(g)]; x=x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()
  print(f'% expressing {g}: {100*(x[m.to_numpy()]>0).mean():.1f}%')

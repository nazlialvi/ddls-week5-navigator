from pathlib import Path
import scanpy as sc
import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
from scipy import sparse

app = FastAPI(title='PBMC Data Navigator API')
DATA = Path('data/pbmc3k.h5ad')
MARKERS = Path('results/marker_table.csv')
CURATED = None
AD = None

def load_data():
    global AD, CURATED
    if AD is None:
        AD = sc.read_h5ad(DATA)
        import pandas as pd
        table = pd.read_csv(MARKERS)
        CURATED = sorted(set(g for s in table.top8_markers for g in s.split(';')) | set('PPBP PF4 GP9 NKG7 GNLY CD3E CD8A MS4A1 CD14 LYZ FCGR3A FCER1A CST3 MKI67 TOP2A'.split()))
    return AD

class RunRequest(BaseModel):
    cluster: str

def median(values):
    return float(np.median(np.asarray(values)))

def counts(ad, gene):
    x = ad.layers['counts'][:, list(ad.var_names).index(gene)]
    return np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()

@app.get('/', response_class=HTMLResponse)
def index():
    return Path('index.html').read_text()

@app.post('/run')
def run(req: RunRequest):
    ad = load_data()
    import pandas as pd
    labels = ad.obs['leiden'].astype(str)
    mask = labels.eq(str(req.cluster)).to_numpy()
    if not mask.any():
        raise HTTPException(404, f'Unknown cluster: {req.cluster}')
    table = pd.read_csv(MARKERS)
    row = table[table.cluster.astype(str).eq(str(req.cluster))].iloc[0]
    n = int(mask.sum())
    qc = {}
    for key in ['n_genes', 'total_counts', 'pct_mito']:
        qc[key] = {'cluster': median(ad.obs.loc[mask, key]), 'other': median(ad.obs.loc[~mask, key])}
    enrichment=[]
    for gene in CURATED:
        x=counts(ad,gene)
        enrichment.append({'gene':gene,'cluster_pct':100*float(((x[mask]>0).sum())/n),'other_pct':100*float(((x[~mask]>0).sum())/(~mask).sum())})
    enrichment.sort(key=lambda z:z['cluster_pct']-z['other_pct'], reverse=True)
    return {'cluster':str(req.cluster),'cell_count':n,'suggested_cell_type':str(row.suggested_cell_type),'top_markers':str(row.top8_markers).split(';'),'qc_medians':qc,'curated_gene_enrichment':enrichment[:10]}

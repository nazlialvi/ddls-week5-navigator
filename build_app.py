import json
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np

ad=sc.read_h5ad('data/pbmc3k.h5ad')
markers=pd.read_csv('results/marker_table.csv')
findings=json.loads(Path('results/findings.json').read_text())
type_map=findings['suggested_types']
Path('results/findings.json').write_text(json.dumps(findings, indent=2))
obs=ad.obs
clusters=obs['leiden'].astype(str).tolist()
curated=[]
for s in markers.top8_markers:
    curated += s.split(';')
curated += 'PPBP PF4 GP9 NKG7 GNLY CD3E CD8A MS4A1 CD14 LYZ FCGR3A FCER1A CST3 MKI67 TOP2A'.split()
curated=sorted(set(g for g in curated if g in ad.var_names))
genes=list(ad.var_names); gi={g:i for i,g in enumerate(genes)}
X=ad.layers['counts']
from scipy import sparse
def col(g):
 x=X[:,gi[g]]
 return (x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()).tolist()
data={
 'umap':np.asarray(ad.obsm['X_umap']).round(6).tolist(),
 'clusters':clusters,
 'obs':{k:obs[k].astype(float).round(6).tolist() for k in ['n_genes','total_counts','pct_mito']},
 'genes':{g:col(g) for g in curated},
 'curated_genes':curated,
 'cluster_table':markers.to_dict('records'),
 'findings':findings,
 'suggested_types':type_map,
}
for r in data['cluster_table']:
 r['suggested_cell_type']=type_map[str(r['cluster'])]
 m=obs['leiden'].astype(str).eq(str(r['cluster']))
 r['cell_count']=int(m.sum())
 r['median_n_genes']=float(obs.loc[m,'n_genes'].median())
 r['median_total_counts']=float(obs.loc[m,'total_counts'].median())
 r['median_pct_mito']=float(obs.loc[m,'pct_mito'].median())
html='''<!doctype html><html><head><meta charset="utf-8"><title>PBMC Data Navigator</title><script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script><style>body{font:14px system-ui;margin:24px;color:#222}h1{margin-bottom:4px}.answer{background:#eef7ee;padding:12px;border-left:4px solid #398439}.controls{margin:14px 0;padding:10px;background:#f4f4f4}select,button{padding:5px;margin-right:8px}#plot{height:600px}table{border-collapse:collapse;width:100%;margin-top:12px}th,td{border:1px solid #ccc;padding:5px;text-align:left;font-size:12px}th{background:#eee}.result{padding:10px;background:#fafafa;border:1px solid #ddd}.muted{color:#666}</style></head><body><h1>PBMC Data Navigator</h1><div class="answer"><b>Answer</b><br>Cluster 6 = platelets — keep: coherent platelet markers and low mitochondrial signal support retaining it.<br>Cluster 7 = proliferating NK/T-like cells — do not fund: high-depth, housekeeping/cell-cycle dominated markers without a clean novel lineage or convincing doublet pattern.</div><div class="controls"><label>Colour by <select id="colour"></select></label><label>Cluster <select id="cluster"></select></label><button id="runBtn">Run</button></div><div id="plot"></div><div id="result" class="result"><span class="muted">Choose a cluster and click Run.</span></div><h2>Cluster summary</h2><table id="summary"><thead><tr><th>Cluster</th><th>Cells</th><th>Top 8 markers</th><th>Suggested type</th><th>Median n_genes</th><th>Median total_counts</th><th>Median % mito</th></tr></thead><tbody></tbody></table><script>
const D=__DATA__;
document.addEventListener('DOMContentLoaded',()=>{colour.onchange=()=>plot();document.getElementById('runBtn').addEventListener('click',runCluster);plot();});
const colour=document.getElementById('colour'), cluster=document.getElementById('cluster');
[['leiden','Leiden cluster'],['n_genes','n_genes'],['total_counts','total_counts'],['pct_mito','pct_mito'],...D.curated_genes.map(g=>[g,g])].forEach(([v,t])=>colour.add(new Option(t,v)));
[...Array(8).keys()].forEach(i=>cluster.add(new Option('Cluster '+i,String(i))));
function values(v){return v==='leiden'?D.clusters:D.obs[v]||D.genes[v]}
const clusterColours=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f'];
function plot(selected=null){let c=colour.value, vals=values(c);if(c==='leiden'){let traces=[];[...new Set(D.clusters)].sort((a,b)=>Number(a)-Number(b)).forEach(cl=>{let ids=D.clusters.map((x,j)=>x===cl?j:-1).filter(j=>j>=0);traces.push({x:ids.map(j=>D.umap[j][0]),y:ids.map(j=>D.umap[j][1]),mode:'markers',type:'scattergl',name:'Cluster '+cl,marker:{size:6,color:clusterColours[Number(cl)]},text:ids.map(j=>'cluster '+D.clusters[j]),hovertemplate:'%{text}<extra></extra>'})});Plotly.react('plot',traces,{title:'UMAP',xaxis:{title:'UMAP 1'},yaxis:{title:'UMAP 2'},margin:{t:45},legend:{title:{text:'Leiden cluster'}}});return}if(selected){let other=D.clusters.map((x,j)=>String(x)!==String(selected)?j:-1).filter(j=>j>=0), chosen=D.clusters.map((x,j)=>String(x)===String(selected)?j:-1).filter(j=>j>=0);Plotly.react('plot',[{x:other.map(j=>D.umap[j][0]),y:other.map(j=>D.umap[j][1]),mode:'markers',type:'scattergl',name:'Other cells',marker:{size:5,color:'#d9d9d9'},text:other.map(j=>'cluster '+D.clusters[j]),hovertemplate:'%{text}<extra></extra>'},{x:chosen.map(j=>D.umap[j][0]),y:chosen.map(j=>D.umap[j][1]),mode:'markers',type:'scattergl',name:'Cluster '+selected,marker:{size:10,color:'#d62728'},text:chosen.map(j=>'cluster '+D.clusters[j]),hovertemplate:'%{text}<extra></extra>'}],{title:'UMAP',xaxis:{title:'UMAP 1'},yaxis:{title:'UMAP 2'},margin:{t:45},legend:{title:{text:'Selection'}}});return}let colors=vals;Plotly.react('plot',[{x:D.umap.map(x=>x[0]),y:D.umap.map(x=>x[1]),mode:'markers',type:'scattergl',marker:{size:6,color:colors,colorscale:'Viridis',showscale:true,colorbar:{title:c}},text:D.clusters.map(x=>'cluster '+x),hovertemplate:'%{text}<extra></extra>'}],{title:'UMAP',xaxis:{title:'UMAP 1'},yaxis:{title:'UMAP 2'},margin:{t:45}})}
function med(a){a=[...a].sort((x,y)=>x-y);return a.length%2?a[(a.length-1)/2]:(a[a.length/2-1]+a[a.length/2])/2}
function runCluster(){let box=document.getElementById('result');let c=String(cluster.value);box.innerHTML='<span>Running cluster '+c+'…</span>';try{let ix=D.clusters.map(x=>String(x)===c), r=D.cluster_table.find(x=>String(x.cluster)===c), n=r.cell_count;let qc=['n_genes','total_counts','pct_mito'].map(k=>{let inx=D.obs[k].filter((_,i)=>ix[i]), out=D.obs[k].filter((_,i)=>!ix[i]);return [k,med(inx),med(out)]});let pct=D.curated_genes.map(g=>{let a=100*D.genes[g].filter((_,i)=>ix[i]&&D.genes[g][i]>0).length/n, other=100*D.genes[g].filter((_,i)=>!ix[i]&&D.genes[g][i]>0).length/(D.clusters.length-n);return [g,a,other,a-other]}).sort((a,b)=>b[3]-a[3]).slice(0,10);let rows=pct.map(x=>'<tr><td>'+x[0]+'</td><td>'+x[1].toFixed(1)+'%</td><td>'+x[2].toFixed(1)+'%</td><td>'+x[3].toFixed(1)+'%</td></tr>').join('');document.getElementById('result').innerHTML='<h3>Cluster '+c+' — '+r.suggested_cell_type+'</h3><b>Cells:</b> '+n+'<br><b>Top markers:</b> '+r.top8_markers.replaceAll(';',', ')+'<br><b>QC medians (cluster vs all other cells):</b> '+qc.map(x=>x[0]+'='+ (x[0]==='pct_mito'?x[1].toFixed(2)+'% vs '+x[2].toFixed(2)+'%':Math.round(x[1])+' vs '+Math.round(x[2]))).join('; ')+'<br><b>Curated gene enrichment:</b><table><tr><th>Gene</th><th>Cluster %</th><th>Other %</th><th>Difference</th></tr>'+rows+'</table>';plot(c)}catch(e){box.innerHTML='<span style="color:red">Error: '+e.message+'</span>';}}

let body=document.querySelector('#summary tbody');D.cluster_table.forEach(r=>body.insertAdjacentHTML('beforeend',`<tr><td>${r.cluster}</td><td>${r.cell_count}</td><td>${r.top8_markers.replaceAll(';',', ')}</td><td>${r.suggested_cell_type}</td><td>${Math.round(r.median_n_genes)}</td><td>${Math.round(r.median_total_counts)}</td><td>${r.median_pct_mito.toFixed(2)}</td></tr>`));
</script></body></html>'''
Path('index.html').write_text(html.replace('__DATA__',json.dumps(data,separators=(',',':'))))
print(f'Wrote index.html: {Path("index.html").stat().st_size} bytes')

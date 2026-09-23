const fs=require('fs');
const html=fs.readFileSync('index.html','utf8');
const m=html.match(/const D=(.*?);\ndocument\.addEventListener/s);
if(!m) throw new Error('inline data not found');
const D=JSON.parse(m[1]);
const rFor=c=>D.cluster_table.find(x=>String(x.cluster)===c);
function med(a){a=[...a].sort((x,y)=>x-y);return a.length%2?a[(a.length-1)/2]:(a[a.length/2-1]+a[a.length/2])/2}
for(const c of ['6','7']){const ix=D.clusters.map(x=>String(x)===c), r=rFor(c), n=r.cell_count;
const qc=k=>med(D.obs[k].filter((_,i)=>ix[i]));
const pct=g=>100*D.genes[g].filter((x,i)=>ix[i]&&x>0).length/n;
const diffs=D.curated_genes.map(g=>{const a=pct(g),o=100*D.genes[g].filter((x,i)=>!ix[i]&&x>0).length/(D.clusters.length-n);return [g,a,o,a-o]}).sort((a,b)=>b[3]-a[3]).slice(0,5);
console.log('Cluster '+c,'cell count:',n); console.log('QC medians:',Math.round(qc('n_genes')),Math.round(qc('total_counts')),qc('pct_mito').toFixed(2)); console.log('top-5 by difference:',JSON.stringify(diffs));}

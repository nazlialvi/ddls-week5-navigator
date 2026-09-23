const fs=require('fs');
const html=fs.readFileSync('index.html','utf8');
const m=html.match(/const D=(.*?);\ndocument\.addEventListener/s);
if(!m) throw new Error('inline data not found');
const D=JSON.parse(m[1]);
const c='7', ix=D.clusters.map(x=>String(x)===c), r=D.cluster_table.find(x=>String(x.cluster)===c);
function med(a){a=[...a].sort((x,y)=>x-y);return a.length%2?a[(a.length-1)/2]:(a[a.length/2-1]+a[a.length/2])/2}
const qc=k=>med(D.obs[k].filter((_,i)=>ix[i]));
const pct=g=>100*D.genes[g].filter((x,i)=>ix[i]&&x>0).length/r.cell_count;
console.log('cell count:',r.cell_count);
console.log('QC medians:',Math.round(qc('n_genes')),Math.round(qc('total_counts')),qc('pct_mito').toFixed(2));
console.log('% expressing MKI67:',pct('MKI67').toFixed(0)+'%');
console.log('% expressing NKG7:',pct('NKG7').toFixed(0)+'%');

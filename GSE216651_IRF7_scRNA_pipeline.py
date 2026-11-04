#!/usr/bin/env python3
import os, gzip, time, gc, shutil, urllib.request
from pathlib import Path
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy.stats import ttest_ind, mannwhitneyu
from itertools import combinations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
SAMPLES=[
('GSM6685302','sc_healthy_1','sc','Healthy'),
('GSM6685303','sc_healthy_2','sc','Healthy'),
('GSM6685304','sc_healthy_3','sc','Healthy'),
('GSM6685305','sc_OA_1','sc','OA'),
('GSM6685306','sc_OA_2','sc','OA'),
('GSM6685307','sc_OA_3','sc','OA'),
('GSM6685308','sn_healthy_1','sn','Healthy'),
('GSM6685309','sn_OA_1','sn','OA'),
('GSM6685310','sn_OA_2','sn','OA')]
O=Path('results/GSE216651'); O.mkdir(parents=True,exist_ok=True)
T=Path('/tmp/GSE216651'); T.mkdir(parents=True,exist_ok=True)
mygenes=['LST1','TYROBP','LYZ','FCER1G','CD14','CSF1R','C1QA','C1QB','MS4A7','FCGR3A']
other=['COL1A1','DCN','LUM','PECAM1','VWF','CD3D','CD3E','NKG7','MS4A1','ADIPOQ']
def fetch(url,p):
 if p.exists() and p.stat().st_size>2000:return
 for attempt in range(5):
  try:
   print('Downloading',url,flush=True)
   req=urllib.request.Request(url,headers={'User-Agent':'NCBI GEO scientific reanalysis'})
   with urllib.request.urlopen(req,timeout=180) as r,open(str(p)+'.part','wb') as f:shutil.copyfileobj(r,f)
   os.replace(str(p)+'.part',p)
   print('Downloaded',p.stat().st_size,'bytes',flush=True)
   return
  except Exception as e:
   print('Download retry',attempt,e,flush=True);time.sleep(3*(attempt+1))
 raise RuntimeError('GEO download failed '+url)
rows=[];pb={};pb_my={}
for gsm,name,mode,group in SAMPLES:
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6685nnn/{gsm}/suppl/{gsm}_{name}_'
 files={t:T/f'{gsm}_{name}_{t}.{ext}.gz' for t,ext in [('features','tsv'),('barcodes','tsv'),('matrix','mtx')]}
 for typ,path in files.items():fetch(base+f'{typ}.'+('mtx' if typ=='matrix' else 'tsv')+'.gz',path)
 features=pd.read_csv(files['features'],sep='\t',header=None)
 genes=features.iloc[:,1].astype(str).to_numpy()
 barcodes=pd.read_csv(files['barcodes'],sep='\t',header=None)
 with gzip.open(files['matrix'],'rb') as f:m=mmread(f).tocsr()
 assert m.shape==(len(genes),len(barcodes)),(name,m.shape,len(genes),len(barcodes))
 n_raw=m.shape[1]; totals=np.asarray(m.sum(axis=0)).ravel()
 n_genes=np.asarray(m.getnnz(axis=0)).ravel()
 mt=np.where(np.char.startswith(genes.astype(str),'MT-'))[0]
 mt_counts=np.asarray(m[mt,:].sum(axis=0)).ravel() if len(mt)>0 else np.zeros(m.shape[1])
 mt_frac=np.divide(mt_counts,totals,out=np.zeros(len(totals)),where=totals>0)
 ok=(totals>=500)&(n_genes>=200)&(mt_frac<=0.2)
 print(name,'raw cells',n_raw,'QC',ok.sum(),flush=True)
 if ok.sum()<10:raise ValueError('Too few cells pass QC for '+name)
 m=m[:,ok];totals=np.asarray(m.sum(axis=0)).ravel(); n=m.shape[1]
 idx=np.where(genes=='IRF7')[0]
 if len(idx)!=1:raise ValueError('IRF7 not unique in '+name+': '+str(idx))
 irf=m[idx[0],:].toarray().ravel()
 norm=np.log1p(1e4*irf/np.maximum(1,totals))
 indices={g:int(i) for i,g in enumerate(genes) if g in set(mygenes+other)}
 def marker(g):
  return m[indices[g],:].toarray().ravel() if g in indices else np.zeros(n)
 A=np.array([marker(g) for g in mygenes])
 B=np.array([marker(g) for g in other])
 mark=A>0
 mys=np.log1p(1e4*A/np.maximum(1,totals)).mean(axis=0)
 oth=np.log1p(1e4*B/np.maximum(1,totals)).mean(axis=0)
 gate=(mark.sum(axis=0)>=3)&(mys>oth)&((marker('LST1')+marker('TYROBP'))>0)
 ngate=int(gate.sum()); t_my=int(totals[gate].sum());i_my=int(irf[gate].sum())
 row=dict(GSM=gsm,sample=name,modality=mode,group=group,n_cells_raw=n_raw,
   n_cells_QC=n,IRF7_detected_cells=int((irf>0).sum()),IRF7_detect_pct=100*(irf>0).mean(),
   IRF7_mean_log1p_10k=float(norm.mean()),IRF7_median_log1p_10k=float(np.median(norm)),
   IRF7_pseudobulk_count=int(irf.sum()),pseudobulk_total_UMI=int(totals.sum()),
   IRF7_CPM=float(1e6*irf.sum()/totals.sum()),
   myeloid_proxy_cells=ngate,myeloid_proxy_IRF7_count=i_my,
   myeloid_proxy_total_UMI=t_my,
   myeloid_proxy_IRF7_CPM=float(1e6*i_my/t_my) if t_my else np.nan,
   myeloid_proxy_IRF7_detect_pct=float(100*(irf[gate]>0).mean()) if ngate else np.nan)
 rows.append(row)
 print('RESULT',name,row['IRF7_CPM'],row['IRF7_detect_pct'],row['myeloid_proxy_cells'],flush=True)
 if mode=='sc':
  whole=np.asarray(m.sum(axis=1)).ravel().astype(np.int64)
  sub=np.asarray(m[:,gate].sum(axis=1)).ravel().astype(np.int64)
  pb[name]=pd.Series(whole,index=genes).groupby(level=0).sum()
  pb_my[name]=pd.Series(sub,index=genes).groupby(level=0).sum()
 pd.DataFrame(rows).to_csv(O/'GSE216651_IRF7_partial.csv',index=False)
 for f in files.values():f.unlink()
 del features,barcodes,m,irf,norm,A,B
 gc.collect()
D=pd.DataFrame(rows);D.to_csv(O/'GSE216651_IRF7_9_donors.csv',index=False)
for tab,what in [(pb,'all_cells'),(pb_my,'myeloid_marker_proxy')]:
 matrix=pd.concat(tab,axis=1).fillna(0).astype(np.int64)
 matrix.index.name='gene'
 matrix.to_csv(O/f'GSE216651_sc_pseudobulk_{what}.csv.gz',compression='gzip')
stats=[]
for mode in ['sc','sn']:
 sub=D[D.modality==mode]
 for variable in ['IRF7_CPM','IRF7_detect_pct','IRF7_mean_log1p_10k',
                  'myeloid_proxy_IRF7_CPM','myeloid_proxy_IRF7_detect_pct']:
  a=sub.loc[sub.group=='Healthy',variable].dropna().to_numpy()
  b=sub.loc[sub.group=='OA',variable].dropna().to_numpy()
  if len(a)==0 or len(b)==0:continue
  x=dict(modality=mode,measure=variable,n_healthy=len(a),n_OA=len(b),
   healthy_mean=float(a.mean()),OA_mean=float(b.mean()),
   OA_div_healthy=float(b.mean()/a.mean()) if a.mean()!=0 else np.nan,
   welch_p=float(ttest_ind(a,b,equal_var=False).pvalue) if len(a)>=2 and len(b)>=2 else np.nan,
   MW_p=float(mannwhitneyu(a,b,alternative='two-sided').pvalue) if len(a)>=2 and len(b)>=2 else np.nan)
  if len(a)==3 and len(b)==3:
   v=np.r_[a,b];delta=abs(a.mean()-b.mean())
   permutations=[]
   for ids in combinations(range(6),3):
    chosen=np.array(ids);op=np.array([i for i in range(6) if i not in ids])
    permutations.append(abs(v[chosen].mean()-v[op].mean()))
   x['exact_permutation_p']=float(np.mean(np.array(permutations)>=delta-1e-12))
  stats.append(x)
P=pd.DataFrame(stats); P.to_csv(O/'GSE216651_IRF7_donor_level_tests.csv',index=False)
for field,label in [('IRF7_CPM','All cells, pseudobulk CPM'),
                    ('IRF7_detect_pct','IRF7-detecting cells (%)'),
                    ('myeloid_proxy_IRF7_CPM','Marker-enriched myeloid proxy CPM')]:
 fig,axes=plt.subplots(1,2,figsize=(10,4.5))
 for ax,mode in zip(axes,['sc','sn']):
  sub=D[D.modality==mode]
  for x,group in enumerate(['Healthy','OA']):
   vals=sub.loc[sub.group==group,field].to_numpy()
   xs=np.linspace(x-.07,x+.07,len(vals)) if len(vals)>1 else np.array([x])
   ax.scatter(xs,vals,s=70,zorder=3)
   ax.hlines(np.nanmean(vals),x-.17,x+.17,colors='black',linewidth=2)
  ax.set_xticks([0,1],['Healthy','OA'])
  ax.set_xlim(-.4,1.4);ax.set_ylabel(label)
  ax.set_title('scRNA-seq (n=3 vs 3)' if mode=='sc' else 'snRNA-seq (n=1 vs 2)')
  ax.spines[['top','right']].set_visible(False)
 fig.suptitle('GSE216651 | IRF7 | each dot = one patient')
 fig.tight_layout()
 fig.savefig(O/('GSE216651_IRF7_'+field+'.png'),dpi=250)
 plt.close(fig)
print('FINAL DONORS\n',D.to_string(index=False),flush=True)
print('FINAL STATS\n',P.to_string(index=False),flush=True)

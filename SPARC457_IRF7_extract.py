import requests,io,zipfile,os,pathlib,pandas as pd, json,traceback
O=pathlib.Path("results/SPARC457_IRF7");O.mkdir(parents=True,exist_ok=True)
url="https://api.pennsieve.io/discover/datasets/457/versions/1/download"
r=requests.get(url,timeout=180)
print("DOWNLOAD HTTP",r.status_code,"size",len(r.content),flush=True);r.raise_for_status()
z=zipfile.ZipFile(io.BytesIO(r.content))
members=z.namelist()
print("FILES IN ZIP",len(members),flush=True)
for f in members:print("FILE",f,flush=True)
res=[]
for name in members:
 if not name.lower().endswith(('.xlsx','.xls','.csv','.tsv','.txt')):continue
 try:
  if name.lower().endswith(('.xlsx','.xls')):
   ex=pd.ExcelFile(io.BytesIO(z.read(name)))
   for sh in ex.sheet_names:
    try:
     d=pd.read_excel(ex,sheet_name=sh,dtype=str)
     mask=d.fillna('').astype(str).apply(lambda s:s.str.contains(r'(?i)(?<![a-z0-9])IRF7(?![a-z0-9])',regex=True)).any(axis=1)
     hits=d[mask]
     print("SHEET",name,sh,"shape",d.shape,"matches",len(hits),flush=True)
     for i,row in hits.iterrows():
      record={"file":name,"sheet":sh,"row_excel":int(i)+2}
      record.update({str(k):str(v) for k,v in row.items() if pd.notna(v)})
      res.append(record);print("IRF7_ROW",json.dumps(record,ensure_ascii=False)[:9000],flush=True)
    except Exception as e:print("SHEET_ERROR",name,sh,str(e)[:300],flush=True)
  else:
   blob=z.read(name)
   if len(blob)>60000000:continue
   for sep in [',','\t']:
    try:
     d=pd.read_csv(io.BytesIO(blob),sep=sep,dtype=str,low_memory=False)
     if d.shape[1]<2:continue
     mask=d.fillna('').astype(str).apply(lambda s:s.str.contains(r'(?i)(?<![a-z0-9])IRF7(?![a-z0-9])',regex=True)).any(axis=1)
     hits=d[mask]
     print("TABULAR",name,"shape",d.shape,"matches",len(hits),flush=True)
     for i,row in hits.iterrows():
      record={"file":name,"sheet":"", "row_excel":int(i)+2}
      record.update({str(k):str(v) for k,v in row.items() if pd.notna(v)})
      res.append(record);print("IRF7_ROW",json.dumps(record,ensure_ascii=False)[:9000],flush=True)
     break
    except Exception:pass
 except Exception as e:print("FILE_ERROR",name,str(e)[:500],flush=True)
pd.DataFrame(res).to_csv(O/"SPARC457_IRF7_all_table_rows.csv",index=False)
(O/"SPARC457_IRF7_rows.json").write_text(json.dumps(res,ensure_ascii=False,indent=2))
print("TOTAL IRF7 rows",len(res),flush=True)

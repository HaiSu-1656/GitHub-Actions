import requests, json, pathlib, sys, time
p=pathlib.Path("results/SPARC457_api_probe")
p.mkdir(parents=True,exist_ok=True)
requests_args=dict(timeout=35,headers={"Accept":"application/json","User-Agent":"research-reanalysis/1.0"})
urls=[
("dataset","https://api.pennsieve.io/discover/datasets/457"),
("versions","https://api.pennsieve.io/discover/datasets/457/versions"),
("browse","https://api.pennsieve.io/discover/datasets/457/versions/1/files/browse?limit=200&offset=0"),
("download","https://api.pennsieve.io/discover/datasets/457/versions/1/download"),
]
for name,url in urls:
 try:
  r=requests.get(url,**requests_args)
  print(f"PROBE {name} HTTP {r.status_code} {len(r.content)} {r.url[:300]}",flush=True)
  (p/(name+".txt")).write_text(r.text[:300000],encoding="utf-8")
  if r.headers.get("Content-Type","").startswith("application/json"):
   try:
    dat=r.json()
    print(f"JSON {name} summary: "+str(dat)[:1800],flush=True)
   except Exception as ex:print("json failure",str(ex),flush=True)
  else:print(f"TEXT {name}: {r.text[:300]}",flush=True)
 except Exception as ex:print(f"ERROR {name}: {str(ex)[:400]}",flush=True)
try:
 url="https://downloads.pennsieve.io/public/manifests"
 for args in [{"datasetId":457,"version":"latest","paths":[""],"limit":1000},{"datasetId":457,"version":1,"paths":[""],"limit":1000}]:
  r=requests.post(url,json=args,timeout=40)
  print("MANIFEST HTTP",r.status_code,len(r.content),r.text[:1400],flush=True)
  if r.ok:(p/"manifest.json").write_text(r.text,encoding="utf-8");break
except Exception as ex:print("MANIFEST ERROR",str(ex),flush=True)

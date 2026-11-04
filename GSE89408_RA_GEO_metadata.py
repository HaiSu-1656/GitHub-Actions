#!/usr/bin/env python3
import csv, gzip, io, re, urllib.request, os
from pathlib import Path
out=Path("results/GSE89408_RA")
out.mkdir(parents=True,exist_ok=True)
url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89408/matrix/GSE89408_series_matrix.txt.gz"
try:
    req=urllib.request.Request(url,headers={"User-Agent":"GSE89408-RA-sample-reanalysis"})
    with urllib.request.urlopen(req,timeout=120) as r:
        data=r.read()
    text=gzip.decompress(data).decode("utf-8",errors="replace")
    metadata={}
    for line in text.splitlines():
        if not line.startswith("!Sample_"):continue
        toks=next(csv.reader([line],delimiter="\t",quotechar='"'))
        if len(toks)<2:continue
        metadata.setdefault(toks[0],[]).append(toks[1:])
    names=metadata.get("!Sample_title",[[]])[0]
    ids=metadata.get("!Sample_geo_accession",[[]])[0]
    chars=[]
    for k,v in metadata.items():
        if k.startswith("!Sample_characteristics_ch1"):chars+=v
    if not names or not ids or len(names)!=len(ids):
        raise ValueError("No reliable title/GSM mapping")
    rows=[]
    for i,(name,gsm) in enumerate(zip(names,ids)):
        char=" | ".join(v[i] for v in chars if len(v)>i)
        group="RA" if "rheumatoid arthritis" in name.lower() else ("Healthy" if "healthy" in name.lower() else ("OA" if "osteoarthritis" in name.lower() else "Other"))
        stage=""
        if group=="RA":
            s=" ".join((name,char)).lower()
            if re.search(r"\b(?:early\s*(?:ra|rheumatoid arthritis)|(?:ra|rheumatoid arthritis)\s*early)\b",s):
                stage="Early"
            if re.search(r"\b(?:established\s*(?:ra|rheumatoid arthritis)|(?:ra|rheumatoid arthritis)\s*established)\b",s):
                stage="Established" if stage=="" else "Ambiguous"
        m=re.search(r"(?:rheumatoid arthritis|healthy|osteoarthritis) tissue\s*(\d+)",name,re.IGNORECASE)
        prefix="RA_tissue_" if group=="RA" else ("normal_tissue_" if group=="Healthy" else ("OA_tissue_" if group=="OA" else ""))
        matrix_col=prefix+m.group(1) if m and prefix else ""
        rows.append(dict(GSM=gsm,title=name,group=group,stage=stage,matrix_col=matrix_col,characteristics=char))
    import pandas as pd
    d=pd.DataFrame(rows)
    d.to_csv(out/"GSE89408_GEO_sample_metadata.csv",index=False)
    print("METADATA rows",len(d),flush=True)
    print("Groups",d.group.value_counts().to_dict(),flush=True)
    print("RA stage distribution",d.loc[d.group=="RA","stage"].value_counts().to_dict(),flush=True)
    for _,r in d[d.group=="RA"].iloc[[0,56,57,151]].iterrows():
        print("RA_CHECK",r['title'],r['GSM'],r['stage'],r['characteristics'][:400],flush=True)
    print("Metadata labels:",list(metadata.keys()),flush=True)
except Exception as ex:
    print("METADATA_UNAVAILABLE",str(ex),flush=True)

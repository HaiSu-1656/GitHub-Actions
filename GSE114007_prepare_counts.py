#!/usr/bin/env python3
from pathlib import Path
from openpyxl import load_workbook
import csv,urllib.request,hashlib
src=Path("GSE114007_raw_counts.xlsx")
if not src.exists():
    urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114007/suppl/GSE114007_raw_counts.xlsx",src)
out=Path("results/GSE114007");out.mkdir(parents=True,exist_ok=True)
wb=load_workbook(src,read_only=True,data_only=True)
assert set(wb.sheetnames)=={"Normal","OA"},wb.sheetnames
heads,sets={},{}
for name in ("Normal","OA"):
    sh=wb[name]; rows=sh.iter_rows(values_only=True);header=list(next(rows))
    assert header[0]=="symbol"
    assert len(header)==(19 if name=="Normal" else 21)
    heads[name]=header[1:]
    genes=[]; vals=[]
    for row in rows:
        if row[0] is None: continue
        genes.append(str(row[0]).strip())
        vals.append([int(v) if isinstance(v,(int,float)) and v==int(v) else v for v in row[1:]])
    assert len(genes)==23710,(name,len(genes))
    assert len(genes)==len(set(genes)),f"duplicate symbols in {name}"
    for v in vals:
        assert all(isinstance(x,int) and x>=0 for x in v),repr(v[:4])
    sets[name]=(genes,vals)
    print(name,len(header)-1,"sample columns",len(genes),"genes",flush=True)
assert sets["Normal"][0]==sets["OA"][0],"Gene symbol lists differ"
normal_cols=heads["Normal"];oa_cols=heads["OA"]
assert sum(s.startswith("Normal_Cart_") for s in normal_cols)==8
assert sum(s.startswith("normal_") for s in normal_cols)==10
assert sum(s.startswith("OA_Cart_") for s in oa_cols)==10
assert sum(s.startswith("OA_") for s in oa_cols)==20
# Paper's two platforms GPL11154, GPL18573; sample mapping from GSM3130531-48 and GSM3130549-68.
metadata=[]
for group,cols in [("Normal",normal_cols),("OA",oa_cols)]:
    for col in cols:
        platform="HiSeq2000" if col.startswith(("Normal_Cart_","OA_Cart_")) else "NextSeq500"
        metadata.append((col,group,platform))
with (out/"sample_metadata.csv").open("w",newline="") as f:
    wr=csv.writer(f); wr.writerow(("sample","condition","platform"));wr.writerows(metadata)
with (out/"all_counts.tsv").open("w",newline="") as f:
    wr=csv.writer(f,delimiter="\t")
    wr.writerow(["symbol"]+normal_cols+oa_cols)
    for i,g in enumerate(sets["Normal"][0]):
        wr.writerow([g]+sets["Normal"][1][i]+sets["OA"][1][i])
i=sets["Normal"][0].index("IRF7")
print("IRF7_NORMAL",sets["Normal"][1][i],flush=True)
print("IRF7_OA",sets["OA"][1][i],flush=True)
print("Counts SHA256",hashlib.sha256((out/"all_counts.tsv").read_bytes()).hexdigest(),flush=True)

import urllib.request, os, hashlib, json, zipfile
from openpyxl import load_workbook
from pathlib import Path
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114007/suppl/GSE114007_raw_counts.xlsx"
p = Path("GSE114007_raw_counts.xlsx")
if not p.exists():
    urllib.request.urlretrieve(url, p)
print("GSE114007 file bytes",p.stat().st_size,"MD5",hashlib.md5(p.read_bytes()).hexdigest(),flush=True)
wb=load_workbook(p,read_only=True,data_only=True)
out=Path("results/GSE114007")
out.mkdir(parents=True,exist_ok=True)
report=[]
for sh in wb.worksheets:
    print("SHEET:",repr(sh.title),"ROWS:",sh.max_row,"COLS:",sh.max_column,flush=True)
    rows=sh.iter_rows(values_only=True)
    for i,r in enumerate(rows):
        if i<6: print(" ROW",i+1,repr(list(r)[:48]),flush=True)
        if i<10: report.append({"sheet":sh.title,"row":i+1,"values":[str(x) if x is not None else "" for x in list(r)[:50]]})
        if any(str(v).strip().lower()=="irf7" for v in r[:4] if v is not None) or any(str(v).strip()=="ENSG00000185507" for v in r[:4] if v is not None):
            print(" IRF7 ROW:",i+1,repr(list(r)[:50]),flush=True)
            report.append({"sheet":sh.title,"IRF7_row":i+1,"values":[str(x) if x is not None else "" for x in list(r)[:50]]})
    print("DONE SHEET",sh.title,flush=True)
(out/"GSE114007_probe.json").write_text(json.dumps(report,indent=2,ensure_ascii=False))

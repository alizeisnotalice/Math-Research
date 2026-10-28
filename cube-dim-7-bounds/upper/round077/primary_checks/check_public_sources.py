from pathlib import Path
import hashlib,json,logging,os
logging.disable(logging.CRITICAL)
from pypdf import PdfReader
root=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/revision_v3/literature/originals')
items=[('stein-stromberg1983.pdf',[0,1,2],['Theorem 1']),('Aldaz2017_1605.05596v2.pdf',[13],['Theorem 5.2','(2 + 3ε)']),('aldaz2020v2.pdf',[10],['Corollary 3.5','uniform bound'])]
out=[]
for name,ix,needles in items:
 p=root/name;r=PdfReader(p);t='\n'.join(r.pages[i].extract_text() for i in ix)
 found={q:q in t for q in needles};assert all(found.values()),(name,found)
 out.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'pdf_pages_read':[i+1 for i in ix],'checked_identifiers':found})
receipt={'status':'PASS','pid':os.getpid(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sources':out,'scope':'Read identified original PDF pages; Aldaz 2020 sharpness applies to the uniform measure-class bound, while the manuscript uses its Lebesgue upper-bound consequence. No mathematical baseline or numerical search rerun.'}
Path(__file__).with_name('check_public_sources.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n');print(json.dumps(receipt,ensure_ascii=False))

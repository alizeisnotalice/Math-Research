from pathlib import Path
from fractions import Fraction as F
import hashlib,json
p=Path(__file__).resolve().parent
out=p/'actual_face_mask_absorption_results_20261007.json';assert not out.exists()
reg=p/'actual_face_mask_absorption_registration_20261007.json'
r=json.loads(reg.read_text());rows=[];checks=0
for n in r['rounds']+r['boundary_dimensions']:
 M=n//65536;a=F(32*M,3*n+5*M);k=F(8192,49)*a
 conditions=[0<=M<=F(n,65536),0<=a<=F(1,6144),k<=F(4,147),1-k>=F(143,147),1/(1-k)<=F(147,143), (M==0)==(n<65536)]
 assert all(conditions);checks+=len(conditions)
 rows.append({'n':n,'M0':M,'alpha':str(a),'absorbed':str(k),'exact_ledger_multiplier':str(1/(1-k))})
res={'status':'PASS_EXACT','checks':checks,'rows':rows,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),'scope':r['scope']}
out.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({'status':res['status'],'checks':checks,'rows':rows}))

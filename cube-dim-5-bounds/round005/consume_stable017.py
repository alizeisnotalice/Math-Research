import sys,json,hashlib,time
from pathlib import Path
from fractions import Fraction
from flint import ctx
p=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(p/'audit_stable017'))
from hermite_pilot import evaluate_certificate
ctx.prec=256
start=time.monotonic(); source=p/'audit_stable017/hermite_results.json'; d=json.loads(source.read_text()); rows=[]
for rec in d['intervals']:
 a,b=Fraction(rec['a']),Fraction(rec['b']);h=b-a
 for l,r in [(Fraction(1,10),Fraction(3,10)),(Fraction(6,10),Fraction(9,10))]:
  x,y=a+l*h,a+r*h; v=evaluate_certificate(rec,x,y)
  assert v.is_finite() and v>0
  rows.append({'x':[str(x),str(y)],'density':v.str(60,more=True)})
 try: evaluate_certificate(rec,a-h,a)
 except ValueError: pass
 else: raise AssertionError('out-of-domain input accepted')
points=[evaluate_certificate(d['intervals'][0],'0.003'),evaluate_certificate(d['intervals'][2],'2.007')]
out={'status':'consumer_checks_passed_not_global_certificate','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rows':rows,'matlab_check_points':[v.str(60,more=True) for v in points],'elapsed_seconds':time.monotonic()-start,'scope':'12 new subinterval API checks plus 6 domain guards; no reintegration'}
(p/'round005/consume_stable017.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))

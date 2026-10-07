from fractions import Fraction as F
from pathlib import Path
import json,hashlib
p=Path(__file__).parent
checks=0;rows=[]
def log_unit(x):
 assert 1<=x<=2
 z=(x-1)/(x+1);s=sum((2*z**(2*j+1)/F(2*j+1) for j in range(48)),F(0))
 tail=2*z**97/(97*(1-z*z))
 return s,s+tail
def log_pos(x):
 if x<=1:return F(0),F(0)
 j=0
 while x>2:x/=2;j+=1
 a,b=log_unit(x);c,d=log_unit(F(2));return a+j*c,b+j*d
for m in [8,32,128]:
 start=checks
 profiles=[[F(2**m)]+[F(0)]*(m-1),[F(m,i+1) for i in range(m)],
 [F(1,4) if i%2 else F(16) for i in range(m)],
 [F(2**(i%12),1+i//12) for i in range(m)]]
 for typ,values in enumerate(profiles):
  weights=[F(1+(i%3),m) for i in range(m)]
  ordered=sorted(zip(values,weights),reverse=True);cum=F(0);A=F(0)
  for b,w in ordered:cum+=w;A=max(A,b*cum)
  measure=sum(weights,F(0))
  for cap in [F(1,16),F(1),F(4),F(64)]:
   exact=sum((w*min(cap,b) for b,w in zip(values,weights)),F(0))
   q=cap*measure/A;lo,hi=log_pos(q)
   assert exact<=A*(1+lo);checks+=1
   # Sharp discrete layercake identity from breakpoints.
   knots=sorted(set([F(0),cap]+[b for b in values if 0<b<cap]))
   lc=sum(((t-s)*sum((w for b,w in zip(values,weights) if b>=(s+t)/2),F(0)) for s,t in zip(knots,knots[1:])),F(0))
   assert lc==exact;checks+=1
   for C in [F(1),F(8192,49),F(147,143)*F(8192,49)]:
    V=C*A;X=measure # lambda=W=1
    # Choosing U with a lower log bound certifies the premise X<=U+V log+(cap X/A).
    U=max(F(0),X-V*lo)
    Klo,Khi=log_pos(2*cap*C)
    assert X<=2*U+2*V*Klo;checks+=1
  rows.append({'round':m,'profile':typ,'weak_A':str(A),'measure':str(measure)})
 rows.append({'round':m,'checks':checks-start})
reg=p/'capped_reference_absorption_registration_20261007.json'
r={'status':'PASS','checks':checks,'rows':rows,'scope':'Scalar budget and finite distribution audit only; no actual reference weak bound proved','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest()}
(p/'capped_reference_absorption_results_20261007.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'status':r['status'],'checks':checks,'rounds':[x for x in rows if 'checks'in x]}))

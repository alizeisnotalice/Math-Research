from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time
import mpmath as mp
import phi_diagnostic as x
p=Path(__file__).resolve().parent;mp.mp.dps=160;start=time.monotonic();rows=[]
for i in range(256):
 q=F(2,3)+(F(917,1000)-F(2,3))*F(2*i+1,512);r=F(3,2)*q;h=r/F(2**48);eps=h/(2*F(3,2)*q);a=(1-eps)*q;R=F(3,2)*q-h;g=R-a
 for fac in (F(4),F(2**256)):
  rho=fac*r;t=rho/q;RR=x.v(R);gg=x.v(g)
  def E(m):
   z=gg-x.v(m)
   return -2*RR*mp.log1p(-z/(2*RR))-z if z>0 else mp.mpf(0)
  for label,w in (('zero',F(0)),('kink',a),('boundary',r)):
   mm=max(F(0),abs(w)-a);B=E(F(0))+E(mm)-x.v(mm);J=1+B/(2*x.v(rho));phi,_=x.phi(q,h,eps,t,w);ratio=phi/J
   assert J>0 and ratio>=0
   rows.append({'i':i,'rho_over_r':str(fac),'point':label,'q':str(q),'rho_scaled_Jphi_over_J':mp.nstr(x.v(rho)*ratio,45)})
out={'scope':'post-run finite H/J diagnostic, not original production run or operator theorem','states':len(rows),'parameters':'256 q midpoints, h/r=2^-48, epsilon factor 1/2, rho/r in {4,2^256}, w in {0,a,r}','rows':rows,'wall_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'phi_source_sha256':hashlib.sha256((p/'phi_diagnostic.py').read_bytes()).hexdigest()};f=p/'ratio_probe.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'states':len(rows),'max_rho_ratio':max(float(x['rho_scaled_Jphi_over_J']) for x in rows),'output_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'source_sha256':out['source_sha256']}))

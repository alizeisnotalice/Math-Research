"""Post-run diagnostic of one initial failed parameter, not an original run receipt."""
import hashlib,json,sys,time
from fractions import Fraction as F
from pathlib import Path
import r057_kinks as m
p=Path(__file__).resolve().parent;rows=[]
for precision in (160,220,256):
 m.mp.mp.dps=precision;d=m.DS[2];q=1-3*d/4+d*F(3,30);r=q+d;h=d/8;gamma=h/2;rho=2**256*r;t=rho/q;a=q*(1-gamma);w=a-F(1,2**96);dq=F(1,2**200)
 z=m.direct(q,d,h,gamma,t,w);ref,der=m.calc(q,d,h,gamma,t,w)
 left=(m.mp.log1p(z)-m.mp.log1p(m.direct(q-dq,d,h,gamma,t,w)))/m.mq(dq)
 right=(m.mp.log1p(m.direct(q+dq,d,h,gamma,t,w))-m.mp.log1p(z))/m.mq(dq)
 e0=abs(m.mq(rho)*(z-ref));e1=max(abs(m.mq(rho)*left-der),abs(m.mq(rho)*right-der))
 rows.append({'decimal_digits':precision,'scaled_integral_error':m.mp.nstr(e0,35),'scaled_derivative_error':m.mp.nstr(e1,35),'integral_error_lt_1e-95':e0<m.mp.mpf('1e-95'),'derivative_error_lt_1e-35':e1<m.mp.mpf('1e-35')})
result={'identity':'post-run one-point precision probe, not original 240-chunk receipt','d':str(d),'q':str(q),'rho_over_r':str(2**256),'h_over_d':'1/8','gamma_over_h':'1/2','w':str(w),'difference_step':'2^-200','rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'executed_unix':time.time(),'limitation':'Original chunk failures did not record assertion line; this does not classify all 240 failures.'}
(p/'postrun_precision_probe.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

"""New analytic error budget for175; no real nodes or old baseline replay."""
import sympy as s,json,os,time
from pathlib import Path
start=time.monotonic();t=s.symbols('t',nonnegative=True);D=1+t+t**2+t**3;R=1+t**4/D
res1=s.expand(D**2-t**3*(4+3*t+2*t**2+t**3))
res2=s.cancel(R-(1+t**4/4)-t**4*(4-D)/(4*D))
assert res1==1+2*t+3*t**2 and res2==0
result={'pid':os.getpid(),'status':'PASS_new_symbolic_precision_budget','seconds':time.monotonic()-start,'one_minus_Rprime_numerator':str(res1),'gap_vs_unrounded164_residual':str(res2),'new_real_roots_gains_H_log':0,'inference':'0<=Rprime<=1 on[0,1],1<=R<=5/4. If a=U^(1/4), y=floor(Ba)/B,z=floor(Bt)/B, then 0<=aR(t)-yR(z)<=(5/4+y)/B.'}
out=Path(__file__).with_name('secant_precision_result.json');assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

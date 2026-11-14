"""Three-point ABSTRACT information spectra. Exact analytic all-t proof separate.
Numerical values are diagnostics; no responsible rectangular realization.
"""
import json,time,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import mpmath as mp
HERE=Path(__file__).resolve().parent

def one(task):
 n,v=task;ans=[]
 for prec in [80,120]:
  with mp.workdps(prec):
   mu=1+mp.lambertw(mp.exp(-1));m=mu*n;delta=mp.log(n)/2+4;h=m-delta;p=mp.exp(-m);a=delta+mp.mpf(v)*mp.log(n) if v!='sqrt' else mp.sqrt(n)/2
   D=1-mp.exp(-h)*(mp.exp(-a)/2+(mp.mpf('.5')-p)*mp.exp(a));b=mp.log(p)-mp.log(D);probs=[mp.mpf('.5'),mp.mpf('.5')-p,p];levels=[h+a,h-a,b];M=mp.fsum(q*x for q,x in zip(probs,levels));var=mp.fsum(q*(x-M)**2 for q,x in zip(probs,levels));V=h+a+mp.log(2)
   assert 0<p<mp.mpf('.5') and h-a>mp.log(2) and a*a<=m/2 and a<=m/4 and b<h-a
   ans.append([mp.nstr(x,55) for x in [V-m,var,a,delta]])
 assert ans[0]==ans[1]
 return dict(n=n,variant=v,precision_pair=[80,120],saved_digits=55,V_minus_mu_n=ans[1][0],variance=ans[1][1],a=ans[1][2],mean_loss_delta=ans[1][3],normal_CDF_distance_lower='1/4 (analytic atom bound)',all_t_conditions_numerically_pass=True,scope='Finite positive three-step rank density, NOT a geometric counterexample')
if __name__=='__main__':
 started=time.monotonic();jobs=[(2**k,v) for k in range(12,35,2) for v in ['1/10','1','4','sqrt']]
 # mpmath parses rational explicitly, avoiding string-to-float conversion.
 jobs=[(n,'0.1' if v=='1/10' else v) for n,v in jobs]
 with ProcessPoolExecutor(max_workers=12) as p:res=list(p.map(one,jobs))
 (HERE/'three_point.json').write_text(json.dumps(dict(status='complete',workers=12,cases=len(res),certified_population=False,formula_scope='Abstract spectra only; all-t bound is analytic',results=res,elapsed=time.monotonic()-started),indent=2)+'\n');print(json.dumps(dict(status='complete',cases=len(res),elapsed=time.monotonic()-started)))

"""Exact old-prior mismatch example plus actual score-tilt prior displacement."""
from pathlib import Path
from fractions import Fraction as F
import math,json,time
import numpy as np
import probe
HERE=Path(__file__).resolve().parent

def main():
 target=HERE/'NORMALIZATION.json';assert not target.exists();start=time.time();rows=[]
 for n in [1,2,4,8,16,32,64]:
  a=sum(F(math.comb(n,k))*F(2,3)**(2*k)*F(1,3)**(2*(n-k))/F(1,2)**n for k in range(n+1))
  expected=F(10,9)**n;assert a==expected
  corrected=sum(F(math.comb(n,k))*F(2,3)**k*F(1,3)**(n-k)for k in range(n+1));assert corrected==1
  rows.append(dict(n=n,old_denominator_trace_exact=str(a),old_denominator_trace_float=float(a),new_denominator_trace_exact=str(corrected)))
 real=[]
 for spec in [dict(K=4,q=1,ell=.1,kind='positive_box',shift=0),dict(K=32,q=1,ell=1/3,kind='positive_box',shift=.5)]:
  for t in [-32,-8,8,32]:
   x,w,a,p=probe.posterior(**spec);z=np.log(w)+(1+t)*np.log(a);ww=np.exp(z-z.max());ww/=ww.sum()
   tr,mu=probe.trace(p,ww);factor=float(np.sum(mu*mu)*spec['K'])
   real.append(dict(**spec,t=t,eta_t=mu.tolist(),new_prior_one_coordinate_trace=tr,
      old_prior_Jensen_factor=factor,
      old_prior_log_lower_bound_at_n1000000=1000000*math.log(factor)))
 out=dict(exact_counter=rows,actual_score_tilt_examples=real,scope='Wrong original-prior denominator creates a spurious exponential normalized trace; not a C_n lower bound',script_SHA256=probe.sha(__file__),seconds=time.time()-start)
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(rows=len(rows),actual_cases=len(real),max_old_log_bound=max(r['old_prior_log_lower_bound_at_n1000000']for r in real))))
if __name__=='__main__':main()

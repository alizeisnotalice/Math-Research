from pathlib import Path
import json,math
from scipy.stats import beta
import probe as p
H=p.HERE;groups=[];a=.05/60
for name,path in [('free_R4',H/'SUMMARY.json'),('balanced_R2',H/'BALANCED_SUMMARY.json')]:
 s=json.loads(path.read_text())
 for rec in s['results']:
  rows=[]
  for rr in rec['rows']:
   k=rr['k_guard'];N=s['config']['test_N'];lam=rec['lambda_'];lo=0. if k==0 else float(beta.ppf(a,k,N-k+1));hi=1. if k==N else float(beta.ppf(1-a,k+1,N-k));rows.append(dict(model=rr['model'],k=k,N=N,lambda_=lam,weak=rr['weak'],CP_weak=[lam*lo,lam*hi]))
  groups.append(dict(group=name,n=rec['n'],rows=rows,gate_minus_bestfixed=rec['gate_minus_bestfixed'],gate_minus_prior=rec['gate_minus_prior'],lower_gate_minus_upper_prior=rows[0]['CP_weak'][0]-rows[1]['CP_weak'][1]))
p.dump(H/'COMBINED_SUMMARY.json',dict(groups=groups,alpha_each=a,simultaneous_directions=60,confidence_scope='Exploration: ideal iid frozen float-classifier probabilities only; not outward interval nor continuous-cell certificate',zero_count_upper=-math.expm1(math.log(a)/8192),full_FP_certificate=False,Cn_growth_proved=False))
print(-math.expm1(math.log(a)/8192))

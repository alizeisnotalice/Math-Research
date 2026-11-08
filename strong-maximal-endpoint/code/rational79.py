import json,sys,math
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from sliver79 import logi,outward
def psi(t):
 l,h=logi(t);return t**3/Q(36)+t**3*l/Q(6),t**3/Q(36)+t**3*h/Q(6)
def phi(t):
 if t<=1:return t,t
 a,b=logi(t);return t*(1+a),t*(1+b)
if __name__=='__main__':
 p=Path(sys.argv[1]);p.mkdir(exist_ok=True,parents=True)
 a,b,c,d=[psi(x) for x in [Q(7,4),Q(3,2),Q(5,4),Q(1)]]
 j=(8*(a[0]-b[1]-c[1]+d[0]),8*(a[1]-b[0]-c[0]+d[1]))
 rows=[]
 for m in [1,2,4,8,16,32,64,128]:
  counts=Counter((3*(8*m+4*i+2*k+3))//8 for i in range(m) for k in range(m))
  lo=sum((Q(v,m*m)*phi(Q(k,3*m))[0] for k,v in counts.items()),Q(0));hi=sum((Q(v,m*m)*phi(Q(k,3*m))[1] for k,v in counts.items()),Q(0))
  mass=sum((Q(v*k,3*m**3) for k,v in counts.items()),Q(0))
  assert hi<=j[0] and sum(counts.values())==m*m
  rows.append(dict(mesh=m,height_denominator=3*m,height_counts=[[k,v] for k,v in sorted(counts.items())],mass=str(mass),modular_interval=[str(lo),str(hi)],L1_error_upper=str(Q(17,24*m)),whole_receiver_square_success=mass>Q(117,100)**2))
 scale=10**40;sq=math.isqrt(2*scale*scale);al=Q(sq,10*scale);au=Q(sq+1,10*scale)
 rr=[]
 for n in [2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576]:
  q=16*n;f=(q*al).__floor__();assert f==(q*au).__floor__()
  l,h=logi(Q(3,2));v,w=logi(1+Q(5,2*q));lm=(l-n*w,h-n*v);assert lm[0]>0
  rr.append(dict(N=n,anchor_denominator=q,lower_anchor_numerator=f,upper_anchor_numerator=f+1,log_rounded_mean_lower_interval=[str(t) for t in lm],source_orthant_budget_ratio='1'))
 (p/'rational79.json').write_text(json.dumps(dict(smooth_source='(1+x/2+y/4) on [0,1]^2',smooth_modular_interval=[str(v) for v in j],density_rows=rows,irrational_anchor_interval=[str(al),str(au)],rounding_rows=rr),indent=2)+'\n')
 print('PASS',len(rows),'density approximants and',len(rr),'zero-budget-loss rational-anchor roundings')

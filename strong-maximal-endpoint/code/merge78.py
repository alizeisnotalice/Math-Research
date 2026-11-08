import math,json,sys
from pathlib import Path
from fractions import Fraction as Q
from scipy.special import lambertw
from allocate77 import phi_interval,exact_dot,rows_at
def inv(x):return x if x<=1 else (x/lambertw(math.e*x).real)
if __name__=='__main__':
 p=Path(sys.argv[1]);z=json.loads((p/'columns78.json').read_text());c=z['cases'][0];assert c['name']=='uniform_single'
 m=len(c['receiver_cells']);n=len(c['source']);full=(1<<m)-1;idx=[j for j,t in enumerate(c['columns']) if t['mask'] in [full^(1<<i) for i in range(m)]]
 assert len(idx)==m
 w=Q(1,m-1);F=list(map(Q,c['source']));Gcols=[list(map(Q,c['columns'][j]['primal'])) for j in idx]
 budget_hi=[sum((w*phi_interval(g[l])[1] for g in Gcols),Q(0)) for l in range(n)]
 out=[]
 for f,b in zip(F,budget_hi):
  if b>=phi_interval(f)[0]:q=f
  else:
   q=Q(math.ceil(inv(float(b))*10**9),10**9)
   for _ in range(8):
    if q>=f:q=f;break
    if phi_interval(q)[0]>=b:break
    q+=Q(1,10**9)
   else:raise ArithmeticError('rounding repair failed')
  out.append(q)
 assert all(0<=g<=f for g,f in zip(out,F))
 rows=list(dict.fromkeys(t for P in c['receiver_cells'] for t in rows_at(c['anchors'][0],P,2)))
 slack=min(exact_dot(a,out)-Q(c['sigma'])*v for a,v in rows);assert slack>=0
 jhi=sum((Q(1,4)*phi_interval(g)[1] for g in out),Q(0));bhi=sum((Q(1,4)*b for b in budget_hi),Q(0))
 assert jhi<=bhi+Q(1,10**6)
 result=dict(case=c['name'],column_indices=idx,weights=[str(w)]*m,rational_dominating_aggregate=[str(q) for q in out],
  minimum_exact_vertex_slack=str(slack),aggregate_modular_upper=str(jhi),weighted_input_modular_upper=str(bhi),rounding_excess_upper=str(max(Q(0),jhi-bhi)),
  scope='Exact inverse-budget aggregate has cost at most weighted actual costs by proof; this rational upward approximation has the displayed certified rounding allowance.')
 (p/'merge78.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS balanced same-anchor cover, exact vertex slack',float(slack),'rounding cost excess <=',float(max(Q(0),jhi-bhi)))
 ns=[2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576]
 cr=[]
 for n in ns:
  e=1/(8*n);d=e/4;h=e*e/(64*n)
  lm=math.log1p(e)+math.log1p(-d)-(n-1)*math.log1p(h)-math.log1p(h/d)
  lq=-math.log(2)+math.log1p(e)+(n-1)*math.log1p(math.log1p(e))
  assert lm>math.log1p(e/3) and lq<0 and h<d
  cr.append(dict(N=n,log_minimum_actual_mean=lm,log_aggregate_height=lq,log_receiver_volume=n*math.log(h)))
 (p/'cross78.json').write_text(json.dumps(dict(rows=cr),indent=2)+'\n')
 print('PASS',len(cr),'same-orientation cross-anchor counterexample dimensions')

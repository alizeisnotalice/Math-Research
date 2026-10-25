#!/usr/bin/env python3
"""Exact six-atom abstract generalized-moment obstruction, no distribution claim."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,time
HERE=Path(__file__).resolve().parent
MOM=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/astra_upper006/verification.json')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def solve(matrix,rhs):
 n=len(rhs);a=[list(row)+[b] for row,b in zip(matrix,rhs)]
 for j in range(n):
  pivot=next(i for i in range(j,n) if a[i][j])
  a[j],a[pivot]=a[pivot],a[j];q=a[j][j];a[j]=[x/q for x in a[j]]
  for i in range(n):
   if i==j:continue
   q=a[i][j]
   if q:a[i]=[x-q*y for x,y in zip(a[i],a[j])]
 return [row[-1] for row in a]
def main():
 started=time.monotonic();reg=json.loads((HERE/'registration.json').read_text());Q=F(reg['Q'])
 moments=json.loads(MOM.read_text())['moments'];m=[F(1)]+[F(moments[str(k)][0]) for k in range(1,5)]
 attempts=[];winner=None
 for phase in range(5):
  master=HERE/f'master_{phase}.json'
  if not master.exists():break
  source=json.loads(master.read_text());support=source.get('dual_support',[])
  if len(support)!=6:
   attempts.append({'phase':phase,'status':'not_six_supports','support_count':len(support)});continue
  zs=[]
  for row in support:
   p=row['source_node'];zs.append(F(p['rational']) if 'rational'in p else F(p['rational_for_exact_witness']))
  assert all(F(0)<z<=F(1) for z in zs)
  matrix=[[z**(20*k) for z in zs] for k in range(5)]+[[z**(-4) for z in zs]]
  rhs=m+[Q];weights=solve(matrix,rhs)
  residuals=[sum(q*w for q,w in zip(row,weights))-b for row,b in zip(matrix,rhs)]
  assert residuals==[F(0)]*6
  positive=all(w>=0 for w in weights);value=sum(w*z**(-3) for w,z in zip(weights,zs))
  attempt={'phase':phase,'status':'exact_nonnegative_witness' if positive else 'negative_weight','z':[str(z) for z in zs],'weights':[str(w) for w in weights],
   'value':str(value),'value_display_only':float(value),'strictly_above_17705_over_12500':value>F(17705,12500),'master_source_sha256':sha(master)}
  attempts.append(attempt)
  if positive and value>F(17705,12500):
   assert value>F(3,2)
   winner={'phase':phase,'atoms':[{'z':str(z),'d':str(z**20),'weight':str(w)} for z,w in zip(zs,weights)],'moments':[str(x) for x in m],
    'negative_moment_Q':str(Q),'moment_residuals':[str(x) for x in residuals],'objective_exact':str(value),'objective_display_only':float(value),
    'objective_minus_3_over_2_exact':str(value-F(3,2)),'objective_minus_threshold_exact':str(value-F(17705,12500)),
    'all_weights_strictly_positive':all(w>0 for w in weights),'all_z_in_open_zero_closed_one':True,'equations_exact':True}
   break
 out={'status':'PASS_exact_abstract_moment_obstruction' if winner else 'NO_CERTIFIED_OBSTRUCTION','pid':os.getpid(),'attempts':attempts,'witness':winner,
  'weak_duality_statement':'For every globally feasible majorant, alpha*Q+A+B*m1+C*m2+D*m3+E*m4 >= sum_i w_i*z_i^(-3) > 3/2 > 17705/12500. This excludes the whole scalar family even without coefficient bounds; nu is only abstract, not real determinant law and not a C5 lower bound.',
  'source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),'moment_source_sha256':sha(MOM),'process_results_sha256':sha(HERE/'process_results.json'),
  'wall_seconds':time.monotonic()-started,'scope':'Conditional use of new Q source only; no actual I or C5 bound, no global majorant certificate'}
 (HERE/'exact_witness_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ('status','pid','wall_seconds')}))
 if winner:print(json.dumps({'phase':winner['phase'],'objective_display_only':winner['objective_display_only'],'all_weights_strictly_positive':winner['all_weights_strictly_positive'],'equations_exact':winner['equations_exact']}))
if __name__=='__main__':main()

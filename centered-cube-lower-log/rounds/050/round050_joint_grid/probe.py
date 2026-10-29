"""Exact scalar budgets for joint GF families and safe finite grids.
No field labels or basis families of size R are allocated.
"""
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
from pathlib import Path
import json,math,hashlib,time
HERE=Path(__file__).resolve().parent

def floor(x):return x.numerator//x.denominator

def ceil(x):return -((-x.numerator)//x.denominator)

def log10_bracket(N=16):
 z=F(9,11);lo=2*sum(z**(2*j+1)/F(2*j+1)for j in range(N));tail=2*z**(2*N+1)/((2*N+1)*(1-z*z))
 return lo,lo+tail

def H_budget(B):
 N=B-1;T=N*(N+1)//2;Dcnt=N*(N-1)
 return N+T*(T-1)+2*Dcnt*Dcnt+4*Dcnt

def dec(x):
 with localcontext() as ctx:
  ctx.prec=45
  return str(D(x.numerator)/D(x.denominator))

def main():
 start=time.perf_counter();loglo,loghi=log10_bracket();assert loglo>F(23,10) and loghi<F(7,3);p=F(99,100);spacing=F(1);rows=[]
 for exponent in [32,42,60]:
  n=10**exponent;threshold_hi=16*exponent*loghi;threshold_lo=16*exponent*loglo
  # Actual exact selection/certification of largest d with v>=16log n.
  d=1
  while 9*p*p*n/(29*((1<<(d+1))-1)**2)>=threshold_hi:d+=1
  K=(1<<d)-1;q=K+1;v=9*p*p*n/(29*K*K);v_next=9*p*p*n/(29*((1<<(d+1))-1)**2)
  assert v>=threshold_hi and v_next<threshold_lo
  B=1<<max(0,ceil(v).bit_length()-1)
  if B<v:B*=2
  assert v<=B<2*v
  base=dict(n_exponent=exponent,n=str(n),d=d,q=str(q),K=str(K),v=str(v),v_decimal=dec(v),B=B,H_B=str(H_budget(B)),K_choice_rationally_certified=True,spacing_constant_a=str(spacing))
  for name,R in [('4',4),('ceil(v)',ceil(v)),('ceil(v^2)',ceil(v*v))]:
   budget=F(R*H_budget(B)+(R*(R-1)//2)*(B-1)**2,q)
   Hcap=1+R*(B-1)*(B-2)+R*(R-1)*(B-1)**2
   Hused=min(K,Hcap);nodes_per_bad_lag=floor(4*v/spacing)+1;Jmax=floor(K*v/(9*spacing))+1;badcap=min(Jmax,Hused*nodes_per_bad_lag)
   sqrt_n=10**(exponent//2);assert sqrt_n*sqrt_n==n
   row=base|dict(R_choice=name,R=str(R),joint_bad_budget=str(budget),joint_bad_budget_decimal=dec(budget),budget_less_than_1=budget<1,existence_conclusion_from_budget='certified existential family'if budget<1 else 'budget inconclusive; not nonexistence',bad_lag_cardinality_polynomial_upper=str(Hcap),bad_lag_cardinality_upper_with_group_cap=str(Hused),grid_cardinality_upper=str(Jmax),nodes_per_bad_lag_upper=str(nodes_per_bad_lag),bad_node_cardinality_upper_conditional_on_joint_family=str(badcap),bad_nodes_over_sqrt_n_decimal=dec(F(badcap,sqrt_n)),scope='Only arithmetic budgets and conditional finite-grid counts; no full R labels generated.')
   rows.append(row)
 # Exact round048 real-bank check of the extra same-endpoint condition at d47.
 old=json.loads((HERE.parent/'round048_banks'/'MANIFEST.json').read_text());bases=[b['basis']for b in old['banks']]
 def rank(xs):
  piv={}
  for x in xs:
   while x:
    j=x.bit_length()-1
    if j in piv:x^=piv[j]
    else:piv[j]=x;break
  return len(piv)
 concrete=[{'banks':[i,j],'same_endpoint_union_rank':rank(bases[i]+bases[j]),'intersection_dimension':22-rank(bases[i]+bases[j])}for i in range(4)for j in range(i+1,4)]
 assert all(r['intersection_dimension']==0 for r in concrete)
 # Adversarial ceil/floor and closed interval counting checks, entirely exact.
 gridchecks=[]
 for vf,af in [(F(3,2),F(1)),(F(7,3),F(2)),(F(4),F(1)),(F(1535,1),F(1))]:
  sep=af/vf;points=[j*sep for j in range(floor(4/sep)+1)];cap=floor(4*vf/af)+1;assert len(points)==cap and points[-1]<=4
  gridchecks.append({'v':str(vf),'a':str(af),'closed_interval_length_in_KL_units':'4','bound':cap,'achieved_by_left_endpoint_grid':len(points)})
 # Small exact combinatorial check: periodic distance, theta boundary floors, bad arcs.
 Ktoy=255;vtoy=F(7,3);atoy=F(1);step=atoy/(Ktoy*vtoy);left=F(10,9);right=F(11,9)
 grid=[left+j*step for j in range(floor((right-left)/step)+1)];Dset={0,31,40,Ktoy-31,Ktoy-40}
 def dist(x,d):
  z=floor((x-d)/Ktoy)
  return min(abs(x-d-z*Ktoy),abs(x-d-(z+1)*Ktoy))
 safe=0;checks=0;perlag={}
 for d0 in Dset:
  count=sum(dist(Ktoy*L,d0)<=2 for L in grid);assert count<=floor(4*vtoy/atoy)+1;perlag[str(d0)]=count
 for L in grid+[F(ceil(Ktoy*left),Ktoy)]:
  x=Ktoy*L;theta_set={F(0),F(1,17),F(16,17),x/2-floor(x/2),-x/2-floor(-x/2)}
  is_safe=min(dist(x,d0)for d0 in Dset)>2
  if is_safe:safe+=1
  for theta in theta_set:
   k=floor(theta+x/2)-floor(theta-x/2);assert k in {floor(x),ceil(x)}
   if is_safe:assert k%Ktoy not in Dset
   checks+=1
 toy={'K':Ktoy,'v':str(vtoy),'a':str(atoy),'D_star':sorted(Dset),'grid_nodes':len(grid),'bad_nodes':sum(min(dist(Ktoy*L,d0)for d0 in Dset)<=2 for L in grid),'per_lag_closed_arc_counts':perlag,'per_lag_floor_upper':floor(4*vtoy/atoy)+1,'safe_test_scales':safe,'floor_ceil_phase_cases_including_boundaries':checks,'passed':True,'scope':'Combinatorial check only; D_star here is a deliberately chosen inverse-closed toy set, not a new GF certificate.'}
 out={'log10_rational_bracket':[str(loglo),str(loghi)],'small_exact_periodic_grid_check':toy,'budget_rows':rows,'round048_four_real_banks_same_endpoint_check':concrete,'closed_interval_floor_count_tests':gridchecks,'elapsed_seconds':time.perf_counter()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Budget<1 proves existence. Budget>=1 is inconclusive, even when a concrete family exists. This is a specified real convex-mixture GF family, not an upper bound for general C_n.'}
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'rows':[{k:r[k]for k in ['n_exponent','d','v_decimal','B','R_choice','R','joint_bad_budget_decimal','budget_less_than_1','bad_nodes_over_sqrt_n_decimal']}for r in rows],'real_bank_same_endpoint_ranks':concrete,'seconds':out['elapsed_seconds']},indent=2))
if __name__=='__main__':main()

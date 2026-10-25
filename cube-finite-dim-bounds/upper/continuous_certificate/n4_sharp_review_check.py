"""Three coefficient sets and one exact heat-kernel anchor; read-only source review."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import prod,comb
import json,hashlib,runpy,time
from flint import fmpq_mpoly_ctx,fmpq,fmpq_poly,arb,ctx as actx
OUT=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')
SOURCE=MAIN/'upper/u11_minimum_global_results.json'
GREEN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-upper/cube-finite-dim-bounds/upper/u11_green.py')
inputs=[MAIN/'rounds/011q-registration.md',MAIN/'src/round011_n4_coefficient_ratio.py',MAIN/'results/round011_n4_coefficient_ratio.json',MAIN/'proofs/round011_n4_minimum_sharpening.md',SOURCE,GREEN]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
assert hashes[str(SOURCE)]=='43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f'
assert hashes[str(GREEN)]=='fdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d'
start=time.perf_counter();C=fmpq_mpoly_ctx.get(('u0','u1','u2','u3'),ordering='lex');u=C.gens()
cones=[[(1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1)]]
def pol(d):return C.from_dict({tuple(map(int,e.split(','))):fmpq(v) for e,v in d.items()})
expected=[F(1,24),F(17,640),F(482191,13088544)];target=F(17,640);rows=[]
for row,rays,want in zip(json.loads(SOURCE.read_text())['physical_cases'],cones,expected):
 sub=[sum((rays[j][i]*u[j] for j in range(4)),C.constant(0)) for i in range(4)]
 A=(u[0]**3*pol(row['numerator_dict'])).compose(*sub);D=pol(row['denominator_dict']).compose(*sub)
 ad={tuple(map(int,e)):F(str(v)) for e,v in A.to_dict().items()};dd={tuple(map(int,e)):F(str(v)) for e,v in D.to_dict().items()};support=set(ad)|set(dd)
 assert dd and all(v>=0 for v in dd.values()) and all(v>=0 for v in ad.values())
 ratios=[(ad.get(e,F(0))/v,e) for e,v in dd.items() if v>0];best=min(v for v,e in ratios);assert best==want
 zeroD=[e for e in support if dd.get(e,F(0))==0];assert all(ad.get(e,F(0))>=0 for e in zeroD)
 diff={e:ad.get(e,F(0))-target*dd.get(e,F(0)) for e in support};assert all(v>=0 for v in diff.values())
 assert all(ad.get(e,F(0))-best*dd.get(e,F(0))>=0 for e in support)
 rows.append({'case':row['case'],'minimum_coefficient_ratio':str(best),'attaining_powers':[list(e) for v,e in sorted(ratios) if v==best],'A_nonzero_terms':len(ad),'D_nonzero_terms':len(dd),'A_only_D_zero_terms':len(zeroD),'D_only_A_zero_terms':len(set(dd)-set(ad)),'minimum_A_on_explicit_D_zero_terms':str(min((ad[e] for e in zeroD),default=F(0))),'A_minus_target_D_all_coefficients_nonnegative':True,'A_minus_own_min_D_all_coefficients_nonnegative':True})
# ONE point, independently sum the exact rational function of the time variable.
g=runpy.run_path(str(GREEN));point=(F(1),F(1),F(1),F(0));rates=g['rates']();nodes=[sum(point[j]*r[j] for j in range(4))/q for r,q in zip(g['RAYS'],rates)]
groups=defaultdict(F)
for ids,det in g['CONES']:groups[tuple(sorted(nodes[i] for i in ids))]+=F(det)/prod(rates[i] for i in ids)
N=fmpq_poly([0]);D=fmpq_poly([1])
for bs,alpha in groups.items():
 den=fmpq_poly([1])
 for b in bs:den*=fmpq_poly([fmpq(-b.numerator,b.denominator),1])
 h=D.gcd(den);dl=D//h;dr=den//h
 N=N*dr+dl*fmpq(alpha.numerator,alpha.denominator);D=dl*den
 h=N.gcd(D);N=N//h;D=D//h
assert D(fmpq(0))!=0 # no zero pole survives in the full time integrand

def taylor(poly,b,count):
 coeff=[F(str(v)) for v in poly.coeffs()]
 return [sum((coeff[j]*comb(j,k)*b**(j-k) for j in range(k,len(coeff))),F(0)) for k in range(count)]
res=[]
for b in sorted(set(nodes)):
 if not b:continue
 rest=D;factor=fmpq_poly([fmpq(-b.numerator,b.denominator),1]);m=0
 while rest(fmpq(b.numerator,b.denominator))==0:
  quo,rem=divmod(rest,factor);assert rem==0;rest=quo;m+=1
 if not m:continue
 nn=taylor(N,b,m);dd=taylor(rest,b,m);assert dd[0]
 f=[]
 for k in range(m):f.append((nn[k]-sum((dd[j]*f[k-j] for j in range(1,k+1)),F(0)))/dd[0])
 res.append({'node':str(b),'pole_order':m,'simple_pole_coefficient':str(f[-1])})
anchor=-sum((F(v['simple_pole_coefficient'])*(1 if F(v['node'])>0 else -1) for v in res),F(0))/32
assert anchor==target
actx.prec=256;gi=(arb(target.numerator)/target.denominator)/arb.pi()**3
record={'status':'PASS narrow coefficient sharpening and exact anchor review; final audit inherited','inputs':hashes,'three_cases':rows,'uniform_a_cubed_pi_cubed_G_lower':str(min(expected)),'anchor':{'point':['1','1','1','0'],'independent_time_integrand_groups':len(groups),'reduced_numerator_degree':N.degree(),'reduced_denominator_degree':D.degree(),'zero_poles_cancelled':True,'nonzero_poles':res,'pi_cubed_G':str(anchor),'method':'sum full cone rational Q(s), t=i*s; exact Laurent coefficients via Taylor recurrence, not author allnodes evaluator'},'minimum_model_G_interval':str(gi),'model_minimum_exact':'17/(640*pi^3)','new_maximal_upper_claim':False,'runtime_seconds':time.perf_counter()-start,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'n4_sharp_review_result.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))

"""Direct exact reconstruction of two tuples and manuscript base constants. No enumeration."""
from fractions import Fraction as F
from pathlib import Path
import math,json,hashlib
f8=math.factorial(8);f16=math.factorial(16);f64=math.factorial(64)
P=F(3579,4096);e0=F(1,2**20);eJ=eS=eT=eD=F(1,1024);eA=F(1,64)

def build(a,s,d):
 sigma=s*d*d;qmin=1-3*d/4;ulo=P/(2*(1+eJ));uhi=1/(P*(1-d/2)*(1-eJ));M=max(F(4),2*uhi/qmin)
 mass=(1-eA)*d/(2*(1+3*d/4));bracket=(1-3/a)*ulo-eS-eT;Gamma=mass*bracket/M;gamma=(1-eD)*Gamma;c=gamma-sigma/(1-sigma)
 return locals()

def check_row(v,L0,pm):
 a,s,d,sigma,qmin,ulo,uhi,M,mass,bracket,Gamma,gamma,c=[v[x] for x in ('a','s','d','sigma','qmin','ulo','uhi','M','mass','bracket','Gamma','gamma','c')]
 h=e0/(L0*d);nlo=F(L0**2,2);eta=d/32-a*d*d;nu=(s/2-a)*d*d
 mean=f8*(24*d**3/F(L0**6)+4*d/F(L0**7));var=48*f8*d**3/F(L0**6)
 cap=24*2**20*d*d/(eta**4*L0**3)
 scalar=F(2,L0**2)+(F(f16,2**176)+F(2**15*f16,2**160) if L0==2**27 else F(f16,2**256)+F(2**17*f16,2**240))
 tests={
 'domain_and_signs':a>3 and s>2*a and 0<d<=F(1,16) and d<1/(32*a) and 0<sigma<F(1,2) and bracket>0,
 'error_domains':0<eJ<1 and 0<eA<1 and 0<eD<1 and eS>0 and eT>0,
 'degree_branch_M4':M==4,'root_budget':0<gamma<Gamma<1,'strict_net_slack':c>F(1,2**16),
 'h_over_d':h/d<F(1,8),'real_epsilon_over_h':2*d/(e0*nlo*nlo)<F(1,2),
 'extended_q_positive':qmin-h>F(1,2) and 1-d/4+h<1,'extended_R_above_one':1+d/4-2*h>1,
 'source_separation':F(L0,2)*(qmin-h)>4*(1+3*d/4),'n_nd_Ld':nlo>=33 and nlo*d>=2 and L0*d>=1,
 'x_t_domains':F(L0,2)>1 and (1-sigma)*L0>1,'target_corner_gap':1/(1+3*d/4)-F(2,L0)>F(1,2),
 'target_mean':mean<=F(1,32),'target_variance':var<=F(3,32)**2*eA/3,'target_count':F(8,9)/(nlo*d)<=eA/3,
 'target_corner_face':F(3,4)**32<eA/3,'individual_J':4*d*f8/F(L0**7)<=eJ,'source_cap_tail':cap<eS,
 'tau_exponent':nu*L0>=64,'tau_prefactor':(2**20*d)**2*L0**3<2**(2*pm),'tau_Taylor64':F(2**pm*f64,2**384)<eT,
 'scalar_Taylor16':scalar<c-F(1,2**16),'P0_product':1-F(517,4096)==P,
 }
 assert all(tests.values()),[k for k,x in tests.items() if not x]
 v.update(tests=tests,L0=L0,scalar_base=scalar,eta=eta,nu=nu)
 return v
simple=check_row(build(F(9),F(19),F(1,1024)),2**27,51)
raw=check_row(build(F(73,8),F(293,16),F(1023,2**20)),2**32,58)
assert simple['gamma']==F(9775559871,281956214374400)
assert raw['gamma']==F(735063024285513,21076775872652902400)
assert simple['c']==F(4893063658542547,295647162275777740800)
assert raw['c']==F(6468305130094850733787980147,370780099533215320770594052505600)
comparison={
 'raw_strictly_below_simple_full_cost':raw['scalar_base']+F(4,2**32)<raw['c']-(simple['gamma']-simple['sigma']),
 'old_two_budget_gap':F(5,2**40)<F(1,2**16),
 'old_logloss_gap':F(1,2**38)+F(14,2**20)<F(1,2**16),
}
assert all(comparison.values())
fields=('a','s','d','sigma','ulo','uhi','M','mass','bracket','Gamma','gamma','c','eta','nu')
receipt={'status':'PASS','method':'stdlib Fraction/integer/factorial only','tuples':{name:{'exact_budget':{k:str(row[k]) for k in fields},'L0':str(row['L0']),'checks':row['tests']} for name,row in [('simple',simple),('raw',raw)]},'comparison':comparison,'scope':'Two designated exact tuples only; no grid traversal, optimality assertion, producer import, old-baseline replay, n=7 or production computation. Infinite tails are proved analytically in main.tex.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
root=Path(__file__).resolve().parent
(root/'verification_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'PASS','tuple_checks_each':len(simple['tests']),'comparison_checks':len(comparison),'thresholds':['log(n)>=2^27','log(n)>=2^32'],'no_search':True},ensure_ascii=False))

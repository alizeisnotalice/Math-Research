"""Exact full-cube certificate, fixed target, no I/tau/cost computation.
Reconstructs the one permitted star/global representation revision, checks null
stars and coefficient positivity after each physical cone substitution.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from functools import reduce
from math import gcd,lcm,prod
import hashlib,json,subprocess,sys,time
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx as arbctx
import u11_minimum_dd as dd
ROOT=Path(__file__).parent;start=time.time()
for name in ['u11_minimum_star_revision.py','u11_minimum_global_revision.py']:
 z=subprocess.run([sys.executable,'-B',str(ROOT/name)],capture_output=True,text=True)
 assert z.returncode==0,z.stderr
ctx=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');xs=ctx.gens();a=xs[0];g=dd.g;q=dd.QBODY

def canonical(z):
 den=lcm(*(x.denominator for x in z));nn=[int(x*den) for x in z];div=reduce(gcd,(abs(x) for x in nn if x));pp=tuple(x//div for x in nn);sc=Q(div,den)
 if next(x for x in pp if x)<0:pp=tuple(-x for x in pp);sc=-sc
 return pp,sc
null=[]
for label in ['A','C','E','F']:
 ray=g.BASE_RAYS[label];rid=g.RAYS.index(ray);terms=[];common=Counter()
 for ids,det in g.CONES:
  if rid not in ids:continue
  ff=Counter();co=Q(det)/prod(q[i] for i in ids)
  for j in ids:
   if j==rid:continue
   vec=tuple(Q(ray[k])/q[rid]-Q(g.RAYS[j][k])/q[j] for k in range(4));fac,sc=canonical(vec);ff[fac]+=1;co/=sc
  common|=ff;terms.append((co,ff))
 assert sum(common.values())<=32
 lin={f:sum((ctx.constant(v)*x for v,x in zip(f,xs)),ctx.constant(0)) for f in common}
 N=ctx.constant(0)
 for co,ff in terms:
  z=ctx.constant(fmpq(co.numerator,co.denominator))
  for f,m in common.items():z*=lin[f]**(m-ff[f])
  N+=z
 assert N==0,(label,N)
 null.append({'ray_type':label,'star_terms':len(terms),'common_degree':sum(common.values()),'numerator_exactly_zero':True})

raw=json.loads((ROOT/'u11_minimum_global_results.json').read_text())
# Same exact physical cones as the proved four-dimensional spectral decomposition.
cones=[[(1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)],
       [(2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)],
       [(2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1)]]
TARGET=Q(21,25000);PIUP=Q(22,7);k=TARGET*PIUP**3
arbctx.prec=128
assert arb.pi()<arb(22)/7
rows=[]
for row,rs in zip(raw['physical_cases'],cones):
 def pol(d):return ctx.from_dict({tuple(map(int,z.split(','))):fmpq(v) for z,v in d.items()})
 N=pol(row['numerator_dict']);D=pol(row['denominator_dict'])
 sub=[sum((rs[j][i]*x for j,x in enumerate(xs)),ctx.constant(0)) for i in range(4)]
 const,factors=D.factor();signproof=[]
 for f,m in factors:
  z=f.compose(*sub);co=list(z.to_dict().values())
  if all(x<=0 for x in co):z=-z;co=[-x for x in co];const*=(-1)**m
  assert co and all(x>=0 for x in co) and any(x>0 for x in co)
  signproof.append({'factor':str(f),'multiplicity':int(m),'cone_coefficients_nonnegative':True})
 assert const>0
 P=(a**3*N-ctx.constant(fmpq(k.numerator,k.denominator))*D).compose(*sub)
 coefficients=P.to_dict();assert coefficients and all(v>=0 for v in coefficients.values())
 rows.append({'case':row['case'],'cone_rays':rs,'N_degree':int(N.total_degree()),'D_degree':int(D.total_degree()),
              'P_degree':int(P.total_degree()),'P_nonzero_coefficients':len(coefficients),
              'P_negative_coefficients':0,'denominator_constant_positive':str(const),
              'denominator_factor_sign_proof':signproof,
              'P_coefficients':[{'powers':list(map(int,e)),'value':str(v)} for e,v in sorted(coefficients.items())]})
 print('CERT case',row['case'],'positive coefficient count',len(coefficients),flush=True)
inputs={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
 ['u11_green.py','u11_representation.md','u11_minimum_dd.py','u11_minimum_star_revision.py','u11_minimum_global_revision.py','u11_minimum_global_results.json']}
out={'status':'EXACT FULL CUBE LOWER CERTIFICATE PASS; independent audit pending',
     'target_g_lower':str(TARGET),'pi_upper':str(PIUP),'pi_upper_verified_by_Arb_bits':128,
     'k_target_pi_upper_cubed':str(k),'null_star_identities':null,'physical_cases':rows,
     'proof_rule':'P=a^3*N-k*D has nonnegative coefficients after exact nonnegative physical cone map; D>0 on interior; extend by G continuity/off-line and infinity/on-line',
     'global_G_lower_certified_by_author':True,'new_weak_type_upper_bound':False,
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inputs':inputs}
(ROOT/'u11_minimum_certificate_results.json').write_text(json.dumps(out,indent=2)+'\n')
print('FULL CUBE G >=',TARGET,'PASS seconds',time.time()-start,flush=True)

"""Read-only fixed-identity and boundary review, independent of author producer scripts."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
from collections import Counter
from math import prod
import runpy,json,time,hashlib
import sympy as sy
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx as actx
OUT=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-upper/cube-finite-dim-bounds/upper')
start=time.perf_counter();files=['u11_green.py','u11_representation.md','u11_minimum_proof.md','u11_minimum_certificate_check.py','u11_minimum_dd.py','u11_minimum_star_results.json','u11_minimum_global_results.json','u11_minimum_certificate_results.json']
hashes={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in files}
assert hashes['u11_minimum_certificate_check.py']=='35489f9fa4c52a09154105ad949a8d3410f61bf226e8a24dd8d57359f0a3e060'
g=runpy.run_path(str(SRC/'u11_green.py'));rates=g['rates']();rays=g['RAYS'];cones=g['CONES']
ctx=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');x=ctx.gens();ZERO=ctx.constant(0);ONE=ctx.constant(1)
def qp(q):q=Q(q);return ctx.constant(fmpq(q.numerator,q.denominator))
def lin(v):return sum((qp(v[i])*x[i] for i in range(4)),ZERO)
def pol(d):return ctx.from_dict({tuple(map(int,e.split(','))):fmpq(v) for e,v in d.items()})
def reducefrac(n,d):
 if not n:return ZERO,ONE
 h=n.gcd(d);return n//h,d//h
def addfrac(n,d,p,q):
 h=d.gcd(q);dl=d//h;ql=q//h
 return reducefrac(n*ql+p*dl,dl*q)
# Verify symmetry using generators of the entire signed-permutation group.
base=Counter((tuple(sorted(rays[i] for i in ids)),det) for ids,det in cones)
symmetry=[]
for tag,perm,sign in [('swap01',(1,0,2,3),(1,1,1,1)),('swap12',(0,2,1,3),(1,1,1,1)),('swap23',(0,1,3,2),(1,1,1,1)),('flip0',(0,1,2,3),(-1,1,1,1))]:
 T=lambda r:tuple(sign[k]*r[perm[k]] for k in range(4))
 image=Counter((tuple(sorted(T(rays[i]) for i in ids)),det) for ids,det in cones)
 assert image==base
 for r,q in zip(rays,rates):assert dict(zip(rays,rates))[T(r)]==q
 symmetry.append(tag)
# Fixed stars: iterative reduced rational addition, no common-denominator star search.
starjson=json.loads((SRC/'u11_minimum_star_results.json').read_text())['stars'];starchecks=[];computed={}
for label in ['A','B','C','D','E','F']:
 r=g['BASE_RAYS'][label];rid=rays.index(r);n,d=ZERO,ONE;count=0
 for ids,det in cones:
  if rid not in ids:continue
  denominator=ONE
  for j in ids:
   if j!=rid:denominator*=lin(tuple(Q(r[k])/rates[rid]-Q(rays[j][k])/rates[j] for k in range(4)))
  n,d=addfrac(n,d,qp(Q(det)/prod(rates[j] for j in ids)),denominator);count+=1
 computed[label]=(n,d)
 if label in ('B','D'):
  row=starjson[0 if label=='B' else 1];assert n*pol(row['denominator_dict'])==d*pol(row['numerator_dict'])
 else:assert n==0
 starchecks.append({'label':label,'terms':count,'identity_pass':True,'zero':not bool(n)})
 print('STAR',label,'PASS',flush=True)
# Verify the compact formulas printed in the complete proof against reconstructed stars.
a,b,c,d=x;u=a+b;v=a*b;ss=c*c+d*d;tt=c*c*d*d
PB=v*u*u*(u*u+v)-v*(5*u*u+4*v)*ss+4*v*ss*ss+(13*u*u-32*v-4*ss)*tt
bn=u**3*PB;bd=2*(a*a-c*c)*(a*a-d*d)*(b*b-c*c)*(b*b-d*d)*(u*u-4*(c+d)**2)*(u*u-4*(c-d)**2)
assert computed['B'][0]*bd==computed['B'][1]*bn
dn=qp(Q(-9,2))*(a+b+c+d)**3
dd=(a-2*b-2*c+d)*(a-2*b+c-2*d)*(a+b-2*c-2*d)*(2*a-b-c+2*d)*(2*a-b+2*c-d)*(2*a+2*b-c-d)
assert computed['D'][0]*dd==computed['D'][1]*dn
physical=[[(1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1)]]
# Construct 20 representatives and their exact pullback substitutions directly.
terms=[]
for i,j in combinations(range(4),2):
 for sign in (-1,1):
  rem=[k for k in range(4) if k not in (i,j)];subs=[x[i],sign*x[j],x[rem[0]],x[rem[1]]]
  ray=[0]*4;ray[i]=1;ray[j]=sign
  n,d=computed['B'];terms.append((tuple(ray),n.compose(*subs),d.compose(*subs),'B'))
for ss in product((-1,1),repeat=3):
 ray=(1,)+ss;subs=[ray[i]*x[i] for i in range(4)];n,d=computed['D'];terms.append((ray,n.compose(*subs),d.compose(*subs),'D'))
assert len(terms)==20 and len(set(t[0] for t in terms))==20
raw=json.loads((SRC/'u11_minimum_global_results.json').read_text())['physical_cases'];cert=json.loads((SRC/'u11_minimum_certificate_results.json').read_text())['physical_cases'];casechecks=[];boundary=[];seen={}
for k,(rs,row,cr) in enumerate(zip(physical,raw,cert)):
 R=sy.Matrix.hstack(*(sy.Matrix(r) for r in rs));assert abs(R.det())==[1,1,2][k]
 sub=[sum((rs[j][i]*x[j] for j in range(4)),ZERO) for i in range(4)]
 n,d=ZERO,ONE;signchecks=[]
 for ray,tn,td,kind in terms:
  coeff=[sum(ray[i]*rs[j][i] for i in range(4)) for j in range(4)]
  assert all(v>=0 for v in coeff) or all(v<=0 for v in coeff)
  sign=1 if any(v>0 for v in coeff) else -1
  assert any(v!=0 for v in coeff)
  if kind=='B':assert sign==1
  signchecks.append({'ray':ray,'coefficients_on_cone_rays':coeff,'interior_sign':sign})
  n,d=addfrac(n,d,qp(Q(-sign,16))*tn,td)
 N,D=pol(row['numerator_dict']),pol(row['denominator_dict'])
 assert n*D==d*N
 const,factors=D.factor();denproof=[]
 for f,m in factors:
  fs=f.compose(*sub);dd=fs.to_dict();sg=-1 if all(v<=0 for v in dd.values()) else 1
  const*=sg**m;coeff={e:sg*v for e,v in dd.items()};assert coeff and all(v>=0 for v in coeff.values()) and any(v>0 for v in coeff.values())
  denproof.append({'factor':str(f),'multiplicity':int(m),'normalized_cone_coefficients':{','.join(map(str,e)):str(v) for e,v in coeff.items()}})
 assert const>0
 P=(x[0]**3*N-qp(Q(3993,153125))*D).compose(*sub)
 assert P and all(v>=0 for v in P.to_dict().values())
 encoded={tuple(v['powers']):fmpq(v['value']) for v in cr['P_coefficients']};assert P.to_dict()==encoded
 casechecks.append({'case':row['case'],'cone_matrix_inverse':str(R.inv()),'all_20_signs_constant_in_strict_interior':signchecks,'global_rational_identity_pass':True,'D_factors_strictly_positive_in_interior':denproof,'P_nonnegative_and_matches_certificate':True,'P_terms':len(P.to_dict())})
 # All 15 nonempty support strata. One exact point diagnoses each zero-factor pattern.
 direction=tuple(sum(r[i] for r in rs) for i in range(4))
 for mask in range(1,16):
  point=tuple(sum(rs[j][i] for j in range(4) if mask>>j&1) for i in range(4));Dvalue=D(*map(lambda v:fmpq(v),point));online=g['frame_line_distance_squared'](point)==0
  key=point
  if key not in seen:
   seen[key]=None if online else g['green_allnodes_pi3_coefficient'](point)
  exact=seen[key]
  # Exact inward path x(t)=point+t*sum(rays); compare leading coefficients.
  path=[qp(point[i])+direction[i]*x[0] for i in range(4)];nt=N.compose(*path).to_dict();dt=D.compose(*path).to_dict()
  nv=min(e[0] for e,v in nt.items() if v);dv=min(e[0] for e,v in dt.items() if v);lead=Q(str(nt[(nv,0,0,0)]))/Q(str(dt[(dv,0,0,0)]))
  if online:assert nv<dv and lead>0
  else:assert nv>=dv and (lead if nv==dv else Q(0))==exact
  boundary.append({'case':row['case'],'support_mask':mask,'point':point,'D_zero':not bool(Dvalue),'body_line':online,'inward_order_numerator':int(nv),'inward_order_denominator':int(dv),'limit_pi3':'+infinity' if online else str(exact),'exact_limit_matches_allnodes':True})
 print('CASE',row['case'],'PASS',flush=True)
# Independent finite basis-completion verification needed in the uniform short-time proof.
pairs=[];V=g['V']
for i,j in combinations(range(8),2):
 dets=[abs(int(sy.Matrix.hstack(*(sy.Matrix(V[k]) for k in (i,j,a,b))).det())) for a,b in combinations([k for k in range(8) if k not in (i,j)],2)]
 assert any(d>=8 for d in dets);pairs.append(min(d for d in dets if d))
actx.prec=128;assert arb.pi()<arb(22)/7
out={'status':'PASS limited exact identity/sign/boundary review; external audit remains required','input_hashes_at_read':hashes,'signed_permutation_generators_checked':symmetry,'star_identities':starchecks,'complete_proof_compact_B_D_formula_checked':True,'case_checks':casechecks,'boundary_checks':boundary,'boundary_unique_points':len(seen),'boundary_D_zero_count':sum(v['D_zero'] for v in boundary),'boundary_line_count':sum(v['body_line'] for v in boundary),'pair_basis_completion_count':len(pairs),'smallest_nonzero_pair_completion_determinant':min(pairs),'runtime_seconds':time.perf_counter()-start,'review_checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'n4_review_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('case_checks','boundary_checks')},indent=2))

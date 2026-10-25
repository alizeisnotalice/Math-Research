"""Second n4 review. No prior checker/producer imports; all writes stay here.
Shared exact arithmetic: SymPy, Fraction, python-flint polynomial compose/Arb.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations,product,combinations
from collections import Counter,defaultdict
from math import prod
import hashlib,json,time,sys
import sympy as s
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
load=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
start=time.monotonic()
lock=load(HERE/'inputs.json')['sha256']
def check_lock():
 for p,h in lock.items(): assert sha(ROOT/p)==h,p
check_lock()
C=fmpq_mpoly_ctx.get(('x0','x1','x2','x3'),ordering='deglex');X=C.gens();zero=C.constant(0);one=C.constant(1)
def cq(q):return C.constant(fmpq(str(q)))
def lin(a):return sum((cq(v)*x for v,x in zip(a,X)),zero)
def pol(d):return C.from_dict({tuple(map(int,k.split(','))):fmpq(v) for k,v in d.items()})
def plus(f,g):
 a,b=f;c,d=g;gg=b.gcd(d);num=a*(d//gg)+c*(b//gg);den=b*(d//gg)
 if not num:return zero,one
 gg=num.gcd(den);return num//gg,den//gg
V=list(product((1,),(-1,1),(-1,1),(-1,1)))
rr=[(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(2,1,1,0),(3,1,1,1)]
base=[tuple(rr[i] for i in j) for j in [(0,1,4,5),(4,1,2,3),(4,1,3,5)]]
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
rate=lambda r:sum(abs(dot(v,r)) for v in V)
det=lambda vs:Q(s.det(s.Matrix(vs)))
cones={}
for p in permutations(range(4)):
 for signs in product((-1,1),repeat=4):
  for r in base:
   t=tuple(tuple(signs[k]*v[p[k]] for k in range(4)) for v in r)
   assert tuple(sorted(t)) not in cones
   cones[tuple(sorted(t))]=(t,abs(det(t)))
assert len(cones)==1152
rays={v for t,d in cones.values() for v in t};assert len(rays)==240
for t,d in cones.values():
 for v in V:
  vals=[dot(v,r) for r in t];assert min(vals)>=0 or max(vals)<=0
assert [rate(r) for r in rr]==[8,8,12,12,16,24]
vol=sum((d/Q(24)/prod(rate(v) for v in t) for t,d in cones.values()),Q(0));assert vol==Q(11,4608)
rank=[s.Matrix([v for v in V if dot(v,r)==0]).rank() for r in rr]
assert rank==[0,3,0,3,2,1]
assert sum(s.Matrix([v for v in V if dot(v,r)==0]).rank()==3 for r in rays)==40
# All pairs extend to a basis. Cap uses the largest determinant of each completion.
maxd=[]
for i,j in combinations(range(8),2):
 ds=[abs(det([V[i],V[j],V[k],V[l]])) for k,l in combinations(set(range(8))-{i,j},2)]
 maxd.append(max(ds));assert max(ds)>0
assert Counter(maxd)=={Q(8):16,Q(16):12}
S=2*sum((1/d for d in maxd),Q(0));assert S==Q(11,2)
assert 4*8**2*4*S==5632 and 5632*8*4==180224
# Stars are exact rational functions summed in reverse cone order.
stars={};starrows=[]
for name,r in zip('ABCDEF',rr):
 f=(zero,one);count=0
 for t,d in reversed(list(cones.values())):
  if r not in t:continue
  denominator=prod((lin([Q(r[k],rate(r))-Q(v[k],rate(v)) for k in range(4)]) for v in t if v!=r),start=one)
  f=plus(f,(cq(d/prod(rate(v) for v in t)),denominator));count+=1
 if name in 'ACEF':assert f[0]==0
 stars[name]=f;starrows.append([name,count,bool(f[0])])
print('new cones and six star identities checked',flush=True)
# Closed formulas independently transcribed from their mathematical statements.
a,b,c,d=X;u=a+b;v=a*b;ss=c*c+d*d;tt=c*c*d*d
PB=v*u*u*(u*u+v)-v*(5*u*u+4*v)*ss+4*v*ss*ss+(13*u*u-32*v-4*ss)*tt
B=(u**3*PB,2*(a*a-c*c)*(a*a-d*d)*(b*b-c*c)*(b*b-d*d)*(u*u-4*(c+d)**2)*(u*u-4*(c-d)**2))
D=(-cq(Q(9,2))*(a+b+c+d)**3,(a-2*b-2*c+d)*(a-2*b+c-2*d)*(a+b-2*c-2*d)*(2*a-b-c+2*d)*(2*a-b+2*c-d)*(2*a+2*b-c-d))
for name,f in [('B',B),('D',D)]:assert f[0]*stars[name][1]==f[1]*stars[name][0]
terms=[]
for i,j in combinations(range(4),2):
 for sg in (-1,1):
  rem=sorted(set(range(4))-{i,j});sub=[X[i],sg*X[j],X[rem[0]],X[rem[1]]]
  r=tuple((1 if k==i else sg if k==j else 0) for k in range(4))
  terms.append((r,B[0].compose(*sub),B[1].compose(*sub)))
for r in V:
 sub=[r[k]*X[k] for k in range(4)];terms.append((r,D[0].compose(*sub),D[1].compose(*sub)))
models=[];minrows=[]
for index,row in enumerate(load(ROOT/'upper/u11_minimum_global_results.json')['physical_cases']):
 f=zero,one
 for r,n,de in terms:
  vals=[dot(r,w) for w in base[index]];assert any(vals) and (min(vals)>=0 or max(vals)<=0)
  sg=1 if max(vals)>0 else -1
  f=plus(f,(-cq(Q(sg,16))*n,de))
 N,De=pol(row['numerator_dict']),pol(row['denominator_dict']);assert f[0]*De==f[1]*N
 assert N.total_degree()==17 and De.total_degree()==20
 sub=[lin([r[k] for r in base[index]]) for k in range(4)]
 factor,factors=De.factor()
 for fp,m in factors:
  assert fp.total_degree()==1
  vals=list(fp.compose(*sub).to_dict().values());sg=-1 if max(vals)<0 or max(vals)==0 else 1
  assert all(sg*z>=0 for z in vals) and any(vals);factor*=sg**int(m)
 assert factor>0
 AA=(a**3*N).compose(*sub).to_dict();DD=De.compose(*sub).to_dict()
 assert AA and DD and all(z>0 for z in AA.values()) and all(z>0 for z in DD.values())
 support=set(AA)|set(DD)
 target=fmpq(17,640)
 assert all(AA.get(e,0)-target*DD.get(e,0)>=0 for e in support)
 best=min(Q(str(AA.get(e,0)/z)) for e,z in DD.items())
 assert best==[Q(1,24),Q(17,640),Q(482191,13088544)][index]
 minrows.append({'case':row['case'],'A_terms':len(AA),'D_terms':len(DD),'A_only':len(set(AA)-set(DD)),'minimum_ratio':str(best)})
 models.append((a**3*N,De))
print('three full sums and sharp minimum coefficients checked',flush=True)
# Different anchor mechanism: combine univariate rational functions before
# partial fractions. SymPy residue differentiates the full simplified function.
z=s.Symbol('z');anchorrows=[]
for point in [(1,0,0,0),(1,1,0,0),(1,1,1,0)]:
 group=defaultdict(Q)
 for t,dt in cones.values():group[tuple(sorted(Q(dot(point,r),rate(r)) for r in t))]+=dt/prod(rate(r) for r in t)
 ff=s.cancel(sum((s.Rational(co.numerator,co.denominator)/prod(z-s.Rational(b.numerator,b.denominator) for b in bs) for bs,co in group.items()),s.S.Zero))
 assert s.denom(ff).subs(z,0)!=0
 nodes=set(b for bs in group for b in bs if b)
 simple={b:Q(s.residue(ff,z,s.Rational(b.numerator,b.denominator))) for b in nodes}
 value=-sum((r*(1 if b>0 else -1) for b,r in simple.items()),Q(0))/32
 assert value=={(1,0,0,0):Q(3,32),(1,1,0,0):Q(1,24),(1,1,1,0):Q(17,640)}[point]
 anchorrows.append({'point':point,'rational_function':str(ff),'G_pi3':str(value),'zero_denominator_at_zero':False})
print('three full univariate rational anchor identities checked',flush=True)
# Geometry independent recursive subdivision; lexical first maximal edge.
def volume(t):return abs(det([tuple(t[j][k]-t[0][k] for k in (1,2,3)) for j in (1,2,3)]))/6
def split(t,depth):
 if not depth:return [t]
 edges=list(combinations(range(4),2));ls=[sum((t[i][k]-t[j][k])**2 for k in (1,2,3)) for i,j in edges]
 i,j=edges[ls.index(max(ls))];mid=tuple((u+v)/2 for u,v in zip(t[i],t[j]));out=[]
 for k in (i,j):
  child=tuple(mid if h==k else t[h] for h in range(4));assert volume(child)==volume(t)/2
  out.extend(split(child,depth-1))
 return out
roots=[tuple(tuple(Q(v,r[0]) for v in r) for r in rows) for rows in base]
assert [volume(t) for t in roots]==[Q(1,36),Q(1,12),Q(1,18)]
cells=[(0,t) for t in split(roots[0],8)];corner=(Q(1),)*4
for j in (1,2):
 other=[r for r in roots[j] if r!=corner];assert len(other)==3
 assert all(sum((x-1)**2 for x in r[1:])<=2 for r in other)
 for k in range(28):
  outer=[tuple(1+(x-1)/2**k for x in r) for r in other]
  inner=[tuple((x+1)/2 for x in r) for r in outer]
  for i in range(3):
   t=tuple(outer[i:]+inner[:i+1]);cells.extend((j,leaf) for leaf in split(t,max(8-2*k,0)))
assert len(cells)==2440
covered=sum((volume(t) for j,t in cells),Q(0));assert covered==Q(1,6)-Q(5,36*2**84)
cert=load(ROOT/'upper/continuous_integral_n4/ci4b_certificate.json');assert len(cert['cells'])==len(cells)
ctx.prec=512
qb=lambda q:arb(Q(q).numerator)/Q(q).denominator
total=arb(0);coeffcount=0;mins=None
for idx,((j,t),row) in enumerate(zip(cells,cert['cells'])):
 assert row['case']==j and tuple(tuple(map(Q,r)) for r in row['vertices'])==t
 p=list(map(Q,row['P']));assert len(p)==4 and min(p)>=0
 sub=[lin([r[k] for r in t]) for k in range(4)]
 NN=models[j][0].compose(*sub).to_dict();DD=models[j][1].compose(*sub).to_dict();assert DD and min(DD.values())>0
 # Assemble all coefficients with an explicit sparse shift, independent of
 # both author's and 035's polynomial product/residual implementation.
 residual=defaultdict(lambda:fmpq(0))
 for dic,weights in [(DD,list(map(lambda q:fmpq(str(q)),p))),(NN,[-1]*4)]:
  for e,val in dic.items():
   for k in range(4):
    ee=list(e);ee[k]+=1;residual[tuple(ee)]+=val*weights[k]
 assert all(val>=0 for val in residual.values()),idx
 nonzero=[val for val in residual.values() if val];coeffcount+=len(nonzero)
 if nonzero:mins=min(nonzero) if mins is None else min(mins,min(nonzero))
 total+=qb(volume(t))*sum((qb(v).root(3)**4 for v in p),arb(0))/4
 if (idx+1)%400==0:print('new residual checker',idx+1,'cells',flush=True)
assert coeffcount==4932474
pi=arb.pi();eps=qb(Q(1,2**27));Dcap=1/(48*pi**4);Bcap=180224/pi**3
cap=arb(2).root(3)*4*pi*(Dcap.root(3)**4*eps**3/3+qb(Q(3,5))*Bcap.root(3)**4*eps.root(3)**5)/8
integral=6*total/pi**4+cap;I0=Q(92072351,500000000000)
assert Q(str(integral.upper().fmpq()))<I0
atan=lambda x,n:sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(n)),Q(0))
pm=16*atan(Q(1,5),13)-4*atan(Q(1,239),4);pup=Q(355,113);assert pm<pup
cost=256*I0**3*Q(640,17)**4*pup**12;gap=Q(14371,1000)**3-cost;assert gap>0
check_lock()
result={'status':'all finite checks passed; analytical interpretation in review.md; shared theory separate','source_lock_sha256':sha(HERE/'inputs.json'),'checker_sha256':sha(Path(__file__)),'cones':1152,'triangulation_rays':240,'true_vertices':40,'spectral_volume':str(vol),'pair_max_det_counts':{str(k):v for k,v in Counter(maxd).items()},'stars':starrows,'minimum':minrows,'anchors':anchorrows,'cells':2440,'residual_nonzero_count':coeffcount,'minimum_positive_residual':str(mins),'covered_volume':str(covered),'cap_expression_512bit':str(cap),'I_expression_512bit':str(integral),'I_expression_exact_endpoints':[str(integral.lower().fmpq()),str(integral.upper().fmpq())],'I_outward':str(I0),'pi_machin_upper':str(pm),'pi_bound_gap':str(pup-pm),'strict_cost_cube_gap':str(gap),'seconds':time.monotonic()-start,'dependencies':{'python':sys.version,'sympy':s.__version__,'flint_shared':'fmpq_mpoly exact compose/gcd/factor and Arb512; Fraction shared; no previous checker or producer imported'}}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','I_expression_512bit','strict_cost_cube_gap','seconds']},indent=2))

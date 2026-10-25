"""11O continuous I_H certificate via exact-point DD, tensor Cheb and convex chord.

Each output patch certifies its entire closed square. Four-patch stage is not
full I certification. All cached intervals include node displacement error.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,prod
import json,hashlib,time,argparse,platform
from flint import arb,arb_mat,ctx
OUT=Path(__file__).resolve().parent
CONES=OUT/'cb11_cones.json'
CONE_SHA='23258e3d0236700ff442d602ccbb738a17c92728f3b40d95ddf9dd5d9ef924f7'
assert hashlib.sha256(CONES.read_bytes()).hexdigest()==CONE_SHA
DEGREE=32;N=33;GRID=128;NODE_BITS=128;PRECISION=384
ctx.prec=PRECISION
EPS=F(9,5000);TAU=F(1,3);FOUR=[(127,57),(127,127),(127,0),(0,0)]
CONFIG={'schema':1,'cone_sha256':CONE_SHA,'epsilon':str(EPS),'tau':str(TAU),'degree':DEGREE,'grid':GRID,'node_bits':NODE_BITS,'precision_bits':PRECISION,'rho':'5/2'}
def aq(x):
 x=F(x);return arb(x.numerator)/x.denominator
def pack(x):
 m,e=x.mid().man_exp();r,f=x.rad().man_exp();return [str(m),int(e),str(r),int(f)]
def unpack(v):return arb((int(v[0]),v[1]),(int(v[2]),v[3]))
def exact(x):
 m,e=x.man_exp();return F(int(m))*F(2)**int(e)
# Full 96 cones, exact opposite pairs; no signed summand omitted.
seen=set();PREP=[]
for c in json.loads(CONES.read_text())['cones']:
 rays=tuple(map(tuple,c['rays']));key=tuple(sorted(rays));opp=tuple(sorted(tuple(-v for v in r) for r in rays));pair=min(key,opp)
 if pair in seen:continue
 seen.add(pair);ks=tuple(F(a)+TAU*b for a,b in c['rates_body_face_coefficients'])
 assert abs(c['det'])==1
 assert all(k>=3*(abs(r[1])+abs(r[2])) for r,k in zip(rays,ks))
 PREP.append((rays,ks,F(2)/prod(ks)))
assert len(PREP)==48
NORMALIZATION=8*arb.pi()**3

def point(y,z):
 """All exact coalescences and zero projections classified with Fraction."""
 X=(F(1),F(y),F(z));E=aq(EPS);ee=E*E;total=arb(0);fcache={};qcache={}
 def ab(b):
  if b not in qcache:qcache[b]=aq(b)
  return qcache[b]
 def f(b):
  b=abs(b)
  if b not in fcache:fcache[b]=(ee+ab(b)**2).log()/2
  return fcache[b]
 for rays,ks,w in PREP:
  bs=tuple(sum((F(r[i])*X[i] for i in range(3)),F(0))/k for r,k in zip(rays,ks))
  if len(set(bs))==1:
   b=ab(bs[0]);dd=(ee-b*b)/(2*(ee+b*b)**2)
  elif len(set(bs))==2:
   a=next(b for b in bs if bs.count(b)==2);b=next(b for b in bs if b!=a)
   df=ab(a)/(ee+ab(a)**2);dd=(f(b)-f(a)-df*ab(b-a))/ab(b-a)**2
  else:dd=sum((f(b)/ab((b-bs[(j+1)%3])*(b-bs[(j+2)%3])) for j,b in enumerate(bs)),arb(0))
  total+=aq(w)*dd
 return total/NORMALIZATION

h=aq(F(1,256));rho=aq(F(5,2));s=h*(rho-1/rho)/2;delta=aq(EPS)-s/3;assert delta>0
ELLIPSE_M=15/(1088*arb.pi()**3*delta**2)
ERROR=(2*DEGREE+2)*4*ELLIPSE_M*rho**(-DEGREE)/(rho-1)
# Uniform real coordinate derivative from Cauchy disk of radius eps.
R=aq(EPS);D=15/(1088*arb.pi()**3*(aq(EPS)-R/3)**2*R)
NODE_ERROR=2*D*h/2**NODE_BITS
ROOTS=[None]*N
for j in range(DEGREE//2):
 x=(arb.pi()*(2*j+1)/(2*N)).cos();a=(x*2**NODE_BITS).floor().unique_fmpz();assert a is not None
 ROOTS[j]=F(int(a),2**NODE_BITS);ROOTS[DEGREE-j]=-ROOTS[j]
 assert abs(x-aq(ROOTS[j]))<aq(F(1,2**NODE_BITS))
ROOTS[DEGREE//2]=F(0)
A=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos()*(1 if k==0 else 2)/N for j in range(N)] for k in range(N)])
# Exact degree-d Bernstein transform of T_k(2u-1).
P=[[F(1)],[F(-1),F(2)]]
for k in range(2,N):
 p=[F(0)]*(k+1)
 for j,a in enumerate(P[-1]):p[j]-=2*a;p[j+1]+=4*a
 for j,a in enumerate(P[-2]):p[j]-=a
 P.append(p)
BF=[[sum((P[k][l]*F(comb(j,l),comb(DEGREE,l)) for l in range(min(j,k)+1)),F(0)) for k in range(N)] for j in range(N)]
for k in range(N):
 assert BF[0][k]==(-1)**k and BF[-1][k]==1
 for t in [F(1,3),F(1,2),F(2,3)]:
  assert sum((BF[j][k]*comb(DEGREE,j)*t**j*(1-t)**(DEGREE-j) for j in range(N)),F(0))==sum((P[k][j]*t**j for j in range(k+1)),F(0))
B=arb_mat([[aq(v) for v in row] for row in BF]);TRANS=B*A
MU=[F(1)]+[F(0) if k%2 else F(-1,k*k-1) for k in range(1,N)]
for k in range(N):assert sum((P[k][j]/(j+1) for j in range(k+1)),F(0))==MU[k]
# Scalar weights for exact mean: mu*A in each coordinate.
MA=arb_mat([[sum((aq(MU[k])*A[k,j] for k in range(N)),arb(0)) for j in range(N)]])
CACHE_FILE=OUT/'cb11_integral_cheb_four_cache.json'

def nodes(i,j,cache=None):
 key=f'{i},{j}'
 if cache is not None and key in cache['patches']:
  return arb_mat([[unpack(v) for v in row] for row in cache['patches'][key]]),True
 ys=[F(2*i+1,256)+t/256 for t in ROOTS];zs=[F(2*j+1,256)+t/256 for t in ROOTS]
 rows=[]
 for y in ys:
  row=[]
  for z in zs:
   v=point(y,z);assert v.is_finite();v+=arb(0,NODE_ERROR.upper());row.append(v)
  rows.append(row)
 return arb_mat(rows),False

def check_patch(i,j,cache=None,save_cache=None):
 assert 0<=j<=i<GRID
 start=time.perf_counter();V,hit=nodes(i,j,cache);node_seconds=time.perf_counter()-start
 if save_cache is not None:save_cache['patches'][f'{i},{j}']=[[pack(V[a,b]) for b in range(N)] for a in range(N)]
 C=TRANS*V*TRANS.transpose();mean=(MA*V*MA.transpose())[0,0]
 lo=min(C[a,b].lower() for a in range(N) for b in range(N));hi=max(C[a,b].upper() for a in range(N) for b in range(N))
 L=max(arb(0),(lo-ERROR).lower());U=(hi+ERROR).upper();assert U>=L and U>0
 mean_upper=min(U,(mean+ERROR).upper());assert mean_upper>=L
 # Stable chord slope, equivalent to (U^1.5-L^1.5)/(U-L).
 slope=(U+(U*L).sqrt()+L)/(U.sqrt()+L.sqrt())
 chord=L*L.sqrt()+slope*(mean_upper-L)
 assert chord>=0
 weight=1 if i==j else 2;contribution=(chord*weight/GRID**2).upper()
 return {'patch':[i,j],'whole_square':[f'{i}/{GRID}',f'{i+1}/{GRID}',f'{j}/{GRID}',f'{j+1}/{GRID}'],'weight':weight,'L':pack(L),'U':pack(U),'mean_H_upper':pack(mean_upper),'mean_P_ball':pack(mean),'chord_mean_upper':pack(chord.upper()),'IH_contribution_upper':pack(contribution),'IH_contribution_exact':str(exact(contribution)),'L_display':str(L),'U_display':str(U),'mean_H_upper_display':str(mean_upper),'contribution_display':str(contribution),'max_node_radius':str(max(V[a,b].rad() for a in range(N) for b in range(N))),'max_Bernstein_radius':str(max(C[a,b].rad() for a in range(N) for b in range(N))),'nodes':N*N,'cache_hit':hit,'node_seconds':node_seconds,'total_seconds':time.perf_counter()-start}

def time_error():
 pi=arb.pi();frac=lambda x,p:(x.log()*aq(p)).exp()
 return frac(aq(2),F(-4,3))*aq(F(3,4))*frac(aq(F(256,3))/pi**3,F(1,3))*frac(aq(F(16,3))/pi,F(2,3))*frac(aq(EPS),F(4,3))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['four','part','aggregate'],default='four');ap.add_argument('--part',type=int,default=0);ap.add_argument('--parts',type=int,default=12);args=ap.parse_args()
 assert sum(1 if i==j else 2 for i in range(GRID) for j in range(i+1))==GRID**2
 start=time.perf_counter();config_sha=hashlib.sha256(json.dumps(CONFIG,sort_keys=True).encode()).hexdigest();records=[]
 if args.stage=='aggregate':
  for k in range(args.parts):
   r=json.loads((OUT/f'cb11_integral_cheb_part_{k:02}.json').read_text());assert r['config']==CONFIG and r['part']==k and r['parts']==args.parts;records+=r['patches']
  assert {tuple(r['patch']) for r in records}=={(i,j) for i in range(GRID) for j in range(i+1)} and len(records)==8256
  Bsum=sum((F(r['IH_contribution_exact']) for r in records),F(0));Etime=time_error();val=((aq(Bsum).log()*aq(F(2,3))).exp()+Etime);I=val*val.sqrt()
  rec={'status':'PASS I target conditional on source/time dependencies' if I<aq(F(29,50000)) else 'INCOMPLETE I target','I_certified_locally':bool(I<aq(F(29,50000))),'config':CONFIG,'patch_count':len(records),'IH_upper_exact':str(Bsum),'IH_upper':str(aq(Bsum)),'Etime':str(Etime),'I_upper':pack(I.upper()),'I_upper_display':str(I.upper()),'source_dependencies_pending_audit':True,'triangle_weighted_area_exact':'1'}
 else:
  patches=FOUR if args.stage=='four' else [(i,j) for i in range(GRID) for j in range(i+1) if (i*(i+1)//2+j)%args.parts==args.part]
  cache=None
  if args.stage=='part' and CACHE_FILE.exists():
   cache=json.loads(CACHE_FILE.read_text());assert cache['config']==CONFIG
  saved={'config':CONFIG,'patches':{}} if args.stage=='four' else None
  for i,j in patches:
   r=check_patch(i,j,cache,saved);records.append(r);print(json.dumps({k:r[k] for k in ['patch','L_display','U_display','contribution_display','cache_hit','total_seconds']}),flush=True)
  if saved is not None:CACHE_FILE.write_text(json.dumps(saved,separators=(',',':'))+'\n')
  rec={'status':'PASS continuous patch chord enclosures, not full I' if args.stage=='four' else 'complete assigned finite patch subset, not full I','config':CONFIG,'config_sha256':config_sha,'part':args.part,'parts':args.parts,'patches':records,'ellipse_M':str(ELLIPSE_M),'interpolation_error':str(ERROR),'node_displacement_error':str(NODE_ERROR),'derivative_bound':str(D),'triangle_weighted_area_exact':'1','I_certified':False}
 rec.update({'python':platform.python_version(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start})
 name='cb11_integral_cheb_'+('four_result' if args.stage=='four' else f'part_{args.part:02}' if args.stage=='part' else 'all_result')+'.json';(OUT/name).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps({k:v for k,v in rec.items() if k!='patches'},indent=2),flush=True)
if __name__=='__main__':main()

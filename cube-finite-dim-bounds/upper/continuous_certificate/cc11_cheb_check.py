"""11J fixed-patch continuous minimum certificate; exact rational nodes + Arb."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import runpy,json,hashlib,time,argparse,platform
from flint import arb,arb_mat,ctx
OUT=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/cb11_frame_green_locked.py')
SOURCE_SHA='ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
OLD=OUT/'cc11_regularized.py'
assert hashlib.sha256(OLD.read_bytes()).hexdigest()==json.loads((OUT/'manifest.json').read_text())['files']['cc11_regularized.py']['sha256']
ns=runpy.run_path(str(OLD));ctx.prec=384
N=21;DEGREE=20;PRECISION=384;NODE_BITS=128;MAX_DEPTH=8
CACHE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/cc11_cheb_cache');CACHE.mkdir(exist_ok=True)
CONFIG={'schema':1,'source_sha256':SOURCE_SHA,'evaluator_sha256':hashlib.sha256(OLD.read_bytes()).hexdigest(),'precision_bits':PRECISION,'node_bits':NODE_BITS,'degree':DEGREE,'epsilon':'1/50','patch_grid':16}
def aq(q):q=F(q);return arb(q.numerator)/q.denominator
def pack(x):
 m,e=x.mid().man_exp();r,f=x.rad().man_exp();return [str(m),int(e),str(r),int(f)]
def unpack(a):return arb((int(a[0]),a[1]),(int(a[2]),a[3]))
def exact_lower(x):return x.lower()
def lower_matrix(A):return min(exact_lower(A[i,j]) for i in range(N) for j in range(N))
# Same signed-cone proof as first batch, rechecked as exact rational arithmetic.
ratios=[]
for cone,ks,w in ns['ns']['_prepared'](F(1,3)):
 for r,k in zip(cone,ks):
  norm=abs(r[1])+abs(r[2])
  if norm:ratios.append(k/norm)
assert min(ratios)==3
EPS=arb(1)/50;R=arb(1)/50
DERIVATIVE_BOUND=15/(1088*arb.pi()**3*(EPS-R/3)**2*R)
NODE_DISPLACEMENT=arb(1)/(32*2**NODE_BITS)
NODE_VALUE_ERROR=2*DERIVATIVE_BOUND*NODE_DISPLACEMENT
rho=arb(3);h=arb(1)/32;s=h*(rho-1/rho)/2;delta=EPS-s/3
M=15/(1088*arb.pi()**3*delta**2)
INTERPOLATION_ERROR=(2*DEGREE+2)*4*M*rho**(-DEGREE)/(rho-1)
TARGET=arb(13)/2500
# Positive roots floor to a dyadic, negative roots are exact negatives.
ROOTS=[None]*N
for j in range(DEGREE//2):
 x=(arb.pi()*(2*j+1)/(2*N)).cos();a=(x*2**NODE_BITS).floor().unique_fmpz();assert a is not None
 ROOTS[j]=F(int(a),2**NODE_BITS);ROOTS[DEGREE-j]=-ROOTS[j]
 assert abs(x-aq(ROOTS[j])) < arb(1)/2**NODE_BITS
ROOTS[DEGREE//2]=F(0)
# A converts values at exact Chebyshev roots to Chebyshev coefficients.
A=arb_mat([[((arb.pi()*k*(2*j+1)/(2*N)).cos()*(1 if k==0 else 2)/N) for j in range(N)] for k in range(N)])
# Exact power coefficients of T_k(2u-1), then exact degree-20 Bernstein transform.
P=[[F(1)],[F(-1),F(2)]]
for k in range(2,N):
 p=[F(0)]*(k+1)
 for j,a in enumerate(P[-1]):p[j]-=2*a;p[j+1]+=4*a
 for j,a in enumerate(P[-2]):p[j]-=a
 P.append(p)
BF=[[sum((P[k][l]*F(comb(j,l),comb(DEGREE,l)) for l in range(min(j,k)+1)),F(0)) for k in range(N)] for j in range(N)]
B=arb_mat([[aq(x) for x in row] for row in BF]);TRANSFORM=B*A
# Verify the rational conversion at u=0,1 and midpoint using exact arithmetic.
for k in range(N):
 assert BF[0][k]==(-1)**k and BF[-1][k]==1
 atmid=sum((BF[j][k]*F(comb(DEGREE,j),2**DEGREE) for j in range(N)),F(0))
 assert atmid==sum((P[k][l]*F(1,2**l) for l in range(k+1)),F(0))

def node_matrix(i,j,use_cache=True):
 path=CACHE/f'patch_{i:02}_{j:02}.json'
 if use_cache and path.exists():
  rec=json.loads(path.read_text());assert rec['config']==CONFIG and rec['patch']==[i,j]
  return arb_mat([[unpack(x) for x in row] for row in rec['values']]),True
 ys=[F(2*i+1,32)+t/32 for t in ROOTS];zs=[F(2*j+1,32)+t/32 for t in ROOTS]
 values=[]
 for y in ys:
  row=[]
  for z in zs:
   val=ns['point_divided_difference'](y,z)
   assert val.is_finite()
   val+=arb(0,NODE_VALUE_ERROR.upper())
   row.append(val)
  values.append(row)
 rec={'config':CONFIG,'patch':[i,j],'values':[[pack(v) for v in row] for row in values]}
 path.write_text(json.dumps(rec,separators=(',',':'))+'\n')
 return arb_mat(values),False

def split(A,axis):
 L=arb_mat(N,N);R=arb_mat(N,N)
 for j in range(N):
  work=[A[k,j] if axis==0 else A[j,k] for k in range(N)]
  left=[work[0]];right=[work[-1]]
  for _ in range(DEGREE):
   work=[(work[k]+work[k+1])/2 for k in range(len(work)-1)]
   left.append(work[0]);right.append(work[-1])
  right.reverse()
  for k in range(N):
   if axis==0:L[k,j]=left[k];R[k,j]=right[k]
   else:L[j,k]=left[k];R[j,k]=right[k]
 return L,R

def certify_matrix(C):
 stack=[('',C)];accepted=[];failed=[];tested=0;min_root=lower_matrix(C)
 while stack:
  path,A=stack.pop();tested+=1;bound=lower_matrix(A)-INTERPOLATION_ERROR
  if bound>TARGET:
   accepted.append({'path':path,'H_lower':pack(bound.lower()),'display':str(bound.lower())});continue
  if len(path)>=MAX_DEPTH:
   failed.append({'path':path,'H_lower':pack(bound.lower()),'display':str(bound.lower())});continue
  L,R=split(A,len(path)%2);stack.append((path+'1',R));stack.append((path+'0',L))
 assert tested<=511 and len(accepted)+len(failed)<=256
 return {'passed':not failed,'nodes_tested':tested,'accepted_leaves':accepted,'failed_leaves':failed,'root_polynomial_lower':str(min_root),'max_depth_used':max(map(lambda r:len(r['path']),accepted+failed))}

def check_patch(i,j,use_cache=True):
 start=time.perf_counter();V,cachehit=node_matrix(i,j,use_cache);node_time=time.perf_counter()-start
 C=TRANSFORM*V*TRANSFORM.transpose();res=certify_matrix(C)
 res.update({'patch':[i,j],'closed_box':[f'{i}/16',f'{i+1}/16',f'{j}/16',f'{j+1}/16'],'cache_hit':cachehit,'node_seconds':node_time,'total_seconds':time.perf_counter()-start,'largest_node_radius':str(max(V[a,b].rad() for a in range(N) for b in range(N))),'largest_Bernstein_radius':str(max(C[a,b].rad() for a in range(N) for b in range(N)))})
 return res

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['four','all'],default='four');parser.add_argument('--fresh',action='store_true');args=parser.parse_args()
 patches=[(15,7),(15,15),(8,3),(0,0)] if args.stage=='four' else [(i,j) for i in range(16) for j in range(i+1)]
 start=time.perf_counter();results=[]
 for i,j in patches:
  res=check_patch(i,j,not args.fresh);results.append(res)
  lo=min((unpack(r['H_lower']) for r in res['accepted_leaves']),default=arb(0))
  print(json.dumps({'patch':[i,j],'passed':res['passed'],'nodes_tested':res['nodes_tested'],'depth':res['max_depth_used'],'accepted_H_lower':str(lo),'seconds':res['total_seconds'],'cache_hit':res['cache_hit']}),flush=True)
 rec={'status':'PASS continuous lower certificate' if all(r['passed'] for r in results) else 'INCOMPLETE, not counterexample','stage':args.stage,'config':CONFIG,'max_binary_depth':MAX_DEPTH,'max_nodes_per_patch':511,'derivative_bound':str(DERIVATIVE_BOUND),'node_value_displacement_error':str(NODE_VALUE_ERROR),'ellipse_M':str(M),'interpolation_error':str(INTERPOLATION_ERROR),'source_sha256':SOURCE_SHA,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':platform.python_version(),'wall_seconds':time.perf_counter()-start,'patches':results,'coverage_complete':args.stage=='all' and len(results)==136 and all(r['passed'] for r in results)}
 (OUT/f'cc11_cheb_{args.stage}_result.json').write_text(json.dumps(rec,indent=2)+'\n')
 print(json.dumps({k:v for k,v in rec.items() if k!='patches'},indent=2),flush=True)
if __name__=='__main__':main()

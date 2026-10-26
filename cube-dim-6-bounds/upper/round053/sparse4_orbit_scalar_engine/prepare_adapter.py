"""Structural conversion only: no Fourier/kernel evaluation or S2 scan."""
import json,gzip,hashlib,itertools,math,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
LIMIT=2**63-1
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def det(M):
 A=[list(map(int,r)) for r in M];n=len(A);prev=1;sign=1
 for k in range(n-1):
  z=next((j for j in range(k,n) if A[j][k]),None)
  if z is None:return 0
  if z!=k:A[k],A[z]=A[z],A[k];sign=-sign
  a=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    v=A[i][j]*a-A[i][k]*A[k][j];assert v%prev==0;A[i][j]=v//prev
  for i in range(k+1,n):A[i][k]=0
  prev=a
 return sign*A[-1][-1]
def frame():
 out=[((1,)+s,Q(1)) for s in itertools.product([-1,1],repeat=5)]
 for support in itertools.combinations(range(6),4):
  for signs in itertools.product([-1,1],repeat=3):
   w=[0]*6
   for j,v in zip(support,(1,)+signs):w[j]=v
   out.append((tuple(w),Q(1,80)))
 assert len(out)==152 and len({w for w,b in out})==152
 return out
def main():
 if (P/'adapter_results.json').exists():raise SystemExit('Refuse structural replay')
 tick=time.monotonic();reg=json.loads((P/'registration.json').read_text())
 for f,h in reg['source_sha256'].items():assert sha(Path(f))==h
 src=next(Path(f) for f in reg['source_sha256'] if f.endswith('global_fan.json.gz'));raw=gzip.decompress(src.read_bytes());assert hashlib.sha256(raw).hexdigest()==reg['fan_uncompressed_sha256']
 d=json.loads(raw);assert len(d['rays'])==478 and len(d['cones'])==13229 and d['beta']=='1/80'
 rays=[];directions=frame();qmax=0;rmax=0;scalemax=0
 for i,v in enumerate(d['rays']):
  assert v['id']==i and len(v['primitive'])==6
  r=v['primitive'];assert all(isinstance(z,int) for z in r) and math.gcd(*r)==1
  q=Q(v['q_true_152']);assert q>0
  assert sum(beta*abs(sum(a*b for a,b in zip(w,r))) for w,beta in directions)==q
  scale=q.denominator;R=[scale*z for z in r];qq=q.numerator
  assert all(abs(z)<=LIMIT for z in R) and qq<=LIMIT and 80%scale==0
  assert all(Q(z,qq)==Q(a)/q for z,a in zip(R,r))
  rays.append(dict(R=R,Q=qq,scale=scale));qmax=max(qmax,qq);rmax=max(rmax,max(map(abs,R)));scalemax=max(scalemax,scale)
 cones=[];maxdet=0;maxdetprim=0;weight_digest=hashlib.sha256()
 for v in d['cones']:
  ids=v['ray_ids'];assert len(ids)==len(set(ids))==6 and all(0<=i<478 for i in ids)
  mat=[[d['rays'][ids[j]]['primitive'][k] for j in range(6)] for k in range(6)]
  dd=abs(det(mat));assert dd==v['abs_primitive_determinant']>0
  scaled=dd*math.prod(rays[i]['scale'] for i in ids);assert scaled<=LIMIT
  w=Q(dd,1)
  for i in ids:w/=Q(d['rays'][i]['q_true_152'])
  ws=Q(scaled,math.prod(rays[i]['Q'] for i in ids));assert w==ws
  maxdet=max(maxdet,scaled);maxdetprim=max(maxdetprim,dd);cones.append([scaled,ids]);weight_digest.update(f'{w}\n'.encode())
 # Independent full finite group action, no fan kernel loop.
 x=(50,50,50,45,27,0);counts=Counter()
 for perm in itertools.permutations(range(6)):
  for signs in itertools.product([-1,1],repeat=6):counts[tuple(signs[j]*x[perm[j]] for j in range(6))]+=1
 assert len(counts)==3840 and set(counts.values())=={12} and counts[x]==12
 pairs={min(z,tuple(-a for a in z)) for z in counts};assert len(pairs)==1920
 assert all(counts[z]==counts[tuple(-a for a in z)]==12 for z in counts)
 canonical=sorted(pairs);assert all(sum(counts[z] for z in [v,tuple(-a for a in v)])==24 for v in canonical)
 packed=dict(rays=rays,cones=cones,point_orbit_representatives=canonical,orbit_weight=24,point_denominator=50)
 digest=hashlib.sha256(json.dumps(packed,separators=(',',':'),sort_keys=True).encode()).hexdigest()
 maxdot=6*50*rmax;maxden=50*qmax;assert maxdot<=LIMIT and maxden<=LIMIT
 out=dict(status='PASS_STRUCTURAL_ADAPTER_ONLY_FAN_COVERAGE_PENDING_F050',all_rays_checked=478,all_simplices_checked=13229,true_rates_independently_checked=478,integer_column_scale_max=scalemax,max_integer_ray_coordinate=rmax,max_integer_rate=qmax,max_primitive_det=maxdetprim,max_scaled_det=maxdet,weight_product_method='arbitrary precision or sequential Arb division; never int64 product of six rates',all_scaled_determinants_int64=True,normalized_input_canonical_sha256=digest,all_weights_digest_sha256=weight_digest.hexdigest(),group_size=46080,stabilizer12_exact=counts[x],distinct_images=3840,negative_orbit_representatives=1920,weight_per_representative=24,actual_terms=13229*1920,logical_terms=13229*46080,max_proposed_point_dot_numerator_bound=maxdot,max_proposed_point_denominator_bound=maxden,frame_direction_count=152,new_true_S2_pairs=152*151//2,S2_preflight_executed=False,positive_T_values=0,full_kernel_calls=0,seconds=time.monotonic()-tick)
 (P/'adapter_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Four old T values; exact G reconstruction and positive-band local boxes."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os,sys,time,math
import flint
from flint import arb,fmpz,fmpq,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
OLD=PROJECT/'upper/round007/g_time_moments/results.json'
GFILE=PROJECT/'upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json'
SHIFT=PROJECT/'audit/round007/sparse_shift_bound.md'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def A(x):return arb(x.numerator)/x.denominator
def ball(x):return {'ball':x.str(100,more=True),'lower':x.lower().str(100,more=True),'upper':x.upper().str(100,more=True)}
def phi(x):return (12*(x/12).asinh()).exp()
def largest(pred):
 lo,hi=0,100001;n=0
 while hi-lo>1:
  mid=(lo+hi)//2;n+=1
  if pred(F(mid,1000000)):lo=mid
  else:hi=mid
 assert n<=30
 return F(lo,1000000),n
start=time.monotonic()
raw=json.loads(OLD.read_text());gr=json.loads(GFILE.read_text())
assert raw['center_numerators']==gr['point_numerators'] and raw['denominator']==gr['denominator']
R=fmpq(fmpz(gr['Q_over_pi5_numerator']),fmpz(gr['Q_over_pi5_denominator']))
G=arb(R)/arb.pi()**5
assert (G-arb(gr['G_arb_384'])).contains(0)
assert sha(SHIFT)=='c386ed307d9df8379c4dcca9983bf1571405b112b093c4bb9c9c0cd1dc04fede'
taus=[F(r['tau']) for r in raw['time_nodes']]; hs=[arb(r['H_arb384'])for r in raw['time_nodes']]
assert len(taus)==8 and all(taus[k]<taus[k+1] for k in range(7))
masses=[(hs[k]-hs[k+1]).lower()for k in range(7)];assert all(v>0 for v in masses)
cc=[F(n,gr['denominator'])for n in gr['point_numerators'][1:]]
rows=[]
for j in [4,5,6,7]:
 T=taus[j];m0=T/4000;threshold=T/4096;KT=G-hs[j]
 def lower(r):return sum((masses[k]/phi(A(r/taus[k]))for k in range(j)),arb(0))
 def cover(r):
  if not r:return F(0),[]
  bs=[(max(F(0),c-r),min(F(1),c+r))for c in cc]
  assert all(bs[i][1]<bs[j][0] or bs[j][1]<bs[i][0] for i in range(5)for j in range(i+1,5))
  vol=F(120)
  for lo,hi in bs:vol*=hi-lo
  assert 0<vol<1
  return vol,[[str(lo),str(hi)]for lo,hi in bs]
 if lower(F(0))>A(m0):
  radius,niter=largest(lambda r:lower(r)>A(m0));lb=lower(radius)
  assert lb>A(m0) and not lower(radius+F(1,1000000))>A(m0)
 else:radius=F(0);niter=0;lb=lower(F(0))
 rbare,nbare=largest(lambda r:lower(r)>A(threshold))if lower(F(0))>A(threshold)else(F(0),0)
 vol,bs=cover(radius)
 rows.append({'T':str(T),'K_T_center':ball(KT),'formal_point_cost':ball(A(T)/(64*KT)),
 'threshold_for_64':str(threshold),'m0_for_cost_125_over_2':str(m0),
 'center_below_threshold_excludes_this_T':bool(KT<A(threshold)),
 'old_time_index':j,'band_start':str(taus[0]),'omitted_nonnegative_time':'[0,1/1000]',
 'band_lower_at_center':ball(lower(F(0))),'local_band_radius_m0':str(radius),
 'local_band_lower_at_radius':ball(lb),'local_strict_pass':bool(radius>0 and lb>A(m0)),
 'radius_for_bare_64_threshold':str(rbare),'bisection_steps':niter,'bare_bisection_steps':nbare,
 'clipped_box_free_bounds':bs,'positive_face_permutation_images':120 if radius else 0,
 'exact_covered_face_volume':str(vol),'exact_uncovered_face_volume':str(1-vol),
 'covered_face_volume_decimal':ball(A(vol)),'global_m_T_certified':False})
inputs=[OLD,GFILE,SHIFT,HERE/'registration.md',Path(__file__)]
out={'status':'PASS_OLD_DATA_POSTPROCESS_LOCAL_ONLY','pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,
 'python':sys.version,'python_executable':sys.executable,'python_flint':flint.__version__,'arb_bits':ctx.prec,
 'new_kernel_passes':0,'T_count':4,'exact_G_reconstructed':ball(G),'center_numerators':gr['point_numerators'],
 'center_denominator':gr['denominator'],'rows':rows,'global_m_T_certified':False,'new_C6_bound':False,
 'inputs':[{'path':str(p.resolve()),'sha256':sha(p)}for p in inputs],'source_sha256':sha(__file__),'exit_code':0}
(HERE/'postprocess_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'elapsed_seconds':out['elapsed_seconds'],
 'rows':[{'T':r['T'],'point_cost':r['formal_point_cost']['ball'],'radius':r['local_band_radius_m0'],
 'volume':r['covered_face_volume_decimal']['ball']}for r in rows]},indent=2))

"""Read saved data; direct Fourier + explicit 2^8 event enumeration checks the C++ interface."""
from pathlib import Path
from hashlib import sha256
from math import pi,comb,isqrt
import json
import numpy as np
H=Path(__file__).resolve().parent
R=json.loads((H/'bernoulli_weak_endpoint_guard_results_20261007.json').read_text())
OUT=H/'bernoulli_weak_endpoint_saved_interface_review_results_20261007.json'
assert not OUT.exists()
checks=[]; differences=[]
for rd in R['rounds']:
 for f in rd['saved_profiles']:
  p=H/f['file']
  checks.append(sha256(p.read_bytes()).hexdigest()==f['sha256'])
for c in (1.,0.5):
 data=np.load(H/f'bernoulli_weak_endpoint_n8_c{c:g}_rep0_20261007.npz')
 x=data['receiver'][[0,7,31]]; rr=data['r_nodes']; n=8; mid=4; sh=isqrt(n)//2
 vtx=((np.arange(2**n)[:,None]>>np.arange(n))&1).astype(float)
 cnt=vtx.sum(1)
 pc=comb(n,mid)/2**n
 ps=(comb(n,mid-sh)+comb(n,mid+sh))/2**n
 pb=sum(comb(n,k) for k in (mid-1,mid,mid+1))/2**n
 shell=(cnt==mid)/pc
 shifted=((cnt==mid-sh)|(cnt==mid+sh))/ps
 band=((cnt>=mid-1)&(cnt<=mid+1))/pb
 direct=[]
 for xx in x:
  allr=np.zeros((4,len(rr)))
  values=[]
  for per in (1.,2.,4.):
   k=np.arange(1,16385,2); vv=(2*pi*k/per)**2
   gh=1/(1+c*(vv/np.log1p(vv)-1))
   gg=0.5+(2/pi)*np.sum((gh/k)[:,None]*np.sin(2*pi*k[:,None]*xx[None,:]/per),axis=0)
   bits=((xx%per)<per/2).astype(float)
   responses=np.zeros((3,len(rr)))
   for ir,r in enumerate(rr):
    prob=(1-r)*bits+r*gg
    atom=np.prod(np.where(vtx==1,prob,1-prob),axis=1)
    checks.append(abs(atom.sum()-1)<1e-12)
    responses[:,ir]=[atom@shell,atom@shifted,atom@band]
   values.append(responses)
  allr[:3]=values[2]
  allr[3]=(values[0][2]+values[1][1]+values[2][0])/3
  direct.append(allr.max(1))
 direct=np.asarray(direct).T
 saved=data['profile'][:,[0,7,31]]
 diff=np.abs(direct-saved)
 reports=next(q for q in R['rounds'] if q['n']==8)['reports']
 errors=np.asarray([v['additive_profile_error_analytic_plus_float_screen'] for v in reports if v['c']==c])
 checks.append((diff<=errors[:,None]+1e-12).all())
 differences.append({'c':c,'max_direct_enumeration_profile_difference':float(diff.max())})
assert all(checks)
res={'status':'PASS_SAVED_INTERFACE_REVIEW','checks':len(checks),
     'direct_source_event_enumeration':'2^8 for 3 saved receiver points, all 129 r nodes, both c',
     'one_dimensional_evaluation':'Direct original Fourier series, independent of FFT interpolation',
     'max_differences':differences,
     'all_24_profile_hashes_match':True,
     'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
     'scope':'No main data rerun, not interval arithmetic, not a weak endpoint proof.'}
OUT.write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res))

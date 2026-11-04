import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import json,itertools,time,multiprocessing as mp
from pathlib import Path
import numpy as np
from scipy.stats import qmc
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'exploration005'))
from lattice_batch import fit
base=Path(__file__).resolve().parent;out=base/'rescramble';out.mkdir(exist_ok=True)
selected=json.loads((out/'selected.json').read_text());N=int(os.environ.get('SOBOL_POWER','16'));CHUNK=1024

def worker(idx):
 start=time.time();u=selected['u'];k=selected['k'];centers=np.array(selected['centers']);dims=np.array(selected['halfwidths']);counts=[];witness=[]
 for site,(center,d) in enumerate(zip(centers,dims)):
  boxes=[]
  for ctr,w in zip(centers,dims):
   for ij in itertools.product((-1,0,1),repeat=2):
    shifted=ctr+k*np.array(ij)
    if np.max(np.abs(shifted-center))<1.5+2*np.max(dims):boxes.append((shifted,w))
  seed=820000+100*idx+site;points=center+(2*qmc.Sobol(2,scramble=True,seed=seed).random_base2(N)-1)*d;count=0
  for i in range(0,len(points),CHUNK):
   for pt in points[i:i+CHUNK]:
    value=fit(pt,boxes)
    if value < u*u*(1-1e-11):
     count+=1
     if len(witness)<16:witness.append({'site':site,'point':pt.tolist(),'best_defect_approx':float(value)})
   partial={'id':idx,'site':site,'points':min(i+CHUNK,len(points)),'count':count,'completed_counts':counts,'elapsed':time.time()-start,'pid':os.getpid()}
   (out/f'replicate_{idx:02d}_progress.json').write_text(json.dumps(partial))
  counts.append(count)
 cap=sum(float(np.prod(2*dims[i]))*counts[i]/len(points) for i in range(len(dims)))/k**2;gain=cap/(1-selected['q']);delta=gain-selected['baseline']
 result={'id':idx,'pid':os.getpid(),'counts':counts,'samples_per_hole':len(points),'gain':gain,'difference':delta,'third_scale':delta/u**4,'elapsed':time.time()-start,'accepted_examples':witness,'note':'independent scrambled QMC estimate; no certified error bound; r<=1.5'}
 (out/f'replicate_{idx:02d}.json').write_text(json.dumps(result,indent=2));return result
if __name__=='__main__':
 print(json.dumps({'master':os.getpid(),'fixed_geometry':'rescramble/selected.json','replicates':15,'samples_per_hole':2**N}),flush=True)
 with mp.Pool(15) as pool:
  for r in pool.imap_unordered(worker,range(15)): print(json.dumps({k:r[k] for k in ('id','pid','third_scale','elapsed')}),flush=True)

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,sys,time,random,math
from pathlib import Path
from fractions import Fraction as Q
from concurrent.futures import ProcessPoolExecutor
import numpy as np
FACES=[0,4,8]
def overlaps(lo,hi):return [max(0,min(hi,4*(k+1))-max(lo,4*k)) for k in [0,1]]
def average(f,lo,hi):
 t=f
 for i in reversed(range(f.ndim)):t=np.tensordot(t,np.array(overlaps(lo[i],hi[i]),dtype=np.int64),axes=([-1],[0]))
 return Q(int(t),8*math.prod(b-a for a,b in zip(lo,hi)))
def vertex(f,a,x):return average(f,[min(u,v) for u,v in zip(a,x)],[max(u,v) for u,v in zip(a,x)])
def gridmax(f,x):
 n=f.ndim;t=f
 for i in reversed(range(n)):
  M=np.array([overlaps(min(a,x[i]),max(a,x[i])) for a in FACES],dtype=np.int64)
  t=np.moveaxis(np.tensordot(M,t,axes=([1],[i])),0,i)
 dens=np.array([1],dtype=np.int64)
 for xi in x:dens=np.multiply.outer(dens,[abs(xi-a) for a in FACES])
 best=Q(-1);ix=None
 for k,(a,b) in enumerate(zip(t.flat,dens.flat)):
  v=Q(int(a),8*int(b))
  if v>best:best=v;ix=k
 anchor=[FACES[i] for i in np.unravel_index(ix,(3,)*n)]
 assert vertex(f,anchor,x)==best
 return best,anchor
def run(job):
 n,seed=job;rng=random.Random(seed);f=np.array([rng.randrange(17) for _ in range(2**n)],dtype=np.int64).reshape((2,)*n)
 x=[rng.choice([-1,1,3,5,7,9]) for _ in range(n)];mx,anc=gridmax(f,x);rects=[]
 for _ in range(128):
  lo=[rng.randint(min(-4,xi-8),xi) for xi in x];hi=[rng.randint(xi,max(12,xi+8)) for xi in x]
  for i in range(n):
   if lo[i]==hi[i]:hi[i]+=1
  v=average(f,lo,hi);assert v<=mx
  rects.append(dict(lower=lo,upper=hi,average=str(v)))
 chains=[]
 for _ in range(8):
  a=[rng.choice([z for z in range(-4,13) if z!=xi]) for xi in x];vs=[vertex(f,a,x)];aa=[a.copy()]
  for i in range(n):
   cands=[]
   for face in FACES:
    b=a.copy();b[i]=face;cands.append((vertex(f,b,x),b))
   v,b=max(cands,key=lambda t:t[0]);assert v>=vs[-1];a=b;vs.append(v);aa.append(a.copy())
  assert all(z in FACES for z in a) and vs[-1]<=mx
  chains.append(dict(anchors=aa,averages=[str(v) for v in vs]))
 return dict(N=n,seed=seed,source_heights_eighths=f.flatten().tolist(),receiver_units_quarters=x,
  grid_maximum=str(mx),maximizing_grid_anchor_units=anc,rectangles=rects,rounding_chains=chains)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[(n,79000+100*n+j) for n in [2,3,4,6,8] for j in range(8)];t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:r=list(pool.map(run,jobs))
 d=dict(workers=12,seconds=time.monotonic()-t,coordinate_scale=4,source_height_scale=8,cases=r)
 (out/'grid79.json').write_text(json.dumps(d,separators=(',',':'))+'\n')
 print('PASS',len(r),'exact grid maxima;',len(r)*128,'arbitrary rectangles;',sum(8*z['N'] for z in r),'endpoint moves;',d['seconds'],'seconds')

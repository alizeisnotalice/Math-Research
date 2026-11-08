import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time,random
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import linprog
import allocate77 as old
old.envelopes=lru_cache(None)(old.envelopes)
CELLS=[[-1,0,3,4],[8,9,4,5],[3,4,-1,0],[4,5,8,9],[0,1,0,1],[7,8,7,8]]
CORNERS=[[0,0],[0,8],[8,0],[8,8]]
def cases():
 out=[]
 for shape,kind in [('uniform','single'),('uniform','corners'),('uniform','five'),('notch','corners'),('checker','corners'),('cross','corners')]:
  f=[Q(2)]*64 if shape=='uniform' else old.f_values(shape,Q(1,2) if shape=='notch' else Q(1),2)
  out.append(dict(name=shape+'_'+kind,source=[str(x) for x in f],anchors=[[4,4]] if kind=='single' else CORNERS+([[4,4]] if kind=='five' else []),seed=None))
 for seed in [78001,78002,78003,78004]:
  rng=random.Random(seed);unit=[Q(5+rng.randrange(4),4) for _ in range(16)]
  f=[unit[(i//2)*4+j//2] for i in range(8) for j in range(8)]
  out.append(dict(name='random_'+str(seed),source=[str(x) for x in f],anchors=CORNERS,seed=seed))
 for c in out:
  F=list(map(Q,c['source']));valid=[]
  for a in c['anchors']:
   v=[i for i,p in enumerate(CELLS) if all(old.exact_dot(row,F)>Q(10011,10000)*vol for row,vol in old.rows_at(a,p,2))]
   valid.append(v)
  allcells=sorted(set(i for v in valid for i in v))
  assert allcells==list(range(len(CELLS))),(c['name'],allcells)
  c.update(receiver_cells=CELLS,eligible=valid,sigma='1001/1000',refinement=2,mesh=32)
 return out
def task(x):
 ci,ai,mask,c=x;F=list(map(Q,c['source']));a=c['anchors'][ai];selected=[i for i in range(len(CELLS)) if mask>>i&1]
 rows=list(dict.fromkeys(z for i in selected for z in old.rows_at(a,CELLS[i],2)))
 ans=old.lp_allocate(F,rows,Q(c['sigma']),Q(1,4),c['mesh'])
 ans.update(case=ci,anchor_index=ai,mask=mask)
 return ans
def integer_dp(cols,m,bound):
 n=1<<m;d=[None]*n;prev=[None]*n;d[0]=Q(0)
 for mask in range(n):
  if d[mask] is None:continue
  for j,c in enumerate(cols):
   nxt=mask|c['mask']
   if nxt==mask:continue
   v=d[mask]+Q(c[bound])
   if d[nxt] is None or v<d[nxt]:d[nxt]=v;prev[nxt]=(mask,j)
 assert d[-1] is not None
 chosen=[];mask=n-1
 while mask:mask,j=prev[mask];chosen.append(j)
 return d[-1],chosen,[None if x is None else str(x) for x in d]
def master(c,cols):
 m=len(c['receiver_cells']);B=np.array([[int(z['mask']>>i&1) for z in cols] for i in range(m)])
 costs=[[float(Q(z[k])) for z in cols] for k in ['lower','upper']]
 opt=dict(method='highs',options={'threads':1,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
 lo=linprog(costs[0],A_ub=-B,b_ub=[-1]*m,bounds=(0,None),**opt)
 up=linprog(costs[1],A_ub=-B,b_ub=[-1]*m,bounds=(0,None),**opt)
 assert lo.success and up.success
 w=[Q(max(0,math.ceil(x*10**9)),10**9) for x in up.x]
 cover=[sum((w[j] for j,z in enumerate(cols) if z['mask']>>i&1),Q(0)) for i in range(m)]
 factor=max([Q(1)]+[1/v for v in cover]);w=[x*factor for x in w]
 y=[Q(max(0,math.floor(-x*10**9)),10**9) for x in lo.ineqlin.marginals]
 factor=min([Q(1)]+[Q(z['lower'])/s for z in cols if (s:=sum((y[i] for i in range(m) if z['mask']>>i&1),Q(0)))>0])
 y=[x*factor for x in y]
 assert all(sum((y[i] for i in range(m) if z['mask']>>i&1),Q(0))<=Q(z['lower']) for z in cols)
 assert all(sum((w[j] for j,z in enumerate(cols) if z['mask']>>i&1),Q(0))>=1 for i in range(m))
 lb=sum(y,Q(0));ub=sum((x*Q(z['upper']) for x,z in zip(w,cols)),Q(0))
 il,_,dl=integer_dp(cols,m,'lower');iu,chosen,du=integer_dp(cols,m,'upper')
 assert lb<=ub and il<=iu
 return dict(fractional_lower=str(lb),fractional_upper=str(ub),weights=[[j,str(x)] for j,x in enumerate(w) if x],dual=[str(x) for x in y],
  integer_lower=str(il),integer_upper=str(iu),integer_columns=chosen,DP_lower=dl,DP_upper=du,
  certified_integer_over_fractional_lower=float(il/ub),certified_integer_over_fractional_upper=float(iu/lb),
  numerical_fractional_lower_model=lo.fun,numerical_fractional_upper_model=up.fun)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);cs=cases();jobs=[]
 for ci,c in enumerate(cs):
  for ai,v in enumerate(c['eligible']):
   em=sum(1<<i for i in v)
   jobs.extend((ci,ai,mask,c) for mask in range(1,1<<len(CELLS)) if mask&em==mask)
 t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:solved=list(pool.map(task,jobs))
 for ci,c in enumerate(cs):
  c['columns']=[z for z in solved if z['case']==ci];c['cover']=master(c,c['columns'])
  print(c['name'],len(c['columns']),'IP/FP',c['cover']['certified_integer_over_fractional_lower'],c['cover']['certified_integer_over_fractional_upper'],flush=True)
 doc=dict(workers=12,seconds=time.monotonic()-t,case_count=len(cs),column_count=len(solved),cases=cs)
 (out/'columns78.json').write_text(json.dumps(doc,separators=(',',':'))+'\n')
 print('PASS',len(cs),'complete finite column families;',len(solved),'certified source optimizations;',doc['seconds'],'seconds')

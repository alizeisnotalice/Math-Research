import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,sys,time,math,warnings,hashlib
from pathlib import Path
from fractions import Fraction as Q
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
warnings.filterwarnings('ignore',message='Unrecognized options detected.*')
DEN=10**12
def log_interval(t,m=20):
 u=(t-1)/(t+1);lo=sum((2*u**(2*k+1)/Q(2*k+1) for k in range(m)),Q(0))
 hi=lo+2*u**(2*m+1)/(Q(2*m+1)*(1-u*u))
 return Q((lo*DEN).__floor__(),DEN),Q((hi*DEN).__ceil__(),DEN)
def phi_interval(t):
 if t<=1:return t,t
 lo,hi=log_interval(t);return t*(1+lo),t*(1+hi)
def envelopes(mesh):
 nodes=[1+Q(k,mesh) for k in range(mesh+1)]
 logs=[log_interval(t) for t in nodes]
 low=[(Q(1),Q(0))]+[(2+l,-t) for t,(l,h) in zip(nodes,logs)]
 upper=[(Q(1),Q(0))]
 for a,b,la,lb in zip(nodes[:-1],nodes[1:],logs[:-1],logs[1:]):
  va=a*(1+la[1]);vb=b*(1+lb[1]);s=(vb-va)/(b-a);upper.append((s,va-s*a))
 return low,upper
EDGES=[-4,-2,0,1,2,3,4,5,6,7,8,9,10,12,16]
CELLS=[(a,b,c,d) for a,b in zip(EDGES[:-1],EDGES[1:]) for c,d in zip(EDGES[:-1],EDGES[1:])]
ANCHORS=[(i,j) for i in range(0,9,2) for j in range(0,9,2)]
def rows_at(anchor,cell,refine):
 step=2//refine
 ans=[]
 for x in [cell[0],cell[1]]:
  for y in [cell[2],cell[3]]:
   area=abs((x-anchor[0])*(y-anchor[1]))
   if area==0:continue
   lx,rx=sorted([x,anchor[0]]);ly,ry=sorted([y,anchor[1]])
   ox=[max(0,min(rx,k+step)-max(lx,k)) for k in range(0,8,step)]
   oy=[max(0,min(ry,k+step)-max(ly,k)) for k in range(0,8,step)]
   vec=tuple(a*b for a in ox for b in oy)
   ans.append((vec,area))
 return ans
def f_values(shape,eps,refine):
 ans=[]
 for i in range(4*refine):
  for j in range(4*refine):
   x=i//refine;y=j//refine
   ans.append(Q(1) if shape=='notch' and x==y==0 else
    1+eps if shape=='notch' else Q(2) if shape=='checker' and (x+y)%2==0 else
    Q(1,2) if shape=='checker' else Q(2) if x in [1,2] or y in [1,2] else Q(0))
 return ans
def exact_dot(a,b):return sum((x*y for x,y in zip(a,b) if x),Q(0))
def lp_allocate(F,rows,sigma,w,mesh):
 n=len(F);m=len(rows);low,upper=envelopes(mesh)
 ai=[];aj=[];av=[];rhs=[]
 for i,(a,v) in enumerate(rows):
  for j,z in enumerate(a):
   if z:ai.append(i);aj.append(j);av.append(-z)
  rhs.append(-float(sigma*v))
 for j in range(n):
  for k,(s,b) in enumerate(low):
   i=m+j*len(low)+k;ai.extend([i,i]);aj.extend([j,n+j]);av.extend([float(s),-1]);rhs.append(-float(b))
 A=coo_matrix((av,(ai,aj)),shape=(len(rhs),2*n)).tocsr()
 start=time.monotonic()
 res=linprog([0.]*n+[float(w)]*n,A_ub=A,b_ub=rhs,bounds=[(0,float(f)) for f in F]+[(0,None)]*n,
  method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9,'threads':1})
 assert res.success,res.message
 raw=[min(f,max(Q(0),Q(round(float(g)*10**9),10**9))) for f,g in zip(F,res.x[:n])]
 eta=Q(0)
 for a,v in rows:
  cur=exact_dot(a,raw);target=sigma*v
  if cur<target:
   den=exact_dot(a,F)-cur;assert den>0
   eta=max(eta,(target-cur)/den)
 assert eta<=1
 G=[(1-eta)*g+eta*f for g,f in zip(raw,F)]
 slacks=[exact_dot(a,G)-sigma*v for a,v in rows];assert min(slacks)>=0
 UB=sum((w*max(s*g+b for s,b in upper) for g in G),Q(0))
 mu=[Q(max(0,math.floor(float(-z)*10**10)),10**10) for z in res.ineqlin.marginals]
 lm=mu[:m];theta=[]
 for j in range(n):
  ts=mu[m+j*len(low):m+(j+1)*len(low)];total=sum(ts,Q(0))
  if total>w:ts=[q*w/total for q in ts]
  assert all(q>=0 for q in ts) and sum(ts,Q(0))<=w
  theta.append(ts)
 d=[sum((t*s for t,(s,b) in zip(ts,low)),Q(0))-sum((l*a[j] for l,(a,v) in zip(lm,rows) if l),Q(0)) for j,ts in enumerate(theta)]
 LB=sum((l*sigma*v for l,(a,v) in zip(lm,rows)),Q(0))+sum((t*b for ts in theta for t,(s,b) in zip(ts,low)),Q(0))+sum((min(Q(0),f*z) for f,z in zip(F,d)),Q(0))
 assert 0<=LB<=UB
 sparse_l=[[i,str(q)] for i,q in enumerate(lm) if q]
 sparse_t=[[j,k,str(q)] for j,ts in enumerate(theta) for k,q in enumerate(ts) if q]
 local=[f if any(a[j] for a,v in rows) else Q(0) for j,f in enumerate(F)]
 jl,jh=map(sum,zip(*(tuple(w*z for z in phi_interval(f)) for f in local)))
 return dict(primal=[str(g) for g in G],success_dual=sparse_l,epigraph_dual=sparse_t,
  lower=str(LB),upper=str(UB),relative_gap=float((UB-LB)/max(Q(1),UB)),repair_mixing=str(eta),
  minimum_vertex_slack=str(min(slacks)),localized_full_source_modular_upper=str(jh),
  LP_status=int(res.status),LP_seconds=time.monotonic()-start,LP_objective=float(res.fun),
  dual_row_count=m,epigraph_pieces=len(low))
def run(job):
 shape=job['shape'];eps=Q(job['epsilon']);refine=job['refinement'];sigma=1+min(eps/20,Q(1,1000));margin=eps/10000
 F=f_values(shape,eps,refine);w=Q(1,refine**2);valid={};score={}
 for ci,c in enumerate(CELLS):
  for ai,a in enumerate(ANCHORS):
   rows=rows_at(a,c,refine)
   if rows and all(exact_dot(z,F)>=(sigma+margin)*v for z,v in rows):
    valid.setdefault(ci,[]).append(ai);score[ci,ai]=min(float(exact_dot(z,F)/v) for z,v in rows)
 groups={}
 if job['assignment']=='margin':
  for ci,candidates in valid.items():
   ai=max(candidates,key=lambda a:(score[ci,a],-a));groups.setdefault(ai,[]).append(ci)
 else:
  left=set(valid)
  while left:
   areas=lambda ai:sum((CELLS[i][1]-CELLS[i][0])*(CELLS[i][3]-CELLS[i][2]) for i in left if ai in valid[i])
   ai=max(range(len(ANCHORS)),key=lambda a:(areas(a),-a));chosen=sorted(i for i in left if ai in valid[i]);assert chosen
   groups[ai]=chosen;left.difference_update(chosen)
 cert=[]
 for ai,cs in sorted(groups.items()):
  rows=list(dict.fromkeys(row for ci in cs for row in rows_at(ANCHORS[ai],CELLS[ci],refine)))
  z=lp_allocate(F,rows,sigma,w,job['mesh']);z.update(anchor=list(ANCHORS[ai]),receiver_cells=cs);cert.append(z)
 jl,jh=map(sum,zip(*(tuple(w*z for z in phi_interval(f)) for f in F)))
 lo=sum((Q(z['lower']) for z in cert),Q(0));hi=sum((Q(z['upper']) for z in cert),Q(0))
 area=sum((Q((CELLS[i][1]-CELLS[i][0])*(CELLS[i][3]-CELLS[i][2]),4) for i in valid),Q(0))
 return dict(job,sigma=str(sigma),source=[str(f) for f in F],source_cell_volume=str(w),
  receiver_cells=CELLS,covered_cells=sorted(valid),groups=cert,
  source_modular_lower=str(jl),source_modular_upper=str(jh),witness_area=str(area),
  allocation_lower=str(lo),allocation_upper=str(hi),reuse_lower=float(lo/jh),reuse_upper=float(hi/jl),
  localized_reuse_upper=float(sum((Q(z['localized_full_source_modular_upper']) for z in cert),Q(0))/jl),
  raw_anchor_reuse=len(cert),relative_primal_dual_gap=float((hi-lo)/max(Q(1),hi)),
  scope='Entire receiver cells certified. Allocation optimum is conditional on this fixed anchor/receiver grouping and source partition.')
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
 jobs=[dict(shape='notch',epsilon=e,refinement=r,assignment=a,mesh=32) for e in ['1/100','1/10','1/2','1'] for r in [1,2] for a in ['greedy','margin']]
 jobs +=[dict(shape=s,epsilon='1',refinement=r,assignment=a,mesh=32) for s in ['checker','cross'] for r in [1,2] for a in ['greedy','margin']]
 mode=sys.argv[2] if len(sys.argv)>2 else 'base'
 if mode=='refine':jobs=[dict(shape=s,epsilon=e,refinement=2,assignment=a,mesh=64) for s,e in [('notch','1/2'),('notch','1'),('checker','1'),('cross','1')] for a in ['greedy','margin']]
 suffix='_refine' if mode=='refine' else ''
 t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 doc=dict(workers=12,seconds=time.monotonic()-t,log_series_terms=20,log_rounding_denominator=DEN,primal_rounding_denominator=10**9,dual_rounding_denominator=10**10,rows=rows)
 (out/('allocation77'+suffix+'.json')).write_text(json.dumps(doc,indent=2)+'\n')
 summary=[{k:v for k,v in r.items() if k not in ['source','receiver_cells','covered_cells','groups']}|dict(groups=len(r['groups']),covered_cell_count=len(r['covered_cells'])) for r in rows]
 (out/('summary77'+suffix+'.json')).write_text(json.dumps(dict(workers=12,seconds=doc['seconds'],rows=summary),indent=2)+'\n')
 print(json.dumps(dict(PASS=True,cases=len(rows),groups=sum(len(r['groups']) for r in rows),seconds=doc['seconds'],max_gap=max(r['relative_primal_dual_gap'] for r in rows),max_reuse=max(r['reuse_upper'] for r in rows)),indent=2))

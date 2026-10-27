#!/usr/bin/env python3
from fractions import Fraction as F
from math import lcm
from pathlib import Path
from collections import defaultdict
import json,time,resource,hashlib,os
start=time.monotonic(); cpu=time.process_time()
p=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round037/compute');s=Path('/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round037/joint_cover/SPEC.md')
w=json.loads((p/'r037_failure_witness.json').read_text());g=json.loads((p/'r037_failure_geometry.json').read_text());nm=json.loads((p/'r037_no_move.json').read_text())
n,a,N=8,4,48;cubes=[(tuple(map(F,x['center'])),F(x['radius'])) for x in g['cube_centers_and_radii']];assert len(cubes)==N
scale=lcm(*(v.denominator for c,r in cubes for v in (*c,r))); boxes=[(int(c[0]*scale),int(c[1]*scale),int(r*scale)) for c,r in cubes]; H=[(2*a*r)**n for c,r in cubes];assert H==list(map(F,w['h']))
def atomize(selected):
 # Label bits 0..N-1 original cubes; N..2N-1 selected enlarged cubes.
 bb=boxes+[(boxes[i][0],boxes[i][1],8*boxes[i][2]) for i in selected];idx=list(range(N))+[N+i for i in selected]; out=defaultdict(int)
 cuts=sorted({0}|{r for _,_,r in bb})
 for lo,hi in zip(cuts,cuts[1:]):
  active=[(x,y,r,idx[k]) for k,(x,y,r) in enumerate(bb) if r>=hi]
  xs=sorted({x+sg*r for x,y,r,i in active for sg in (-1,1)});ys=sorted({y+sg*r for x,y,r,i in active for sg in (-1,1)})
  xx=[]; yy=[]
  for vals,axis,outaxis in [(xs,0,xx),(ys,1,yy)]:
   for l,u in zip(vals,vals[1:]):
    mid=l+u;mask=0
    for x,y,r,i in active:
     c=(x,y)[axis]
     if 2*(c-r)<mid<2*(c+r):mask|=1<<i
    outaxis.append((u-l,mask))
  W=(2*hi)**6-(2*lo)**6
  for dx,mx in xx:
   for dy,my in yy:
    m=mx&my
    if m&((1<<N)-1):out[m]+=dx*dy*W
 return {m:F(v,scale**n) for m,v in out.items()}
atoms=atomize([]); U=sum(atoms.values(),F()); assert U==F(w['union_volume'])
order=w['tie_order'];assert sorted(order)==list(range(N));rank={i:k for k,i in enumerate(order)};tau=F(w['tau']); sel=[]; psel={}; vs={}; remaining=set(range(N))
while True:
 eligible=[i for i in remaining if sum((psel[j] for j in sel if max(abs(cubes[i][0][k]-cubes[j][0][k]) for k in range(2))<7*cubes[j][1]),F())<=tau]
 if not eligible:break
 i=min(eligible,key=lambda j:(-cubes[j][1],rank[j]));mask=sum(1<<j for j in sel)
 v=sum((vol for m,vol in atoms.items() if m>>i&1 and not m&mask),F()); vs[i]=v;psel[i]=v/H[i];sel.append(i);remaining.remove(i)
assert all(max(abs(cubes[i][0][k]-cubes[j][0][k]) for k in range(2))!=7*cubes[j][1] for i in range(N) for j in sel)
assert sel==w['selected_indices'] and {str(i):str(v) for i,v in psel.items()}==g['selected_p'] and {str(i):str(v) for i,v in vs.items()}==g['selected_increment_v']
S=sum(1<<j for j in sel); ea=atomize(sel); rej=defaultdict(F)
for mask,vol in ea.items():
 qm=mask&((1<<N)-1)
 if not qm&S:
  em=mask>>N;rej[em&S]+=vol
assert sum(atoms.values(),F())==sum(ea.values(),F());flow=defaultdict(F)
for em,vol in rej.items():
 P=sum((psel[j] for j in sel if em>>j&1),F());assert P>tau
 for j in sel:
  if em>>j&1:flow[j]+=vol*psel[j]/P
q=[(vs[i]+flow[i])/H[i] if i in sel else F(0) for i in range(N)];assert q==list(map(F,w['q'])) and sum(H[i]*q[i] for i in range(N))==U
assert set(sel)=={i for i in range(N) if q[i]>0};assert all(H[i]==max(H) for i in sel)
# Strict witness box intersection, all n dimensions; inactive centers zero.
adj=[set() for _ in range(N)]
for i,(ci,ri) in enumerate(cubes):
 for j in range(i):
  cj,rj=cubes[j]
  if all(abs(ci[k]-cj[k])<a*(ri+rj) for k in range(2)):
   adj[i].add(j);adj[j].add(i)
cliques=[];nodes=0
def bk(R,P,X):
 global nodes
 nodes+=1
 if not P and not X:cliques.append(tuple(sorted(R)));return
 pivot=max(P|X,key=lambda v:len(P&adj[v])) if P|X else None
 for v in sorted(P-(adj[pivot] if pivot is not None else set())):
  bk(R|{v},P&adj[v],X&adj[v]);P.remove(v);X.add(v)
bk(set(),set(range(N)),set()); cliques.sort();assert len(cliques)==40 and cliques==sorted(tuple(x) for x in w['maximal_cliques'])
assert any(len(C)<N for C in cliques)
load=max(sum((q[i] for i in C),F()) for C in cliques);assert load==F(w['L0_exact'])==F(w['Lalg_exact'])
z=list(map(F,w['M_dual_feasible_vector']));assert all(x>=0 for x in z) and all(sum((z[i] for i in C),F())<=1 for C in cliques);Mlo=sum((H[i]*z[i] for i in range(N)),F()); assert Mlo==F(w['M_lower_exact']);optup=U/Mlo;assert optup==F(w['Lopt_upper_exact'])
assert 3 not in adj[2] and H[2]==H[3] and Mlo==H[2]+H[3]
q_alt=[F(0)]*N;q_alt[2]=q_alt[3]=U/(H[2]+H[3])
assert sum(H[i]*q_alt[i] for i in range(N))==U
assert max(sum((q_alt[i] for i in C),F()) for C in cliques)==optup
assert load-optup>=load/F(10000) and 2*load>load+optup
# For any clique anchor has maximal h within C: nonzero q has globally maximal h.
# Thus sum q - sum hq/max_C h == 0 for every clique, including all subcliques.
assert all(H[i]==max(H) for i in range(N) if q[i]>0)
for C in cliques:
 hmax=max(H[i] for i in C);assert sum((q[i] for i in C),F())==sum((H[i]*q[i] for i in C),F())/hmax
files=[p/x for x in ['r037_failure_witness.json','r037_failure_geometry.json','r037_no_move.json','r037_failure_witness.py','r037_failure_geometry.py','r037_no_move.py']]+[s,Path(__file__)]
out={'scope':'one frozen rational geometry and one listed legal tie order; not registered batch','checks':{'exact_union':True,'selection':True,'increments':True,'rejected_incidence_and_flow':True,'payment_identity':True,'complete_maximal_cliques':True,'feasible_dual_packing':True,'strict_gap':True,'explicit_alternative_attains_dual_bound':True,'no_positive_s_for_any_clique':True},'n':n,'N':N,'scale':scale,'tau':str(tau),'selected':sel,'union_volume':str(U),'selected_increment':{str(i):str(vs[i]) for i in sel},'selected_p':{str(i):str(psel[i]) for i in sel},'rejected_classes':{str(k):str(v) for k,v in rej.items()},'flow':{str(i):str(flow[i]) for i in sel},'q_positive':{str(i):str(q[i]) for i in sel},'clique_count':len(cliques),'bronkerbosch_nodes':nodes,'M_lower':str(Mlo),'explicit_alternative_q_2_q_3':str(optup),'L0':str(load),'Lopt_upper':str(optup),'gap':str(load-optup),'approx_gap':float(load-optup),'gap_relative':float((load-optup)/load),'sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in files},'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

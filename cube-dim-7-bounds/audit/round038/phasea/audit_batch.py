#!/usr/bin/env python3
import json,os,time,resource,hashlib,sys,traceback,re
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
sys.set_int_max_str_digits(0)
SRC=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round038/compute')
DEST=Path(__file__).resolve().parent/'results';DEST.mkdir(exist_ok=True)
NODECAP=2**20; CLIQUECAP=65536
class Censored(Exception):pass

def run(pathstr):
 st=time.monotonic();cpu=time.process_time();path=Path(pathstr); raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest();d=json.loads(raw);out={'input':str(path),'input_sha256':digest,'pid':os.getpid()}
 try:
  n=int(d['n']);a=F(d['a']);N=d['N'];boxes=[(list(map(F,x['center'])),F(x['radius'])) for x in d['geometry']];assert len(boxes)==N
  h=[(2*a*r)**n for _,r in boxes];assert h==list(map(F,d['witness_budgets'])) and all(x>0 for x in h)
  adj=[set() for _ in range(N)]
  for i,(ci,ri) in enumerate(boxes):
   for j in range(i):
    cj,rj=boxes[j];dims=max(len(ci),len(cj));assert dims<=n
    if all(abs((ci[k] if k<len(ci) else 0)-(cj[k] if k<len(cj) else 0))<a*(ri+rj) for k in range(dims)):
     adj[i].add(j);adj[j].add(i)
  cs=[];nodes=0
  def search(R,P,X):
   nonlocal nodes
   nodes+=1
   if nodes>NODECAP or len(cs)>=CLIQUECAP:raise Censored('node_or_clique_cap')
   if not P and not X:cs.append(tuple(sorted(R)));return
   u=max(P|X,key=lambda z:len(P&adj[z])) if P|X else None
   for v in sorted(P-(adj[u] if u is not None else set())):
    search(R|{v},P&adj[v],X&adj[v]);P.remove(v);X.add(v)
  try:search(set(),set(range(N)),set())
  except Censored as e:
   out.update(status='graph_censored',reason=str(e),graph_nodes=nodes,partial_cliques=len(cs));return finish(out,st,cpu)
  cs.sort();assert cs==sorted(tuple(x) for x in d['maximal_cliques'])
  U=F(d['U']);q=list(map(F,d['q_before']));assert len(q)==N and min(q)>=0 and sum((h[i]*q[i] for i in range(N)),F())==U
  def stats(q):
   w=[sum((q[i] for i in c),F()) for c in cs];return w,max(w,default=F()),sum((h[i]*q[i]**2 for i in range(N)),F())
  weights,L,E=stats(q);assert L==F(d['L_before'])
  steps=d['steps'];assert len(steps)<=128
  for step in steps:
   i,j=step['i'],step['j'];assert 0<=i<N and 0<=j<N and i!=j and j not in adj[i]
   assert q[i]==F(step['old_i']) and q[j]==F(step['old_j']); assert L==F(step['load_before']) and E==F(step['energy_before'])
   pay=h[i]*q[i]+h[j]*q[j];assert pay==F(step['pair_payment'])
   q[i]=F(step['new_i']);q[j]=F(step['new_j']);assert q[i]>=0 and q[j]>=0 and h[i]*q[i]+h[j]*q[j]==pay
   weights,newL,newE=stats(q);assert newL==F(step['load_after']) and newE==F(step['energy_after']) and (newL,newE)<(L,E)
   L,E=newL,newE
  assert q==list(map(F,d['q_after'])) and L==F(d['L_after']) and sum((h[i]*q[i] for i in range(N)),F())==U
  stop=d['stop'];assert stop in ('iteration_cap','stop_pair_coordinate')
  paircount=0
  if stop=='iteration_cap':assert len(steps)==128
  else:
   # For each nonedge, exact 1D minimization: maximum of constant, increasing and decreasing affine functions.
   contains=[set(c) for c in cs]
   for i in range(N):
    for j in range(i+1,N):
     if j in adj[i]:continue
     paircount+=1;pay=h[i]*q[i]+h[j]*q[j]
     if not pay:continue
     X=max((weights[k]-q[i] for k,C in enumerate(contains) if i in C),default=F())
     Y=max((weights[k]-q[j] for k,C in enumerate(contains) if j in C),default=F())
     R=max((weights[k] for k,C in enumerate(contains) if i not in C and j not in C),default=F())
     upper=pay/h[i]; cross=(pay+h[j]*(Y-X))/(h[i]+h[j]);x=max(F(),min(upper,cross))
     optimal=max(R,X+x,Y+(pay-h[i]*x)/h[j]); lo=max(F(),(pay-h[j]*(optimal-Y))/h[i]);hi=min(upper,optimal-X);assert lo<=hi
     x=max(lo,min(hi,pay/(h[i]+h[j])));y=(pay-h[i]*x)/h[j]
     e=E+h[i]*(x*x-q[i]*q[i])+h[j]*(y*y-q[j]*q[j])
     assert (optimal,e)>=(L,E),f'positive pair {i},{j}'
  out.update(status='pass',graph_nodes=nodes,cliques=len(cs),steps=len(steps),stop=stop,stop_nonedges_checked=paircount)
 except Exception as e:out.update(status='failure',error=repr(e),traceback=traceback.format_exc())
 return finish(out,st,cpu)

def finish(out,st,cpu):
 out.update(wall_seconds=time.monotonic()-st,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 target=DEST/(Path(out['input']).stem+'.audit.json');tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2)+'\n');os.replace(tmp,target);return out
if __name__=='__main__':
 files=[]
 for path in sorted(SRC.glob('phasea_*.json')):
  if not re.fullmatch(r'phasea_[ABC]_[0-4]_\d{2}_(?:partial|full_disjoint)\.json',path.name):continue
  target=DEST/(path.stem+'.audit.json')
  if target.exists() and json.loads(target.read_text())['status'] != 'failure':continue
  try:json.loads(path.read_bytes())
  except (OSError,ValueError):continue
  files.append(str(path))
 print('submitted',len(files),'workers',10,flush=True)
 with ProcessPoolExecutor(max_workers=10) as ex:
  for future in as_completed([ex.submit(run,p) for p in files]):
   r=future.result(); print(Path(r['input']).name,r['status'],r.get('steps'),flush=True)

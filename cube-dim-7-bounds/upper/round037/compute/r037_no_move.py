#!/usr/bin/env python3
"""Exact verification of final algorithm-D pool for frozen R037 witness."""
import json,hashlib,os,time,resource
from fractions import Fraction as Q
from pathlib import Path
from joint_core import load,maskids
p=Path(__file__).resolve();d=json.loads(p.with_name('r037_failure_witness.json').read_text());q=[Q(v) for v in d['q']];h=[Q(v) for v in d['h']];n=len(q);M=[sum((1<<i for i in c),0) for c in d['maximal_cliques']];adj=[0]*n
for C in M:
 for i in maskids(C):adj[i]|=C & ~(1<<i)
allmask=(1<<n)-1;start=time.monotonic();cpu=time.process_time();groups={}
for i in range(n):groups[adj[i]|(1<<i)]=groups.get(adj[i]|(1<<i),0)|(1<<i)
twins={v for v in groups.values() if v.bit_count()>1}
def calc(C):
 ids=tuple(maskids(C));T=sum((q[i] for i in ids),Q());k=min(ids,key=lambda i:(-h[i],i));after=sum((h[i]*q[i] for i in ids),Q())/h[k]
 return T,T-after,k,after,ids
allmetrics=[(C,calc(C)) for C in M];basis=sorted(allmetrics,key=lambda z:(-z[1][0],z[1][4]))[:32]+sorted(allmetrics,key=lambda z:(-z[1][1],-z[1][0],z[1][4]))[:32]
candidates=set()
for C,_ in basis:
 ids=tuple(maskids(C));candidates.add(C)
 for th in set(h[i] for i in ids):
  for C2 in (sum((1<<i for i in ids if h[i]<=th),0),sum((1<<i for i in ids if h[i]>=th),0)):
   if C2.bit_count()>1:candidates.add(C2)
ordinary=sorted(((C,calc(C)) for C in candidates if C.bit_count()>1),key=lambda z:(-z[1][1],-z[1][0],z[1][4]))[:64];pool={C for C,_ in ordinary}|twins;L=load(q,M);positive=0;safe=0;hidden=[];bad=[]
for C in sorted(pool):
 T,s,k,qk,ids=calc(C)
 if s<=0:continue
 positive+=1;common=allmask
 for i in ids:common&=adj[i]
 boundary=adj[k]&~(C|common)&allmask
 delta=max((sum((q[i] for i in maskids(D&boundary)),Q()) for D in M),default=Q())
 if s>=delta:safe+=1
 q2=q.copy()
 for i in ids:q2[i]=Q()
 q2[k]=qk;Lp=load(q2,M)
 if s<delta and Lp<L:hidden.append({'C':list(ids),'anchor':k,'saving':str(s),'boundary':str(delta),'new_load':str(Lp)})
 bad.append({'C':list(ids),'anchor':k,'saving':str(s),'boundary':str(delta),'passes_safety':s>=delta})
assert safe==0 and d['accepted_moves']==[]
assert all(h[i]==max(h) for i in range(n) if q[i]>0)
out={'case':'A/0/6/partial','candidate_pool_count':len(pool),'positive_savings_count':positive,'all_positive_payments_on_global_max_h':True,'safety_accepted_count':safe,'rejected_candidates':bad,'conservative_boundary_improving_moves':hidden,'status':'zero accepted moves independently reproduced against frozen clique list','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'witness_sha256':hashlib.sha256(p.with_name('r037_failure_witness.json').read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
p.with_name('r037_no_move.json').write_text(json.dumps(out,indent=2)+'\n');print('pool',len(pool),'positive',positive,'safe',safe,'boundary_conservative',len(hidden))

"""R037 algorithm D, deterministic boundary-safe clique compression."""
from fractions import Fraction as Q
from joint_core import maskids,load

def indices(C):return tuple(maskids(C))
def compress(q,h,adj,cliques,steps=24):
 q=list(q);n=len(q);full=(1<<n)-1;groups={}
 for i in range(n):groups.setdefault(adj[i]|(1<<i),0);groups[adj[i]|(1<<i)]|=1<<i
 twins={C for C in groups.values() if C.bit_count()>1};moves=[]
 def metrics(C):
  ids=indices(C);T=sum((q[i] for i in ids),Q());k=min(ids,key=lambda i:(-h[i],i));after=sum((h[i]*q[i] for i in ids),Q())/h[k]
  return T,T-after,k,after,ids
 for it in range(steps):
  L=load(q,cliques);mdata=[(C,metrics(C)) for C in cliques]
  byT=sorted(mdata,key=lambda z:(-z[1][0],z[1][4]))[:32]
  byS=sorted(mdata,key=lambda z:(-z[1][1],-z[1][0],z[1][4]))[:32]
  candidates=set()
  for C,_ in byT+byS:
   ids=indices(C);candidates.add(C)
   for threshold in set(h[i] for i in ids):
    lo=sum((1<<i for i in ids if h[i]<=threshold),0);hi=sum((1<<i for i in ids if h[i]>=threshold),0)
    if lo.bit_count()>1:candidates.add(lo)
    if hi.bit_count()>1:candidates.add(hi)
  ordered=sorted(((C,metrics(C)) for C in candidates if C.bit_count()>1),key=lambda z:(-z[1][1],-z[1][0],z[1][4]))[:64]
  pool={C for C,_ in ordered}|twins;valid=[]
  for C in pool:
   T,s,k,after,ids=metrics(C)
   if s<=0:continue
   common=full
   for i in ids:common &=adj[i]
   boundary=adj[k]&~(C|common)&full
   delta=max((sum((q[i] for i in maskids(D&boundary)),Q()) for D in cliques),default=Q())
   if s<delta:continue
   qp=q.copy()
   for i in ids:qp[i]=Q()
   qp[k]=after;Lp=load(qp,cliques)
   assert Lp<=L and sum((h[i]*q[i] for i in range(n)),Q())==sum((h[i]*qp[i] for i in range(n)),Q())
   valid.append(((Lp,sum(qp,Q()),-s,ids,k),qp,C,s,delta,L))
  if not valid:return q,moves,'stop_candidate_pool'
  _,new,C,s,delta,L=min(valid,key=lambda x:x[0]);Lp=load(new,cliques)
  moves.append({'C':list(indices(C)),'anchor':metrics(C)[2],'saving':str(s),'boundary':str(delta),'before':str(L),'after':str(Lp),'payment_conserved':True});q=new
 return q,moves,'iteration_cap'

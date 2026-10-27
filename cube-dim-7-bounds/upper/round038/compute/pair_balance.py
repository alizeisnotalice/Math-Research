"""R038 Phase A: exact nonedge pair step, minimum-load plateau tie-broken by energy."""
from fractions import Fraction as Q
from joint_core import maskids,load

def pair_step(q,h,adj,cliques,i,j,clique_weights=None,precomputed=None):
 if i==j or (adj[i]>>j)&1:raise ValueError('requires a nonedge')
 P=h[i]*q[i]+h[j]*q[j]
 if not P:return q[i],q[j],load(q,cliques),Q()
 if clique_weights is None:clique_weights=[sum((q[k] for k in maskids(C)),Q()) for C in cliques]
 if precomputed is None:
  Ai=max((w-q[i] for C,w in zip(cliques,clique_weights) if C>>i&1),default=Q());Aj=max((w-q[j] for C,w in zip(cliques,clique_weights) if C>>j&1),default=Q());rest=max((w for C,w in zip(cliques,clique_weights) if not C & ((1<<i)|(1<<j))),default=Q())
 else:
  outside,ranked=precomputed;Ai=outside[i];Aj=outside[j];rest=next((w for w,C in ranked if not C&((1<<i)|(1<<j))),Q())
 cross=max(Q(),min(P/h[i],(P+h[j]*(Aj-Ai))/(h[i]+h[j])))
 m=max(rest,Ai+cross,Aj+(P-h[i]*cross)/h[j]);lo=max(Q(),(P-h[j]*(m-Aj))/h[i]);hi=min(P/h[i],m-Ai)
 assert lo<=hi
 xi=max(lo,min(hi,P/(h[i]+h[j])));xj=(P-h[i]*xi)/h[j]
 energy_delta=h[i]*(xi*xi-q[i]*q[i])+h[j]*(xj*xj-q[j]*q[j])
 assert xi>=0 and xj>=0 and h[i]*xi+h[j]*xj==P
 return xi,xj,m,energy_delta

def balance(q,h,adj,cliques,limit=128):
 q=list(q);n=len(q);pairs=[(i,j) for i in range(n) for j in range(i+1,n) if not (adj[i]>>j)&1];steps=[];evaluations=0
 for z in range(limit):
  weights=[sum((q[k] for k in maskids(C)),Q()) for C in cliques];L=max(weights,default=Q());E=sum((h[i]*q[i]*q[i] for i in range(n)),Q());best=None
  ranked=sorted(((weights[k],C) for k,C in enumerate(cliques)),reverse=True)
  outside=[max((weights[k]-q[i] for k,C in enumerate(cliques) if C>>i&1),default=Q()) for i in range(n)]
  # Skip pairs with no payment: their state is exactly unchanged.
  for i,j in pairs:
   if q[i]==0 and q[j]==0:continue
   xi,xj,newL,dE=pair_step(q,h,adj,cliques,i,j,weights,(outside,ranked));evaluations+=1
   if (newL,dE)<(L,Q()):
    key=(newL,dE,i,j)
    if best is None or key<best[0]:best=(key,xi,xj)
  if best is None:return q,steps,'stop_pair_coordinate',evaluations
  (newL,dE,i,j),xi,xj=best;P=h[i]*q[i]+h[j]*q[j];qi=q[i];qj=q[j];q[i]=xi;q[j]=xj
  assert h[i]*q[i]+h[j]*q[j]==P and load(q,cliques)==newL and sum((h[k]*q[k]*q[k] for k in range(n)),Q())==E+dE
  steps.append({'i':i,'j':j,'old_i':str(qi),'old_j':str(qj),'new_i':str(xi),'new_j':str(xj),'load_before':str(L),'load_after':str(newL),'energy_before':str(E),'energy_after':str(E+dE),'pair_payment':str(P)})
 return q,steps,'iteration_cap',evaluations

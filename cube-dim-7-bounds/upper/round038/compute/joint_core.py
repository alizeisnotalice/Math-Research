"""Exact clique and rational LP bound primitives; no tasks are launched here."""
from fractions import Fraction as Q
import math
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

class Censored(Exception):pass

def maximal_cliques(adj,node_limit=2**20,clique_limit=65536):
 n=len(adj); nodes=0; out=[]
 def visit(R,P,X):
  nonlocal nodes
  nodes+=1
  if nodes>=node_limit:raise Censored('graph_node_limit')
  if not P and not X:
   out.append(R)
   if len(out)>=clique_limit:raise Censored('graph_clique_limit')
   return
  union=P|X
  pivot=max((i for i in range(n) if union>>i&1),key=lambda i:(P&adj[i]).bit_count(),default=-1)
  candidates=P & ~(adj[pivot] if pivot>=0 else 0)
  while candidates:
   low=candidates & -candidates; v=low.bit_length()-1;candidates-=low
   visit(R|low,P&adj[v],X&adj[v]);P&=~low;X|=low
 visit(0,(1<<n)-1,0)
 return out,nodes

def maskids(mask):
 while mask:
  bit=mask&-mask;mask-=bit;yield bit.bit_length()-1

def load(q,cliques):return max((sum((q[i] for i in maskids(C)),Q()) for C in cliques),default=Q())

def rational_bounds(h,cliques):
 """M=min sum z_C with incidence*z >= h. Return exact feasible primal/dual bounds."""
 n=len(h);k=len(cliques);rows=[];cols=[]
 for j,C in enumerate(cliques):
  for i in maskids(C):rows.append(i);cols.append(j)
 mat=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(n,k)); hf=np.array([float(v) for v in h]); scale=max(hf); hf=hf/scale
 primal=linprog(np.ones(k),A_ub=-mat,b_ub=-hf,bounds=(0,None),method='highs')
 dual=linprog(-hf,A_ub=mat.T,b_ub=np.ones(k),bounds=(0,None),method='highs')
 if not (primal.success and dual.success):raise RuntimeError('LP non-success')
 # Treat floating solutions only as proposals: convert to exact rationals, then repair.
 den=10**12; scl=max(h); z=[Q(max(0,math.ceil(float(v)*den)),den)*scl for v in primal.x]
 for i in range(n):
  cover=sum((z[j] for j,C in enumerate(cliques) if C>>i&1),Q())
  if cover<h[i]:z[next(j for j,C in enumerate(cliques) if C>>i&1)]+=h[i]-cover
 assert all(sum((z[j] for j,C in enumerate(cliques) if C>>i&1),Q())>=h[i] for i in range(n))
 y=[Q(max(0,math.floor(float(v)*den)),den) for v in dual.x]
 worst=max((sum((y[i] for i in maskids(C)),Q()) for C in cliques),default=Q())
 if worst>1:y=[v/worst for v in y]
 assert all(sum((y[i] for i in maskids(C)),Q())<=1 for C in cliques)
 lower=sum((h[i]*y[i] for i in range(n)),Q());upper=sum(z,Q());assert 0<lower<=upper
 return lower,upper,z,y

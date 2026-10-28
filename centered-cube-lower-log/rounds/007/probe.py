from pathlib import Path
import json,numpy as np
P=Path(__file__).resolve().parent
rows=[]
for r in json.loads((P.parent/'005/refinement.json').read_text()):
 b=np.array(r['coefficients']);m=np.arange(1,len(b)+1);z=np.sqrt(b)
 stats=[]
 for l,u in [(10/9,11/9),(1.3,1.4)]:
  L=np.linspace(l,u,2001)
  A=np.sin(np.pi*L[:,None]*m)*z
  D=np.pi*m*np.cos(np.pi*L[:,None]*m)*z
  norm=np.linalg.norm(A,axis=1);U=A/norm[:,None]
  deriv=(D-U*np.sum(U*D,axis=1)[:,None])/norm[:,None]
  length=float(np.trapezoid(np.linalg.norm(deriv,axis=1),L))
  counts={}
  for sep in [.1,.25,.5]:
   picked=[]
   for i in range(len(L)):
    if not picked or np.min(np.linalg.norm(U[picked]-U[i],axis=1))>=sep:picked.append(i)
   counts[str(sep)]=len(picked)
  stats.append(dict(window=[l,u],curve_length=length,greedy_packing=counts))
 inc=[]
 L=np.linspace(10/9,11/9,1001)
 for K in [1000,10000,100000]:
  J=L+1/K
  diff=(np.sin(np.pi*J[:,None]*m)/J[:,None]-np.sin(np.pi*L[:,None]*m)/L[:,None])*np.sqrt(b/3)
  actual=float(np.max(np.sum(diff*diff,axis=1)))
  a=10/9;bound=(np.pi*np.sqrt(np.sum(m*m*b)/3)/a+np.sqrt(b.sum()/3)/a**2)**2/K**2
  inc.append(dict(K=K,increment_squared=actual,analytic_bound=float(bound)))
 rows.append(dict(M=len(b),S0=float(b.sum()),S2=float(np.sum(m*m*b)),source_sup_bound=float(np.pi*np.sum(m*np.sqrt(2*b/3))),windows=stats,increments=inc))
(P/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))

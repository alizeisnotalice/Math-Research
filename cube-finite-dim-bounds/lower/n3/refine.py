"""Continuous full-joint-union optimization of a finite two-orbit witness family."""
import json,time
import numpy as np
from scipy.optimize import minimize
from search import ROOT,decode
from witness import union
inp=json.loads((ROOT/'witness_results.json').read_text())
keys=json.loads((ROOT/'witness_family.json').read_text())
kl=np.array([k[:6] for k in keys]).reshape(-1,2,3);kh=np.array([k[6:] for k in keys]).reshape(-1,2,3)
cnt=np.prod(np.maximum(0,kh-kl+1),axis=2)

def value(z):
 if np.max(np.abs(z[:3]))>.05 or np.max(np.abs(z[3:6]))>.1 or abs(z[6])>.22:return 0.
 p,phase,w=decode(z);ph=np.array([np.zeros(3),phase]);mass=cnt@w;r=mass**(1/3)/2
 high=(kh+ph[None,:,:])*p;low=(kl+ph[None,:,:])*p
 high[cnt==0]=-np.inf;low[cnt==0]=np.inf
 lo=np.maximum(np.max(high,axis=1)-r[:,None],0);hi=np.minimum(np.min(low,axis=1)+r[:,None],p)
 b=np.stack([lo,hi],axis=2);b=b[np.all(hi>lo+1e-14,axis=1)]
 return union(b)[0]/2

def main():
 start=time.monotonic();history=[];runs=[]
 def obj(z):
  v=value(z);history.append((v,z.tolist()));return -v
 starts=[c['z'] for c in sorted(inp,key=lambda c:c['joint_union_value'],reverse=True) if np.linalg.norm(c['z'])>1e-5][:4]
 starts += [[0,0,0,0,.025,.025,.06],[0,0,0,.02,.025,-.025,-.06]]
 for z in starts:
  res=minimize(obj,z,method='Nelder-Mead',options={'maxfev':300,'xatol':1e-8,'fatol':1e-10})
  runs.append({'start':z,'z':res.x.tolist(),'joint_union_value':-res.fun,'nfev':res.nfev,'message':str(res.message)})
  print(runs[-1],flush=True)
 best=max(history);p,ph,w=decode(best[1]);out={'status':'numerical_only','evaluation_count':len(history),'elapsed_seconds':time.monotonic()-start,'witness_count':len(keys),'runs':runs,'best':{'joint_union_value':best[0],'z':best[1],'p':p.tolist(),'phase':ph.tolist(),'weights':w.tolist()}}
 (ROOT/'refine_results.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()

"""Numerical pilot; sampled g and quadrature I are NOT a certified bound."""
from pathlib import Path
import importlib.util,hashlib,json,time,math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
P=Path(__file__).resolve().parents[1]
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/cb11_frame_green_locked.py')
HASH='ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==HASH
spec=importlib.util.spec_from_file_location('frame11',SOURCE);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
records=[];started=time.monotonic()
for tau_s in ('1/10','1/3','1','3'):
 tau=F(tau_s);r={'tau':tau_s,'status':'numerical only'};seen=[]
 try:
  for i in range(9):
   for j in range(9):
    x=(F(i,8),F(j,8));g=mod.green((1,*x),tau)
    if math.isfinite(g):seen.append((g,list(map(float,x))))
  gridbest=min(seen);calls=[0]
  def objective(x):
   calls[0]+=1
   if calls[0]>256:raise RuntimeError('finite evaluation budget')
   g=mod.green((1,*map(float,x)),tau)
   if math.isfinite(g):seen.append((g,list(map(float,x))));return g
   return 1e3
  try:
   opt=minimize(objective,gridbest[1],method='SLSQP',bounds=((0,1),(0,1)),options={'maxiter':60,'ftol':1e-13})
   r['optimizer']={'success':bool(opt.success),'message':str(opt.message),'calls':calls[0]}
  except Exception as e:r['optimizer']={'exception':str(e),'calls':calls[0]}
  g,x=min(seen);r.update(g_sample=g,g_location=x,grid_g=gridbest[0])
  r['integrals']=[]
  for order in (8,16):
   z,w=np.polynomial.legendre.leggauss(order);z=(z+1)/2;w=w/2
   val=math.fsum(float(w[i]*w[j])*mod.green((1,float(z[i]),float(z[j])),tau)**1.5 for i in range(order) for j in range(order))
   r['integrals'].append({'order':order,'I_estimate':val,'optimistic_cost_estimate':3**1.5*val/g**1.5})
 except Exception as e:r.update(status='calculation error, no mathematical conclusion',exception=repr(e))
 records.append(r);print(json.dumps(r),flush=True)
 out={'status':'finite numerical pilot only; no certified g/J or C3 upper bound','source_sha256':HASH,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'records':records,'elapsed_seconds':time.monotonic()-started}
 (P/'results/round011_mixed_cost.json').write_text(json.dumps(out,indent=2)+'\n')

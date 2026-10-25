from pathlib import Path
import importlib.util,hashlib,json,time,math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc
P=Path(__file__).resolve().parents[1];SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/u11_green_locked.py');HASH='fdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==HASH
spec=importlib.util.spec_from_file_location('u11locked',SOURCE);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
started=time.monotonic();fallback=[];seen=[]
def green(y):
 x=[1.,*map(float,y)]
 try:g=mod.green_float_generic(x)
 except ValueError as e:
  fallback.append({'x':x,'reason':str(e)})
  g=float(mod.green_allnodes_pi3_coefficient(tuple(F(str(a)) for a in x)))/math.pi**3
 if not(g>0 and math.isfinite(g)):raise ValueError('nonpositive/nonfinite pilot Green')
 return g
out={'status':'numerical pilot only, no certified g/I','source_sha256':HASH,'min_samples':64,'integrals':[]}
try:
 for y in 1e-5+(1-2e-5)*qmc.Sobol(3,scramble=True,seed=20261025).random_base2(6):seen.append((green(y),list(map(float,y))))
 best=min(seen);calls=[0]
 def obj(y):
  calls[0]+=1
  if calls[0]>192:raise RuntimeError('finite min budget')
  g=green(y);seen.append((g,list(map(float,y))));return g
 try:
  opt=minimize(obj,best[1],method='SLSQP',bounds=[(1e-5,1-1e-5)]*3,options={'maxiter':40,'ftol':1e-13});out['optimizer']={'success':bool(opt.success),'message':str(opt.message),'calls':calls[0]}
 except Exception as e:out['optimizer']={'exception':repr(e),'calls':calls[0]}
 g,y=min(seen);out.update(g_sample=g,g_location=y)
 for orders in [(6,7,8),(10,11,12)]:
  nodes=[];weights=[]
  for order in orders:
   z,w=np.polynomial.legendre.leggauss(order);nodes.append((z+1)/2);weights.append(w/2)
  val=math.fsum(float(weights[0][i]*weights[1][j]*weights[2][k])*green((nodes[0][i],nodes[1][j],nodes[2][k]))**(4/3) for i in range(orders[0]) for j in range(orders[1]) for k in range(orders[2]))
  row={'orders':orders,'I_estimate':val,'optimistic_cost_estimate':4**(4/3)*val/g**(4/3)};out['integrals'].append(row);print(json.dumps(row),flush=True)
except Exception as e:out.update(status='calculation incomplete, not mathematical obstacle',exception=repr(e))
out['fallbacks']=fallback;out['elapsed_seconds']=time.monotonic()-started;out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(P/'results/round011_n4_cost.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='fallbacks'}),flush=True)

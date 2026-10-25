from pathlib import Path
import json,math,time,hashlib
import numpy as np
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent
start=time.perf_counter()
M=np.arange(24);O=2*M+1

def kernel(t,a):
 if t<.1*a*a:
  k=np.arange(-24,25)
  return np.sum(np.exp(-(1+4*a*k)**2/(4*t))-np.exp(-(1+2*a+4*a*k)**2/(4*t)))/math.sqrt(4*math.pi*t)
 return np.sum(np.exp(-O**2*math.pi**2*t/(4*a*a))*np.cos(O*math.pi/(2*a)))/a

def survival(t,a):
 if t<.004*a*a:return 1.0  # exp(-1/(4*.004)) boundary loss negligible for this pilot only
 return 4/math.pi*np.sum((-1.)**M/O*np.exp(-O**2*math.pi**2*t/(4*a*a)))
rows=[]
for a in [1.6,1.8,2.,2.2,2.4,2.6]:
 mass,em=quad(lambda t:survival(t,a)**6,0,20*a*a,epsabs=1e-11,limit=200,points=[.004*a*a,.1*a*a,1*a*a])
 g,eg=quad(lambda t:kernel(t,a)**6,0,20*a*a,epsabs=1e-12,limit=200,points=[.1*a*a,.2,.5,1.,2.])
 value=mass/(64*g)
 assert g>0 and mass>0
 rows.append(dict(a=a,mass=mass,corner_green=g,formal_cost=value,quad_error_mass=em,quad_error_corner=eg))
 if value<64:break
result=dict(status='numeric_candidate_only',rows=rows,seconds=time.perf_counter()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),registration_sha256=hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),new_C6_bound=False)
(HERE/'pilot_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

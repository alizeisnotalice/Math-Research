"""Finite exploratory batch: two 4096-point Sobol sets, 128/256 t nodes; not a certificate."""
import json,time,hashlib
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.special import roots_legendre
from scipy.stats import qmc
from flint import arb,ctx
start=time.time(); results=[]
q=6/5
g_scaled=3*np.pi/512 # integral t^6/(t^2+1)^6 dt, omit pi^-6
for seed in [601,602]:
 y=qmc.Sobol(5,scramble=True,seed=seed).random_base2(12)
 x=np.column_stack((np.ones(len(y)),y))
 for nodes in [128,256]:
  u,w=roots_legendre(nodes); u=(u+1)/2;w=w/2;t=u/(1-u)
  G=np.zeros(len(x))
  for j in range(nodes): G+=w[j]/(1-u[j])**2 * t[j]**6/np.prod(t[j]**2+x*x,axis=1)
  vals=6**q*(G/g_scaled)**q
  results.append({'seed':seed,'points':len(x),'time_nodes':nodes,'model_cost_estimate':float(vals.mean()),'max_integrand':float(vals.max())})
t=sp.symbols('t',positive=True)
g=sp.beta(sp.Rational(7,2),sp.Rational(5,2))/2
assert sp.simplify(sp.expand_func(g)-3*sp.pi/512)==0
ctx.prec=256
radial=arb(6)**(arb(21)/5)*arb.pi()**3/384
out={'status':'numerical_diagnostic_only','coordinate_g_scaled_exact':'3*pi/512','coordinate_runs':results,'radial_cost_exact':'6^(21/5)*pi^3/384','radial_cost_arb':str(radial),'radial_cost_gt_64':bool(radial>64),'elapsed_seconds':time.time()-start,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('simple_models_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

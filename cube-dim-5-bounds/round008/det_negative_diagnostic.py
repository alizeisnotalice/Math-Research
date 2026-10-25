"""New determinant statistic only; no rigorous bound."""
from pathlib import Path
from itertools import product
import numpy as np
from scipy.stats import qmc
from math import pi,gamma
import time,json,os,hashlib,resource
P=Path(__file__).resolve().parent;began=time.monotonic();beta=.75
V=np.array([(1,)+s for s in product((-1,1),repeat=4)],dtype=float)
uv=qmc.Sobol(d=32,scramble=True,seed=20261026083).random_base2(15).reshape(-1,16,2)
assert np.all((uv>0)&(uv<1));records=[]
for start in range(0,len(uv),4096):
 t=pi*uv[start:start+4096,:,0];e=-np.log(uv[start:start+4096,:,1])
 S=np.exp(np.log(np.sin(beta*t))-np.log(np.sin(t))/beta+(1-beta)/beta*(np.log(np.sin((1-beta)*t))-np.log(e)))
 T=S.sum(axis=1);M=np.einsum('bi,ij,ik->bjk',S/T[:,None],V,V);L=np.linalg.cholesky(M)
 logdetA=5*np.log(T)+2*np.log(np.diagonal(L,axis1=1,axis2=2)).sum(axis=1)
 f=16*gamma(1.75)*np.exp(-.15*logdetA);assert np.isfinite(f).all() and np.all(f>0)
 records.append({'n':len(f),'mean':float(f.mean()),'maximum':float(f.max()),'minimum':float(f.min())})
 assert time.monotonic()-began<120
out={'status':'NUMERICAL_DIAGNOSTIC_not_J_or_I_bound','formula':'J_true=16Gamma(7/4) E[detA^(-3/20)]','seed':20261026083,'n':32768,'mean':sum(x['mean'] for x in records)/8,'blocks_descriptive_not_CI':records,'limits':'No QMC/floating error bound; no confidence interval or global C5 claim. New statistic, no old g or spatial/time baseline replay.','source_sha256':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [Path(__file__),P/'det_negative_registration.md']},'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began,'peak_RSS_mac_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}}
(P/'det_negative_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='blocks_descriptive_not_CI'},indent=2))

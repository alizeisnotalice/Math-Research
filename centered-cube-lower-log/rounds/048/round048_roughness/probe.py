from pathlib import Path
import numpy as np,json,hashlib
from scipy.special import logsumexp
H=Path(__file__).resolve().parent
rng=np.random.default_rng(202610294801);rows=[];v=16.;c=.2;N=1024;rep=32
for R in [1,4,16,64]:
 u=rng.standard_normal((N,R));s=np.sqrt((1-c)*v);logbase=logsumexp(s*u,axis=1)-np.log(R);pi=np.exp(s*u-logsumexp(s*u,axis=1)[:,None]);hh=(pi*pi).sum(1);target=2*(1-c)*hh.mean()
 for dt in [.01,.001,.0001]:
  rho=np.exp(-dt/v);z=rng.standard_normal((N,rep,R));un=rho*u[:,None,:]+np.sqrt(-np.expm1(-2*dt/v))*z
  d=logsumexp(s*un,axis=2)-np.log(R)-logbase[:,None];vv=d.var(axis=1,ddof=1)/dt
  rows.append(dict(v=v,c=c,R=R,dt=dt,base_samples=N,replicates=rep,mean_participation=float(hh.mean()),predicted_variance_rate=float(target),observed_variance_rate=float(vv.mean()),Monte_Carlo_SE=float(vv.std(ddof=1)/np.sqrt(N))))
out=dict(rows=rows,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Actual finite OU increment experiment, Gaussian proxy only, not an HL bound or an outward certificate.')
(H/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

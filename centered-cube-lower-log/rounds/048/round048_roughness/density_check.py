from pathlib import Path
import numpy as np,json,hashlib
from scipy.special import logsumexp
H=Path(__file__).resolve().parent;rng=np.random.default_rng(480199);v=16.;R=3;J=5;dt=.25;t=np.arange(J)*dt;C=np.exp(-abs(t[:,None]-t[None,:])/v);D=1+R*J;Sigma=np.eye(D)
for b in range(R):Sigma[1+b*J:1+(b+1)*J,1+b*J:1+(b+1)*J]=C
x=rng.multivariate_normal(np.zeros(D),Sigma,size=100);errs=[]
for c in [0.,.2,.5,1.]:
 scores=[]
 for b in range(R):
  for j in range(J):
   mu=np.zeros(D);mu[0]=np.sqrt(c*v);mu[1+b*J:1+(b+1)*J]=np.sqrt((1-c)*v)*C[:,j]
   solved=np.linalg.solve(Sigma,mu);direct=x@solved-.5*mu@solved
   expected=-v/2+np.sqrt(c*v)*x[:,0]+np.sqrt((1-c)*v)*x[:,1+b*J+j];errs.append(float(abs(direct-expected).max()));scores.append(direct)
 q=np.logaddexp(np.log(.3),np.log(.7)+logsumexp(np.array(scores),axis=0)-np.log(R*J));lm=logsumexp(np.array(scores).reshape(R,J,len(x)),axis=0)-np.log(R);mx=lm.max(axis=0);assert np.exp(mx-q).max()<=J/.7
assert max(errs)<1e-11
out=dict(status='PASS',full_covariance_cases=4*R*J*len(x),max_log_density_error=max(errs),path_target_cap=J/.7,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Independent dense Gaussian path density check; no HL transfer.')
(H/'DENSITY_CHECK.json').write_text(json.dumps(out,indent=2)+'\n');print(out)

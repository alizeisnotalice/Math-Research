"""One new all-time matrix-mixture diagnostic, never a rigorous Green bound."""
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from math import pi,gamma
from itertools import product
import time,os,json,hashlib,resource,sys
P=Path(__file__).resolve().parent
BETA=.75;R=1.75
V=np.array([(1,)+s for s in product((-1,1),repeat=4)],dtype=float)
X=np.array([1,1,1,1,0.],dtype=float)
K=3*gamma(1.75)/(2**3.5*pi**2.5)
EXPECTED_TNEG=1/(16*gamma(1.75))
def main():
 began=time.monotonic();records=[]
 for seed in (20261026081,20261026082):
  uniforms=qmc.Sobol(d=32,scramble=True,seed=seed).random_base2(15)
  assert np.all((uniforms>0)&(uniforms<1))
  blocks=[]
  for start in range(0,2**15,4096):
   uv=uniforms[start:start+4096].reshape(-1,16,2)
   angle=pi*uv[:,:,0];exponential=-np.log(uv[:,:,1])
   logS=np.log(np.sin(BETA*angle))-np.log(np.sin(angle))/BETA+(1-BETA)/BETA*(np.log(np.sin((1-BETA)*angle))-np.log(exponential))
   S=np.exp(logS);T=S.sum(axis=1);weights=S/T[:,None]
   M=np.einsum('bi,ij,ik->bjk',weights,V,V)
   chol=np.linalg.cholesky(M)
   y=np.linalg.solve(chol,np.broadcast_to(X,(len(T),5))[:,:,None])[:,:,0]
   quadratic=np.sum(y*y,axis=1)
   logdet=2*np.log(np.diagonal(chol,axis1=1,axis2=2)).sum(axis=1)
   value=K*np.exp(-BETA*np.log(T)-logdet/2-R*np.log(quadratic))
   assert np.isfinite(value).all() and (value>0).all()
   reconstruct=chol@chol.swapaxes(-1,-2)
   residual=float(np.max(np.abs(reconstruct-M)))
   eig=np.linalg.eigvalsh(M)
   assert np.all(eig>0)
   blocks.append({'n':len(T),'G_mean':float(value.mean()),'G_max':float(value.max()),'laplace_S1_at1_mean':float(np.exp(-S[:,0]).mean()),'T_negative_beta_mean':float(np.exp(-BETA*np.log(T)).mean()),'min_M_eigenvalue':float(eig[:,0].min()),'max_M_condition':float((eig[:,-1]/eig[:,0]).max()),'cholesky_max_absolute_residual':residual})
   if time.monotonic()-began>120:raise TimeoutError('registered checkpoint; completed blocks only, no automatic restart')
  records.append({'seed':seed,'n':2**15,'blocks_descriptive_not_CI':blocks,'G_estimate':sum(x['G_mean'] for x in blocks)/8,'laplace_S1_at1_estimate':sum(x['laplace_S1_at1_mean'] for x in blocks)/8,'T_negative_beta_estimate':sum(x['T_negative_beta_mean'] for x in blocks)/8})
 out={'status':'NUMERICAL_DIAGNOSTIC_no_rigorous_bounds','model':'alpha=3/2, 16 unnormalized body directions, independent unit positive-stable variables','point':[1,1,1,1,0],'formula':'K E[T^-beta det(M)^-1/2 (x^T M^-1 x)^(-7/4)]','K_double':K,'expected_checks':{'E_exp_minus_S1':float(np.exp(-1)),'E_T_negative_beta':EXPECTED_TNEG},'records':records,'limitations':['No certified random/QMC or floating-point error.','Block ranges are not confidence intervals.','No assertion about global g, I, C5, or even a rigorous point interval.','No density interpolation or finite-time truncation is used; Kanter and matrix arithmetic are double precision.'],'source_hashes':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in (Path(__file__),P/'direct_shape_registration.md')},'runtime':{'pid':os.getpid(),'python':sys.executable,'numpy':np.__version__,'nice':os.getpriority(os.PRIO_PROCESS,0),'threads_env':{k:os.getenv(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},'seconds':time.monotonic()-began,'peak_RSS_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'platform':sys.platform}}
 (P/'direct_shape.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'records':[{k:v for k,v in r.items() if k!='blocks_descriptive_not_CI'} for r in records],'runtime':out['runtime']},indent=2))
if __name__=='__main__':main()

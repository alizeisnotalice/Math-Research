"""One all-time alpha=1 matrix diagnostic, never a rigorous Green bound."""
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from math import pi,gamma
from itertools import product
import time,os,json,hashlib,resource,sys

P=Path(__file__).resolve().parent
BETA=.5;ALPHA=1.;R=2.;SEED=202610260141;N=32768
V=np.array([(1,)+s for s in product((-1,1),repeat=4)],dtype=float)
X=np.array([1,7/8,3/4,5/8,1/2],dtype=float)
K=BETA*gamma(R)/(2**ALPHA*pi**2.5)
EXPECTED_TNEG=1/(16*gamma(1+BETA))

def main():
    began=time.monotonic();blocks=[];failure=None
    checkpoint=Path('/tmp/astra_upper014_diagnostic_checkpoint.json')
    uniforms=qmc.Sobol(d=32,scramble=True,seed=SEED).random_base2(15)
    assert np.all((uniforms>0)&(uniforms<1))
    try:
        for start in range(0,N,4096):
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
            eig=np.linalg.eigvalsh(M);assert np.all(eig>0)
            blocks.append({'n':len(T),'G_mean':float(value.mean()),'G_max':float(value.max()),
                           'laplace_S1_at1_mean':float(np.exp(-S[:,0]).mean()),
                           'T_negative_beta_mean':float(np.exp(-BETA*np.log(T)).mean()),
                           'min_M_eigenvalue':float(eig[:,0].min()),
                           'max_M_condition':float((eig[:,-1]/eig[:,0]).max()),
                           'cholesky_max_absolute_residual':float(np.max(np.abs(reconstruct-M)))})
            checkpoint.write_text(json.dumps({'seed':SEED,'completed_n':sum(b['n'] for b in blocks),'blocks':blocks},indent=2)+'\n')
            if time.monotonic()-began>120:raise TimeoutError('registered checkpoint after complete block')
    except Exception as e:
        failure=repr(e)
    completed=sum(b['n'] for b in blocks)
    summary={k:sum(b[k]*b['n'] for b in blocks)/completed for k in ['G_mean','laplace_S1_at1_mean','T_negative_beta_mean']} if completed else {}
    budget=json.loads((P/'budget.json').read_text())['cases'][0]
    # Only a numerical diagnostic cost, substituting a finite one-point sample mean.
    Iup=float(Fraction(budget['I_upper_expression']['rational_upper']))
    if completed==N and not failure:
        summary['cost_proxy_using_point_sample_mean_NOT_global_g']=5**1.25*Iup/summary['G_mean']**1.25
    out={'status':'NUMERICAL_DIAGNOSTIC_no_rigorous_bounds' if not failure else 'STOPPED_WITH_RECORDED_FAILURE',
         'model':'alpha=1, 16 unnormalized body directions, independent unit positive-stable variables',
         'point':[str(v) for v in [Fraction(1),Fraction(7,8),Fraction(3,4),Fraction(5,8),Fraction(1,2)]],
         'known_infinite_11110_point_was_not_run':True,'seed':SEED,'requested_n':N,'completed_n':completed,
         'formula':'K E[T^-beta det(M)^-1/2 (x^T M^-1 x)^(-2)]','K_double':K,
         'expected_checks':{'E_exp_minus_S1':float(np.exp(-1)),'E_T_negative_beta':EXPECTED_TNEG},
         'summary':summary,'blocks_descriptive_not_CI':blocks,'failure':failure,
         'limitations':['No certified QMC or floating-point error.','Block ranges are not confidence intervals.',
                        'A single-point sample mean is not a global g and not a rigorous bound in either direction.',
                        'No time truncation, density interpolation, sample rejection, matrix repair or condition clipping was used.'],
         'source_hashes':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [Path(__file__),P/'registration.md',P/'budget.json']},
         'runtime':{'pid':os.getpid(),'python':sys.executable,'numpy':np.__version__,'nice':os.getpriority(os.PRIO_PROCESS,0),
                    'threads_env':{k:os.getenv(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},
                    'seconds':time.monotonic()-began,'peak_RSS_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'platform':sys.platform}}
    (P/'diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','point','seed','completed_n','summary','expected_checks','failure','runtime']},indent=2))

from fractions import Fraction
if __name__=='__main__':main()

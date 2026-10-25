"""Two-orbit common-radius periodic discovery; samples are NOT certificates."""
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import differential_evolution, minimize
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parents[1]/'baseline/original_certificates.json'
cert=json.loads(BASE.read_text())[0]
u=2**(1/3);v=3**(1/3)
a=np.array([1+u/2,(2+u*u*v)/4,(1+u)/2])
pbase=a*np.array([2,1,1])

def decode(z):
 p=pbase*np.exp(z[:3]); phase=np.array([.5,0,0])+z[3:6]
 return p,phase,np.array([1+z[6],1-z[6]])

def eval_receivers(z,t,details=False):
 p,phase,w=decode(np.asarray(z));q=np.mod(t[:,None,:]-np.array([np.zeros(3),phase])[None,:,:],1)
 # All point-entry radii <=2.8. k=-3..3 suffices when min(P)>.94.
 rr=np.abs((q[:,:,:,None]-np.arange(-3,4))*p[None,None,:,None]).reshape(len(t),-1)
 best=np.zeros(len(t));br=np.zeros(len(t))
 # chunk keeps large validation memory bounded.
 for j in range(rr.shape[1]):
  r=rr[:,j];N=np.zeros(len(t))
  for b in range(2):
   cnt=np.ones(len(t))
   for d in range(3):
    h=r/p[d]
    cnt*=np.floor(q[:,b,d]+h+1e-11)+np.floor(h-q[:,b,d]+1e-11)+1
   N+=w[b]*cnt
  avg=np.divide(N,8*r**3,out=np.zeros_like(r),where=r>1e-10)
  avg[r>2.8]=0
  improve=avg>best;best[improve]=avg[improve];br[improve]=r[improve]
 bad=best>1+1e-10
 value=float(np.mean(bad)*np.prod(p)/2)
 return (value,best,br) if details else value

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--budget',type=int,default=1120);args=ap.parse_args()
 start=time.monotonic();t=qmc.Sobol(3,scramble=True,seed=202610253).random_base2(12)
 history=[]
 def obj(z):
  val=eval_receivers(z,t);history.append((val,np.array(z).tolist()))
  if len(history)%100==0: print('eval',len(history),'best',max(x[0] for x in history),'elapsed',time.monotonic()-start,flush=True)
  return -val
 bounds=[(-.035,.035)]*3+[(-.08,.08),(-.08,.08),(-.08,.08),(-.18,.18)]
 # Baseline and genuinely different alternating/staggered starts precede population search.
 for z in [np.zeros(7),[0,0,0,0,.06,.06,0],[0,0,0,0,.04,-.04,.12],[0,0,0,.06,0,0,-.12]]:obj(z)
 res=differential_evolution(obj,bounds,popsize=10,maxiter=max(1,args.budget//70-1),tol=0,polish=False,seed=202610253,x0=np.zeros(7),updating='immediate',workers=1)
 # Validate best separated parameter candidates on common independent Sobol points.
 candidates=[]
 for val,z in sorted(history,reverse=True):
  if all(np.linalg.norm(np.array(z)-np.array(c['z']))>.015 for c in candidates):candidates.append({'sample_value':val,'z':z})
  if len(candidates)==8:break
 candidates.append({'sample_value':eval_receivers(np.zeros(7),t),'z':[0.]*7,'baseline':True})
 tv=qmc.Sobol(3,scramble=True,seed=202610254).random_base2(17)
 for c in candidates:
  vals=[];best=[];radii=[]
  for chunk in np.array_split(tv,32):
   vv,bb,rr=eval_receivers(c['z'],chunk,True);best.extend(bb.tolist());radii.extend(rr.tolist())
  c['validation_value']=float(np.mean(np.array(best)>1+1e-10)*np.prod(decode(c['z'])[0])/2)
  c['p'],c['phase'],c['weights']=[x.tolist() for x in decode(c['z'])]
  c['validation_bad_count']=int(np.sum(np.array(best)>1+1e-10));c['validation_size']=len(tv)
  c['maximum_selected_radius']=float(max(radii))
  print('validation',c,flush=True)
 out={'status':'numerical_only','mechanism':'two_orbit_periodic_weighted_atoms','parameterization':'P=pbase exp(z[0:3]); phase=(1/2,0,0)+z[3:6]; masses=(1+z6,1-z6)','discovery_seed':202610253,'validation_seed':202610254,'evaluation_count':len(history),'discovery_size':len(t),'radius_cap':2.8,'baseline_formula':cert['parameter_formula'],'baseline_cert_sha256':hashlib.sha256(BASE.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'candidates':candidates,'bounds':bounds}
 (ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()

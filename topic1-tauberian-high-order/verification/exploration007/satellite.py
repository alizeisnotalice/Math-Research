import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import json,itertools,time,math,multiprocessing as mp,sys
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from numpy.polynomial import polynomial as P
base=Path(__file__).resolve().parent;out=base/'satellite';out.mkdir(exist_ok=True)
POWER=int(os.environ.get('SOBOL_POWER','6'))

def direct(x,r,boxes,n):
 volume=0.
 for c,d in boxes:
  v=np.maximum(0.,np.minimum(x+r,c+d)-np.maximum(x-r,c-d));volume+=float(np.prod(v))
 return volume/(2*r)**n

def min_defect(x,boxes,n,rlo,rhi=1.):
 events=[rlo,rhi]
 for ctr,w in boxes:
  for z in np.abs(x-ctr):
   for v in (abs(z-w),z+w):
    if rlo<v<rhi:events.append(float(v))
 ev=np.unique(events); best=(float('inf'),None)
 for lo,hi in zip(ev[:-1],ev[1:]):
  mid=(lo+hi)/2;poly=np.zeros(n+1)
  for ctr,w in boxes:
   leng=np.maximum(0.,np.minimum(x+mid,ctr+w)-np.maximum(x-mid,ctr-w))
   if np.any(leng<=0):continue
   slope=(x+mid<ctr+w).astype(float)+(x-mid>ctr-w).astype(float);coeff=np.array([1.])
   for l,s in zip(leng,slope):coeff=P.polymul(coeff,[l-s*mid,s])
   poly[:len(coeff)]+=coeff
  if not np.any(poly):continue
  deriv=(np.arange(n+1)-n)*poly
  roots=P.polyroots(deriv) if np.any(np.abs(deriv[1:])>1e-25) else []
  cand=[lo,hi]+[float(z.real) for z in roots if abs(z.imag)<1e-7 and lo<z.real<hi]
  for rad in cand:
   v=direct(x,rad,boxes,n)
   if v<best[0]:best=(v,rad)
 return best

def eval_one(n,lam,beta,c,u,seed):
 delta=u/(1+2*u)-c*u**3;eps=(lam/2**n)**(1/n)*u**(1+2/n);z=beta*u**(2/n)
 gap=z-delta-eps; rhosat=(eps/u+gap)/2
 assert delta>0 and gap>eps/u and rhosat < .5-2*eps and (2*delta<1)
 centers=[]
 for ij in itertools.product((-1,0,1),repeat=n):
  ctr=np.array(ij,float);centers.append((ctr,delta));b=ctr.copy();b[0]+=z;centers.append((b,eps))
 # For x in the central main pore and r<=1, shifts outside {-1,0,1} cannot intersect.
 w=min((4+lam+abs(c))*u*u,delta*.49); cuts=[-delta,-delta+w,delta-w,delta];
 maincap=0.; strata=[];accepted=[];satpoint=np.zeros(n);satpoint[0]=z
 sv=direct(satpoint,rhosat,centers,n);assert sv<u**n
 for idx,zone in enumerate(itertools.product(range(3),repeat=n)):
  low=np.array([cuts[i] for i in zone]);high=np.array([cuts[i+1] for i in zone]);vol=float(np.prod(high-low))
  pts=low+(high-low)*qmc.Sobol(n,scramble=True,seed=seed+idx).random_base2(POWER);count=0
  for x in pts:
   v,r=min_defect(x,centers,n,delta/u)
   if v<u**n*(1-1e-11):
    count+=1
    if len(accepted)<16:accepted.append({'x':x.tolist(),'radius':r,'defect':v,'margin':u**n-v})
  maincap+=vol*count/len(pts);strata.append({'zone':zone,'volume':vol,'count':count,'samples':len(pts)})
 p=(2*eps)**n;q=(2*delta)**n+p;gain=(maincap+p)/(1-q);q0=(2*u/(1+2*u))**n;baseline=q0/(1-q0)
 return {'n':n,'lambda':lam,'beta':beta,'c':c,'u':u,'delta':delta,'epsilon':eps,'z':z,'gap':gap,'satellite_radius':rhosat,'satellite_center_defect':sv,'satellite_center_margin':u**n-sv,'strata':strata,'accepted_examples':accepted,'maincap_estimate':maincap,'main_missed_scaled':((2*delta)**n-maincap)/u**(n+1),'gamma_estimate':(gain-baseline)/u**(n+2),'note':'floating-point stratified Sobol estimate; numerical radii roots; no QMC error bound'}

def worker(job):
 idx,n,lam,beta=job;start=time.time();cc=2*lam/(n*2**n);done=[]
 global POWER
 POWER = 8 if n==2 else 5
 for ui,u in enumerate((1/40,1/80,1/160)):
  for ci,c in enumerate((cc-1/8,cc,cc+1/8)):

   try: result=eval_one(n,lam,beta,c,u,330000+idx*1000+ui*100+ci*30)
   except AssertionError: result={'n':n,'lambda':lam,'beta':beta,'c':c,'u':u,'status':'geometry_or_satellite_witness_preconditions_failed'}
   done.append(result)
   (out/f'branch_{idx:02d}.json').write_text(json.dumps({'id':idx,'n':n,'lambda':lam,'beta':beta,'pid':os.getpid(),'completed':done,'elapsed':time.time()-start},indent=2))
 print(json.dumps({'id':idx,'pid':os.getpid(),'seconds':time.time()-start}),flush=True);return idx

if __name__=='__main__':
 tasks=[]
 for lam in (1/16,1/4,1):
  for extra in (1/16,1/2,2):tasks.append((len(tasks),2,lam,1+math.sqrt(lam)/2+extra))
 for lam in (1/16,1/4,1):
  for extra in (1/4,2):tasks.append((len(tasks),3,lam,(lam/8)**(1/3)+extra))
 print(json.dumps({'master':os.getpid(),'tasks':tasks,'sobol_per_stratum':2**POWER}),flush=True)
 if len(sys.argv)>1 and sys.argv[1]=='pilot':
  t=time.time();r=eval_one(2,1/4,1.3125,1/16,1/40,5481);print(json.dumps({'pilot_seconds':time.time()-t,'gap':r['gap'],'satmargin':r['satellite_center_margin'],'strata':len(r['strata']),'gamma':r['gamma_estimate']}),flush=True)
 else:
  with mp.Pool(15) as pool:
   for _ in pool.imap_unordered(worker,tasks):pass

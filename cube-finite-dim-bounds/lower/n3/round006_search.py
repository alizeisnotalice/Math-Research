"""Primary plus two independent small satellites, complete finite joint unions."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib,json,time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
from witness import union
from round006_union import volume as mask_union
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'round005c_certificate_input.json').read_text());P0=np.array(list(map(float,SEED['periods'])));PH1=np.array(list(map(float,SEED['phases'][1])));EPS1=float(F(SEED['weights'][1]));TARGET=json.loads((ROOT/'round006_target.json').read_text())['ratio_interval'][1]
OLD=sorted(set(tuple(k[:6]+[0,0,0]+k[6:]+[-1,-1,-1]) for k in SEED['witness_keys']))

def decode(z):
 z=np.asarray(z);return P0*np.exp(z[:3]),np.array([[0.,0.,0.],z[3:6],z[6:9]]),np.array([1.,np.exp(z[9]),np.exp(z[10])])

def boxes(z,keys):
 p,ph,w=decode(z);keys=np.asarray(keys,dtype=int);lo=keys[:,:9].reshape(-1,3,3);hi=keys[:,9:].reshape(-1,3,3);cnt=np.prod(np.maximum(0,hi-lo+1),axis=2);mass=cnt@w;r=mass**(1/3)/2
 low=(lo+ph[None,:,:])*p;high=(hi+ph[None,:,:])*p;low[cnt==0]=np.inf;high[cnt==0]=-np.inf;lower=np.maximum(np.max(high,axis=1)-r[:,None],0);upper=np.minimum(np.min(low,axis=1)+r[:,None],p);good=np.all(upper>lower+1e-13,axis=1)
 return np.stack([lower[good],upper[good]],axis=2),keys[good].tolist()

def value(z,keys):
 z=np.asarray(z)
 if max(abs(z[:3]))>.04 or min(z[3:9])<0 or max(z[3:9])>1 or not np.log(2**-16)<=z[9]<=np.log(2**-3) or not np.log(2**-18)<=z[10]<=np.log(2**-3):return 0.
 b,kk=boxes(z,keys)
 if not len(b):return 0.
 cells=int(np.prod([len(np.unique(b[:,j,:])) for j in range(3)]))
 if cells>20_000_000:return mask_union(b)/sum(decode(z)[2])
 return float(union(b)[0]/sum(decode(z)[2]))

def clusters():
 points=list(product(range(-1,2),repeat=3));empty=(0,0,0,-1,-1,-1);keys=set()
 def add(bs):
  keys.add(tuple(c for b in bs for c in b[:3])+tuple(c for b in bs for c in b[3:]))
 for b in range(3):
  for p in points:
   bs=[empty]*3;bs[b]=p+p;add(bs)
 for pair in [(0,1),(0,2),(1,2)]:
  for p,q in product(points,repeat=2):
   bs=[empty]*3;bs[pair[0]]=p+p;bs[pair[1]]=q+q;add(bs)
 for p,q,s in product(points,repeat=3):add([p+p,q+q,s+s])
 return sorted(keys)

CLUSTERS=clusters()

def capture(z,t):
 p,ph,w=decode(z);q=np.mod(t[:,None,:]-ph[None,:,:],1);rr=np.abs((q[:,:,:,None]-np.arange(-4,5))*p[None,None,:,None]).reshape(len(t),-1);keys=set()
 for j in range(rr.shape[1]):
  r=rr[:,j];lo=np.ceil(t[:,None,:]-ph[None,:,:]-r[:,None,None]/p-1e-11).astype(int);hi=np.floor(t[:,None,:]-ph[None,:,:]+r[:,None,None]/p+1e-11).astype(int);cnt=np.prod(np.maximum(0,hi-lo+1),axis=2);mass=cnt@w;good=(mass>8*r**3+1e-10)&(r<=3.2)
  for l,h in zip(lo[good],hi[good]):
   for b in range(3):
    if np.any(h[b]<l[b]):l[b]=[0,0,0];h[b]=[-1,-1,-1]
   keys.add(tuple(map(int,np.concatenate([l.ravel(),h.ravel()]))))
 return sorted(keys)

def calibration():
 keys=[]
 for k in SEED['witness_keys']:keys.append(tuple(k[:6]+k[3:6]+k[6:]+k[9:12]))
 z=np.r_[np.zeros(3),PH1,PH1,np.log(EPS1/2),np.log(EPS1/2)];v=value(z,keys)
 assert abs(v-float(F(TARGET)))<2e-12
 # Exact selected masses and extrema are unchanged when two coincident blocks split weight.
 ep=F(SEED['weights'][1]);assert ep/2+ep/2==ep
 for k in SEED['witness_keys']:
  c0=np.prod(np.maximum(0,np.array(k[6:9])-k[:3]+1));c1=np.prod(np.maximum(0,np.array(k[9:12])-k[3:6]+1));assert F(int(c0))+ep*int(c1)==F(int(c0))+(ep/2)*int(c1)+(ep/2)*int(c1)
 return {'status':'duplicate-satellite interface regression PASS','joint_union':v,'seed_target_display':float(F(TARGET)),'exact_mass_identity_checked_for_keys':len(keys),'phases_coincident':True,'mass_split':['1',str(ep/2),str(ep/2)]}

class Significant(Exception):pass

def save(out,z,keys):
 p,ph,w=decode(z);b,active=boxes(z,keys);saved=sorted(set(map(tuple,active))|set(OLD));assert set(OLD)<=set(saved)
 out['final']={'z':np.asarray(z).tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'joint_union':value(z,keys),'active_keys':len(active),'saved_keys':len(saved),'old170_preserved':set(OLD)<=set(saved)}
 (ROOT/'round006_results.json').write_text(json.dumps(out,indent=2)+'\n');inp={'periods':[format(x,'.15f') for x in p],'phases':[['0','0','0']]+[[format(x,'.15f') for x in q] for q in ph[1:]],'weights':['1',format(w[1],'.15f'),format(w[2],'.15f')],'keys':[list(k) for k in saved]};(ROOT/'round006_certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n');print('FINAL',out['final'],flush=True)

def main():
 start=time.monotonic();cal=calibration();(ROOT/'round006_calibration.json').write_text(json.dumps(cal,indent=2)+'\n');base=sorted(set(OLD)|set(CLUSTERS));pool=set(base);t=qmc.Sobol(3,scramble=True,seed=20261025602).random_base2(12)
 struct=[1-PH1,np.roll(PH1,1),np.roll(PH1,2),PH1[::-1],np.mod(PH1+[.5,.25,.125],1),np.mod(PH1+[.25,.5,.375],1),PH1/2,np.mod(2*PH1,1)];phases=qmc.Sobol(3,scramble=True,seed=20261025601).random_base2(5).tolist()+[p.tolist() for p in struct];cases=[]
 for i,ph in enumerate(phases):
  for factor in [.25,1,4]:
   z=np.r_[np.zeros(3),PH1,ph,np.log(EPS1),np.log(EPS1*factor)];rec=capture(z,t);kk=sorted(set(base)|set(rec));pool.update(rec);val=value(z,kk);cases.append({'phase_id':i,'theta2':ph,'eps2_factor':factor,'z':z.tolist(),'joint_union':val,'captured_keys':len(rec)});(ROOT/'round006_scan_receipt.json').write_text(json.dumps(cases,indent=2)+'\n')
 print('scan',len(cases),'raw clusters',len(CLUSTERS),'pool',len(pool),'best',max(c['joint_union'] for c in cases),flush=True)
 out={'target_exact_upper':TARGET,'seed_input_sha256':hashlib.sha256((ROOT/'round005c_certificate_input.json').read_bytes()).hexdigest(),'initial_cases':cases,'phase_seed':20261025601,'capture_seed':20261025602,'raw_cluster_keys':len(CLUSTERS),'pool_keys':len(pool),'calibration':cal,'runs':[],'optimization_evaluations':0}
 bestcase=max(cases,key=lambda c:c['joint_union']);pool=sorted(pool)
 if bestcase['joint_union']>float(F(TARGET))+1e-4:out.update({'stop':'significant structural scan excess; certify before refinement','elapsed_seconds':time.monotonic()-start});save(out,bestcase['z'],pool);return
 starts=[]
 for c in sorted(cases,key=lambda c:c['joint_union'],reverse=True):
  if all(np.linalg.norm(np.array(c['theta2'])-np.array(s['theta2']))>.15 for s in starts):starts.append(c)
  if len(starts)==3:break
 count=0;runs=[];best=[-1.,None]
 for stage,(c,budget) in enumerate(zip(starts,[1200,600,600])):
  zstart=np.array(c['z']);nfev=0;current=[-1.,None]
  def obj(z):
   nonlocal count,nfev
   val=value(z,pool);count+=1;nfev+=1
   if val>current[0]:current[:]=[val,np.asarray(z).copy()]
   if val>float(F(TARGET))+1e-4:raise Significant()
   return -val
  simplex=np.tile(zstart,(12,1));steps=[.0003,.0001,.0001,.01,.01,.01,.02,.02,.02,.2,.3]
  for j in range(11):simplex[j+1,j]+=steps[j]
  triggered=False
  try:res=minimize(obj,zstart,method='Nelder-Mead',options={'maxfev':budget,'initial_simplex':simplex,'xatol':1e-9,'fatol':1e-11});stop=str(res.message)
  except Significant:triggered=True;stop='significant gain>1e-4; stop and certify'
  z=current[1];before=value(z,pool);rec=capture(z,qmc.Sobol(3,scramble=True,seed=20261025603).random_base2(15));newpool=sorted(set(pool)|set(rec));after=value(z,newpool);row={'stage':stage,'start':c['z'],'z':z.tolist(),'nfev':nfev,'stop':stop,'before_recapture':before,'after_recapture':after,'extra_keys':len(set(rec)-set(pool))};runs.append(row);print(row,flush=True);pool=newpool
  if after>best[0]:best[:]=[after,z.copy()]
  if triggered or after>float(F(TARGET))+1e-4 or after-before>1e-8:break
 z=best[1];rec=capture(z,qmc.Sobol(3,scramble=True,seed=20261025604).random_base2(16));oldpool=pool;pool=sorted(set(pool)|set(rec));out.update({'runs':runs,'optimization_evaluations':count,'final_recapture_seed':20261025604,'final_extra_keys':len(set(rec)-set(oldpool)),'elapsed_seconds':time.monotonic()-start,'stop':'finite batch stop; exact certificate follows'});save(out,z,pool)
if __name__=='__main__':main()

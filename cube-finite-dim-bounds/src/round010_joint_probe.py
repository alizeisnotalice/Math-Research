"""Two spatial correction bases, finite numerical LP only; no global certificate."""
from pathlib import Path
import sys,hashlib,json,math,itertools,time,warnings,os
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.integrate import quad,IntegrationWarning
from scipy.optimize import linprog
from scipy.stats import qmc
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'upper/continuous_bridge'))
from cb9_frame_green import green
SHA=hashlib.sha256((P/'upper/continuous_bridge/cb9_frame_green.py').read_bytes()).hexdigest()
assert SHA=='d23501eeaf417a736730a6e61398e27b9364c231d35ec6ee8bf573e357058bf5'
V=np.array([(1,)+v for v in itertools.product((-1,1),repeat=2)],float)
CENTERS=[np.eye(3)[i]*sign*.75 for i in range(3) for sign in (-1,1)]
GE=(4*math.log(2)-1)/(8*math.pi**3);C=GE/3;D=8*(GE-C)
VOL=np.array([32*math.pi/105,6*(.5**3)*32*math.pi/105])
def unit_b(x):return max(1-float(x@x),0)**2
def unit_lb(x):
 total=0.;n=3
 for v in V:
  dot=float(x@v);a=-dot/n;A=1-float(x@x)+dot*dot/n
  if A<=0:continue
  h=math.sqrt(A/n)
  if abs(a)>2*h:
   z=(h/a)**2;sm=0.;power=z
   for j in range(1,41):
    term=power/((2*j+1)*(2*j+3));sm+=term;power*=z
   val=16*n*n*h**3/math.pi*sm
  else:
   term=0. if a==0 or abs(abs(a)-h)<1e-13 else a*(h*h-a*a)*math.log(abs((a+h)/(a-h)))
   val=4*n*n/math.pi*(2*a*a*h-4*h**3/3+term)
  total+=val
 return total
def bases(x):
 x=np.array(x,float)
 return np.array([unit_b(x),sum(unit_b((x-a)/.5) for a in CENTERS)])
def lbs(x):
 x=np.array(x,float)
 return np.array([unit_lb(x),sum(unit_lb((x-a)/.5)/.5 for a in CENTERS)])
def source(job):
 x,tol=job;x=np.array(x,float);g=float(green(x));row={'x':x.tolist(),'tolerance':tol,'G':g if math.isfinite(g) else 'infinity','Lb':lbs(x).tolist()}
 if g<=C:return dict(row,status='outside_Omega_formula_not_used')
 ts=time.time();messages=[];value=error=0.
 with warnings.catch_warnings(record=True) as ww:
  warnings.simplefilter('always',IntegrationWarning)
  for v in V:
   def f(z):
    if z==0:return 0.
    s=z/(1-z)
    return (max(C-float(green(x+s*v)),0)+max(C-float(green(x-s*v)),0))/(z*z)
   for a,b in zip([0,.1,1/3,.5,.75],[.1,1/3,.5,.75,1]):
    val,err=quad(f,a,b,epsabs=tol,epsrel=tol,limit=100);value+=val;error+=err
  messages=sorted(set(str(w.message) for w in ww))
 row.update(mu=value/math.pi,error_estimate_only=error/math.pi,warnings=messages,seconds=time.time()-ts,status='numerical_only' if not messages and math.isfinite(value) else 'quadrature_unresolved')
 return row

def run_jobs(jobs,cache):
 results={};pending=[]
 for x,t in jobs:
  key=json.dumps([list(x),t])
  if key in cache:results[key]=cache[key]
  else:pending.append((key,(x,t)))
 with ProcessPoolExecutor(max_workers=min(12,os.cpu_count() or 1)) as ex:
  fut={ex.submit(source,job):key for key,job in pending}
  for ff in as_completed(fut):
   key=fut[ff];cache[key]=ff.result();results[key]=cache[key]
   (P/'results/round010_source_cache.json').write_text(json.dumps({'module_sha256':SHA,'rows':cache},indent=2)+'\n')
   print('source',len(results),'/',len(jobs),cache[key]['status'],flush=True)
 return list(results.values())
def solve(rows,dom):
 usable=[r for r in rows if r['status']=='numerical_only']
 AA=[r['Lb'] for r in usable]+[r['b'] for r in dom];BB=[r['mu'] for r in usable]+[r['rhs'] for r in dom]
 res=linprog(-VOL,A_ub=np.array(AA),b_ub=np.array(BB),bounds=[(0,None),(0,None)],method='highs')
 if not res.success:return {'success':False,'message':res.message},None
 w=res.x;slack=np.array(BB)-np.array(AA)@w
 return {'success':True,'weights':w.tolist(),'mass_reduction':float(VOL@w),'optimistic_cost_estimate':8.171963453262865-float(VOL@w)/D,'max_finite_constraint_violation':float(max(0,-min(slack))),'source_rows':len(usable),'domain_rows':len(dom),'message':res.message},[(r,float(r['mu']-np.dot(r['Lb'],w))) for r in usable]
if __name__=='__main__':
 start=time.time();cache={};cp=P/'results/round010_source_cache.json'
 if cp.exists():
  old=json.loads(cp.read_text());assert old['module_sha256']==SHA;cache=old['rows']
 shapes=[(0,0),(.5,0),(1,0),(.5,.5),(1,.5),(1,1),(.95,.9)]
 points=[(r,r*y,r*z) for r in (.4,.8,1,1.4,2,3,6) for y,z in shapes]
 rows=run_jobs([(x,1e-6) for x in points],cache)
 dom=[]
 for x in itertools.combinations_with_replacement(np.linspace(0,1.25,11),3):
  x=np.array(x[::-1]);b=bases(x)
  if not np.any(b):continue
  g=float(green(x))
  if not math.isfinite(g):continue
  rhs=g-(GE if max(x)<=1 else C)
  if rhs < -1e-10:raise RuntimeError('observed G below assumed normalizer; candidate baseline invalid')
  dom.append({'x':x.tolist(),'b':b.tolist(),'rhs':max(rhs,0),'unclipped_rhs':rhs})
 one,sl=solve(rows,dom);out={'status':'numerical finite LP only; no full g/J, source or dominance certification','module_sha256':SHA,'fixed_c':C,'normalizer_unproved_global_g':GE,'basis_volumes':VOL.tolist(),'first_LP':one,'domain_rows':dom,'source_rows':rows}
 (P/'results/round010_joint_probe.json').write_text(json.dumps(out,indent=2)+'\n')
 if sl is not None:
  active=[r['x'] for r,s in sorted(sl,key=lambda t:t[1])[:8]];fine=run_jobs([(x,1e-8) for x in active],cache)
  fineby={tuple(r['x']):r for r in fine};rows2=[fineby.get(tuple(r['x']),r) for r in rows];two,_=solve(rows2,dom);out['second_LP']=two
  if two['success']:
   vv=qmc.Sobol(3,scramble=False).random_base2(5)*1.5;validation=run_jobs([(tuple(sorted(x,reverse=True)),1e-7) for x in vv if np.any(x)],cache);w=np.array(two['weights'])
   for r in validation:
    r['source_slack']=float(r['mu']-np.dot(r['Lb'],w)) if r['status']=='numerical_only' else None
    x=np.array(r['x']);g=float(green(x));r['dominance_slack']=float(g-(GE if max(x)<=1 else C)-bases(x)@w) if math.isfinite(g) else 'infinity'
   out['validation']=validation;out['validation_has_negative_source_slack']=any(r['source_slack'] is not None and r['source_slack'] < -1e-7 for r in validation)
 out['elapsed_seconds']=time.time()-start;out['scope']='two fixed bases, at most2 LPs/8 refined rows/32 Sobol points; excluded quadrature rows are unresolved, not passed constraints'
 (P/'results/round010_joint_probe.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('domain_rows','source_rows','validation')},indent=2))

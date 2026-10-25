"""Scaled bounded Powell blocks through contact switches, full finite box union."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
import round006_search as eng
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'round006_certificate_input.json').read_text());eng.P0=np.array(list(map(float,SEED['periods'])));PH=np.array([list(map(float,q)) for q in SEED['phases'][1:]]).ravel();LOGW=np.log(np.array(list(map(float,SEED['weights'][1:]))));ORIGIN=np.r_[np.zeros(3),PH,LOGW];SCALE=np.array([.0005]*3+[.01]*6+[.2]*2);OLD=sorted(set(map(tuple,SEED['keys'])));TARGET=json.loads((ROOT/'round007_target.json').read_text())['ratio_interval'][1]
BOUNDS=[(-3.,3.)]*3+[(max(-4.,-PH[j]/.01),min(4.,(1-PH[j])/.01)) for j in range(6)]+[(-5.,5.)]*2

def zdecode(x):return ORIGIN+np.asarray(x)*SCALE

def value(x,keys):return float(eng.value(zdecode(x),keys))
class Stop(Exception):pass

def main():
 start=time.monotonic();pool=sorted(set(OLD)|set(eng.CLUSTERS));x=np.zeros(11);seedval=value(x,pool);best=[seedval,x.copy()];trace=[];runs=[];total=0;trigger=False;(ROOT/'round007_best_trace.json').write_text('[]\n')
 blocks=[('period_mass',[0,1,2,9,10],300),('phase',list(range(3,9)),400),('joint',list(range(11)),1000)]
 for stage,(name,indices,budget) in enumerate(blocks):
  xstart=best[1].copy();base=xstart.copy();nfev=0
  def obj(y):
   nonlocal total,nfev,trigger
   trial=base.copy();trial[indices]=y;v=value(trial,pool);total+=1;nfev+=1
   if not np.isfinite(v):trigger=True;raise Stop('nonfinite objective anomaly')
   if v>best[0]:
    best[:]=[v,trial.copy()];trace.append({'nfev':total,'stage':name,'joint_union':v,'x':trial.tolist()});(ROOT/'round007_best_trace.json').write_text(json.dumps(trace,separators=(',',':'))+'\n')
   (ROOT/'round007_progress.json').write_text(json.dumps({'stage':name,'stage_evaluations':nfev,'total_evaluations':total,'best_joint_union':best[0],'best_x':best[1].tolist(),'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
   if v>float(F(TARGET))+1e-4:trigger=True;raise Stop('significant gain>1e-4; prioritize exact certification')
   return -v
  try:
   result=minimize(obj,xstart[indices],method='Powell',bounds=[BOUNDS[j] for j in indices],options={'maxfev':budget,'xtol':1e-7,'ftol':1e-10});message=str(result.message)
  except Stop as e:message=str(e)
  x=best[1].copy();before=value(x,pool);size=16 if name=='joint' or trigger else 15;cap_seed=20261025701+stage;prior_pool=set(pool);extra=eng.capture(zdecode(x),qmc.Sobol(3,scramble=True,seed=cap_seed).random_base2(size));newpool=sorted(set(pool)|set(extra));after=value(x,newpool);best[0]=after;pool=newpool
  if after>before:
   trace.append({'nfev':total,'stage':name+'_recapture','event':'post-block recapture','joint_union':after,'x':x.tolist(),'new_keys':len(set(extra)-prior_pool)});(ROOT/'round007_best_trace.json').write_text(json.dumps(trace,separators=(',',':'))+'\n')
  row={'stage':name,'indices':indices,'start_x':xstart.tolist(),'best_x':x.tolist(),'nfev':nfev,'stop':message,'before_recapture':before,'after_recapture':after,'extra_keys':len(set(extra)-prior_pool),'capture_seed':cap_seed,'capture_size':2**size};runs.append(row);print(row,flush=True)
  if trigger or after-before>1e-6:break
 z=zdecode(best[1]);p,ph,w=eng.decode(z);b,active=eng.boxes(z,pool);keys=sorted(set(map(tuple,active))|set(OLD));assert set(OLD)<=set(keys)
 out={'status':'finite Powell blocks finished; exact certification required','target_exact_upper':TARGET,'seed_input_sha256':hashlib.sha256((ROOT/'round006_certificate_input.json').read_bytes()).hexdigest(),'seed_joint_union':seedval,'blocks':runs,'optimization_evaluations':total,'new_best_count':len(trace),'elapsed_seconds':time.monotonic()-start,'final':{'scaled_x':best[1].tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'joint_union':value(best[1],pool),'active_keys':len(active),'saved_keys':len(keys),'old619_preserved':set(OLD)<=set(keys)}}
 (ROOT/'round007_results.json').write_text(json.dumps(out,indent=2)+'\n');inp={'periods':[format(v,'.15f') for v in p],'phases':[['0','0','0']]+[[format(v,'.15f') for v in q] for q in ph[1:]],'weights':['1',format(w[1],'.15f'),format(w[2],'.15f')],'keys':[list(k) for k in keys]};(ROOT/'round007_certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n');print('FINAL',out['final'],flush=True)
if __name__=='__main__':main()

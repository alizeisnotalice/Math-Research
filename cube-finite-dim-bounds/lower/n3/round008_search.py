"""Single joint11-variable Powell refinement from immutable015 accurate seed."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
import round006_search as eng
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'round007_certificate_input.json').read_text());eng.P0=np.array(list(map(float,SEED['periods'])));PH=np.array([list(map(float,q)) for q in SEED['phases'][1:]]).ravel();LOGW=np.log(np.array(list(map(float,SEED['weights'][1:]))));ORIGIN=np.r_[np.zeros(3),PH,LOGW];SCALE=np.array([.0005]*3+[.01]*6+[.2]*2);OLD=sorted(set(map(tuple,SEED['keys'])));TARGET=json.loads((ROOT/'round008_target.json').read_text())['ratio_interval'][1]
BOUNDS=[(-6.,6.)]*3+[(max(-8.,-PH[j]/.01),min(8.,(1-PH[j])/.01)) for j in range(6)]+[(-7.5,7.5)]*2

def zdecode(x):return ORIGIN+np.asarray(x)*SCALE

def value(x,keys):return float(eng.value(zdecode(x),keys))
class Stop(Exception):pass

def main():
 start=time.monotonic();pool=sorted(set(OLD)|set(eng.CLUSTERS));best=[value(np.zeros(11),pool),np.zeros(11)];seedvalue=best[0];trace=[];nfev=0;(ROOT/'round008_best_trace.json').write_text('[]\n')
 def obj(x):
  nonlocal nfev
  val=value(x,pool);nfev+=1
  if not np.isfinite(val):raise Stop('nonfinite objective anomaly')
  if val>best[0]:
   best[:]=[val,np.asarray(x).copy()];trace.append({'nfev':nfev,'joint_union':val,'scaled_x':np.asarray(x).tolist()});(ROOT/'round008_best_trace.json').write_text(json.dumps(trace,separators=(',',':'))+'\n')
  (ROOT/'round008_progress.json').write_text(json.dumps({'evaluations':nfev,'best_joint_union':best[0],'best_scaled_x':best[1].tolist(),'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
  if val>float(F(TARGET))+1e-4:raise Stop('significant finite-union gain>1e-4; prioritize exact certification')
  return -val
 try:
  res=minimize(obj,np.zeros(11),method='Powell',bounds=BOUNDS,options={'maxfev':1000,'xtol':1e-7,'ftol':1e-10});message=str(res.message)
 except Stop as e:message=str(e)
 x=best[1];before=value(x,pool);rec=eng.capture(zdecode(x),qmc.Sobol(3,scramble=True,seed=20261025801).random_base2(15));newpool=sorted(set(pool)|set(rec));after=value(x,newpool)
 if after>before:
  trace.append({'nfev':nfev,'event':'single terminal recapture; no additional optimizer evaluation','joint_union':after,'scaled_x':x.tolist(),'extra_keys':len(set(rec)-set(pool))});(ROOT/'round008_best_trace.json').write_text(json.dumps(trace,separators=(',',':'))+'\n')
 p,ph,w=eng.decode(zdecode(x));b,active=eng.boxes(zdecode(x),newpool);keys=sorted(set(map(tuple,active))|set(OLD));assert set(OLD)<=set(keys)
 out={'status':'single joint batch completed; exact certification required','target_exact_upper':TARGET,'seed_input_sha256':hashlib.sha256((ROOT/'round007_certificate_input.json').read_bytes()).hexdigest(),'seed_joint_union':seedvalue,'optimizer':'bounded Powell all11 simultaneously','optimization_evaluations':nfev,'stop':message,'new_best_count':len(trace),'capture_seed':20261025801,'capture_size':32768,'before_terminal_recapture':before,'after_terminal_recapture':after,'extra_keys':len(set(rec)-set(pool)),'elapsed_seconds':time.monotonic()-start,'final':{'scaled_x':x.tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'joint_union':after,'active_keys':len(active),'saved_keys':len(keys),'old637_preserved':set(OLD)<=set(keys)}}
 (ROOT/'round008_results.json').write_text(json.dumps(out,indent=2)+'\n');inp={'periods':[format(v,'.15f') for v in p],'phases':[['0','0','0']]+[[format(v,'.15f') for v in q] for q in ph[1:]],'weights':['1',format(w[1],'.15f'),format(w[2],'.15f')],'keys':[list(k) for k in keys]};(ROOT/'round008_certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n');print(out,flush=True)
if __name__=='__main__':main()

"""Final 1000-evaluation refinement, preserving every new best and recapturing."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
import round005_search as eng
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'round005b_certificate_input.json').read_text());eng.P0=np.array(list(map(float,SEED['periods'])));PH=np.array(list(map(float,SEED['phases'][1])));EPS=float(F(SEED['weights'][1]));OLD=sorted(set(tuple(k) for k in SEED['witness_keys']));TARGET=json.loads((ROOT/'round005c_target.json').read_text())['ratio_interval'][1]
class Stop(Exception):pass

def main():
 start=time.monotonic();z=np.r_[np.zeros(3),PH,np.log(EPS)];rec=eng.capture(z,qmc.Sobol(3,scramble=True,seed=20261025521).random_base2(15));pool=sorted(set(OLD)|set(rec));initialextra=len(set(rec)-set(OLD));seedvalue=eng.value(z,pool);nfev=0;best=[seedvalue,z.copy()];trace=[];trigger=False
 def obj(x):
  nonlocal nfev,trigger
  val=eng.value(x,pool);nfev+=1
  if not np.isfinite(val):trigger=True;raise Stop('nonfinite union anomaly')
  if val>best[0]:
   best[:]=[val,np.asarray(x).copy()];trace.append({'nfev':nfev,'joint_union':val,'z':np.asarray(x).tolist()});(ROOT/'round005c_best_trace.json').write_text(json.dumps(trace,separators=(',',':'))+'\n')
  if val>float(F(TARGET))+1e-4:trigger=True;raise Stop('gain above independently audited checkpoint >1e-4')
  return -val
 simplex=np.tile(z,(8,1));steps=[.0005,.0001,.0001,.01,.01,.01,.25]
 for j in range(7):simplex[j+1,j]+=steps[j]
 if seedvalue>float(F(TARGET))+1e-8 and initialextra:
  stop='initial recapture produced material structural union gain; certify before optimization'
 else:
  try:
   res=minimize(obj,z,method='Nelder-Mead',options={'maxfev':1000,'initial_simplex':simplex,'xatol':1e-10,'fatol':1e-12});stop=str(res.message)
  except Stop as e:stop=str(e)
 z=best[1];before=eng.value(z,pool);extra=eng.capture(z,qmc.Sobol(3,scramble=True,seed=20261025522).random_base2(16));finalkeys=sorted(set(pool)|set(extra));after=eng.value(z,finalkeys);p,ph,w=eng.decode(z);bb,active=eng.make_boxes(z,finalkeys);keys=sorted(set(active)|set(OLD));assert set(OLD)<=set(keys)
 out={'status':'final finite batch finished; exact certification required','target_checkpoint_upper':TARGET,'seed_input_sha256':hashlib.sha256((ROOT/'round005b_certificate_input.json').read_bytes()).hexdigest(),'seed_value':seedvalue,'capture_seeds':[20261025521,20261025522],'capture_sizes':[32768,65536],'initial_extra_keys':initialextra,'optimization_evaluations':nfev,'new_best_count':len(trace),'stop':stop,'material_threshold':1e-4,'before_final_recapture':before,'after_final_recapture':after,'final_extra_keys':len(set(extra)-set(pool)),'elapsed_seconds':time.monotonic()-start,'final':{'z':z.tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'joint_union':after,'active_keys':len(active),'saved_keys':len(keys),'original_170_preserved':set(OLD)<=set(keys)}}
 (ROOT/'round005c_results.json').write_text(json.dumps(out,indent=2)+'\n');inp={'periods':[format(x,'.15f') for x in p],'phases':[['0','0','0'],[format(x,'.15f') for x in ph[1]]],'weights':['1',format(w[1],'.15f')],'witness_keys':[list(k) for k in keys]};(ROOT/'round005c_certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n');print(out,flush=True)
if __name__=='__main__':main()

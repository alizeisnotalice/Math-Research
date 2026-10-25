"""Audited-checkpoint refinement; retains old evidence and all 170 keys."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
import round005_search as eng
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'round005_certificate_input.json').read_text());eng.P0=np.array(list(map(float,SEED['periods'])));PH=np.array(list(map(float,SEED['phases'][1])));EPS=float(F(SEED['weights'][1]));OLD=sorted(set(tuple(k) for k in SEED['witness_keys']));TARGET=json.loads((ROOT/'round005_certificate_result.json').read_text())['ratio_interval'][1]
class Enough(Exception):pass

def main():
 start=time.monotonic();z=np.r_[np.zeros(3),PH,np.log(EPS)];extra=eng.capture(z,qmc.Sobol(3,scramble=True,seed=20261025511).random_base2(15));pool=sorted(set(OLD)|set(extra));initial=len(set(extra)-set(OLD));seedval=eng.value(z,pool);runs=[];count=0
 for stage in range(2):
  before_start=z.copy();best=[-1.,None];nfev=0;trigger=False
  def obj(x):
   nonlocal count,nfev,trigger
   value=eng.value(x,pool);count+=1;nfev+=1
   if value>best[0]:best[:]=[value,np.asarray(x).copy()]
   if value>float(F(TARGET))+1e-6:trigger=True;raise Enough()
   return -value
  simplex=np.tile(z,(8,1));steps=np.array([.0003,.00005,.00005,.005,.005,.005,.2])/(2**stage)
  for j in range(7):simplex[j+1,j]+=steps[j]
  try:
   res=minimize(obj,z,method='Nelder-Mead',options={'maxfev':500,'initial_simplex':simplex,'xatol':1e-10,'fatol':1e-12});z=best[1];stop=str(res.message)
  except Enough:z=best[1];stop='gain above checkpoint >1e-6; immediately stop optimization'
  before=eng.value(z,pool);rec=eng.capture(z,qmc.Sobol(3,scramble=True,seed=20261025512+stage).random_base2(16));newpool=sorted(set(pool)|set(rec));after=eng.value(z,newpool);row={'stage':stage,'start':before_start.tolist(),'z':z.tolist(),'nfev':nfev,'stop':stop,'before_recapture':before,'after_recapture':after,'extra_keys':len(set(rec)-set(pool))};runs.append(row);print(row,flush=True);pool=newpool
  if trigger or after>float(F(TARGET))+1e-6:break
 p,ph,w=eng.decode(z);bb,active=eng.make_boxes(z,pool);keys=sorted(set(active)|set(OLD));assert set(OLD)<=set(keys)
 out={'status':'finite-union candidate; exact comparison required','target_checkpoint_upper':TARGET,'seed_input_sha256':hashlib.sha256((ROOT/'round005_certificate_input.json').read_bytes()).hexdigest(),'seed_value':seedval,'initial_extra_keys':initial,'capture_seeds':[20261025511,20261025512,20261025513],'optimization_evaluations':count,'runs':runs,'elapsed_seconds':time.monotonic()-start,'final':{'z':z.tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'joint_union':eng.value(z,pool),'active_keys':len(active),'saved_keys':len(keys),'original_170_preserved':set(OLD)<=set(keys)}}
 (ROOT/'round005b_results.json').write_text(json.dumps(out,indent=2)+'\n');inp={'periods':[format(x,'.15f') for x in p],'phases':[['0','0','0'],[format(x,'.15f') for x in ph[1]]],'weights':['1',format(w[1],'.15f')],'witness_keys':[list(k) for k in keys]};(ROOT/'round005b_certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n');print('FINAL',out['final'],flush=True)
if __name__=='__main__':main()

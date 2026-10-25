"""Reproduce structured diagnostic and deterministic final-path recovery from saved completed starts."""
import sys,json,ast,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
from round003_search import capture,value,boxes,decode

def rectangle(key):
 cnt=0;low=[];high=[]
 for b in range(3):
  l=key[3*b:3*b+3];h=key[9+3*b:12+3*b]
  if not all(x<=y for x,y in zip(l,h)):continue
  cnt+=int(np.prod(np.array(h)-l+1));low.append([3*l[0]+b,l[1],l[2]]);high.append([3*h[0]+b,h[1],h[2]])
 bbox=int(np.prod(np.max(high,axis=0)-np.min(low,axis=0)+1))
 return cnt==bbox
start=time.monotonic();t=qmc.Sobol(3,scramble=True,seed=20261025303).random_base2(12);cases=[];pool=set();allcases=[]
zs=[('baseline',np.zeros(9))]
for mode in ['monotone','opposite','rotated']:
 for f in [.01,.04,.1,.2]:
  for mass in [[0,0],[.1,-.05],[-.1,.1]]:
   z=np.zeros(9);z[3:7]={'monotone':[f,0,2*f,0],'opposite':[f,0,-f,0],'rotated':[f,0,0,f]}[mode];z[7:]=mass;zs.append((f'{mode}_{f}_{mass}',z))
for name,z in zs:
 keys=capture(z,t);rect=[k for k in keys if rectangle(k)];val=value(z,keys);rv=value(z,rect);pool.update(keys)
 c={'name':name,'z':z.tolist(),'joint_union':val,'rectangular_union':rv,'staircase_gain':val-rv,'key_count':len(keys),'rectangle_key_count':len(rect)};cases.append(c);allcases.append((c,keys,rect))
pool=sorted(pool);runs=json.loads((ROOT/'round003_completed_runs.json').read_text());best=max(runs,key=lambda x:x['after_recapture']);keys=sorted(set(pool)|set(capture(best['z'],qmc.Sobol(3,scramble=True,seed=20261025304).random_base2(14))));count=0

def obj(z):
 global count;count+=1;return -value(z,keys)
res=minimize(obj,best['z'],method='Nelder-Mead',options={'maxfev':400,'xatol':1e-8,'fatol':1e-10});z=res.x;before=value(z,keys);rec=capture(z,qmc.Sobol(3,scramble=True,seed=20261025305).random_base2(14));finalkeys=sorted(set(keys)|set(rec));after=value(z,finalkeys);p,ph,w=decode(z);bb,active=boxes(z,finalkeys)
final={'z':z.tolist(),'periods':p.tolist(),'phases':ph.tolist(),'weights':w.tolist(),'before_recapture':before,'joint_union':after,'extra_keys':len(set(rec)-set(keys)),'active_key_count':len(active),'nfev_final':res.nfev,'stop':str(res.message)}
out={'status':'numerical finite-witness unions only','initial_seed':20261025303,'recapture_seeds':[20261025304,20261025305],'initial_cases':cases,'initial_pool_size':len(pool),'runs':runs,'unique_optimization_evaluations':1800,'actual_optimization_evaluations':2200,'recovery_optimization_evaluations':count,'diagnostic_replay_seconds':time.monotonic()-start,'serialization_recovery':'1800 original optimizations completed; final result serialization failed on numpy.int64; same final 400 deterministic evaluations replayed after integer conversion fix','final':final}
(ROOT/'round003_keys.json').write_text(json.dumps([list(k) for k in active],separators=(',',':'))+'\n');(ROOT/'round003_results.json').write_text(json.dumps(out,indent=2)+'\n')
# Rationalize positions and masses without deriving anything from R3 decimals.
def exact_input(z,keys):
 p,ph,w=decode(z);ww=[format(w[0],'.12f'),format(w[1],'.12f')]
 from fractions import Fraction as F
 ww.append(str(F(3)-sum(map(F,ww),F(0))))
 return {'periods':[format(x,'.12f') for x in p],'phases':[[format(x,'.12f') for x in q] for q in ph],'weights':ww,'keys':[list(map(int,k)) for k in keys]}
c,full,rect=max(allcases[1:],key=lambda item:item[0]['staircase_gain'])
inputs={'optimized':exact_input(z,active),'staircase':exact_input(c['z'],full),'rectangular':exact_input(c['z'],rect)}
(ROOT/'round003_certificate_input.json').write_text(json.dumps(inputs,indent=2)+'\n');(ROOT/'round003_staircase_diagnostic.json').write_text(json.dumps({'selected':c,'nonrectangular_keys':[list(k) for k in full if not rectangle(k)],'criterion':'Map orbit index k to original layer coordinate (3*k1+b,k2,k3); selected count smaller than cardinality of its minimum enclosing integer rectangle means a nonrectangular staircase subset.'},indent=2)+'\n')
print('FINAL',final);print('DIAGNOSTIC',c)

"""One targeted eight-basis LP; numerical only."""
import json,itertools,math,time
import numpy as np
from scipy.optimize import linprog
from scipy.stats import qmc
import round010_joint_probe as old
P=old.P
SPECS=[((0,0,0),1.),((.75,0,0),.5),((0,0,0),1.3),((.4,.4,.4),.5),((.8,.4,.4),.35),((.4,.4,0),.4),((.4,0,0),.4),((0,0,0),.6)]
def orbit(seed):return sorted(set(tuple(float(a*b) for a,b in zip(t,sg)) for t in set(itertools.permutations(seed)) for sg in itertools.product((-1,1),repeat=3)))
BASIS=[(np.array(orbit(s)),R) for s,R in SPECS]
VOL=np.array([len(a)*R**3*32*math.pi/105 for a,R in BASIS])
def bs(x):return np.array([sum(old.unit_b((x-a)/R) for a in centers) for centers,R in BASIS])
def ls(x):return np.array([sum(old.unit_lb((x-a)/R)/R for a in centers) for centers,R in BASIS])
if __name__=='__main__':
 start=time.time();cache=json.loads((P/'results/round010_source_cache.json').read_text())['rows'];best={}
 for rr in cache.values():
  key=tuple(rr['x'])
  if key not in best or rr['tolerance']<best[key]['tolerance']:best[key]=rr
 rows=[r for r in best.values() if r['status']=='numerical_only'];A=[ls(np.array(r['x'])).tolist() for r in rows];rhs=[r['mu'] for r in rows];ids=[{'kind':'source','x':r['x']} for r in rows]
 for xx in itertools.combinations_with_replacement(np.linspace(0,1.5,13),3):
  x=np.array(xx[::-1]);b=bs(x)
  if not np.any(b):continue
  g=float(old.green(x))
  if not math.isfinite(g):continue
  cap=g-(old.GE if max(x)<=1 else old.C)
  if cap < -1e-10:raise RuntimeError('G below normalizer')
  A.append(b.tolist());rhs.append(max(0,cap));ids.append({'kind':'domain','x':x.tolist(),'unclipped_rhs':cap})
 A=np.array(A);rhs=np.array(rhs);res=linprog(-VOL,A_ub=A,b_ub=rhs,bounds=[(0,None)]*8,method='highs');report={'status':'one finite numerical LP, not continuous bound','specs':[{'seed':s,'radius':R,'centers':a.tolist(),'mass':float(v)} for (s,R),(a,_),v in zip(SPECS,BASIS,VOL)],'source_rows':len(rows),'domain_rows':len(A)-len(rows),'success':bool(res.success),'message':res.message,'module_sha256':old.SHA}
 if res.success:
  w=res.x;sl=rhs-A@w;report.update(weights=w.tolist(),mass_reduction=float(VOL@w),optimistic_cost_estimate=8.171963453262865-float(VOL@w)/old.D,active_rows=[dict(ids[i],slack=float(sl[i])) for i in np.argsort(sl)[:12]])
  (P/'results/round010b_joint_probe.json').write_text(json.dumps(report,indent=2)+'\n');print('LP cost',report['optimistic_cost_estimate'],flush=True)
  points=np.sort(qmc.Sobol(3,scramble=True,seed=20261025).random_base2(6)*1.6,axis=1)[:,::-1];val=old.run_jobs([(tuple(x),1e-7) for x in points],cache)
  for r in val:
   x=np.array(r['x']);r['new_source_slack']=float(r['mu']-ls(x)@w) if r['status']=='numerical_only' else None
   g=float(old.green(x));r['new_dominance_slack']=float(max(g-old.C,0)-bs(x)@w-(old.GE-old.C if max(x)<=1 else 0)) if math.isfinite(g) else 'infinity'
  report['validation']=val;report['negative_source_rows']=[r['x'] for r in val if r['new_source_slack'] is not None and r['new_source_slack'] < -1e-7];report['negative_dominance_rows']=[r['x'] for r in val if isinstance(r['new_dominance_slack'],float) and r['new_dominance_slack'] < -1e-7]
  (P/'results/round010b_matlab_input.json').write_text(json.dumps({'A':A.tolist(),'rhs':rhs.tolist(),'objective':(-VOL).tolist(),'weights':w.tolist(),'gain':float(VOL@w),'row_ids':ids},indent=2)+'\n')
 report['elapsed_seconds']=time.time()-start;report['scope']='one8basisLP and64fixedvalidation; unresolved quadrature excluded, no global certificate';(P/'results/round010b_joint_probe.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('validation','specs')},indent=2))

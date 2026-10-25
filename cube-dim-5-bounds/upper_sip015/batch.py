#!/usr/bin/env python3
"""One registered finite SIP batch; numerical masters are diagnostic only."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,time
HERE=Path(__file__).resolve().parent
PY='/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python'
MOM=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/astra_upper006/verification.json')
CAPDIR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper-audit-cont/cube-dim-5-bounds/upper_I027')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def worker(index):
 import numpy as np
 import scipy
 from scipy.optimize import linprog,minimize_scalar
 started=time.monotonic();data=json.loads((HERE/'master_input.json').read_text());z=np.array(data['nodes']);Q=float(F(data['Q']));m=[float(F(x)) for x in data['moments']]
 exps=np.array([0,4,24,44,64,84]);mat=z[:,None]**exps
 cost=np.array([Q,*m]);bounds=[(0,10)]+[(-1e6,1e6)]*5
 sol=linprog(cost,A_ub=-mat,b_ub=-z,bounds=bounds,method='highs',options={'time_limit':20.,'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
 out={'phase':index,'pid':os.getpid(),'scipy_version':scipy.__version__,'success':sol.success,'status':int(sol.status),'message':sol.message,'source_sha256':sha(Path(__file__))}
 if not sol.success:
  (HERE/f'master_{index}.json').write_text(json.dumps(out,indent=2)+'\n');return
 coeff=sol.x
 def poly(x):return sum(float(a)*x**int(e) for a,e in zip(coeff,exps))-x
 grid=np.linspace(0,1,4097);values=np.array([poly(x) for x in grid]);tests=[(float(grid[i]),float(values[i]),'sample') for i in (0,int(np.argmin(values)),4096)]
 for i in range(1,4096):
  if values[i]<=values[i-1] and values[i]<=values[i+1]:
   opt=minimize_scalar(poly,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-14,'maxiter':200});tests.append((float(opt.x),float(opt.fun),'bounded_local'))
 p=np.zeros(85);p[exps]=coeff;p[1]=-1
 roots=np.polynomial.polynomial.polyroots(np.polynomial.polynomial.polyder(p))
 for root in roots:
  if abs(root.imag)<1e-7 and 0<float(root.real)<1:tests.append((float(root.real),float(poly(root.real)),'floating_derivative_root'))
 best=min(tests,key=lambda x:x[1]);duals=-sol.ineqlin.marginals
 support=[{'node_index':int(i),'z':float(z[i]),'lambda':float(q),'w':float(q*z[i]**4),'source_node':data['node_provenance'][i]} for i,q in enumerate(duals) if q>1e-9]
 out.update(objective=float(sol.fun),coefficients=coeff.tolist(),separation={'z':best[0],'P_display_only':best[1],'method':best[2],'tested_count':len(tests)},dual_support=support,
  dual_objective=float(np.dot(duals,z)),dual_equation_residual=(mat.T@duals-cost).tolist(),lower_bound_marginals=sol.lower.marginals.tolist(),upper_bound_marginals=sol.upper.marginals.tolist(),wall_seconds=time.monotonic()-started)
 (HERE/f'master_{index}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('phase','pid','success','objective','separation','dual_support','wall_seconds')}),flush=True)

def main():
 moments=json.loads(MOM.read_text())['moments'];m=['1']+[moments[str(k)][0] for k in range(1,5)]
 registration=HERE/'registration.json';reg=json.loads(registration.read_text())
 nodes=[i/160 for i in range(161)];provenance=[{'rational':str(F(i,160))} for i in range(161)]
 out={'status':'running','parent_pid':os.getpid(),'source_sha256':sha(Path(__file__)),'registration_sha256':sha(registration),'moment_source_sha256':sha(MOM),'I027_sources':{p.name:sha(p) for p in [CAPDIR/'manifest.json',CAPDIR/'results.json',CAPDIR/'report.md']},'Q':reg['Q'],'moments':m,'workers':[]}
 for index in range(5):
  (HERE/'master_input.json').write_text(json.dumps({'Q':reg['Q'],'moments':m,'nodes':nodes,'node_provenance':provenance},indent=2)+'\n')
  begun=time.monotonic();p=subprocess.Popen([PY,'-B',str(Path(__file__)),'--worker',str(index)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
  timedout=False
  try:stdout,stderr=p.communicate(timeout=30)
  except subprocess.TimeoutExpired:
   timedout=True;p.kill();stdout,stderr=p.communicate()
  rec={'phase':index,'actual_pid':p.pid,'exit_code':p.returncode,'timeout':timedout,'wall_seconds':time.monotonic()-begun,'stdout':stdout,'stderr':stderr,'input_sha256':sha(HERE/'master_input.json')}
  out['workers'].append(rec);(HERE/'process_results.json').write_text(json.dumps(out,indent=2)+'\n');print(stdout,end='',flush=True)
  if p.returncode!=0 or timedout:out['status']='stopped_process_failure';break
  result=json.loads((HERE/f'master_{index}.json').read_text())
  if not result['success']:out['status']='stopped_solver_status';break
  if result['separation']['P_display_only']>=-1e-9:out['status']='completed_numerical_separation_no_violation';break
  if index==4:out['status']='completed_four_exchange_cap';break
  z=result['separation']['z'];nodes.append(z);provenance.append({'float_exchange_value':z,'rational_for_exact_witness':str(F(round(z*10**8),10**8))})
 out['completed_at_unix']=time.time();(HERE/'process_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'workers':len(out['workers'])}),flush=True)
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='--worker':worker(int(sys.argv[2]))
 else:main()

"""One MATLAB discovery LP; saved exact local rows checked after rationalization."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json,os,sys,time,subprocess
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import savemat
from scipy import sparse
ROOT=Path(__file__).resolve().parent
start=time.monotonic();OUT=ROOT/'solve_result.json'
if OUT.exists():raise SystemExit('Existing solver run; no retry.')
r={'status':'started','LP_requested':0,'LP_started':0,'LP_completed':0,'exact_rows_checked':0,'candidate_exported':False}
def save():
 r['elapsed_seconds']=time.monotonic()-start
 OUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def guard():
 if time.monotonic()-start>=175:raise TimeoutError('175s internal deadline.')
def read(path):return json.loads(path.read_text())
def checked(path,wanted):
 assert sha256(path.read_bytes()).hexdigest()==wanted,path.name
 return path
save()
try:
 assert os.getpgrp()==os.getpid(),'Must use owned wrapper process group.'
 b=read(ROOT/'build_result.json');assert b['complete']
 directions=read(checked(ROOT/'directions.json',b['directions_sha256']))['directions'];assert len(directions)==12
 assert [d['old_column'] for d in directions]==[2,3,4,5,198,199,201,202,222,244,267,268]
 rows=read(checked(ROOT/'rows.json',b['rows_sha256']))['rows'];assert len(rows)==16272
 baseline_path=ROOT.parent/'round01/scaled_certificate.json'
 baseline=read(checked(baseline_path,'21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'))
 c1=Q(baseline['c']);alpha=Q(baseline['alpha']);R=Q(baseline['R']);h=Q(baseline['h'])
 old=ROOT.parent/'round07';cache=ROOT.parent/'round06/cache.local'
 manifest=read(old/'local_artifact_manifest.json')['files']
 base_manifest=read(ROOT.parent/'round06/result.json')['local_cache_files']
 for item in base_manifest:checked(ROOT.parent/'round06'/item['file'],item['sha256']);guard()
 point_path=ROOT.parent/'round06/point_manifest.json'
 points=read(checked(point_path,'513c8d953695d9ddd1c7a097ee9af2ab624f79b875330d02e0d42ad2dcc70ebd'))
 allmat=[];allrhs=[];targets=[]
 for kind,n in [('K',26604),('A',24970)]:
  name=kind+'_directions_288.local.npy';checked(old/name,manifest[name]['sha256'])
  m=np.load(old/name,mmap_mode='r');assert m.shape==(n,288)
  sing=np.load(cache/(kind+'_singular.npy'));base=np.load(cache/(kind+'_baseline.npy'))
  target=np.array([float(Q(u)+Q(v)<=1) for u,v in points[kind]]) if kind=='K' else np.zeros(n)
  slack=base-target;assert np.min(slack)>=0
  mat=np.column_stack([sing,m[:,1]]+[m[:,int(d['old_column'])] for d in directions]+[np.zeros(n)]*12)
  allmat.append(mat);allrhs.append(target+slack/2);targets.append(target)
 local=np.array([[float(Q(row['c_column'])),float(Q(row['s_column']))]+[float(Q(v)) for v in row['t_midpoints']]+[-float(Q(v)) for v in row['t_penalties']] for row in rows])
 lrhs=np.array([float(Q(row['baseline_slack'])/2) for row in rows])
 assert local.shape==(16272,26)
 absrows=np.zeros((24,26))
 for j in range(12):absrows[2*j,2+j]=1;absrows[2*j+1,2+j]=-1;absrows[2*j:2*j+2,14+j]=-1
 A=np.vstack([-allmat[0],-allmat[1],-local,absrows]);rhs=np.r_[-allrhs[0],-allrhs[1],-lrhs,np.zeros(24)]
 assert A.shape==(67870,26) and np.isfinite(A).all() and np.isfinite(rhs).all()
 def mult(i,j):return 1 if i==j==0 else 4 if j==0 or i==j else 8
 smass=sum((h*h*mult(*map(int,k.split(',')))*Q(v)/2 for k,v in baseline['a'].items()),Q(0))
 f=np.r_[float(alpha)*float(R)**(2-float(alpha))/(2-float(alpha)),float(smass),[float(Q(d['half_mass'])) for d in directions],np.zeros(12)]
 bound=[Q(d['bound']) for d in directions]
 lb=np.r_[float(c1/2),.5,[-float(v) for v in bound],np.zeros(12)]
 ub=np.r_[float(3*c1/2),1.5,[float(v) for v in bound],[float(v) for v in bound]]
 basevec=np.r_[float(c1),1.,np.zeros(24)]
 r['floating_baseline_min_slack']=float(np.min(rhs-A@basevec))
 assert r['floating_baseline_min_slack']>=-1e-8
 savemat(ROOT/'model.local.mat',{'A':sparse.csr_matrix(A),'b':rhs[:,None],'f':f[:,None],'lb':lb[:,None],'ub':ub[:,None]})
 r['matrix_sha256']=sha256((ROOT/'model.local.mat').read_bytes()).hexdigest()
 r['solver_source_sha256']=sha256((ROOT/'solver.m').read_bytes()).hexdigest();save();guard()
 script=str(ROOT/'solver.m').replace("'","''")
 with (ROOT/'matlab.log').open('x') as log:
  r['LP_requested']=1;save()
  proc=subprocess.Popen(['/Applications/MATLAB_R2026a.app/bin/matlab','-batch',"run('"+script+"')"],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
  r['matlab_pid']=proc.pid;save()
  code=proc.wait(timeout=max(.1,175-(time.monotonic()-start)))
 r['matlab_exit_code']=code;guard();sol=read(ROOT/'solver_result.json')
 assert sol['input_sha256']==r['matrix_sha256']
 r['LP_started']=sol.get('calls_started',0)
 r['LP_completed']=sol.get('calls_completed',0);r['solver_status']=sol['status'];r['exitflag']=sol.get('exitflag');save()
 assert code==0 and sol['status']=='returned' and sol['exitflag']==1,'Solver status rejected, no retry.'
 x=np.array(sol['x']).reshape(-1);assert x.shape==(26,) and np.isfinite(x).all()
 r['float_cost']=float(f@x);r['float_worst_slack']=float(np.min(rhs-A@x));r['float_bound_slack']=float(min(np.min(x-lb),np.min(ub-x)))
 q=[Q(round(float(v)*2**50),2**50) for v in x[:14]];q+=list(map(abs,q[2:14]))
 assert c1/2<=q[0]<=3*c1/2 and Q(1,2)<=q[1]<=Q(3,2)
 assert all(abs(q[2+j])<=bound[j] for j in range(12))
 bad=0;min_slack=None;worst=None
 for row in rows:
  guard();v=q[0]*Q(row['c_column'])+q[1]*Q(row['s_column'])+sum((q[2+j]*Q(row['t_midpoints'][j])-q[14+j]*Q(row['t_penalties'][j]) for j in range(12)),Q(0))-Q(row['baseline_slack'])/2
  bad+=int(v<0)
  if min_slack is None or v<min_slack:min_slack=v;worst={'node_id':row['leaf_node_id'],'bernstein_index':row['bernstein_index']}
  r['exact_rows_checked']+=1
 r['exact_bad_count']=bad;r['exact_min_slack']=str(min_slack);r['worst_exact_row']=worst
 qf=np.array(list(map(float,q)));r['rational_float_cost']=float(f@qf)
 r['point_checks']={}
 for k,kind,tol in [(0,'K',1e-9),(1,'A',1e-8)]:
  vals=allmat[k]@qf;assert np.isfinite(vals).all()
  r['point_checks'][kind]={'strengthened_bad':int(np.count_nonzero(vals-allrhs[k]<-tol)),'original_bad':int(np.count_nonzero(vals-targets[k]<-tol)),'tolerance':tol}
 r['point_bad_count']=sum(v['strengthened_bad']+v['original_bad'] for v in r['point_checks'].values())
 r['point_rows_checked']=51574
 (ROOT/'parameters.json').write_text(json.dumps({'q':[str(v) for v in q],'direction_names':[d['name'] for d in directions]},indent=2)+'\n')
 if bad==0 and r['point_bad_count']==0 and float(f@qf)<3.590901146403276-5e-7:
  coeff={k:q[1]*Q(v) for k,v in baseline['a'].items()}
  for t,d in zip(q[2:14],directions):
   for k,v in d['coefficients'].items():coeff[k]=coeff.get(k,Q(0))+t*Q(v)
  candidate={k:baseline[k] for k in ('alpha','R','h','N')};candidate.update(c=str(q[0]),a={k:str(v) for k,v in coeff.items() if v})
  (ROOT/'candidate.json').write_text(json.dumps(candidate,indent=2)+'\n');r['candidate_exported']=True
 r['status']='one_finite_lp_and_exact_local_check_closed'
except Exception as e:r.update(status='failed_or_rejected',reason=f'{type(e).__name__}: {e}')
r['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();save();print(json.dumps(r,indent=2))

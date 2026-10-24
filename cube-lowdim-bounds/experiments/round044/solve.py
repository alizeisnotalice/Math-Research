"""One MATLAB discovery LP; fixed contraction and structural rationalization; no retry."""
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
r={'status':'started','LP_requested':0,'LP_started':0,'LP_completed':0,'exact_rows_checked':0,'candidate_exported':False,'cost_enclosures_started':0,'cost_enclosures_completed':0,'point_rows_checked':0,'inputs':[]}
def save():
 r['elapsed_seconds']=time.monotonic()-start
 OUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def guard():
 if time.monotonic()-start>=175:raise TimeoutError('175s internal deadline.')
def read(path):return json.loads(path.read_text())
def checked(path,wanted):
 actual=sha256(path.read_bytes()).hexdigest();assert actual==wanted,path.name
 r['inputs'].append({'file':os.path.relpath(path,ROOT.parent),'sha256':actual})
 return path
save()
try:
 assert os.getpgrp()==os.getpid(),'Must use owned wrapper process group.'
 b=read(ROOT/'build_result.json');assert b['complete']
 directions=read(checked(ROOT/b['directions_file'],b['directions_sha256']))['directions'];assert len(directions)==12
 assert [d['old_column'] for d in directions]==[2,3,4,5,198,199,201,202,222,244,267,268]
 old_rows=read(checked(ROOT/b['legacy_rows_file'],b['legacy_rows_sha256']))['rows'];assert len(old_rows)==16272
 new_rows=read(checked(ROOT/b['new_rows_file'],b['new_rows_sha256']))['rows'];assert len(new_rows)==20448
 rows=old_rows+new_rows;assert len(rows)==36720
 r['build_result_sha256']=sha256((ROOT/'build_result.json').read_bytes()).hexdigest()
 r['new_rows_sha256']=b['new_rows_sha256'];r['legacy_rows_sha256']=b['legacy_rows_sha256']
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
 assert local.shape==(36720,26)
 absrows=np.zeros((24,26))
 for j in range(12):absrows[2*j,2+j]=1;absrows[2*j+1,2+j]=-1;absrows[2*j:2*j+2,14+j]=-1
 structure=np.zeros((1,26));structure[0,0]=-1;structure[0,1]=float(c1)
 A=np.vstack([-allmat[0],-allmat[1],-local,absrows,structure]);rhs=np.r_[-allrhs[0],-allrhs[1],-lrhs,np.zeros(25)]
 assert A.shape==(88319,26) and np.isfinite(A).all() and np.isfinite(rhs).all()
 def mult(i,j):return 1 if i==j==0 else 4 if j==0 or i==j else 8
 smass=sum((h*h*mult(*map(int,k.split(',')))*Q(v)/2 for k,v in baseline['a'].items()),Q(0))
 f=np.r_[float(alpha)*float(R)**(2-float(alpha))/(2-float(alpha)),float(smass),[float(Q(d['half_mass'])) for d in directions],np.zeros(12)]
 bound=[Q(d['bound']) for d in directions]
 lb=np.r_[float(c1/2),.5,[-float(v) for v in bound],np.zeros(12)]
 ub=np.r_[float(3*c1/2),1.5,[float(v) for v in bound],[float(v) for v in bound]]
 basevec=np.r_[float(c1),1.,np.zeros(24)]
 r['floating_baseline_min_local_slack']=float(np.min(local@basevec-lrhs))
 r['floating_baseline_min_slack']=min(r['floating_baseline_min_local_slack'],0.)
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
 r['solver_result_sha256']=sha256((ROOT/'solver_result.json').read_bytes()).hexdigest()
 r['matlab_log_sha256']=sha256((ROOT/'matlab.log').read_bytes()).hexdigest()
 assert sha256((ROOT/'model.local.mat').read_bytes()).hexdigest()==r['matrix_sha256']
 assert sha256((ROOT/'solver.m').read_bytes()).hexdigest()==r['solver_source_sha256']
 assert sol['input_sha256']==r['matrix_sha256']
 r['LP_started']=sol.get('calls_started',0)
 r['LP_completed']=sol.get('calls_completed',0);r['solver_status']=sol['status'];r['exitflag']=sol.get('exitflag');save()
 assert code==0 and sol['status']=='returned' and sol['exitflag']==1,'Solver status rejected, no retry.'
 x=np.array(sol['x']).reshape(-1);assert x.shape==(26,) and np.isfinite(x).all()
 r['float_cost']=float(f@x);float_slacks=rhs-A@x
 assert np.isfinite(float_slacks).all() and np.isfinite(r['float_cost'])
 r['float_worst_slack']=float(np.min(float_slacks));r['float_bound_slack']=float(min(np.min(x-lb),np.min(ub-x)))
 r['float_point_checks']={};offset=0
 for k,kind,tol in [(0,'K',1e-9),(1,'A',1e-8)]:
  strengthened=float_slacks[offset:offset+len(targets[k])];original=strengthened+allrhs[k]-targets[k];offset+=len(targets[k])
  r['float_point_checks'][kind]={'strengthened_bad':int(np.count_nonzero(strengthened < -tol)),'original_bad':int(np.count_nonzero(original < -tol)),'tolerance':tol}
 r['point_rows_checked']=51574
 lam=Q(1023,1024)
 xc=[Q(float(v)) for v in x[:14]]
 contracted=[(1-lam)*c1+lam*xc[0],(1-lam)+lam*xc[1]]+[lam*v for v in xc[2:14]]
 g0=contracted[0]-contracted[1]*c1
 grid=2**50;g=Q(round(g0*grid),grid);scale=Q(round(contracted[1]*grid),grid)
 q=[scale*c1+g,scale]+[Q(round(v*grid),grid) for v in contracted[2:14]];q+=list(map(abs,q[2:14]))
 r['fixed_lambda']=str(lam);r['structural_g_before_rounding']=str(g0);r['structural_g']=str(g)
 (ROOT/'parameters.json').write_text(json.dumps({'q':[str(v) for v in q],'g':str(g),'lambda':str(lam),'direction_names':[d['name'] for d in directions]},indent=2)+'\n')
 assert g>=0 and q[0]-q[1]*c1==g,'Structural quantization rejected; no clamping or retry.'
 assert c1/2<=q[0]<=3*c1/2 and Q(1,2)<=q[1]<=Q(3,2)
 assert all(abs(q[2+j])<=bound[j] for j in range(12))
 bad=0;min_slack=None;worst=None
 for row in rows:
  guard();v=q[0]*Q(row['c_column'])+q[1]*Q(row['s_column'])+sum((q[2+j]*Q(row['t_midpoints'][j])-q[14+j]*Q(row['t_penalties'][j]) for j in range(12)),Q(0))-Q(row['baseline_slack'])/2
  bad+=int(v<0)
  if min_slack is None or v<min_slack:min_slack=v;worst={'node_id':row['leaf_node_id'],'region_key':row.get('region_key'),'parent':row.get('parent'),'bernstein_index':row['bernstein_index']}
  r['exact_rows_checked']+=1
 r['exact_bad_count']=bad;r['exact_min_slack']=str(min_slack);r['worst_exact_row']=worst
 qf=np.array(list(map(float,q)));r['rational_float_cost']=float(f@qf)
 assert np.isfinite(qf).all() and np.isfinite(r['rational_float_cost'])
 r['point_checks']={}
 for k,kind,tol in [(0,'K',1e-9),(1,'A',1e-8)]:
  vals=allmat[k]@qf;assert np.isfinite(vals).all()
  r['point_checks'][kind]={'strengthened_bad':int(np.count_nonzero(vals-allrhs[k]<-tol)),'original_bad':int(np.count_nonzero(vals-targets[k]<-tol)),'tolerance':tol}
 r['point_bad_count']=sum(v['strengthened_bad']+v['original_bad'] for batch in (r['float_point_checks'],r['point_checks']) for v in batch.values())
 r['point_rows_checked']+=51574
 if bad==0 and r['point_bad_count']==0 and float(f@qf)<3.590901146403276-5e-7:
  guard()
  coeff={k:q[1]*Q(v) for k,v in baseline['a'].items()};T={}
  for t,d in zip(q[2:14],directions):
   dm=sum((h*h*mult(*map(int,k.split(',')))*Q(v)/2 for k,v in d['coefficients'].items()),Q(0))
   assert dm==Q(d['half_mass'])
   assert all(h*(sum(map(int,k.split(',')))+4)<=Q(3,4) for k,v in d['coefficients'].items() if Q(v))
   for k,v in d['coefficients'].items():T[k]=T.get(k,Q(0))+t*Q(v)
  for k,v in T.items():coeff[k]=coeff.get(k,Q(0))+v
  assert all(0<=int(k.split(',')[1])<=int(k.split(',')[0])<=baseline['N'] for k in coeff)
  offset_exact=h*h*sum((mult(*map(int,k.split(',')))*v for k,v in coeff.items()),Q(0))/2
  assert offset_exact==q[1]*smass+sum((q[2+j]*Q(d['half_mass']) for j,d in enumerate(directions)),Q(0))
  structured={k:baseline[k] for k in ('alpha','R','h','N')}
  structured.update(beta=baseline['alpha'],c=str(q[0]),a={k:str(v) for k,v in coeff.items() if v},
    structure={'baseline_sha256':'21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172',
      's':str(q[1]),'g':str(g),'T':{k:str(v) for k,v in T.items() if v},'T_support_bound':'3/4',
      'scope':'Outside T support, K=s*K1+g*Psi>=0; dominance and global generator still need full certification.'})
  sp=ROOT/'structured_kernel.json';sp.write_text(json.dumps(structured,indent=2)+'\n')
  import importlib.util
  cp=checked(ROOT.parent/'round10/enclose_cost.py','f5de232467e0060ffd3f0495b028cc6328372ef67c086a0a7693d5f38e762b67')
  spec=importlib.util.spec_from_file_location('round25_exact_cost',cp);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
  r['cost_enclosures_started']=1;save();exact=mod.enclose(sp,digits=60);r['cost_enclosures_completed']=1;guard()
  ap=checked(ROOT.parent/'round01/scaling_arithmetic.json','1f7b1a0184101a93f5d64e463377ee4e86374660fee3e3c2de7676ff46237159')
  baseline_lower=Q(read(ap)['new_cost_lower']);upper=Q(exact['upper'])
  exact.update(baseline_lower=str(baseline_lower),strict_gain_lower=str(baseline_lower-upper),
    improves_baseline=upper<baseline_lower,below_user_target=upper<Q(17,5),user_target='17/5',
    cost_source_sha256=sha256(cp.read_bytes()).hexdigest(),
    scope='Exact mass and structure only; no new general weak bound until full-domain checks and analytic bridge.')
  (ROOT/'exact_cost.json').write_text(json.dumps(exact,indent=2)+'\n')
  r['exact_cost_sha256']=sha256((ROOT/'exact_cost.json').read_bytes()).hexdigest()
  if upper<baseline_lower:
   candidate={k:structured[k] for k in ('alpha','R','h','N','c','a')}
   (ROOT/'candidate.json').write_text(json.dumps(candidate,indent=2)+'\n');r['candidate_exported']=True
   r['candidate_sha256']=sha256((ROOT/'candidate.json').read_bytes()).hexdigest()
 r['status']='one_finite_lp_and_exact_local_check_closed'
except Exception as e:r.update(status='failed_or_rejected',reason=f'{type(e).__name__}: {e}')
if (ROOT/'parameters.json').exists():r['parameters_sha256']=sha256((ROOT/'parameters.json').read_bytes()).hexdigest()
r['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();save();print(json.dumps(r,indent=2))

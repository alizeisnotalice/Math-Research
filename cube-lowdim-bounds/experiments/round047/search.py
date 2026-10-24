#!/usr/bin/env python3
"""One necessary-point LP discovery, then <=26 strictly enclosed support rows.

No mathematical work at import. Lead-owned wrapper supplies the 300s hard limit.
Unenclosed floating cache rows have zero weight in the final exact lower bound.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,MAX_COLUMNS,PRECISION=300.,392,192
COLUMNS=[2,3,4,5,198,199,201,202,222,244,267,268]
PINS={
 '../round01/scaled_certificate.json':'21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172',
 '../round14/directions.json':'98c0e2fcc7f5d448aef33acc24f2b7fe22977f7a2361b20dc7c2ade926a0fd6e',
 '../round25/solve_result.json':'50bed7523f577382aee952be871ad248f3efb38a46441d4a50d7596c3592b239',
 '../round25/exact_cost.json':'0120e7a55961ae7858541fb95319a1e4454c3c4e953c603478e51e488b7f8b33',
 '../round06/point_manifest.json':'513c8d953695d9ddd1c7a097ee9af2ab624f79b875330d02e0d42ad2dcc70ebd',
 '../upstream-v4/checker/independent/hl_core.py':'c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf',
}
POINTS=((Q(23497,65536),Q(10819,65536)),(Q(47105,131072),Q(21905,131072)))


def digest(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def dump(path,value):
 temp=path.with_suffix(path.suffix+'.tmp')
 temp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n');temp.replace(path)
def mult(key):
 i,j=map(int,key.split(','));return 1 if i==j==0 else 4 if j==0 or i==j else 8
def endpoint(value):
 if not value.is_finite() or not value.is_exact():raise ValueError('Nonexact/nonfinite endpoint.')
 m,e=value.man_exp();m,e=int(m),int(e)
 return Q(m*2**e) if e>=0 else Q(m,2**(-e))
def ball(value):
 if not value.is_finite():raise ValueError('Nonfinite coefficient.')
 lo,hi=endpoint(value.lower()),endpoint(value.upper())
 if lo>hi:raise ValueError('Reversed interval.')
 return {'lower':str(lo),'upper':str(hi),'ball':value.str(70)}
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
 sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
 if any((ROOT/n).exists() for n in ('result.json','calls.jsonl','model.local.mat','certificate.json','solver_result.json')):
  raise SystemExit('Existing evidence: no repeated optimization or coefficient evaluation.')
 start=monotonic();inputs={};proc=None
 counts={'columns_started':0,'columns_completed':0,'columns_failed':0,'columns_interrupted':0,
         'LP_requested':0,'LP_started':0,'LP_completed':0,'integral_calls':0,'box_model_calls':0}
 result={'status':'started','counts':counts,'limits':{'seconds':300,'column_values':392,'LP':1,'LP_soft_seconds':20,'bits':192},
         'call_unit':'One complete K or AK scalar for one unit basis column at one exact point',
         'source_sha256_before':digest(Path(__file__)),'in_flight':None,'candidate_exported':False}
 dump(ROOT/'result.json',result);trace=(ROOT/'calls.jsonl').open('x')
 def guard():
  if monotonic()-start>=SECONDS:raise TimeoutError('300s internal deadline.')
 def save():
  result['elapsed_seconds']=monotonic()-start;dump(ROOT/'result.json',result)
 def checked(relative,wanted):
  guard();path=ROOT/relative;actual=digest(path)
  if actual!=wanted:raise ValueError('Pinned input changed: '+relative)
  inputs[relative]=actual;return path
 def emit(value):trace.write(json.dumps(value,allow_nan=False)+'\n');trace.flush()
 try:
  if os.getpgrp()!=os.getpid():raise ValueError('Requires lead-owned isolated process group.')
  for rel,sha in PINS.items():checked(rel,sha)
  baseline=read(ROOT/'../round01/scaled_certificate.json');directions=read(ROOT/'../round14/directions.json')['directions']
  if [d['old_column'] for d in directions]!=COLUMNS:raise ValueError('Frozen direction order changed.')
  alpha,R,h=map(Q,(baseline['alpha'],baseline['R'],baseline['h']));c1=Q(baseline['c'])
  if (alpha,R,h,baseline['N'])!=(Q(47,40),Q(33,20),Q(1,32),80):raise ValueError('Wrong fixed family.')
  masses=[]
  for d in directions:
   mass=sum((h*h*mult(k)*Q(v)/2 for k,v in d['coefficients'].items()),Q(0))
   if mass!=Q(d['half_mass']):raise ValueError('Exact direction mass mismatch.')
   if any(Q(v) and (sum(map(int,k.split(',')))+4)*h>Q(3,4) for k,v in d['coefficients'].items()):
    raise ValueError('Direction support exceeds fixed 3/4.')
   masses.append(mass)
  smass=sum((h*h*mult(k)*Q(v)/2 for k,v in baseline['a'].items()),Q(0))
  cost=read(ROOT/'../round25/exact_cost.json');scale,root=int(cost['root_scale']),int(cost['root_floor'])
  if cost['root_degree']!=40 or cost['root_power']!=33 or scale!=10**60:raise ValueError('Wrong saved radical.')
  if not root**40*R.denominator**33<=R.numerator**33*scale**40<(root+1)**40*R.denominator**33:
   raise ValueError('Saved root bracket invalid.')
  f_lower=[alpha/(2-alpha)*Q(root,scale),smass]+masses+[Q(0)]*12
  f_upper=[alpha/(2-alpha)*Q(root+1,scale)]+f_lower[1:]
  bounds=[Q(d['bound']) for d in directions]
  lb=[c1/2,Q(1,2)]+[-v for v in bounds]+[Q(0)]*12
  ub=[3*c1/2,Q(3,2)]+bounds+bounds
  if not all(l<=u for l,u in zip(lb,ub)) or lb[0]<0:raise ValueError('Invalid exact variable box.')
  result['exact_family']={'alpha':str(alpha),'R':str(R),'h':str(h),'N':80,
   'c1':str(c1),'lower_bounds':list(map(str,lb)),'upper_bounds':list(map(str,ub)),
   'cost_lower_coefficients':list(map(str,f_lower)),'cost_upper_coefficients':list(map(str,f_upper))}
  save();guard()
  for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
  import numpy as np
  from scipy import sparse
  from scipy.io import savemat
  previous=read(ROOT/'../round25/solve_result.json')
  cache_pins={item['file']:item['sha256'] for item in previous['inputs']}
  point_manifest=read(ROOT/'../round06/point_manifest.json')
  matrices=[];descriptions=[];targets=[]
  for kind,n in [('K',26604),('A',24970)]:
   rel='round07/'+kind+'_directions_288.local.npy';checked('../'+rel,cache_pins[rel])
   matrix=np.load(ROOT/'..'/rel,mmap_mode='r')
   srel='round06/cache.local/'+kind+'_singular.npy';checked('../'+srel,cache_pins[srel])
   singular=np.load(ROOT/'..'/srel,mmap_mode='r')
   if matrix.shape!=(n,288) or singular.shape!=(n,) or len(point_manifest[kind])!=n:raise ValueError('Wrong cache dimensions.')
   selected=np.column_stack([singular,matrix[:,1]]+[matrix[:,d['old_column']] for d in directions]+[np.zeros(n)]*12)
   if not np.isfinite(selected).all():raise ValueError('Nonfinite frozen point matrix.')
   matrices.append(selected)
   for index,pair in enumerate(point_manifest[kind]):
    guard();u,v=map(Q,pair)
    target=Q(int(u+v<=1)) if kind=='K' else Q(0)
    targets.append(target);descriptions.append({'kind':kind,'u':str(u),'v':str(v),'target':str(target),'old_index':index})
  if len(descriptions)!=51574:raise ValueError('Wrong original necessary pool size.')
  core=load('round28_independent_core',ROOT/'../upstream-v4/checker/independent/hl_core.py')
  from flint import ctx
  import flint
  ctx.prec=PRECISION
  result['flint_version']=flint.__version__
  specs_path=ROOT/'basis.local';specs_path.mkdir(exist_ok=False)
  specs=[]
  for j in range(14):
   guard()
   # Spec requires nonempty acls. One baseline Spec supplies unscaled Psi
   # and S1 separately; no dummy coefficients are introduced.
   data=baseline if j<2 else {'alpha':baseline['alpha'],'R':baseline['R'],'h':baseline['h'],'N':80,
      'c':'0','a':directions[j-2]['coefficients']}
   path=specs_path/(str(j)+'.json');dump(path,data)
   specs.append(specs[0] if j==1 else core.Spec(path))
  result['basis_sha256']={str(j):digest(specs_path/(str(j)+'.json')) for j in range(14)}
  enclosed={}
  def strict_row(desc):
   kind,u,v=desc['kind'],Q(desc['u']),Q(desc['v']);key=(kind,u,v)
   if key in enclosed:return enclosed[key]
   if u<0 or v<0 or u+v<=0:raise ValueError('Point outside nonnegative off-origin chart.')
   if kind=='A':
    if not (u>0 and v>0 and u+v!=R):raise ValueError('Unsupported active AK point; no replacement support.')
    if u+v>R:
     # Each active J(x,y) has 0<y<R and x+y>R, hence x-(R-y)>0.
     for x,y in ((u,v),(v,u)):
      if y<R and not (0<y/(x+y)<R/(x+y)<1 and x-(R-y)>0):raise ValueError('Exterior J branch guard failed.')
   coefficients=[]
   for j,spec in enumerate(specs):
    guard()
    if counts['columns_started']>=MAX_COLUMNS:raise TimeoutError('392 coefficient call cap.')
    event={'call':counts['columns_started']+1,'kind':kind,'u':str(u),'v':str(v),'basis_column':j,'status':'started'}
    counts['columns_started']+=1;result['in_flight']=event.copy();emit(event);save()
    try:
     if j==0:value=spec.Psi(core.A(u+v)) if kind=='K' else spec.APsi_point(u,v)
     elif kind=='K':value=core.A(spec.S_point(u,v))
     else:
      one,two=spec.AS_point(u,v);value=one+two
     guard();packed=ball(value);coefficients.append(packed)
     counts['columns_completed']+=1;event.update(status='completed',**packed);emit(event)
     result['in_flight']=None;save()
    except TimeoutError:
     counts['columns_interrupted']+=1;event['status']='interrupted';emit(event);result['in_flight']=None;save();raise
    except Exception as exc:
     counts['columns_failed']+=1;event.update(status='failed',reason=str(exc));emit(event);result['in_flight']=None;save();raise
   record={**desc,'coefficients':coefficients}
   enclosed[key]=record;dump(ROOT/'point_enclosures.json',{'rows':list(enclosed.values())})
   return record
  new_matrix=[]
  for i,(u,v) in enumerate(POINTS):
   desc={'kind':'A','u':str(u),'v':str(v),'target':'0','new_point':'P'+str(i+1)}
   record=strict_row(desc)
   new_matrix.append([float((Q(b['lower'])+Q(b['upper']))/2) for b in record['coefficients']]+[0.]*12)
   descriptions.append(desc);targets.append(Q(0))
  points=np.vstack(matrices+[np.array(new_matrix)])
  absrows=np.zeros((24,26))
  for j in range(12):absrows[2*j,2+j]=1;absrows[2*j+1,2+j]=-1;absrows[2*j:2*j+2,14+j]=-1
  structure=np.zeros((1,26));structure[0,0]=-1;structure[0,1]=float(c1)
  A=np.vstack([-points,absrows,structure]);b=np.r_[-np.array(list(map(float,targets))),np.zeros(25)]
  if A.shape!=(51601,26) or not np.isfinite(A).all() or not np.isfinite(b).all():raise ValueError('Wrong necessary model.')
  savemat(ROOT/'model.local.mat',{'A':sparse.csr_matrix(A),'b':b[:,None],
     'f':np.array(list(map(float,f_lower)))[:,None],'lb':np.array(list(map(float,lb)))[:,None],
     'ub':np.array(list(map(float,ub)))[:,None]})
  dump(ROOT/'model_metadata.json',{'point_rows':51576,'structural_rows':25,'rows':51601,'variables':26,
      'necessary_targets':'Original 1_D for K, zero for AK; no halfslack or Bernstein rows.',
      'model_sha256':digest(ROOT/'model.local.mat'),'row_descriptions':descriptions})
  model_sha=digest(ROOT/'model.local.mat');solver_sha=digest(ROOT/'solver.m')
  result.update(model_sha256=model_sha,solver_source_sha256=solver_sha,point_row_count=51576);save();guard()
  counts['LP_requested']=1;save()
  with (ROOT/'matlab.log').open('x') as log:
   script=str(ROOT/'solver.m').replace("'","''")
   proc=subprocess.Popen(['/Applications/MATLAB_R2026a.app/bin/matlab','-batch',"run('"+script+"')"],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
   result['matlab_pid']=proc.pid;result['matlab_pgid']=os.getpgid(proc.pid);result['coordinator_pgid']=os.getpgrp()
   if result['matlab_pgid']!=os.getpgrp():raise ValueError('Batch escaped owned process group.')
   save();code=proc.wait(timeout=max(.01,SECONDS-(monotonic()-start)))
  guard();sol=read(ROOT/'solver_result.json');counts['LP_started']=sol.get('calls_started',0);counts['LP_completed']=sol.get('calls_completed',0)
  result.update(matlab_exit_code=code,solver_result_sha256=digest(ROOT/'solver_result.json'),matlab_log_sha256=digest(ROOT/'matlab.log'),
      solver_status=sol['status'],exitflag=sol.get('exitflag'),float_discrete_objective=sol.get('objective'))
  if digest(ROOT/'model.local.mat')!=model_sha or digest(ROOT/'solver.m')!=solver_sha or sol['input_sha256']!=model_sha:
   raise ValueError('Model/solver provenance changed.')
  if code!=0 or sol['status']!='returned' or sol['exitflag']!=1:raise ValueError('Single LP not optimal; no retry and no infeasibility conclusion.')
  raw=sol['lambda']['ineqlin']
  if len(raw)!=51601 or not all(math.isfinite(float(x)) for x in raw):raise ValueError('Missing/nonfinite inequality dual.')
  quantized=[max(Q(0),Q(round(Q(float(v))*2**40),2**40)) for v in raw]
  support=sorted((i for i,y in enumerate(quantized) if y>0),key=lambda i:(-quantized[i],i))[:26]
  support.sort();result['support_row_indices']=support;save()
  if not support:raise ValueError('No positive quantized dual support; no second support selection.')
  M=[max(abs(l),abs(u)) for l,u in zip(lb,ub)];certificate_rows=[]
  for index in support:
   guard();y=quantized[index]
   if index<51576:
    desc=descriptions[index];record=strict_row(desc)
    alo=[-Q(z['upper']) for z in record['coefficients']]+[Q(0)]*12
    ahi=[-Q(z['lower']) for z in record['coefficients']]+[Q(0)]*12
    target=Q(desc['target']);rhs=-target
    origin={'type':'necessary_point','point':desc,'enclosure_key':[desc['kind'],desc['u'],desc['v']]}
   else:
    row=index-51576;alo=[Q(0)]*26
    if row<24:
     j=row//2;alo[2+j]=Q(1 if row%2==0 else -1);alo[14+j]=Q(-1)
    elif row==24:alo[0]=Q(-1);alo[1]=c1
    else:raise ValueError('Unexpected structural row.')
    ahi=alo.copy();rhs=Q(0);origin={'type':'exact_structure','structural_row':row}
   mid=[(l+u)/2 for l,u in zip(alo,ahi)];radius=[(u-l)/2 for l,u in zip(alo,ahi)]
   error=sum((d*m for d,m in zip(radius,M)),Q(0))
   if any(d<0 for d in radius):raise ValueError('Negative interval radius.')
   certificate_rows.append({'index':index,'multiplier':str(y),'origin':origin,
       'A_lower':list(map(str,alo)),'A_upper':list(map(str,ahi)),
       'A_midpoint':list(map(str,mid)),'A_radius':list(map(str,radius)),
       'b':str(rhs),'box_error':str(error),'relaxed_b':str(rhs+error)})
  residual=f_lower.copy();weighted_rhs=Q(0)
  for row in certificate_rows:
   guard();y=Q(row['multiplier']);weighted_rhs+=y*Q(row['relaxed_b'])
   residual=[r+y*Q(a) for r,a in zip(residual,row['A_midpoint'])]
  box_terms=[min(r*l,r*u) for r,l,u in zip(residual,lb,ub)]
  lower=-weighted_rhs+sum(box_terms,Q(0))
  certificate={'status':'exact_sparse_necessary_lower_bound','rows':certificate_rows,'maximum_support':26,
      'lower_bounds':list(map(str,lb)),'upper_bounds':list(map(str,ub)),
      'cost_lower_coefficients':list(map(str,f_lower)),'cost_upper_coefficients':list(map(str,f_upper)),
      'residual':list(map(str,residual)),'box_terms':list(map(str,box_terms)),
      'weighted_relaxed_rhs':str(weighted_rhs),'lower_bound':str(lower),'target':'17/5','excludes_cost_at_most_target':lower>Q(17,5),
      'root_certificate':{k:cost[k] for k in ('root_degree','root_power','root_scale','root_floor','R')},
      'baseline_sha256':PINS['../round01/scaled_certificate.json'],'directions_sha256':PINS['../round14/directions.json'],
      'model_sha256':model_sha,'solver_result_sha256':digest(ROOT/'solver_result.json'),
      'scope':'Only fixed alpha/R, 12 T directions, exact parameter box and c>=s*c1. Not all kernels or true C2; no new upper bound.'}
  dump(ROOT/'certificate.json',certificate)
  result.update(status='one_lp_and_sparse_exact_lower_bound_closed',certificate_sha256=digest(ROOT/'certificate.json'),
      exact_lower_bound=str(lower),exact_lower_bound_decimal=float(lower),excludes_cost_at_most_target=lower>Q(17,5),
      point_enclosures_sha256=digest(ROOT/'point_enclosures.json'),
      qualification='Only recomputed support necessary rows enter exact lower bound. Remaining floating rows have zero multipliers.')
 except Exception as exc:
  result.update(status='failed_or_inconclusive',reason=f'{type(exc).__name__}: {exc}')
 finally:
  if proc is not None and proc.poll() is None:
   proc.terminate()
   try:proc.wait(timeout=.5)
   except subprocess.TimeoutExpired:proc.kill();proc.wait()
  trace.close()
  after={rel:digest(ROOT/rel) for rel in inputs if (ROOT/rel).is_file()}
  result.update(input_sha256_before=inputs,input_sha256_after=after,inputs_unchanged=after==inputs,
      source_sha256_after=digest(Path(__file__)),trace_sha256=digest(ROOT/'calls.jsonl'))
  if after!=inputs or result['source_sha256_before']!=result['source_sha256_after']:
   result['status']='provenance_changed_invalid'
  if (ROOT/'solver_result.json').exists():
   saved=read(ROOT/'solver_result.json');counts['LP_started']=saved.get('calls_started',0);counts['LP_completed']=saved.get('calls_completed',0)
  save()
 assert counts['columns_started']<=MAX_COLUMNS
 assert counts['columns_started']==counts['columns_completed']+counts['columns_failed']+counts['columns_interrupted']
 print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds']}),flush=True)
 return 0 if result['status']=='one_lp_and_sparse_exact_lower_bound_closed' else 1


if __name__=='__main__':raise SystemExit(main())

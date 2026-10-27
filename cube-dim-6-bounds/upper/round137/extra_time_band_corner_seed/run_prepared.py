#!/usr/bin/env python3
"""UNRUN C137:one boundary point,five disjoint positive windows,max15 attempts."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from fractions import Fraction as Q
import concurrent.futures as cf,hashlib,importlib.util,json,resource,subprocess,threading,time,traceback
P=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(p):return json.loads(Path(p).read_text())
def req(x,s):
 if not x:raise RuntimeError(s)
def write(path,d,exclusive=False):
 path=Path(path)
 if exclusive:
  with path.open('x')as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 else:
  t=path.with_name(path.name+'.own-tmp')
  with t.open('w')as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
  os.replace(t,path)
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def strints(data):return[[str(a),str(b)]for a,b in[data[0],*data[1],data[2],data[3],data[4]]]
def overlap(a,b):return all(max(Q(x[0]),y[0])<=min(Q(x[1]),y[1])for x,y in zip(a,[b[0],*b[1],b[2],b[3],b[4]]))

def main():
 start=time.monotonic();cpu=time.process_time();inp=rd(P/'inputs.json');reg=rd(P/'registration.json');rel=rd(P/'parent_release.json');qual=rd(P/'cache_qualification.json')
 req(rel.get('authorized')is True and str(rel.get('round'))=='137','primary137 release')
 req(rel['script_sha256']==sha(__file__)and rel['registration_sha256']==sha(P/'registration.json'),'released source/registration')
 req(rel.get('exclusive_owner_for_registered_semantic_keys')is True,'main must reserve all5 semantic keys against other native dispatch')
 lanes=int(rel['parallel_windows']);req(1<=lanes<=5,'registered max5 parallel windows')
 for n,h in reg['prepared_pins'].items():req(sha(P/n)==h,'prepared pin '+n)
 for row in rd(P/'source_manifest.json')['sources']:req(sha(row['path'])==row['sha256'],'source '+row['role'])
 req(os.getpriority(os.PRIO_PROCESS,0)>=10,'nice10')
 req(len(qual['windows'])==5 and all(r['action']=='NEW'for r in qual['windows']),'this frozen batch has5 genuinely missing queries, no old reissues')
 U=load('c137_original_family_helpers',inp['native_family_helpers_path']);V=load('c137_original_native_validator',inp['native_family_validator_path']);base=load('c137_original_base',inp['base_path'])
 base._c016_validate_jet_512=base.validate_jet
 family=rd(inp['native_family_registration_path']);pins=reg['native_pins'];G=rd(inp['fullG135_results_path']);GL,GU=map(Q,G['G_interval'])
 req(G['status']=='COMPLETE_ONE_FULL_GREEN_POINT'and G['full_positive_time_G']is True and G['point']==inp['center'],'same original complete Green point')
 enc=base.encode_center(inp['center'],max_ray_abs=5,max_rate=160)
 req(enc['denominator']==8 and enc['numerators']==[8,8,8,3,2,1],'exact boundary point')
 write(P/'started.json',{'pid':os.getpid(),'release_sha256':sha(P/'parent_release.json'),'unix_time':time.time(),'semantic_keys':[r['semantic_key']for r in qual['windows']]},True)
 stop=threading.Event();lock=threading.RLock();state={'status':'RUNNING','windows':{},'actual_attempts_started':0,'active_native':0,'max_active_native':0};nativeCPU=0.;accepted={};failures=[]
 (P/'windows').mkdir(exist_ok=False)
 def publish():write(P/'progress.json',state)
 def worker(row):
  nonlocal nativeCPU
  wi=row['window_index'];w=row['window'];key=row['semantic_key'];e={k:inp[k]for k in ['model','model_definition','center','quadratic']};e['window']=w
  req(U.semantic_key(e,reg['fan_sha256'])==key,'physical semantic key')
  wd=P/'windows'/str(wi);wd.mkdir(exist_ok=False);prior=[];levels=[];owned=None
  with lock:state['windows'][str(wi)]={'semantic_key':key,'window':w,'status':'READY','levels':levels};publish()
  try:
   for bits in (128,256,512):
    with lock:
     if stop.is_set():state['windows'][str(wi)]['status']='STOPPED_BEFORE_NEW_ATTEMPT';publish();return
     req(state['actual_attempts_started']<15,'global15 actual attempt cap')
    ad=wd/str(bits);ad.mkdir(exist_ok=False);b=family['binaries'][str(bits)];cmd=[b['path'],'full',inp['native_fan_path'],str(Q(w['a']).numerator),str(Q(w['a']).denominator),str(Q(w['T']).numerator),str(Q(w['T']).denominator),'8',*map(str,enc['numerators'])]
    apins={**pins,b['path']:b['sha256'],b['source_path']:b['source_sha256']}
    intent={'schema':'C137-single-window-intent','semantic_key':key,'expected':e,'precision_bits':bits,'command':cmd,'formula_family_id':U.FORMULA_FAMILY_ID,'pins':apins,'release_sha256':sha(P/'parent_release.json'),'time_or_CPU_ceiling':None}
    write(ad/'intent.json',intent,True);t0=time.monotonic()
    with(ad/'native.json').open('xb')as out,(ad/'stderr.txt').open('xb')as err:
     with lock:
      if stop.is_set():state['windows'][str(wi)]['status']='STOPPED_BEFORE_LAUNCH';publish();return
      proc=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=out,stderr=err,cwd=ad);owned=proc
      state['actual_attempts_started']+=1;state['active_native']+=1;state['max_active_native']=max(state['max_active_native'],state['active_native'])
      level={'precision_bits':bits,'status':'NATIVE_RUNNING','pid':proc.pid,'token':'C137.'+key+'.p'+str(bits),'attempt_directory':str(ad)};levels.append(level)
      state['windows'][str(wi)]['status']='NATIVE_RUNNING';publish();write(ad/'original_process.json',{'pid':proc.pid,'wrapper_pid':os.getpid(),'command':cmd,'intent_sha256':sha(ad/'intent.json')},True)
     # wait4 is scoped to this exact PID; each real child CPU is recorded separately.
     pid,status,usage=os.wait4(proc.pid,0);req(pid==proc.pid,'reap exact original native PID');code=os.waitstatus_to_exitcode(status);proc.returncode=code;owned=None
    childcpu=usage.ru_utime+usage.ru_stime
    with lock:nativeCPU+=childcpu;state['active_native']-=1;level.update(status='REAPED',exit_code=code);publish()
    raw=(ad/'native.json').read_bytes();err=(ad/'stderr.txt').read_bytes()
    receipt={'semantic_key':key,'precision_bits':bits,'native_pid':proc.pid,'wrapper_pid':os.getpid(),'exit_code':code,'wall_seconds':time.monotonic()-t0,'native_CPU_seconds':childcpu,'child_reaped':True,'command':cmd,'pins':apins,'stdout_sha256':sha(ad/'native.json'),'stderr_sha256':sha(ad/'stderr.txt'),'intent_sha256':sha(ad/'intent.json'),'process_sha256':sha(ad/'original_process.json'),'formula_family_id':U.FORMULA_FAMILY_ID,'one_window_evaluations':1,'time_or_CPU_ceiling':None}
    write(ad/'receipt.json',receipt,True)
    with lock:level.update(receipt_path=str(ad/'receipt.json'));publish()
    if U.valid_precision_guard(code,err,raw):
     with lock:level['status']='CLASSIFIED_DISTINCT_NODE_PRECISION_GUARD';publish()
     continue
    req(code==0 and err==b'','fatal native error/stderr, no retry')
    data=V.validate_native(U.strict_json_bytes(raw),e,bits);req(all(overlap(p,data)for p in prior),'disjoint trusted same-semantic precision intervals')
    ok,width=V.width_acceptable(data);prior.append(strints(data));receipt.update(width_acceptable=ok,width_report=width,validated_intervals=strints(data));write(ad/'receipt.json',receipt)
    with lock:level.update(status='VALID_COMPLETE_WIDTH_PASS'if ok else'VALID_COMPLETE_WIDTH_FAIL',width_report=width,validated_intervals=strints(data));publish()
    if not ok:continue
    prov={'kind':'C016_NATIVE_SAME_WINDOW_POSITIVE_BAND','precision_bits':bits,'full_group_terms':967680,'computed_sign_pairs':483840,'synthetic_fixture':False,'native_sha256':sha(ad/'native.json'),'receipt_sha256':sha(ad/'receipt.json'),'gradient_gate_sha256':inp['gradient_gate_sha256'],'formula_family_id':U.FORMULA_FAMILY_ID,'binary_sha256':b['sha256'],'source_sha256':b['source_sha256'],'independent_full_gradient_crosschecked_for_this_new_attempt':False}
    jet={**e,'F_interval':list(map(str,data[0])),'face_gradient_intervals':[list(map(str,g))for g in data[1][1:]],'source_proof':'C016 native same-window all-real gradient; exact original attempt '+str(ad/'receipt.json'),'provenance':prov}
    U.validate_jet_family(base,jet,e,pinned_binary_sha=b['sha256'],pinned_source_sha=b['source_sha256']);write(ad/'accepted_jet.json',jet,True)
    record={'status':'AVAILABLE','semantic_key':key,'window':w,'accepted_precision_bits':bits,'attempt_directory':str(ad),'raw_path':str(ad/'native.json'),'raw_sha256':sha(ad/'native.json'),'receipt_path':str(ad/'receipt.json'),'receipt_sha256':sha(ad/'receipt.json'),'ten_intervals':strints(data),'jet_path':str(ad/'accepted_jet.json'),'jet_sha256':sha(ad/'accepted_jet.json')}
    write(wd/'accepted.json',record,True)
    with lock:
     accepted[wi]=(data,record);state['windows'][str(wi)].update(status='AVAILABLE',accepted_precision_bits=bits)
     for left in (wi-1,wi):
      if left in accepted and left+1 in accepted:
       pa=accepted[left][0][2];pt=accepted[left+1][0][3];req(max(pa[0],pt[0])<=min(pa[1],pt[1]),'fatal adjacent endpoint density mismatch in partial batch')
     partial=sum((v[0][0][0]for v in accepted.values()),Q(0));req(partial<=GU,'CONTRADICTION: disjoint positive-band lower sum exceeds certified full G upper');publish()
    return
   with lock:state['windows'][str(wi)]['status']='FINITE_PRECISION_LADDER_EXHAUSTED';publish()
  except BaseException:
   stop.set();error=traceback.format_exc();emergency=None
   if owned is not None:
    pid,st,ru=os.wait4(owned.pid,0);owned.returncode=os.waitstatus_to_exitcode(st)
    emergency={'native_pid':pid,'child_reaped':True,'exit_code':owned.returncode,'native_CPU_seconds':ru.ru_utime+ru.ru_stime}
    with lock:nativeCPU+=ru.ru_utime+ru.ru_stime;state['active_native']-=1
   write(wd/'failure.json',{'semantic_key':key,'error':error,'emergency_original_child_receipt':emergency,'automatic_retry':False})
   with lock:failures.append({'window_index':wi,'error':error});state['windows'][str(wi)]['status']='FATAL_NO_RETRY';publish()
 def combined_prefix(n):
  items=[accepted[j][0]for j in range(n)];F=tuple(sum(d[0][k]for d in items)for k in(0,1));g=[tuple(sum(d[1][j][k]for d in items)for k in(0,1))for j in range(6)];inter=[]
  for j in range(1,n):
   # Previous p_a and next p_T refer to the same shared exact time.
   lo=max(items[j-1][2][0],items[j][3][0]);hi=min(items[j-1][2][1],items[j][3][1]);req(lo<=hi,'fatal adjacent endpoint density mismatch');inter.append([str(lo),str(hi)])
  a=Q(inp['windows'][n-1]['a']);T=Q(inp['windows'][0]['T']);terms=[(5*F[0],5*F[1]),(T*items[0][3][0],T*items[0][3][1]),(-a*items[-1][2][1],-a*items[-1][2][0])]
  terms +=[tuple(sorted(Q(x)*v for v in it))for x,it in zip(inp['center'],g)];eu=tuple(sum(it[k]for it in terms)for k in(0,1));req(eu[0]<=0<=eu[1],'merged exact Euler')
  req(F[0]<=GU,'positive prefix/fullG consistency')
  return{'length':n,'outer_window':{'kind':'B','a':str(a),'T':str(T)},'F_interval':list(map(str,F)),'six_gradients':[list(map(str,z))for z in g],'exact_Euler_interval':list(map(str,eu)),'adjacent_density_intersections':inter,'lower_at_least_g0':F[0]>=Q(inp['comparison_g0']),'lower_at_least_gstar':F[0]>=Q(inp['target']),'fullG_upper_consistency':True,'sum_upper_at_most_fullG_lower_strong_check':F[1]<=GL}
 failure=None;status='PARTIAL_FINITE_BATCH'
 try:
  with cf.ThreadPoolExecutor(max_workers=lanes)as pool:
   futures=[pool.submit(worker,r)for r in qual['windows']]
   for f in cf.as_completed(futures):f.result()
  prefixes=[]
  for n in range(1,6):
   if all(j in accepted for j in range(n)):prefixes.append(combined_prefix(n))
   else:break
  status='COMPLETE_FIVE_WINDOW_POINT_SEED'if len(accepted)==5 and not failures else'PARTIAL_FINITE_BATCH'
  write(P/'results.json',{'status':status,'point':inp['center'],'accepted_windows':[accepted[j][1]for j in sorted(accepted)],'prefixes':prefixes,'fullG135_interval':list(map(str,(GL,GU))),'actual_native_attempts':state['actual_attempts_started'],'native_CPU_seconds':nativeCPU,'max_active_native':state['max_active_native'],'failures':failures,'new_point_window_semantics':5,'precision_ladder':[128,256,512],'scope':'same single point only; no box or global floor claimed','no_short_time_tail_needed_for_lower':True,'global_g_proved':False,'new_C6_bound':False})
  req(not failures,'batch has fatal window failure, no restart')
 except BaseException as ex:
  failure=repr(ex);status='FAILED_NO_RETRY';raise
 finally:
  write(P/'run_receipt.json',{'status':status,'wrapper_pid':os.getpid(),'wall_before_receipt':time.monotonic()-start,'controller_CPU_seconds':time.process_time()-cpu,'native_CPU_seconds':nativeCPU,'actual_attempts_started':state['actual_attempts_started'],'active_native_at_exit':state['active_native'],'max_active_native':state['max_active_native'],'parallel_windows':lanes,'priority':os.getpriority(os.PRIO_PROCESS,0),'time_or_CPU_ceiling':None,'automatic_restarts':0,'failure':failure,'fatal_windows':len(failures)})
if __name__=='__main__':main()

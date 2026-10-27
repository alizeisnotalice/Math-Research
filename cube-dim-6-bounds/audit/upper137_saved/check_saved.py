#!/usr/bin/env python3
"""Independent C137 saved-only exact-rational audit; zero native/Arb calls."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,itertools,json,os,sys,time,traceback
sys.dont_write_bytecode=True
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round137/extra_time_band_corner_seed')
A=Path(__file__).resolve().parent
G=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round137/primary_gate/terminal_handoff.json')
PRE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-cover090-audit/cube-dim-6-bounds/audit/upper137_preflight')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def write(n,d):
 with (A/n).open('x')as f:json.dump(d,f,indent=2);f.write('\n')
def load(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def strings(v):return [list(map(str,x))for x in v]
def independent_dyadic(d):
 assert set(d)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
 a,b,e=[int(d[k])for k in ['lower_mantissa','upper_mantissa','common_binary_exponent']]
 assert all(str(int(v))==v for v in d.values());assert a<=b
 factor=F(2)**e;return a*factor,b*factor
reg=read(P/'registration.json');inp=read(P/'inputs.json');r=read(P/'results.json');receipt=read(P/'run_receipt.json');rel=read(P/'parent_release.json');qual=read(P/'cache_qualification.json')
start=time.monotonic();cpu=time.process_time()
files=[P/n for n in ['registration.json','inputs.json','source_manifest.json','cache_qualification.json','run_prepared.py','observe_once.py','parent_release.json','started.json','progress.json','results.json','run_receipt.json','observed_tool_receipt.json']]+[G]
for i in range(5):files += [P/'windows'/str(i)/n for n in ['accepted.json','128/native.json','128/stderr.txt','128/intent.json','128/original_process.json','128/receipt.json','128/accepted_jet.json']]
pins={str(p):sha(p)for p in files}
write('finite_saved_registration.json',{'status':'REGISTERED_SAVED_ONLY_EXACT_RATIONAL','payload_hashes':pins,'windows':5,'raw_native_files':5,'ten_intervals_per_window':10,'physical_gradients_per_window':6,'prefixes':5,'partial_subsets':32,'new_native_calls':0,'new_Arb_calls':0,'new_exp_or_log_calls':0,'time_limit':None,'reviewer_script_sha256':sha(__file__)})
try:
 assert sha(P/'run_prepared.py')==reg['script_sha256']==rel['script_sha256']=='dc85a6eb06bcaeb9b937a79b08d8c5ad0eb32e7bcca7500a162d31b559da83c2'
 assert sha(P/'registration.json')==rel['registration_sha256'] and sha(P/'observe_once.py')==reg['observer_sha256']==rel['observer_sha256']
 assert rel['authorized'] and rel['exclusive_owner_for_registered_semantic_keys'] and rel['parallel_windows']==5
 assert sha(PRE/'static_review.json')==rel['independent_static_review_sha256']
 for n,h in reg['prepared_pins'].items():assert sha(P/n)==h
 sources=read(P/'source_manifest.json')['sources']
 for v in sources:assert sha(v['path'])==v['sha256'] and Path(v['path']).stat().st_size==v['bytes']
 V=load('c137_independent_original_validator',inp['native_family_validator_path']);U=load('c137_independent_original_helpers',inp['native_family_helpers_path']);B=load('c137_independent_original_base',inp['base_path']);B._c016_validate_jet_512=B.validate_jet
 family=read(inp['native_family_registration_path']);binary=family['binaries']['128']
 assert r['status']==receipt['status']=='COMPLETE_FIVE_WINDOW_POINT_SEED' and r['point']==inp['center']==['1','1','1','3/8','1/4','1/8']
 assert r['actual_native_attempts']==receipt['actual_attempts_started']==5 and r['max_active_native']==receipt['max_active_native']==5 and receipt['active_native_at_exit']==0
 assert receipt['automatic_restarts']==0 and receipt['failure']is None and receipt['fatal_windows']==0 and receipt['priority']>=10 and receipt['time_or_CPU_ceiling']is None
 assert not r['global_g_proved'] and not r['new_C6_bound'] and not r['failures']
 observed=read(P/'observed_tool_receipt.json');handoff=read(G)
 assert handoff['receipt_sha256']==sha(P/'observed_tool_receipt.json') and handoff['actual_receipt']==observed
 assert observed['exit_code']==0 and observed['child_reaped'] and observed['automatic_restarts']==0 and observed['time_or_CPU_ceiling']is None and observed['release_sha256']==sha(P/'parent_release.json')
 assert sha(P/'science_stdout.txt')==observed['stdout_sha256'] and sha(P/'science_stderr.txt')==observed['stderr_sha256']
 assert read(P/'started.json')['release_sha256']==sha(P/'parent_release.json')
 assert len(r['accepted_windows'])==len(qual['windows'])==5 and len({v['semantic_key']for v in qual['windows']})==5 and all(v['action']=='NEW'for v in qual['windows'])
 datas=[];window_audits=[];pids=[];childcpu=0.;progress=read(P/'progress.json')
 for i,row in enumerate(r['accepted_windows']):
  wd=P/'windows'/str(i);ad=wd/'128';expected={k:inp[k]for k in ['model','model_definition','center','quadratic']};expected['window']=inp['windows'][i]
  assert row==read(wd/'accepted.json') and row['window']==expected['window']==qual['windows'][i]['window']
  semantic=U.semantic_key(expected,reg['fan_sha256']);assert semantic==row['semantic_key']==qual['windows'][i]['semantic_key']
  assert row['accepted_precision_bits']==128 and row['status']=='AVAILABLE' and sorted(p.name for p in wd.iterdir())==['128','accepted.json']
  for f,h in [('raw_path','raw_sha256'),('receipt_path','receipt_sha256'),('jet_path','jet_sha256')]:assert sha(row[f])==row[h]
  raw=U.strict_json_bytes((ad/'native.json').read_bytes());d=V.validate_native(raw,expected,128)
  own=[independent_dyadic(raw['band']),*[independent_dyadic(z)for z in raw['gradient']],independent_dyadic(raw['p_a']),independent_dyadic(raw['p_T']),independent_dyadic(raw['Euler_xgrad_plus5B_plusTpT_minusapa'])]
  assert own==[d[0],*d[1],d[2],d[3],d[4]] and raw['full_group_terms']==967680 and raw['computed_sign_pairs']==483840 and len(d[1])==6
  assert own[0][0]>0 and own[7][0]>0 and own[8][0]>0 and own[9][0]<=0<=own[9][1]
  own_width=[]
  for j,(lo,hi)in enumerate(own):
   limit=F(1,65000*2**40)*max(F(1),abs(lo),abs(hi));assert hi-lo<=limit;own_width.append({'index':j,'width':str(hi-lo),'limit':str(limit)})
  ok,width=V.width_acceptable(d);assert ok and width=={'widths':own_width,'failed_indices':[]}
  a,T=F(expected['window']['a']),F(expected['window']['T']);assert 0<a<T
  terms=[(5*d[0][0],5*d[0][1]),(T*d[3][0],T*d[3][1]),(-a*d[2][1],-a*d[2][0])]+[tuple(sorted(F(c)*v for v in g))for c,g in zip(inp['center'],d[1])]
  eu=tuple(sum(t[k]for t in terms)for k in (0,1));assert eu==d[5] and eu[0]<=0<=eu[1]
  rec=read(ad/'receipt.json');intent=read(ad/'intent.json');proc=read(ad/'original_process.json');jet=read(row['jet_path'])
  assert rec['validated_intervals']==row['ten_intervals']==strings(own) and rec['width_acceptable'] and rec['width_report']==width
  assert rec['exit_code']==0 and rec['child_reaped'] and rec['precision_bits']==128 and rec['semantic_key']==semantic and rec['one_window_evaluations']==1 and rec['time_or_CPU_ceiling']is None
  assert sha(ad/'native.json')==rec['stdout_sha256'] and sha(ad/'stderr.txt')==rec['stderr_sha256'] and (ad/'stderr.txt').read_bytes()==b''
  assert sha(ad/'intent.json')==rec['intent_sha256']==proc['intent_sha256'] and sha(ad/'original_process.json')==rec['process_sha256']
  assert rec['native_pid']==proc['pid'] and rec['wrapper_pid']==proc['wrapper_pid']==receipt['wrapper_pid'];pids.append(rec['native_pid']);childcpu+=rec['native_CPU_seconds']
  cmd=[binary['path'],'full',inp['native_fan_path'],str(a.numerator),str(a.denominator),str(T.numerator),str(T.denominator),'8','8','8','8','3','2','1']
  assert intent['expected']==expected and intent['command']==rec['command']==proc['command']==cmd and intent['precision_bits']==128 and intent['semantic_key']==semantic and intent['release_sha256']==sha(P/'parent_release.json')
  for path,h in rec['pins'].items():assert sha(path)==h
  assert jet['F_interval']==list(map(str,d[0])) and jet['face_gradient_intervals']==strings(d[1][1:]) and jet['quadratic']==expected['quadratic']
  assert jet['provenance']['native_sha256']==sha(ad/'native.json') and jet['provenance']['receipt_sha256']==sha(ad/'receipt.json') and jet['provenance']['gradient_gate_sha256']==inp['gradient_gate_sha256']
  assert not jet['provenance']['independent_full_gradient_crosschecked_for_this_new_attempt']
  U.validate_jet_family(B,jet,expected,pinned_binary_sha=binary['sha256'],pinned_source_sha=binary['source_sha256'])
  ps=progress['windows'][str(i)];assert ps['status']=='AVAILABLE' and ps['accepted_precision_bits']==128 and len(ps['levels'])==1 and ps['levels'][0]['status']=='VALID_COMPLETE_WIDTH_PASS' and ps['levels'][0]['pid']==rec['native_pid']
  datas.append(d);window_audits.append({'index':i,'semantic_key':semantic,'window':expected['window'],'accepted_bits':128,'ten_interval_identity':True,'six_gradient_identity':True,'all_ten_widths_pass':True,'exact_Euler_interval':list(map(str,eu)),'native_Euler_interval':list(map(str,d[4])),'F_interval':list(map(str,d[0]))})
 assert len(set(pids))==5 and abs(childcpu-receipt['native_CPU_seconds'])<1e-12
 G135=read(inp['fullG135_results_path']);GL,GU=map(F,G135['G_interval']);assert G135['status']=='COMPLETE_ONE_FULL_GREEN_POINT' and G135['full_positive_time_G'] and G135['point']==inp['center'] and r['fullG135_interval']==G135['G_interval']
 intersections=[]
 for i in range(4):
  assert F(inp['windows'][i]['a'])==F(inp['windows'][i+1]['T']);inter=(max(datas[i][2][0],datas[i+1][3][0]),min(datas[i][2][1],datas[i+1][3][1]));assert inter[0]<=inter[1];intersections.append(list(map(str,inter)))
 # Every partial accepted subset satisfies its available adjacency and fullG lower-sum gate.
 for mask in range(32):
  ids=[i for i in range(5)if mask>>i&1];assert sum((datas[i][0][0]for i in ids),F(0))<=GU
  for i in range(4):
   if i in ids and i+1 in ids:assert max(datas[i][2][0],datas[i+1][3][0])<=min(datas[i][2][1],datas[i+1][3][1])
 prefix_audits=[];assert len(r['prefixes'])==5
 for n,saved in enumerate(r['prefixes'],1):
  items=datas[:n];value=tuple(sum(d[0][k]for d in items)for k in (0,1));grad=[tuple(sum(d[1][j][k]for d in items)for k in (0,1))for j in range(6)]
  a=F(inp['windows'][n-1]['a']);T=F(inp['windows'][0]['T']);terms=[(5*value[0],5*value[1]),(T*items[0][3][0],T*items[0][3][1]),(-a*items[-1][2][1],-a*items[-1][2][0])]+[tuple(sorted(F(c)*v for v in g))for c,g in zip(inp['center'],grad)]
  eu=tuple(sum(t[k]for t in terms)for k in (0,1));assert eu[0]<=0<=eu[1]
  want={'length':n,'outer_window':{'kind':'B','a':str(a),'T':str(T)},'F_interval':list(map(str,value)),'six_gradients':strings(grad),'exact_Euler_interval':list(map(str,eu)),'adjacent_density_intersections':intersections[:n-1],'lower_at_least_g0':value[0]>=F(inp['comparison_g0']),'lower_at_least_gstar':value[0]>=F(inp['target']),'fullG_upper_consistency':value[0]<=GU,'sum_upper_at_most_fullG_lower_strong_check':value[1]<=GL}
  assert saved==want and value[1]<GL and value[0]>F(inp['target']) and value[0]>F(inp['comparison_g0'])
  prefix_audits.append({**want,'lower_minus_gstar_exact':str(value[0]-F(inp['target'])),'fullG_lower_minus_sum_upper_exact':str(GL-value[1])})
 SL,SU=map(F,r['prefixes'][-1]['F_interval']);remainder=(GL-SU,GU-SL);capture=(SL/GU,SU/GL);assert 0<remainder[0]<=remainder[1] and 0<capture[0]<=capture[1]<1
 for p,h in pins.items():assert sha(p)==h
 out={'verdict':'PASS_C137_SAVED_EXACT_RATIONAL_REVIEW','source_sha256':sha(P/'run_prepared.py'),'saved_result_sha256':sha(P/'results.json'),'reviewer_registration_sha256':sha(A/'finite_saved_registration.json'),'source_manifest_hashes_verified':len(sources),'new_windows':5,'accepted_first_attempt_bits':[128]*5,'raw_ten_field_intervals_checked':50,'physical_gradient_intervals_checked':30,'all_partial_subsets_checked':32,'window_audits':window_audits,'adjacent_density_intersections':intersections,'prefix_audits':prefix_audits,'fullG135_interval':list(map(str,(GL,GU))),'positive_omitted_time_remainder_interval':list(map(str,remainder)),'captured_fraction_interval':list(map(str,capture)),'new_native_calls':0,'new_Arb_calls':0,'new_exp_or_log_calls':0,'single_point_only':True,'global_g_certificate':False,'new_C6_bound':False,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
 write('verification_result.json',out);write('exit_receipt.json',{'status':out['verdict'],'exit_code':0,'result_sha256':sha(A/'verification_result.json'),'wall_seconds':out['wall_seconds'],'CPU_seconds':out['CPU_seconds'],'new_native_calls':0,'new_Arb_calls':0,'new_exp_or_log_calls':0,'time_limit':None});print(json.dumps({k:out[k]for k in ['verdict','raw_ten_field_intervals_checked','physical_gradient_intervals_checked','all_partial_subsets_checked','new_native_calls','new_Arb_calls','wall_seconds','CPU_seconds']}))
except BaseException:
 write('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'new_native_calls':0,'new_Arb_calls':0});raise

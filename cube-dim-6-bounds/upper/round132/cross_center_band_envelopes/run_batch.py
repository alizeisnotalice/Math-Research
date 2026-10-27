#!/usr/bin/env python3
"""Gated local cross-center positive-band envelopes, exact rational + Arb512."""
import sys,os,json,hashlib,datetime,importlib.util,itertools,traceback,time
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def write_new(p,x):
 with Path(p).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
def atomic(p,x):
 p=Path(p);t=p.with_name(p.name+'.pending')
 with t.open('wb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
 os.replace(t,p);fd=os.open(str(p.parent),os.O_RDONLY)
 try:os.fsync(fd)
 finally:os.close(fd)
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def load(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def release_requirements():
 return {'round':'132','authorized':True,'scope':'LOCAL_CROSS_CENTER_BAND_ENVELOPES_ONLY','precision_bits':512,'workers':1,'nice_minimum':10,'time_limit':None,'native_calls_permitted':0,'max_corner_exp':327680,'run_batch_sha256':sha(P/'run_batch.py'),'input_hashes_sha256':sha(P/'input_hashes.json'),'selected_boxes_sha256':sha(P/'selected_boxes_and_centers.json'),'selected_source_entries_sha256':sha(P/'selected_source_cache_entries.json'),'conditional_M6_sha256':sha(P/'conditional_M6_check.json'),'proof_sha256':sha(P/'PROOF.md'),'prepared_packet_sha256':sha(P/'prepared_packet.json')}
def release_gate():
 rp=P/'parent_release.json';require(rp.is_file(),'MAIN_AND_INDEPENDENT_REVIEW_RELEASE_REQUIRED')
 r=read(rp)
 for k,v in release_requirements().items():require(r.get(k)==v,'RELEASE_MISMATCH:'+k)
 require(r.get('main_approved')is True and r.get('independent_review_approved')is True,'BOTH_GATES_REQUIRED')
 return rp
def authenticate_sources():
 hashes=read(P/'input_hashes.json')
 for group in ['frozen_sources','selected_original_artifacts']:
  for p,h in hashes[group].items():require(sha(p)==h,'INPUT_HASH:'+p)
 cfg=read(P/'selected_boxes_and_centers.json');entries=read(P/'selected_source_cache_entries.json')['entries']
 cube=P.parents[2];src=cube/'upper/round127/continuation_v121';reg=read(src/'registration_127.json');inp=read(reg['inputs']['117_inputs_path'])
 U=load('c132_saved_helpers',src/'precision_family_helpers.py');V=load('c132_saved_validator',src/'native_validator_family.py');B=load('c132_original_base',inp['base_path']);B._c016_validate_jet_512=B.validate_jet
 jets={}
 for k,m in entries.items():
  e=m['expected'];require(U.semantic_key(e,reg['runtime_and_fan_pins'][inp['native_paths']['fan']])==k,'SEMANTIC_IDENTITY')
  # Authenticate the original expected center/window, never a forged target-box expected input.
  if'jet_path'in m:jet=read(m['jet_path'])
  else:
   jet=read(m['jet_ref']['path'])
   for sel in m['jet_ref'].get('selector',[]):jet=jet[sel]
  bits=m['precision_bits'];U.validate_jet_family(B,jet,e,allowed_bits=(128,256,512),pinned_binary_sha=reg['binaries'][str(bits)]['sha256'],pinned_source_sha=reg['binaries'][str(bits)]['source_sha256'])
  require(jet['provenance']['native_sha256']==m['native_sha256'],'SOURCE_JET_RAW_HASH')
  raw=V.validate_native(read(m['native_path']),e,bits);require(jet['F_interval']==list(map(str,raw[0]))and jet['face_gradient_intervals']==[list(map(str,g))for g in raw[1][1:]],'SOURCE_JET_RAW_VALUES')
  nr=read(m['receipt_path']);require(nr['exit_code']==0,'SOURCE_RECEIPT_EXIT')
  require(list(map(Q,jet['center']))[0]==1,'CENTER_FACE')
  jets[k]=jet
 return cfg,entries,jets

def certify_center(box,jet,window,budget,*,arb,ctx):
 """Independent implementation of A057 (12) allowing a true center outside the new box."""
 require(jet['window']==window,'SAME_BAND_REQUIRED');center=list(map(Q,jet['center']));bounds=[tuple(map(Q,b))for b in box];lo,hi=map(Q,jet['F_interval']);grad=[tuple(map(Q,b))for b in jet['face_gradient_intervals']]
 require(lo>0 and len(center)==6 and center[0]==1 and len(bounds)==len(grad)==5,'TRUE_JET_SCHEMA')
 require(all(a<=b for a,b in bounds),'BOX_ORDER');a=Q(window['a']);require(a>0 and a<Q(window['T']),'POSITIVE_BAND')
 def ball(q):return arb(q.numerator)/arb(q.denominator)
 ctx.prec=512;rows=[];lows=[];highs=[]
 for bits in itertools.product((0,1),repeat=5):
  d=[bounds[i][bit]-center[i+1]for i,bit in enumerate(bits)];b=sum(min(gl*di,gu*di)for(gl,gu),di in zip(grad,d));f=min(hi,max(lo,b));q=sum(x*x for x in d)/32;exponent=b/f-q/a**2
  require(budget['used']<budget['cap'],'EXP_BUDGET_EXHAUSTED');budget['used']+=1
  enclosure=ball(f)*ball(exponent).exp();require(enclosure.is_finite(),'NONFINITE_ARB_ENCLOSURE');lower=Q(str(enclosure.lower().fmpq()));upper=Q(str(enclosure.upper().fmpq()));lows.append(lower);highs.append(upper)
  rows.append({'corner_bits':list(bits),'displacement':list(map(str,d)),'gradient_dot_lower':str(b),'clamped_true_value':str(f),'frame32_quadratic':str(q),'exponent':str(exponent),'enclosure_lower':str(lower),'enclosure_upper':str(upper)})
 return {'true_center':jet['center'],'window':window,'precision_bits':512,'complete_box_lower':str(min(lows)),'robust_tangent_min_enclosure_upper':str(min(highs)),'source_proof':jet['source_proof'],'corners':rows,'scope':'one reliable lower bound on this complete new box; no center-in-box assumption'}
def run():
 rp=release_gate() # first, before any scientific work or start marker
 require(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART');os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 cfg,entries,jets=authenticate_sources()
 import flint
 from flint import arb,ctx
 require(flint.__FLINT_VERSION__=='3.6.0','PINNED_FLINT_VERSION')
 require(cfg['selected_box_count']<=128 and len(cfg['windows'])==5,'FINITE_SELECTION')
 require(cfg['actual_planned_corner_exp']<=cfg['max_corner_exp']==327680,'FINITE_BUDGET')
 require(all(Q(cfg['windows'][j+1]['T'])==Q(cfg['windows'][j]['a'])for j in range(4)),'ADJACENT_DISJOINT_BANDS')
 write_new(P/'started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'release_sha256':sha(rp),'precision_bits':512,'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'native_calls':0})
 budget={'used':0,'cap':cfg['actual_planned_corner_exp']};results=[];(P/'box_results').mkdir(exist_ok=False);t0=time.monotonic()
 try:
  for item in cfg['selected_boxes']:
   lower=Q(0);bandrows=[]
   for band in item['bands']:
    centerrows=[]
    for ref in band['candidate_keys']:
     k=ref['semantic_key'];jet=jets[k];require(jet['center']==ref['true_center'],'TRUE_CENTER_PRESERVED')
     cert=certify_center(item['box'],jet,band['window'],budget,arb=arb,ctx=ctx);centerrows.append({'semantic_key':k,'source_artifact_key':entries[k]['artifact_key'],**cert})
    best=max((Q(r['complete_box_lower'])for r in centerrows),default=Q(0));lower+=best
    bandrows.append({'window':band['window'],'same_band_aggregation':'max complete-box lower, never sum centers','band_complete_box_lower':str(best),'maximizing_keys':[r['semantic_key']for r in centerrows if Q(r['complete_box_lower'])==best],'center_bounds':centerrows,'missing_band_zero':not centerrows})
   result={'node_id':item['node_id'],'box':item['box'],'band_envelopes':bandrows,'summed_disjoint_band_lower':str(lower),'passes_new_local_target':lower>=Q(cfg['target']),'passes_reference_local_target':lower>=Q(cfg['reference_target']),'target':cfg['target'],'reference_target':cfg['reference_target'],'global_g_certificate':False,'new_C6_bound':False,'scope':'only this originally unresolved box'}
   write_new(P/'box_results'/(item['node_id']+'.json'),result);results.append({k:v for k,v in result.items()if k!='band_envelopes'})
   atomic(P/'progress.json',{'completed_boxes':len(results),'corner_exp_used':budget['used'],'native_calls':0,'global_g_certificate':False})
  require(budget['used']==budget['cap'],'CORNER_COUNT')
  write_new(P/'result.json',{'status':'LOCAL_BATCH_COMPLETE','root_id':cfg['root_id'],'source_unresolved_leaf_count':cfg['source_unresolved_leaf_count'],'selected_boxes':len(results),'new_target_local_passes':sum(r['passes_new_local_target']for r in results),'reference_target_local_passes':sum(r['passes_reference_local_target']for r in results),'all_selected_pass_new_target':all(r['passes_new_local_target']for r in results),'root_new_target_closed':False,'reason_no_root_or_global_claim':'Existing passed leaves were only certified at1/65000 and were not rerun at1/63500; selection remains a local experiment.','corner_exp_used':budget['used'],'native_calls':0,'wall_seconds':time.monotonic()-t0,'boxes':results,'global_g_certificate':False,'new_C6_bound':False})
 except BaseException:
  atomic(P/'exit_receipt.json',{'status':'FAILED_LOCAL_PARTIAL','traceback':traceback.format_exc(),'corner_exp_used':budget['used'],'completed_boxes':len(results),'native_calls':0,'global_g_certificate':False,'new_C6_bound':False});raise
 else:atomic(P/'exit_receipt.json',{'status':'LOCAL_BATCH_COMPLETE','corner_exp_used':budget['used'],'native_calls':0,'global_g_certificate':False,'new_C6_bound':False})
if __name__=='__main__':
 if sys.argv[1:]==['--source-check-only']:
  cfg,entries,jets=authenticate_sources();print(json.dumps({'status':'ORIGINAL_SOURCE_CHECK_ONLY_PASS','jets':len(jets),'native_calls':0,'exp_calls':0}))
 elif not sys.argv[1:]:run()
 else:raise SystemExit('usage: run_batch.py [--source-check-only]')

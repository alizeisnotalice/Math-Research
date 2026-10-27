#!/usr/bin/env python3
"""Dual-gated C141: exact Fraction odd-Taylor lower probes, zero exp/native/checker calls."""
import json,hashlib,sys,os,itertools,datetime,time,traceback
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
ODD=(1,3,5,7,9,11,13,15)
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def atomic(p,x):
 p=Path(p);t=p.with_name(p.name+'.pending')
 with t.open('wb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def write_new(p,x):
 with Path(p).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
def release_requirements():
 return {'round':'141','authorized':True,'scope':'128_FIXED_TERMINAL_UNRESOLVED_BOXES_WITH_FIVE_EXISTING_R137_JETS_ONLY','native_calls_permitted':0,'exp_calls_permitted':0,'checker_calls_permitted':0,'subdivisions_permitted':0,'workers':1,'nice_minimum':10,'time_or_CPU_ceiling':None,'odd_degrees':list(ODD),'max_box_band_corner_cases':20480,'max_Taylor_recurrence_updates':307200,'g0':'1/65000','gstar':'1/63500','run_batch_sha256':sha(P/'run_batch.py'),'selection_sha256':sha(P/'selected_boxes_and_existing_jets.json'),'input_hashes_sha256':sha(P/'input_hashes.json'),'proof_sha256':sha(P/'PROOF.md'),'prepared_packet_sha256':sha(P/'prepared_packet.json')}
def release_gate():
 rp=P/'parent_release.json';require(rp.is_file(),'MAIN_AND_INDEPENDENT_REVIEW_RELEASE_REQUIRED');r=read(rp)
 for k,v in release_requirements().items():require(r.get(k)==v,'RELEASE_MISMATCH:'+k)
 require(r.get('main_approved')is True and r.get('independent_review_approved')is True,'DUAL_GATE_REQUIRED');return rp

def source_preflight():
 for path,h in read(P/'input_hashes.json')['frozen_inputs'].items():require(sha(path)==h,'FROZEN_INPUT_HASH:'+path)
 cfg=read(P/'selected_boxes_and_existing_jets.json');require(len(cfg['boxes'])==128 and cfg['odd_degrees']==list(ODD),'FIXED_SELECTION_AND_ODD_SET');require(len(cfg['jets'])==5,'FIVE_JETS_ONLY');require(len({(b['root_id'],b['node_id'])for b in cfg['boxes']})==128,'UNIQUE_ORIGINAL_BOXES');jets=[]
 for box in cfg['boxes']:require(len(box['box'])==5 and box['status']in('queued','blocked')and all(Q(lo)<=Q(hi)for lo,hi in box['box']),'ORIGINAL_UNRESOLVED_RECTANGLE')
 for meta in cfg['jets']:
  jet=read(meta['jet_path']);intent=read(meta['intent_path']);e=intent['expected'];nr=read(meta['receipt_path']);require(jet['center']==e['center']==cfg['anchor']and jet['window']==e['window']==meta['window'],'ORIGINAL_TRUE_CENTER_WINDOW');require(jet['model']==e['model']and jet['model_definition']==e['model_definition'],'ORIGINAL_MODEL');require(intent['semantic_key']==meta['semantic_key']==nr['semantic_key'],'ORIGINAL_SEMANTIC_IDENTITY');require(jet['provenance']['native_sha256']==meta['raw_sha256']and nr['stdout_sha256']==meta['raw_sha256'],'ORIGINAL_RAW_LINK');require(nr['exit_code']==0 and nr['width_acceptable'],'ORIGINAL_ACCEPTED_JET');require(len(jet['face_gradient_intervals'])==5 and Q(jet['F_interval'][0])>0,'POSITIVE_VALUE_FIVE_GRADIENTS');require(Q(jet['F_interval'][0])<=Q(jet['F_interval'][1])and all(Q(lo)<=Q(hi)for lo,hi in jet['face_gradient_intervals']),'SOURCE_INTERVAL_ORDER');require(Q(meta['window']['a'])>0 and Q(meta['window']['T'])>Q(meta['window']['a']),'POSITIVE_OWN_BAND');jets.append(jet)
 require(all(Q(cfg['bands'][j+1]['T'])==Q(cfg['bands'][j]['a'])for j in range(4)),'DISJOINT_ADJACENT_TIME_BANDS')
 return cfg,jets

def odd_polynomial_lower(r):
 """P_odd(r)<=exp(r) for every real r; this function computes no exponential."""
 term=Q(1);partial=Q(1);best=Q(0);degree=0
 for k in range(1,16):
  term=term*r/k;partial+=term
  if k in ODD and partial>best:best=partial;degree=k
 return best,degree

def robust_corner_lower(jet,point):
 c=list(map(Q,jet['center']));lo,hi=map(Q,jet['F_interval']);require(0<lo<=hi,'POSITIVE_CENTER_INTERVAL');grad=[tuple(map(Q,g))for g in jet['face_gradient_intervals']];d=[Q(x)-y for x,y in zip(point[1:],c[1:])];b=sum(min(gl*di,gu*di)for(gl,gu),di in zip(grad,d));f=min(hi,max(lo,b));q=sum(x*x for x in d)/32;a=Q(jet['window']['a']);r=b/f-q/a**2;poly,degree=odd_polynomial_lower(r);lower=f*poly
 return {'true_center':jet['center'],'window':jet['window'],'point':point,'displacement':list(map(str,d)),'gradient_dot_lower':str(b),'clamped_true_value':str(f),'frame32_quadratic':str(q),'own_time_a':str(a),'exponent_exact':str(r),'chosen_odd_degree':degree,'nonnegative_exp_lower_polynomial_exact':str(poly),'band_corner_lower_exact':str(lower),'strictly_positive_lower':lower>0,'zero_lower_is_not_strict_positivity':lower==0}

def run():
 rp=release_gate();require(not(P/'started.json').exists(),'SOLE_BATCH_NO_AUTOMATIC_RESTART');os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));cfg,jets=source_preflight();write_new(P/'started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'release_sha256':sha(rp),'nice':os.getpriority(os.PRIO_PROCESS,0),'native_calls':0,'exp_calls':0,'checker_calls':0});t0=time.monotonic();results=[];cache={};cases=0;updates=0;zero_band_count=0
 try:
  for choice in cfg['boxes']:
   perband=[];total=Q(0);prefix=[]
   for wi,(meta,jet)in enumerate(zip(cfg['jets'],jets)):
    refs=[];bounds=[]
    for bits in itertools.product((0,1),repeat=5):
     point=['1']+[choice['box'][i][bit]for i,bit in enumerate(bits)];ck=hashlib.sha256(canon({'semantic_key':meta['semantic_key'],'point':point})).hexdigest();cases+=1;require(cases<=20480,'CORNER_BAND_CASE_CAP')
     if ck not in cache:cache[ck]={'semantic_key':meta['semantic_key'],**robust_corner_lower(jet,point)};updates+=15;require(updates<=307200,'FINITE_POLYNOMIAL_CAP')
     refs.append({'corner_bits':list(bits),'cache_key':ck});bounds.append(Q(cache[ck]['band_corner_lower_exact']))
    lower=min(bounds);zero_band_count+=lower==0;total+=lower;perband.append({'window':meta['window'],'semantic_key':meta['semantic_key'],'whole_box_rational_lower_exact':str(lower),'strictly_positive_lower':lower>0,'corner_references':refs});prefix.append({'consecutive_band_count':wi+1,'sum_of_bandwise_box_lowers_exact':str(total),'passes_g0':total>=Q(cfg['target_g0']),'passes_gstar':total>=Q(cfg['target_gstar'])})
   g0=total>=Q(cfg['target_g0']);gs=total>=Q(cfg['target_gstar']);require(not gs or g0,'TARGET_ORDER');classification='GSTAR_AND_G0_LOCAL_BOX_PASS'if gs else'G0_ONLY_LOCAL_BOX_PASS'if g0 else'LOWER_CERTIFICATE_UNRESOLVED_BOTH';results.append({'root_id':choice['root_id'],'node_id':choice['node_id'],'box':choice['box'],'sum_of_bandwise_min32_lowers_exact':str(total),'passes_g0':g0,'passes_gstar':gs,'classification':classification,'per_band':perband,'prefix_box_lower_comparisons':prefix,'global_g0_complete':False,'global_gstar_complete':False,'true_G_upper_bound':False,'scope':'this originally unresolved terminal box only; not inherited gstar coverage of prior g0-paid leaves'})
  require(cases==20480,'REGISTERED_FINITE_BATCH_SIZE');g0=sum(r['passes_g0']for r in results);gs=sum(r['passes_gstar']for r in results);output={'status':'SOLE_PURE_RATIONAL_BATCH_COMPLETE','selected_boxes':128,'terminal_residual_leaf_universe':cfg['source_residual_leaf_count'],'selected_root_count':len(set(r['root_id']for r in results)),'g0':'1/65000','gstar':'1/63500','g0_local_box_pass_count':g0,'gstar_local_box_pass_count':gs,'g0_only_count':g0-gs,'unresolved_both_count':128-g0,'point_upper_obstruction_count':0,'true_g_counterexample_count':0,'actual_G_upper_bound_count':0,'reason_no_counterexample_claim':'Only lower certificates are computed; lower certificate failure is unresolved and cannot exclude actual G or even the true robust family.','corner_band_case_count':cases,'unique_rational_source_point_evaluations':len(cache),'Taylor_recurrence_updates':updates,'zero_band_lower_count':zero_band_count,'finite_odd_degrees':list(ODD),'new_native_calls':0,'new_exp_calls':0,'checker_calls':0,'new_subdivisions':0,'wall_seconds':time.monotonic()-t0,'boxes':results,'global_g_certificate':False,'global_gstar_certificate':False,'new_C6_bound':False,'old_g0_paid_leaves_not_promoted_to_gstar':True}
  write_new(P/'rational_point_cache.json',{'entries':cache,'unique_entries':len(cache),'native_calls':0,'exp_calls':0});write_new(P/'result.json',output);write_new(P/'exit_receipt.json',{'status':output['status'],'exit_code':0,'result_sha256':sha(P/'result.json'),'rational_cache_sha256':sha(P/'rational_point_cache.json'),'corner_band_cases':cases,'Taylor_updates':updates,'native_calls':0,'exp_calls':0,'checker_calls':0});print(json.dumps({k:v for k,v in output.items()if k!='boxes'},sort_keys=True))
 except BaseException:
  atomic(P/'failure_receipt.json',{'status':'FAILED_STOP_NO_RETRY','traceback':traceback.format_exc(),'completed_boxes':len(results),'cases':cases,'Taylor_updates':updates,'native_calls':0,'exp_calls':0,'checker_calls':0});raise
if __name__=='__main__':
 if sys.argv[1:]==['--source-check-only']:
  cfg,jets=source_preflight();print(json.dumps({'status':'FROZEN_SOURCE_HASH_AND_METADATA_ONLY_PASS','selected_boxes':len(cfg['boxes']),'jets':len(jets),'corner_lower_evaluations':0,'native_calls':0,'exp_calls':0,'checker_calls':0}))
 elif not sys.argv[1:]:run()
 else:raise SystemExit('usage: run_batch.py [--source-check-only]')

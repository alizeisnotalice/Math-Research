#!/usr/bin/env python3
"""The sole registered C133 batch: existing C132 enclosures, Fraction/max/sum only."""
import json,hashlib,os,sys,itertools,datetime,traceback,time
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def write_new(n,x):
 with(P/n).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
def require(ok,msg):
 if not ok:raise ValueError(msg)
def main():
 require(not(P/'started.json').exists(),'SOLE_BATCH_ALREADY_STARTED_NO_RERUN');reg=read(P/'registration.json');require(sha(P/'diagnose.py')==reg['program_sha256'],'PROGRAM_HASH')
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 for row in reg['inputs']:require(sha(row['path'])==row['sha256'],'INPUT_HASH:'+row['path'])
 write_new('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json'),'nice':os.getpriority(os.PRIO_PROCESS,0),'native_calls':0,'exp_calls':0})
 t0=time.monotonic();cfg=read(reg['selected_config_path']);summary=read(reg['source_result_path']);bitslist=list(itertools.product((0,1),repeat=5));results=[];corner_count=0
 try:
  require(len(cfg['selected_boxes'])==93 and len(cfg['windows'])==5,'FIXED_93_5')
  for choice in cfg['selected_boxes']:
   node=choice['node_id'];src=read(Path(reg['box_result_directory'])/(node+'.json'));require(src['node_id']==node and src['box']==choice['box'],'BOX_IDENTITY');require(len(src['band_envelopes'])==5,'BAND_COUNT')
   totals={bits:[Q(0),Q(0),[]]for bits in bitslist}
   for bi,(band,selected)in enumerate(zip(src['band_envelopes'],choice['bands'])):
    require(band['window']==selected['window']==cfg['windows'][bi],'WINDOW_IDENTITY');centers=band['center_bounds'];keys=[x['semantic_key']for x in centers];require(keys==[r['semantic_key']for r in selected['candidate_keys']]and len(keys)==len(set(keys)),'FIXED_CANDIDATE_KEYS')
    maps=[]
    for center,ref in zip(centers,selected['candidate_keys']):
     require(center['true_center']==ref['true_center']and center['window']==band['window'],'UNCHANGED_TRUE_CENTER');cm={tuple(c['corner_bits']):c for c in center['corners']};require(set(cm)==set(bitslist)and len(center['corners'])==32,'FIXED32_CORNERS');maps.append(cm)
     for bits,row in cm.items():
      lower,upper=Q(row['enclosure_lower']),Q(row['enclosure_upper']);require(0<lower<=upper,'POSITIVE_ORDERED_SAVED_ENCLOSURE')
      exactd=[Q(choice['box'][i][bit])-Q(center['true_center'][i+1])for i,bit in enumerate(bits)];require(list(map(str,exactd))==row['displacement'],'CORNER_SOURCE_BINDING')
    require(bool(centers)!=band['missing_band_zero'],'MISSING_BAND_ZERO')
    for bits in bitslist:
     lower=max((Q(cm[bits]['enclosure_lower'])for cm in maps),default=Q(0));upper=max((Q(cm[bits]['enclosure_upper'])for cm in maps),default=Q(0));require(lower<=upper,'MAX_ENCLOSURE');totals[bits][0]+=lower;totals[bits][1]+=upper;totals[bits][2].append({'window':band['window'],'max_center_robust_lower':str(lower),'max_center_robust_upper':str(upper),'lower_max_keys':[center['semantic_key']for center,cm in zip(centers,maps)if Q(cm[bits]['enclosure_lower'])==lower],'upper_max_keys':[center['semantic_key']for center,cm in zip(centers,maps)if Q(cm[bits]['enclosure_upper'])==upper],'missing_band_zero':not centers})
   corners=[]
   for bits in bitslist:
    lo,hi,bands=totals[bits];corners.append({'corner_bits':list(bits),'face_point':[choice['box'][i][bit]for i,bit in enumerate(bits)],'fixed_robust_family_lower':str(lo),'fixed_robust_family_upper':str(hi),'band_envelopes':bands});corner_count+=1
   classes={}
   for target in reg['targets']:
    T=Q(target);failed=[c for c in corners if Q(c['fixed_robust_family_upper'])<T];alllower=all(Q(c['fixed_robust_family_lower'])>=T for c in corners);prior=Q(src['summed_disjoint_band_lower'])>=T
    require(not prior or alllower,'WHOLEBOX_PASS_MUST_PAY_CORNERS')
    worst=min(corners,key=lambda c:(Q(c['fixed_robust_family_upper']),c['corner_bits']));lowest=min(corners,key=lambda c:(Q(c['fixed_robust_family_lower']),c['corner_bits']))
    classification='FIXED_CERTIFICATE_FAMILY_OBSTRUCTED_AT_CORNER'if failed else'WHOLEBOX_ALREADY_PASSED'if prior else'ALL_CORNERS_PAID_LOCAL_SUBDIVISION_CANDIDATE_ONLY'if alllower else'INTERVAL_UNRESOLVED'
    classes[target]={'classification':classification,'original_wholebox_pass':prior,'obstruction_corner_count':len(failed),'all32_corner_lowers_ge_target':alllower,'interval_unresolved_corner_count':sum(Q(c['fixed_robust_family_lower'])<T<=Q(c['fixed_robust_family_upper'])for c in corners),'worst_upper_corner':{'corner_bits':worst['corner_bits'],'face_point':worst['face_point'],'upper_exact':worst['fixed_robust_family_upper'],'upper_minus_target_exact':str(Q(worst['fixed_robust_family_upper'])-T),'target_minus_upper_exact':str(T-Q(worst['fixed_robust_family_upper']))},'lowest_lower_corner':{'corner_bits':lowest['corner_bits'],'lower_exact':lowest['fixed_robust_family_lower'],'lower_minus_target_exact':str(Q(lowest['fixed_robust_family_lower'])-T)}}
   results.append({'node_id':node,'box':choice['box'],'targets':classes,'corners':corners,'scope':'only fixed selected robust certificate family; never a bound from above on G or a true-g counterexample'})
  require(corner_count==93*32,'FIXED_CORNER_BUDGET');aggregates={}
  for target in reg['targets']:
   counts={name:sum(r['targets'][target]['classification']==name for r in results)for name in ['FIXED_CERTIFICATE_FAMILY_OBSTRUCTED_AT_CORNER','WHOLEBOX_ALREADY_PASSED','ALL_CORNERS_PAID_LOCAL_SUBDIVISION_CANDIDATE_ONLY','INTERVAL_UNRESOLVED']};worst=min(results,key=lambda r:(Q(r['targets'][target]['worst_upper_corner']['upper_exact']),r['node_id']));aggregates[target]={'counts':counts,'obstruction_corners_total':sum(r['targets'][target]['obstruction_corner_count']for r in results),'local_subdivision_candidate_node_ids':[r['node_id']for r in results if r['targets'][target]['classification']=='ALL_CORNERS_PAID_LOCAL_SUBDIVISION_CANDIDATE_ONLY'],'global_worst_saved_family_upper_corner':{'node_id':worst['node_id'],**worst['targets'][target]['worst_upper_corner']}}
  output={'status':'SOLE_FIXED_BATCH_COMPLETE','boxes':93,'corners':corner_count,'targets':aggregates,'results':results,'native_calls':0,'exp_calls':0,'checker_calls':0,'new_subdivisions':0,'source_result_sha256':sha(reg['source_result_path']),'registration_sha256':sha(P/'registration.json'),'wall_seconds':time.monotonic()-t0,'global_g_certificate':False,'new_C6_bound':False,'true_g_counterexample':False,'scope':'obstructions exclude only same candidate jets plus same robust formula, even after axis box subdivisions containing the exhibited corner; paid corners alone do not prove box interior coverage'}
  write_new('diagnostic.json',output);write_new('exit_receipt.json',{'status':'SOLE_FIXED_BATCH_COMPLETE','exit_code':0,'boxes':93,'corners':corner_count,'native_calls':0,'exp_calls':0,'checker_calls':0,'diagnostic_sha256':sha(P/'diagnostic.json')});print(json.dumps({'status':output['status'],'boxes':93,'corners':corner_count,'targets':aggregates,'native_calls':0,'exp_calls':0},sort_keys=True))
 except BaseException:
  write_new('failure_receipt.json',{'status':'FAILED_STOP_NO_RETRY','traceback':traceback.format_exc(),'completed_boxes':len(results),'native_calls':0,'exp_calls':0});raise
if __name__=='__main__':main()

#!/usr/bin/env python3
"""One registered independent saved141 audit: Fraction/Horner only, no producer imports."""
import sys
sys.set_int_max_str_digits(0)
sys.dont_write_bytecode=True
import hashlib,itertools,json,math,os,resource,time,traceback
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).resolve().parent
ODD=(1,3,5,7,9,11,13,15)
checks=0

def check(x,msg):
 global checks
 checks+=1
 if not x:raise ArithmeticError(msg)
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(v):return(json.dumps(v,sort_keys=True,separators=(',',':'))+'\n').encode()
def save(n,d):
 p=R/n;t=p.with_name(p.name+'.own-tmp')
 with t.open('w')as h:json.dump(d,h,indent=2);h.write('\n');h.flush();os.fsync(h.fileno())
 os.replace(t,p)
def interval(v):
 check(isinstance(v,list)and len(v)==2,'two endpoints');lo,hi=map(F,v);check(lo<=hi,'interval order');return lo,hi
def dyadic(v):
 check(set(v)=={'lower_mantissa','upper_mantissa','common_binary_exponent'},'raw dyadic schema')
 m=[int(v[n])for n in('lower_mantissa','upper_mantissa','common_binary_exponent')]
 check(all(str(z)==v[n]for z,n in zip(m,('lower_mantissa','upper_mantissa','common_binary_exponent'))),'canonical raw integers')
 check(m[0]<=m[1],'raw order');return(F(m[0])*F(2)**m[2],F(m[1])*F(2)**m[2])
def horner(r,n):
 # Independent of the producer's shared Taylor recurrence.
 value=F(1,math.factorial(n))
 for k in range(n-1,-1,-1):value=value*r+F(1,math.factorial(k))
 return value

def main():
 start=time.monotonic();c0=time.process_time();reg=rd(R/'registration.json');sources=rd(R/'source_manifest.json')['sources'];spec=rd(R/'inputs.json')
 check(reg['checker_sha256']==sha(__file__),'registered independent source')
 check(reg['source_manifest_sha256']==sha(R/'source_manifest.json')and reg['inputs_sha256']==sha(R/'inputs.json'),'registered fixed inputs')
 for row in sources:check(sha(row['path'])==row['sha256'],'pre source '+row['path'])
 check(not(R/'started.json').exists(),'sole audit no restart')
 n=os.getpriority(os.PRIO_PROCESS,0)
 if n<10:os.nice(10-n)
 with(R/'started.json').open('x')as h:json.dump({'pid':os.getpid(),'unix_time':time.time(),'checker_sha256':sha(__file__),'registration_sha256':sha(R/'registration.json')},h);h.write('\n');h.flush();os.fsync(h.fileno())
 status='FAILED_STOP_NO_RETRY';failure=None;pointcount=boxcount=0
 try:
  cfg=rd(spec['selection_path']);author=rd(spec['result_path']);saved=rd(spec['cache_path']);terminal=rd(spec['main_terminal_receipt_path']);ex=rd(spec['author_exit_receipt_path']);forest=rd(spec['terminal_forest_path'])
  check(terminal['exit_code']==0 and terminal['restarts']==0,'actual primary sole exit0')
  check(terminal['source_sha256']==spec['author_script_sha256']and terminal['release_sha256']==sha(spec['author_release_path']),'primary actual source/release')
  check(ex['exit_code']==0 and ex['result_sha256']==sha(spec['result_path'])and ex['rational_cache_sha256']==sha(spec['cache_path']),'actual author output links')
  check(len(cfg['boxes'])==128 and len(cfg['jets'])==5 and cfg['odd_degrees']==list(ODD),'registered finite data size')
  check(cfg['anchor']==['1','1','1','3/8','1/4','1/8'],'actual boundary anchor')
  check(cfg['bands']==[v['window']for v in cfg['jets']],'exact five window list')
  for j in range(4):check(F(cfg['bands'][j]['a'])==F(cfg['bands'][j+1]['T']),'disjoint adjacent bands')
  check(author['status']=='SOLE_PURE_RATIONAL_BATCH_COMPLETE','author terminal')
  check(len(saved['entries'])==saved['unique_entries']==5970,'exact saved rational point universe')
  sources_by_key={};coeffs={k:F(1,math.factorial(k))for k in range(16)}
  for row in cfg['jets']:
   jet=rd(row['jet_path']);intent=rd(row['intent_path']);receipt=rd(row['receipt_path']);raw=rd(row['raw_path']);key=row['semantic_key']
   for p,h in[(row['jet_path'],row['jet_sha256']),(row['intent_path'],row['intent_sha256']),(row['receipt_path'],row['receipt_sha256']),(row['raw_path'],row['raw_sha256'])]:check(sha(p)==h,'authentic unchanged saved137 source')
   check(intent['semantic_key']==receipt['semantic_key']==key,'semantic source identity')
   check(intent['expected']['center']==jet['center']==cfg['anchor']and intent['expected']['window']==jet['window']==row['window'],'physical same center/window')
   check(jet['model']==intent['expected']['model']=='body32_unit_cauchy_v1'and jet['model_definition']==intent['expected']['model_definition'],'original model')
   check(receipt['exit_code']==0 and receipt['child_reaped']is True and receipt['precision_bits']==128 and receipt['width_acceptable']is True,'real137 accepted128 receipt')
   check(raw['status']=='all_real_nodes_supported'and raw['precision_bits']==128 and raw['full_group_terms']==967680 and raw['computed_sign_pairs']==483840,'saved native complete group')
   val=dyadic(raw['band']);grad=[dyadic(v)for v in raw['gradient']];pa=dyadic(raw['p_a']);pt=dyadic(raw['p_T']);ne=dyadic(raw['Euler_xgrad_plus5B_plusTpT_minusapa'])
   check(len(grad)==6 and val[0]>0 and pa[0]>0 and pt[0]>0 and ne[0]<=0<=ne[1],'saved native positivity/Euler')
   fields=[val,*grad,pa,pt,ne]
   check(len(fields)==10,'all ten fields')
   for lo,hi in fields:check(hi-lo<=F(1,65000*(1<<40))*max(F(1),abs(lo),abs(hi)),'saved complete ten-width gate')
   a,T=F(row['window']['a']),F(row['window']['T']);check(0<a<T,'positive own band')
   eu=[(5*val[0],5*val[1]),(T*pt[0],T*pt[1]),(-a*pa[1],-a*pa[0])]
   eu +=[tuple(sorted(F(c)*x for x in v))for c,v in zip(jet['center'],grad)]
   check(sum(v[0]for v in eu)<=0<=sum(v[1]for v in eu),'source rational same-band Euler')
   check(interval(jet['F_interval'])==val and[interval(v)for v in jet['face_gradient_intervals']]==grad[1:],'jet/raw exact physical fields')
   check(key not in sources_by_key,'unique band source')
   sources_by_key[key]={'center':list(map(F,jet['center'])),'window':row['window'],'L':val[0],'U':val[1],'gradient':grad[1:]}
  fresh={};degrees={};clampbranches={'L':0,'interior':0,'U':0};zero_points=0
  for key,row in saved['entries'].items():
   skey=row['semantic_key'];check(skey in sources_by_key,'one of five actual windows');src=sources_by_key[skey];point=row['point'];check(len(point)==6 and F(point[0])==1,'actual face point')
   check(key==hashlib.sha256(canon({'semantic_key':skey,'point':point})).hexdigest(),'saved sourcepoint key')
   d=[F(point[j])-src['center'][j]for j in range(1,6)]
   b=sum(min(lo*x,hi*x)for(lo,hi),x in zip(src['gradient'],d));L,U=src['L'],src['U']
   if b<L:f=L;branch='L'
   elif b>U:f=U;branch='U'
   else:f=b;branch='interior'
   check(f>0,'positive common source value')
   q=sum(x*x for x in d)/32;a=F(src['window']['a']);r=b/f-q/(a*a)
   polys={0:F(0)};polys.update({n:horner(r,n)for n in ODD});best=max(polys.values());degree=min(n for n,v in polys.items()if v==best);lower=f*best
   check(row['true_center']==cfg['anchor']and row['window']==src['window'],'cache physical origin')
   for name,actual in [('gradient_dot_lower',b),('clamped_true_value',f),('frame32_quadratic',q),('own_time_a',a),('exponent_exact',r),('nonnegative_exp_lower_polynomial_exact',best),('band_corner_lower_exact',lower)]:check(F(row[name])==actual,'independent exact '+name)
   check(list(map(F,row['displacement']))==d and row['chosen_odd_degree']==degree,'exact displacement and odd maximizer')
   check(row['strictly_positive_lower']==(lower>0)and row['zero_lower_is_not_strict_positivity']==(lower==0),'safe zero positivity distinction')
   fresh[key]=lower;degrees[degree]=degrees.get(degree,0)+1;clampbranches[branch]+=1;zero_points+=lower==0;pointcount+=1
   if pointcount%1000==0:save('progress.json',{'phase':'independent_Horner_saved_points','points_checked':pointcount,'assertions':checks})
  check(pointcount==5970,'all registered points recomputed once')
  check(len(author['boxes'])==128,'all saved boxes')
  g0,gs=F(1,65000),F(1,63500);check(cfg['target_g0']==author['g0']=='1/65000'and cfg['target_gstar']==author['gstar']=='1/63500','separate exact targets')
  used=set();passrows=[];unresolved=[];allids=set();cases=zero_bands=g0only=0
  for choice,result in zip(cfg['boxes'],author['boxes']):
   rid,nid=choice['root_id'],choice['node_id'];identity=(rid,nid);check(identity not in allids,'no duplicate original partition leaf');allids.add(identity)
   root=forest['roots'][rid];node=root['nodes'][nid]
   check(root['root_complete']is False and node['children']==[] and nid in root['queue']+root['blocked_leaves'],'original terminal unresolved partition leaf')
   check(node['status']==choice['status']in('queued','blocked')and node['box']==choice['box']==result['box'],'original partition source coordinates')
   check(result['root_id']==rid and result['node_id']==nid,'result/source leaf identity')
   check(len(result['per_band'])==len(result['prefix_box_lower_comparisons'])==5,'five band minima/prefixes')
   total=F(0);bandbounds=[]
   for wi,meta in enumerate(cfg['jets']):
    band=result['per_band'][wi];check(band['semantic_key']==meta['semantic_key']and band['window']==meta['window'],'same-band lower identity')
    check(len(band['corner_references'])==32,'all32 target rectangle corners')
    vals=[];seenbits=set()
    for ref,bits in zip(band['corner_references'],itertools.product((0,1),repeat=5)):
     check(ref['corner_bits']==list(bits)and bits not in seenbits,'complete deterministic original32 corners');seenbits.add(bits)
     point=['1']+[choice['box'][j][bits[j]]for j in range(5)];k=hashlib.sha256(canon({'semantic_key':meta['semantic_key'],'point':point})).hexdigest()
     check(ref['cache_key']==k and k in fresh,'proper sourcepoint corner reference');check(saved['entries'][k]['point']==point,'actual original target corner')
     used.add(k);vals.append(fresh[k]);cases+=1
    minimum=min(vals);check(F(band['whole_box_rational_lower_exact'])==minimum and band['strictly_positive_lower']==(minimum>0),'per-band min32 exact and positivity')
    total+=minimum;bandbounds.append(str(minimum));zero_bands+=minimum==0;prefix=result['prefix_box_lower_comparisons'][wi]
    check(prefix['consecutive_band_count']==wi+1 and F(prefix['sum_of_bandwise_box_lowers_exact'])==total,'sum after per-band minima')
    check(prefix['passes_g0']==(total>=g0)and prefix['passes_gstar']==(total>=gs),'independent dual prefix target flags')
   p0,pstar=total>=g0,total>=gs;check(F(result['sum_of_bandwise_min32_lowers_exact'])==total,'exact full bandwise sum')
   check(result['passes_g0']==p0 and result['passes_gstar']==pstar,'exact final local target flags')
   cls='GSTAR_AND_G0_LOCAL_BOX_PASS'if pstar else'G0_ONLY_LOCAL_BOX_PASS'if p0 else'LOWER_CERTIFICATE_UNRESOLVED_BOTH'
   check(result['classification']==cls and result['global_g0_complete']is False and result['global_gstar_complete']is False and result['true_G_upper_bound']is False,'local-only no false upper/scope')
   record={'root_id':rid,'node_id':nid,'box':choice['box'],'verified_sum_lower':str(total),'verified_bandwise_min32_lowers':bandbounds,'passes_g0':p0,'passes_gstar':pstar}
   if pstar:passrows.append(record)
   elif not p0:unresolved.append({'root_id':rid,'node_id':nid})
   else:g0only+=1
   boxcount+=1
  check(cases==20480 and used==set(fresh),'no omitted/unused saved rational points or corners')
  check(len(passrows)==26 and len(unresolved)==102 and g0only==0,'26 both-pass102 unresolved exact IDs')
  check(author['g0_local_box_pass_count']==author['gstar_local_box_pass_count']==26 and author['g0_only_count']==0 and author['unresolved_both_count']==102,'saved result local counts')
  check(author['corner_band_case_count']==cases and author['unique_rational_source_point_evaluations']==5970 and author['Taylor_recurrence_updates']==5970*15==89550,'saved finite counters')
  check(author['zero_band_lower_count']==zero_bands==349,'zero band bounds count')
  check(all(author[k]==0 for k in['point_upper_obstruction_count','true_g_counterexample_count','actual_G_upper_bound_count','new_native_calls','new_exp_calls','checker_calls','new_subdivisions']),'no extra calls or unjustified exclusions')
  check(author['old_g0_paid_leaves_not_promoted_to_gstar']is True and author['global_g_certificate']is False and author['global_gstar_certificate']is False and author['new_C6_bound']is False,'no global shortcut or old-g0 promotion')
  for row in sources:check(sha(row['path'])==row['sha256'],'post source '+row['path'])
  save('verified_local_passes.json',{'status':'INDEPENDENT_EXACT_26_LOCAL_CERTIFICATES','author_result_sha256':sha(spec['result_path']),'terminal_forest_sha256':sha(spec['terminal_forest_path']),'pass_count':26,'passes':passrows,'unresolved_count':102,'unresolved_selected_ids':unresolved,'global_residual_merge_or_count_certified':False})
  status='PASS_ALL5970_POINTS_128BOXES_26_LOCAL_CERTIFICATES'
  save('review_results.json',{'status':status,'assertions':checks,'independent_saved_points':pointcount,'boxes_checked':boxcount,'box_band_corner_references_checked':cases,'saved137_jets_authenticated':5,'zero_rational_points':zero_points,'zero_band_lower_count':zero_bands,'g0_local_passes':26,'gstar_local_passes':26,'unresolved_selected':102,'g0_only':0,'odd_degree_histogram':degrees,'clamp_branch_histogram':clampbranches,'independent_polynomial_method':'separate factorial-coefficient Horner evaluation at all eight odd degrees, not producer shared recurrence','sources_hashed_before_and_after':len(sources),'new_native_exp_or_old_baseline_replays':0,'global_residual_merge_or_count_certified':False,'global_g_certificate':False,'new_C6_bound':False,'verified_passes_sha256':sha(R/'verified_local_passes.json'),'wall_before_receipt':time.monotonic()-start,'CPU_before_receipt':time.process_time()-c0})
 except BaseException as e:
  failure=repr(e);save('failure.json',{'failure':failure,'traceback':traceback.format_exc(),'saved_points_completed':pointcount,'boxes_completed':boxcount,'automatic_retry':False});raise
 finally:save('run_receipt.json',{'status':status,'pid':os.getpid(),'failure':failure,'wall_before_receipt':time.monotonic()-start,'CPU_before_receipt':time.process_time()-c0,'assertions':checks,'points_checked':pointcount,'boxes_checked':boxcount,'priority':os.getpriority(os.PRIO_PROCESS,0),'cores':1,'children':0,'time_or_CPU_ceiling':None,'automatic_retries':0,'new_native_exp_checker_or_old_baseline_replays':0})
 print(json.dumps({'status':status,'assertions':checks,'points':pointcount,'boxes':boxcount,'local_passes':26}),flush=True)
if __name__=='__main__':main()

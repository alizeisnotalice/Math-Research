#!/usr/bin/env python3
"""Dual-gated, finite closed-rectangle bisection with fixed source atlas; Arb512 only for new midfaces."""
import os,sys,json,hashlib,datetime,importlib.util,itertools,collections,time,traceback
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def atomic(p,x):
 p=Path(p);t=p.with_name(p.name+'.pending')
 with t.open('wb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
 os.replace(t,p);fd=os.open(str(p.parent),os.O_RDONLY)
 try:os.fsync(fd)
 finally:os.close(fd)
def write_new(p,x):
 with Path(p).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def load(name,p):
 sp=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def sourcepointkey(key,point):return hashlib.sha256(canon({'semantic_key':key,'point':point})).hexdigest()
def points(box):return [['1']+[box[i][b]for i,b in enumerate(bits)]for bits in itertools.product((0,1),repeat=5)]
def volume(box):
 v=Q(1)
 for lo,hi in box:v*=Q(hi)-Q(lo)
 return v
def split_geometry(box):
 widths=[Q(hi)-Q(lo)for lo,hi in box];axis=max(range(5),key=lambda i:(widths[i],-i));require(widths[axis]>0,'NO_POSITIVE_SPLIT_AXIS');mid=(Q(box[axis][0])+Q(box[axis][1]))/2
 lower=[b.copy()for b in box];upper=[b.copy()for b in box];lower[axis][1]=str(mid);upper[axis][0]=str(mid)
 others=[i for i in range(5)if i!=axis];face=[]
 for bits in itertools.product((0,1),repeat=4):
  coords=[None]*5;coords[axis]=str(mid)
  for i,b in zip(others,bits):coords[i]=box[i][b]
  face.append(['1']+coords)
 require(volume(lower)+volume(upper)==volume(box),'EXACT_SPLIT_VOLUME')
 return {'axis_free_index':axis,'axis_full_index':axis+1,'midpoint':str(mid),'lower_box':lower,'upper_box':upper,'midface_points':face,'closed_union':'parent rectangle equals two closed children; common boundary is exact midface; no ordered tightening'}
def ordered(point):
 x=list(map(Q,point));return len(x)==6 and x[0]==1 and all(x[i]>=x[i+1]for i in range(5))and x[-1]>=0

def release_requirements():
 return {'round':'134','authorized':True,'scope':'ONLY_C133_SEVEN_FIXED_ATLAS_LOCAL_RECTANGLES','target':'1/63500','precision_bits':512,'workers':1,'nice_minimum':10,'time_limit':None,'native_calls_permitted':0,'max_depth':6,'max_leaves_per_original_box':64,'max_nodes_total':889,'max_new_exp':225792,'run_batch_sha256':sha(P/'run_batch.py'),'detached_runner_sha256':sha(P/'detached_runner.py'),'detached_launcher_sha256':sha(P/'launch_detached.py'),'fixed_roots_and_atlas_sha256':sha(P/'fixed_roots_and_atlas.json'),'fixed_source_entries_sha256':sha(P/'fixed_source_entries.json'),'initial_point_cache_sha256':sha(P/'initial_point_cache.json'),'input_hashes_sha256':sha(P/'input_hashes.json'),'proof_sha256':sha(P/'PROOF.md'),'prepared_packet_sha256':sha(P/'prepared_packet.json')}
def release_gate():
 rp=P/'parent_release.json';require(rp.is_file(),'MAIN_AND_INDEPENDENT_REVIEW_RELEASE_REQUIRED');r=read(rp)
 for k,v in release_requirements().items():require(r.get(k)==v,'RELEASE_MISMATCH:'+k)
 require(r.get('main_approved')is True and r.get('independent_review_approved')is True,'DUAL_GATE_REQUIRED')
 return rp

def preflight():
 pins=read(P/'input_hashes.json')
 for group in ['frozen_inputs','fixed_selected_source_artifacts']:
  for path,h in pins[group].items():require(sha(path)==h,'FROZEN_SOURCE_HASH:'+path)
 cfg=read(P/'fixed_roots_and_atlas.json');entries=read(P/'fixed_source_entries.json')['entries'];initial=read(P/'initial_point_cache.json')['entries'];CUBE=P.parents[2];src=CUBE/'upper/round127/continuation_v121';reg=read(src/'registration_127.json');inp=read(reg['inputs']['117_inputs_path']);U=load('c134_orig_helpers',src/'precision_family_helpers.py');V=load('c134_orig_validator',src/'native_validator_family.py');B=load('c134_orig_base',inp['base_path']);B._c016_validate_jet_512=B.validate_jet;jets={}
 require(len(cfg['roots'])==7 and cfg['target']=='1/63500'and cfg['max_depth']==6,'FIXED_SELECTION')
 for k,m in entries.items():
  if'jet_path'in m:jet=read(m['jet_path'])
  else:
   jet=read(m['jet_ref']['path'])
   for sel in m['jet_ref'].get('selector',[]):jet=jet[sel]
  e=m['expected'];bits=m['precision_bits'];require(U.semantic_key(e,reg['runtime_and_fan_pins'][inp['native_paths']['fan']])==k,'SEMANTIC_IDENTITY');U.validate_jet_family(B,jet,e,allowed_bits=(128,256,512),pinned_binary_sha=reg['binaries'][str(bits)]['sha256'],pinned_source_sha=reg['binaries'][str(bits)]['source_sha256']);raw=V.validate_native(read(m['native_path']),e,bits);require(jet['F_interval']==list(map(str,raw[0]))and jet['face_gradient_intervals']==[list(map(str,g))for g in raw[1][1:]],'ORIGINAL_TRUE_JET_VALUES');require(read(m['receipt_path'])['exit_code']==0,'ORIGINAL_EXIT');require(len(jet['center'])==6 and Q(jet['center'][0])==1,'TRUE_FACE_CENTER');jets[k]=jet
 for ck,row in initial.items():require(ck==sourcepointkey(row['semantic_key'],row['point'])and Q(0)<Q(row['lower'])<=Q(row['upper']),'INITIAL_CACHE_BINDING')
 for root in cfg['roots']:
  require([len(b['candidate_keys'])for b in root['bands']]==[16,16,0,0,0],'NO_ATLAS_EXPANSION')
  for band in root['bands']:
   require(len(set(band['candidate_keys']))==len(band['candidate_keys']),'UNIQUE_FIXED_SOURCE_KEYS')
   for k in band['candidate_keys']:
    require(jets[k]['window']==band['window'],'MATCHING_FIXED_BAND')
    for point in points(root['box']):require(sourcepointkey(k,point)in initial,'ALL_OLD_CORNERS_REUSED')
 return cfg,entries,jets,initial

class CachedFamily:
 def __init__(self,cfg,jets,initial):
  self.cfg=cfg;self.jets=jets;self.cache=dict(initial);self.initial_keys=set(initial);self.math={k:([Q(x)for x in j['center']],tuple(map(Q,j['F_interval'])),[tuple(map(Q,g))for g in j['face_gradient_intervals']],Q(j['window']['a']))for k,j in jets.items()};self.bounds={k:(Q(r['lower']),Q(r['upper']))for k,r in initial.items()};self.counts={'new_exp':0,'cache_hits_initial':0,'cache_hits_new':0,'old_corner_exp_recomputed':0};self.arb=None;self.ctx=None
 def get(self,k,point,*,allow_new=False):
  ck=sourcepointkey(k,point)
  if ck in self.cache:
   self.counts['cache_hits_initial'if ck in self.initial_keys else'cache_hits_new']+=1;return self.bounds[ck]
  require(allow_new,'MISSING_EXPECTED_CACHED_VERTEX');require(ck not in self.initial_keys,'OLD_CORNER_RECOMPUTE_FORBIDDEN');require(self.arb is not None,'EXP_ONLY_AFTER_GATE');require(self.counts['new_exp']<self.cfg['max_new_exp'],'NEW_EXP_CAP')
  c,(lo,hi),grad,a=self.math[k];d=[Q(x)-y for x,y in zip(point[1:],c[1:])];b=sum(min(gl*di,gu*di)for(gl,gu),di in zip(grad,d));f=min(hi,max(lo,b));quad=sum(x*x for x in d)/32;exponent=b/f-quad/a**2
  def ball(q):return self.arb(q.numerator)/self.arb(q.denominator)
  self.ctx.prec=512;self.counts['new_exp']+=1;enc=ball(f)*ball(exponent).exp();require(enc.is_finite(),'FINITE_ARB_POINT');lower=Q(str(enc.lower().fmpq()));upper=Q(str(enc.upper().fmpq()));require(0<lower<=upper,'OUTWARD_POINT_ENCLOSURE');row={'semantic_key':k,'point':point,'lower':str(lower),'upper':str(upper),'origin':'C134_NEW_MIDFACE_POINT_ONLY','exact_formula_fields':{'displacement':list(map(str,d)),'gradient_dot_lower':str(b),'clamped_true_value':str(f),'frame32_quadratic':str(quad),'exponent':str(exponent)},'precision_bits':512};self.cache[ck]=row;self.bounds[ck]=(lower,upper);return lower,upper
 def point_family(self,root,point,*,allow_new=False):
  lo=Q(0);hi=Q(0);bands=[]
  for band in root['bands']:
   values=[(k,self.get(k,point,allow_new=allow_new))for k in band['candidate_keys']];bl=max((v[0]for k,v in values),default=Q(0));bu=max((v[1]for k,v in values),default=Q(0));lo+=bl;hi+=bu;bands.append({'window':band['window'],'max_robust_lower':str(bl),'max_robust_upper':str(bu),'lower_max_keys':[k for k,v in values if v[0]==bl],'upper_max_keys':[k for k,v in values if v[1]==bu],'missing_band_zero':not values})
  return {'point':point,'fixed_family_lower':str(lo),'fixed_family_upper':str(hi),'point_is_ordered':ordered(point),'bands':bands}
 def box_family(self,root,box):
  vertices=points(box);total=Q(0);bands=[]
  for band in root['bands']:
   values=[]
   for k in band['candidate_keys']:
    data=[self.get(k,p)for p in vertices];values.append((k,min(x[0]for x in data)))
   best=max((v for k,v in values),default=Q(0));total+=best;bands.append({'window':band['window'],'max_center_complete_box_lower':str(best),'maximizing_keys':[k for k,v in values if v==best],'missing_band_zero':not values})
  return {'summed_complete_box_lower':str(total),'bands':bands,'corner_points':vertices}

def run():
 rp=release_gate();require(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART');os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));cfg,entries,jets,initial=preflight()
 import flint
 from flint import arb,ctx
 require(flint.__FLINT_VERSION__=='3.6.0','PINNED_FLINT');family=CachedFamily(cfg,jets,initial);family.arb=arb;family.ctx=ctx;target=Q(cfg['target']);root_results=[];global_counts={'nodes_evaluated':0,'nodes_created':0,'split_attempts':0,'geometry_splits':0};(P/'roots').mkdir(exist_ok=False)
 write_new(P/'started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'release_sha256':sha(rp),'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'native_calls':0,'max_new_exp':225792,'max_nodes_total':889});t0=time.monotonic()
 try:
  for root in cfg['roots']:
   rid=root['original_node_id'];nodes={'r':{'id':'r','parent':None,'depth':0,'box':root['box'],'status':'QUEUED'}};queue=collections.deque(['r']);frontier={'r'};counts={'nodes_evaluated':0,'split_attempts':0,'geometry_splits':0};obstruction=None;before=family.counts.copy();global_counts['nodes_created']+=1
   while queue:
    nid=queue.popleft();node=nodes[nid];counts['nodes_evaluated']+=1;global_counts['nodes_evaluated']+=1;require(global_counts['nodes_evaluated']<=889,'NODE_EVALUATION_CAP');bound=family.box_family(root,node['box']);node['whole_box_certificate']=bound
    if Q(bound['summed_complete_box_lower'])>=target:node['status']='LOCAL_LEAF_PASS'
    elif node['depth']>=6:node['status']='DEPTH6_UNRESOLVED'
    else:
     counts['split_attempts']+=1;global_counts['split_attempts']+=1;require(counts['split_attempts']<=63 and global_counts['split_attempts']<=441,'SPLIT_ATTEMPT_CAP');split=split_geometry(node['box']);midcerts=[family.point_family(root,p,allow_new=True)for p in split['midface_points']];bad=[r for r in midcerts if Q(r['fixed_family_upper'])<target];node['split_attempt']={'geometry':split,'midface_point_certificates':midcerts}
     if bad:
      witness=min(bad,key=lambda r:(Q(r['fixed_family_upper']),r['point']));obstruction={'at_node':nid,'point':witness['point'],'point_is_ordered':witness['point_is_ordered'],'fixed_family_upper':witness['fixed_family_upper'],'target_minus_fixed_family_upper_exact':str(target-Q(witness['fixed_family_upper'])),'midface_certificate':witness,'scope':'only this original fixed-atlas full rectangle; if point outside wedge, no obstruction claim for ordered-clipped methods; never actual G upper'};node['status']='FIXED_FAMILY_FULL_RECTANGLE_OBSTRUCTION'
      for pending in queue:nodes[pending]['status']='NOT_EVALUATED_STOP_AFTER_OBSTRUCTION'
      queue.clear()
     else:
      node['status']='SPLIT_CLOSED_RECTANGLE';node['children']=[nid+'0',nid+'1'];frontier.remove(nid);counts['geometry_splits']+=1;global_counts['geometry_splits']+=1
      for suffix,box in [('0',split['lower_box']),('1',split['upper_box'])]:
       cid=nid+suffix;nodes[cid]={'id':cid,'parent':nid,'depth':node['depth']+1,'box':box,'status':'QUEUED'};queue.append(cid);frontier.add(cid);global_counts['nodes_created']+=1
      require(len(nodes)<=127 and len(frontier)<=64 and global_counts['nodes_created']<=889,'TREE_CAP')
    atomic(P/'roots'/(rid+'.json'),{'original_node_id':rid,'original_box':root['box'],'nodes':nodes,'closed_partition_leaf_ids':sorted(frontier),'queue':list(queue),'counts':counts,'obstruction':obstruction,'running':True,'global_g_certificate':False,'new_C6_bound':False})
   require(sum(volume(nodes[n]['box'])for n in frontier)==volume(root['box']),'EXACT_PARTITION_VOLUME');allpass=all(nodes[n]['status']=='LOCAL_LEAF_PASS'for n in frontier);result={'original_node_id':rid,'original_box':root['box'],'status':'ONLY_THIS_ORIGINAL_BOX_LOCAL_PASS'if allpass else'FIXED_FULL_RECTANGLE_FAMILY_OBSTRUCTED'if obstruction else'DEPTH6_UNRESOLVED','local_pass':allpass,'nodes':nodes,'closed_partition_leaf_ids':sorted(frontier),'closed_partition_volume_exact':str(sum(volume(nodes[n]['box'])for n in frontier)),'original_volume_exact':str(volume(root['box'])),'leaf_count':len(frontier),'counts':counts,'new_exp':family.counts['new_exp']-before['new_exp'],'cache_hits_initial':family.counts['cache_hits_initial']-before['cache_hits_initial'],'cache_hits_new':family.counts['cache_hits_new']-before['cache_hits_new'],'obstruction':obstruction,'global_g_certificate':False,'new_C6_bound':False};atomic(P/'roots'/(rid+'.json'),result);root_results.append({k:v for k,v in result.items()if k!='nodes'});atomic(P/'point_cache.json',{'entries':family.cache,'initial_entry_count':len(initial),'total_entries':len(family.cache),'actual_counts':family.counts,'old_corner_exp_recomputed':0})
  require(family.counts['new_exp']==len(family.cache)-len(initial),'UNIQUE_NEW_POINT_EXP_COUNT');write_new(P/'result.json',{'status':'FINITE_FIXED_BATCH_COMPLETE','target':'1/63500','selected_original_boxes':7,'locally_passed_original_boxes':sum(r['local_pass']for r in root_results),'original_box_results':root_results,'actual_counts':{**global_counts,**family.counts},'max_depth':6,'max_leaves_per_original_box':64,'max_nodes_total':889,'new_exp_cap':225792,'native_calls':0,'wall_seconds':time.monotonic()-t0,'global_g_certificate':False,'new_C6_bound':False,'scope':'only these7 original unresolved rectangles under unchanged fixed source atlas; no transfer from prior passed leaves and no actual-G upper claim'})
 except BaseException:
  atomic(P/'point_cache.json',{'entries':family.cache,'initial_entry_count':len(initial),'actual_counts':family.counts,'old_corner_exp_recomputed':0});atomic(P/'exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'actual_counts':{**global_counts,**family.counts},'native_calls':0,'global_g_certificate':False,'new_C6_bound':False});raise
 else:atomic(P/'exit_receipt.json',{'status':'FINITE_FIXED_BATCH_COMPLETE','actual_counts':{**global_counts,**family.counts},'result_sha256':sha(P/'result.json'),'point_cache_sha256':sha(P/'point_cache.json'),'native_calls':0,'global_g_certificate':False,'new_C6_bound':False})
if __name__=='__main__':
 if sys.argv[1:]==['--saved-only-preflight']:
  cfg,entries,jets,initial=preflight();family=CachedFamily(cfg,jets,initial)
  for root in cfg['roots']:
   source=read(root['original_corner_certificate_path']);require(family.box_family(root,root['box'])['summed_complete_box_lower']==source['summed_disjoint_band_lower'],'OLD_CORNER_REUSE_REPRODUCES_C132');require(all(Q(family.point_family(root,p)['fixed_family_lower'])>=Q(cfg['target'])for p in points(root['box'])),'C133_ALL_ORIGINAL_CORNERS_PAID')
  require(family.counts['new_exp']==0,'NO_PREFLIGHT_EXP');print(json.dumps({'status':'SAVED_ONLY_PREFLIGHT_PASS','boxes':7,'source_jets':len(jets),'initial_source_point_entries':len(initial),'original_wholebox_bounds_reproduced_without_exp':7,'all_original_corners_paid':224,'native_calls':0,'new_exp':0}))
 elif not sys.argv[1:]:run()
 else:raise SystemExit('usage: run_batch.py [--saved-only-preflight]')

import os,sys,json,hashlib,itertools,collections,time,math
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,fmpq,ctx
ctx.prec=768
P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();canon=lambda x:(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
checks=0;expchecks=0;t0=time.monotonic()
def ck(x,label):
 global checks
 if not x:raise AssertionError(label)
 checks+=1
def progress(phase,**values):
 p=P/'verification_checkpoint.json';t=p.with_name(p.name+'.pending');t.write_text(json.dumps({'phase':phase,'pid':os.getpid(),'assertions':checks,'elapsed_seconds':time.monotonic()-t0,**values},indent=2)+'\n');os.replace(t,p)
def key(k,p):return hashlib.sha256(canon({'semantic_key':k,'point':p})).hexdigest()
def points(box):return[['1']+[box[j][bit]for j,bit in enumerate(bits)]for bits in itertools.product((0,1),repeat=5)]
def ordered(p):return all(a>=b for a,b in zip(map(Q,p),list(map(Q,p))[1:]))and Q(p[0])==1 and Q(p[-1])>=0

def split(box):
 axis=max(range(5),key=lambda j:(Q(box[j][1])-Q(box[j][0]),-j));m=(Q(box[axis][0])+Q(box[axis][1]))/2;l=[b.copy()for b in box];u=[b.copy()for b in box];l[axis][1]=str(m);u[axis][0]=str(m);face=[]
 for bits in itertools.product((0,1),repeat=4):
  p=['1']+[None]*5;p[axis+1]=str(m)
  for j,b in zip([j for j in range(5)if j!=axis],bits):p[j+1]=box[j][b]
  face.append(p)
 return{'axis_free_index':axis,'axis_full_index':axis+1,'midpoint':str(m),'lower_box':l,'upper_box':u,'midface_points':face,'closed_union':'parent rectangle equals two closed children; common boundary is exact midface; no ordered tightening'}
def volume(box):return math.prod(Q(u)-Q(l)for l,u in box)
reg=rd(P/'registration.json');ck(sha(__file__)==reg['checker_sha256'],'registered checker');ck(os.getpriority(os.PRIO_PROCESS,0)>=10,'nice>=10')
A=Path(reg['author_directory']);prior=rd(Path(reg['preflight_review_path']));ck(prior['status']=='PASS_C134_STATIC_FIXED_ATLAS_CACHE_AND_CLOSED_GEOMETRY_GATE','sealed preflight PASS')
for n,h in prior['frozen_pins'].items():ck(sha(A/n)==h,'same frozen source/prepared data')
cfg=rd(A/'fixed_roots_and_atlas.json');entries=rd(A/'fixed_source_entries.json')['entries'];initial=rd(A/'initial_point_cache.json')['entries'];pc=rd(A/'point_cache.json');cache=pc['entries'];result=rd(A/'result.json');target=Q(cfg['target'])
ck(target==Q(1,63500)and len(cfg['roots'])==7 and len(entries)==82,'frozen seven local atlases')
ck(len(initial)==16929 and all(cache[k]==r for k,r in initial.items()),'all old cache values/formula references retained byte-equivalent; no old exp')
newkeys=set(cache)-set(initial);ck(len(newkeys)==21126==result['actual_counts']['new_exp'],'actual finite new endpoint count')
jets={}
for k,m in entries.items():
 ck(sha(m['jet_path'])==m['jet_sha256'],'original fixed source jet bytes')
 j=rd(m['jet_path']);jets[k]=(list(map(Q,j['center'])),tuple(map(Q,j['F_interval'])),[tuple(map(Q,g))for g in j['face_gradient_intervals']],Q(j['window']['a']))
progress('NEW_SOURCEPOINT_ARB768',done=0,total=len(newkeys))
for ni,kp in enumerate(sorted(newkeys)):
 r=cache[kp];k=r['semantic_key'];point=r['point'];ck(k in entries and kp==key(k,point)and all(str(Q(v))==v for v in point),'new cache exact sourcepoint identity')
 c,(L,U),grad,a=jets[k];d=[Q(v)-w for v,w in zip(point[1:],c[1:])];b=sum(min(g[0]*z,g[1]*z)for g,z in zip(grad,d));f=min(U,max(L,b));q=sum(z*z for z in d)/32;e=b/f-q/a**2
 fields={'displacement':list(map(str,d)),'gradient_dot_lower':str(b),'clamped_true_value':str(f),'frame32_quadratic':str(q),'exponent':str(e)}
 ck(r['origin']=='C134_NEW_MIDFACE_POINT_ONLY'and r['precision_bits']==512 and r['exact_formula_fields']==fields,'independent new point formula')
 enc=arb(fmpq(f.numerator,f.denominator))*arb(fmpq(e.numerator,e.denominator)).exp();ck(enc.is_finite(),'finite fresh768 exp');ends=[]
 for z in(enc.lower(),enc.upper()):m,ex=z.man_exp();ends.append(Q(int(m))*Q(2)**int(ex))
 ck(0<Q(r['lower'])<=ends[0]<=ends[1]<=Q(r['upper']),'new saved512 contains fresh768');expchecks+=1
 if(ni+1)%5000==0:progress('NEW_SOURCEPOINT_ARB768',done=ni+1,total=len(newkeys))
# Independently simulate only finite saved-cache lookups and BFS decisions, with no further exp.
bounds={k:(Q(r['lower']),Q(r['upper']))for k,r in cache.items()};seen=set(initial);counts={'new_exp':0,'cache_hits_initial':0,'cache_hits_new':0,'old_corner_exp_recomputed':0};globalcounts={'nodes_evaluated':0,'nodes_created':0,'split_attempts':0,'geometry_splits':0};mid_ordered=mid_outside=0;mid_unique=set();mid_ranges={tag:{name:None for name in ['min_lower','max_lower','min_upper','max_upper']}for tag in ['ordered','outside_ordered']};ordered_obstructions=outside_obstructions=0;all_leaves=passed_leaves=unresolved_leaves=0;summary=[];root_core=[]
def get(k,p,allow_new=False):
 kp=key(k,p);ck(kp in cache,'saved cache complete for every request')
 if kp in seen:counts['cache_hits_initial'if kp in initial else'cache_hits_new']+=1
 else:
  ck(allow_new and kp in newkeys,'only new proposed midface may create identity');seen.add(kp);counts['new_exp']+=1
 return bounds[kp]
def point_bound(root,p):
 global mid_ordered,mid_outside
 lo=hi=Q(0);bands=[]
 for band in root['bands']:
  vals=[(k,get(k,p,True))for k in band['candidate_keys']];l=max((v[0]for k,v in vals),default=Q(0));u=max((v[1]for k,v in vals),default=Q(0));lo+=l;hi+=u
  bands.append({'window':band['window'],'max_robust_lower':str(l),'max_robust_upper':str(u),'lower_max_keys':[k for k,v in vals if v[0]==l],'upper_max_keys':[k for k,v in vals if v[1]==u],'missing_band_zero':not vals})
 isordered=ordered(p);mid_ordered+=isordered;mid_outside+=not isordered;mid_unique.add(tuple(p))
 stats=mid_ranges['ordered'if isordered else'outside_ordered']
 for name,value in [('min_lower',lo),('max_lower',lo),('min_upper',hi),('max_upper',hi)]:stats[name]=value if stats[name]is None else min(stats[name],value)if name.startswith('min')else max(stats[name],value)
 return{'point':p,'fixed_family_lower':str(lo),'fixed_family_upper':str(hi),'point_is_ordered':isordered,'bands':bands}
def box_bound(root,box):
 vertices=points(box);total=Q(0);bands=[]
 for band in root['bands']:
  vals=[(k,min(get(k,p)[0]for p in vertices))for k in band['candidate_keys']];v=max((v for k,v in vals),default=Q(0));total+=v
  bands.append({'window':band['window'],'max_center_complete_box_lower':str(v),'maximizing_keys':[k for k,value in vals if value==v],'missing_band_zero':not vals})
 return{'summed_complete_box_lower':str(total),'bands':bands,'corner_points':vertices}
for ri,root in enumerate(cfg['roots']):
 rid=root['original_node_id'];saved=rd(A/'roots'/(rid+'.json'));nodes=saved['nodes'];expected={'r':{'id':'r','parent':None,'depth':0,'box':root['box'],'status':'QUEUED'}};queue=collections.deque(['r']);frontier={'r'};rc={'nodes_evaluated':0,'split_attempts':0,'geometry_splits':0};before=counts.copy();obstruction=None;globalcounts['nodes_created']+=1
 for band in root['bands']:ck(len(band['candidate_keys'])==len(set(band['candidate_keys']))<=16 and all(entries[k]['expected']['window']==band['window']for k in band['candidate_keys']),'same original per-box atlas/window')
 while queue:
  nid=queue.popleft();n=expected[nid];rc['nodes_evaluated']+=1;globalcounts['nodes_evaluated']+=1;bound=box_bound(root,n['box']);n['whole_box_certificate']=bound
  if Q(bound['summed_complete_box_lower'])>=target:n['status']='LOCAL_LEAF_PASS'
  elif n['depth']>=6:n['status']='DEPTH6_UNRESOLVED'
  else:
   rc['split_attempts']+=1;globalcounts['split_attempts']+=1;geom=split(n['box']);midcerts=[point_bound(root,p)for p in geom['midface_points']];n['split_attempt']={'geometry':geom,'midface_point_certificates':midcerts};bad=[r for r in midcerts if Q(r['fixed_family_upper'])<target]
   ck(len(geom['midface_points'])==16 and {tuple(p)for p in points(geom['lower_box'])+points(geom['upper_box'])}=={tuple(p)for p in points(n['box'])+geom['midface_points']},'closed split geometry/inherited and midface corners')
   if bad:
    witness=min(bad,key=lambda r:(Q(r['fixed_family_upper']),r['point']));obstruction={'at_node':nid,'point':witness['point'],'point_is_ordered':witness['point_is_ordered'],'fixed_family_upper':witness['fixed_family_upper'],'target_minus_fixed_family_upper_exact':str(target-Q(witness['fixed_family_upper'])),'midface_certificate':witness,'scope':'only this original fixed-atlas full rectangle; if point outside wedge, no obstruction claim for ordered-clipped methods; never actual G upper'};n['status']='FIXED_FAMILY_FULL_RECTANGLE_OBSTRUCTION';ordered_obstructions+=witness['point_is_ordered'];outside_obstructions+=not witness['point_is_ordered']
    for pending in queue:expected[pending]['status']='NOT_EVALUATED_STOP_AFTER_OBSTRUCTION'
    queue.clear()
   else:
    n['status']='SPLIT_CLOSED_RECTANGLE';n['children']=[nid+'0',nid+'1'];frontier.remove(nid);rc['geometry_splits']+=1;globalcounts['geometry_splits']+=1
    for suffix,box in [('0',geom['lower_box']),('1',geom['upper_box'])]:
     child=nid+suffix;expected[child]={'id':child,'parent':nid,'depth':n['depth']+1,'box':box,'status':'QUEUED'};queue.append(child);frontier.add(child);globalcounts['nodes_created']+=1
  ck(nodes[nid]==n,'every saved BFS node bound/decision/newpoint envelope exactly matches')
 ck(nodes==expected,'all stopped/queued/created node records accounted')
 ck(saved['closed_partition_leaf_ids']==sorted(frontier)and sum((Q(1,2**nodes[n]['depth'])for n in frontier),Q(0))==1 and sum(volume(nodes[n]['box'])for n in frontier)==volume(root['box']),'exact complete closed recursive partition/frontier')
 ck(not any(b.startswith(a)for a,b in zip(sorted(frontier),sorted(frontier)[1:])),'prefix-free leaf partition')
 allpass=all(nodes[n]['status']=='LOCAL_LEAF_PASS'for n in frontier);status='ONLY_THIS_ORIGINAL_BOX_LOCAL_PASS'if allpass else'FIXED_FULL_RECTANGLE_FAMILY_OBSTRUCTED'if obstruction else'DEPTH6_UNRESOLVED'
 ck(saved['local_pass']==allpass and saved['status']==status and saved['obstruction']==obstruction and saved['counts']==rc,'strict original-box pass/stop policy')
 ck(saved['leaf_count']==len(frontier)<=64 and max(n['depth']for n in nodes.values())<=6 and len(nodes)<=127 and saved['closed_partition_volume_exact']==saved['original_volume_exact']==str(volume(root['box'])),'finite exact volume/depth caps')
 ck(saved['new_exp']==counts['new_exp']-before['new_exp']and saved['cache_hits_initial']==counts['cache_hits_initial']-before['cache_hits_initial']and saved['cache_hits_new']==counts['cache_hits_new']-before['cache_hits_new'],'per-original-box exact cache operation counts')
 ck(saved['global_g_certificate']is False and saved['new_C6_bound']is False,'local scope flags')
 all_leaves+=len(frontier);paid=sum(nodes[n]['status']=='LOCAL_LEAF_PASS'for n in frontier);passed_leaves+=paid;unresolved_leaves+=len(frontier)-paid
 summary.append({k:v for k,v in saved.items()if k!='nodes'});root_core.append({'original_node_id':rid,'minimum_leaf_lower_minus_target_exact':str(min(Q(nodes[n]['whole_box_certificate']['summed_complete_box_lower'])for n in frontier)-target),'local_pass':allpass,'status':status,'leaves':len(frontier),'paid_leaves':paid,'unresolved_leaves':len(frontier)-paid,'counts':rc,'new_exp':saved['new_exp']});progress('SAVED_BFS_PARTITION_AND_COUNT_REVIEW',roots_complete=ri+1,new768_exp=expchecks)
ck(seen==set(cache)and counts['new_exp']==expchecks==len(cache)-len(initial)==21126,'all and only reached unique new identities were evaluated')
ck(result['original_box_results']==summary and result['actual_counts']=={**globalcounts,**counts}==pc['actual_counts'],'all consolidated saved counts/results')
ck(result['locally_passed_original_boxes']==sum(r['local_pass']for r in root_core)==3 and result['selected_original_boxes']==7,'exact three complete original rectangle passes')
ck(globalcounts=={'nodes_evaluated':139,'nodes_created':139,'split_attempts':66,'geometry_splits':66}and counts=={'new_exp':21126,'cache_hits_initial':28062,'cache_hits_new':126940,'old_corner_exp_recomputed':0},'exact reported finite operation accounting')
ck(sum(counts[k]for k in ['new_exp','cache_hits_initial','cache_hits_new'])==139*32*32+66*16*32 and expchecks<=225792 and globalcounts['nodes_created']<=889,'new-exp cap and exact lookup budget identity')
ck(result['native_calls']==0 and result['global_g_certificate']is False and result['new_C6_bound']is False,'no global/native/C6 claims')
out={'status':'PASS_ONE_COMPLETE_C134_SAVED_NEW_EXP_ARB768_AND_CLOSED_PARTITIONS','assertions':checks,'new_exp_768_verified':expchecks,'old16929_exp_reevaluated':0,'old_cache_entries_retained':len(initial),'actual_counts':{**globalcounts,**counts},'original_rectangle_passes':3,'total_partition_leaves':all_leaves,'paid_leaves':passed_leaves,'unresolved_leaves':unresolved_leaves,'midface_point_requests_ordered':mid_ordered,'midface_point_requests_outside_ordered':mid_outside,'midface_fixed_family_enclosure_ranges':{tag:{name:None if value is None else str(value)for name,value in ranges.items()}for tag,ranges in mid_ranges.items()},'unique_midface_physical_points':len(mid_unique),'ordered_obstructions':ordered_obstructions,'outside_ordered_full_rectangle_obstructions':outside_obstructions,'roots':root_core,'native_producer_old_baseline_evaluator_calls':0,'global_g_or_C6_proved':False,'precision_bits':768}
(P/'review_results.json').write_text(json.dumps(out,indent=2)+'\n');progress('COMPLETE',new768_exp=expchecks,roots=7);print(json.dumps({k:v for k,v in out.items()if k!='roots'}),flush=True)

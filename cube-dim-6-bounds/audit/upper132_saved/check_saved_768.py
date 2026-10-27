import sys,os,time,json,hashlib,itertools
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,fmpq,ctx
ctx.prec=768
P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks=0;corners=0
def ck(x,label):
 global checks
 if not x:raise AssertionError(label)
 checks+=1
def progress(data):
 t=P/'verification_checkpoint.own-tmp';t.write_text(json.dumps(data,indent=2)+'\n');os.replace(t,P/'verification_checkpoint.json')
def ends(value):
 ck(value.is_finite(),'finite fresh768 enclosure');result=[]
 for z in (value.lower(),value.upper()):
  m,e=z.man_exp();result.append(Q(int(m))*Q(2)**int(e))
 return result
reg=rd(P/'registration.json');ck(sha(__file__)==reg['checker_sha256'],'registered independent checker')
ck(os.getpriority(os.PRIO_PROCESS,0)>=10,'nice>=10')
A=Path(reg['author_directory']);cfg=rd(A/'selected_boxes_and_centers.json');entries=rd(A/'selected_source_cache_entries.json')['entries'];saved=rd(A/'result.json')
prior=rd(Path(reg['preflight_review_path']));ck(prior['status']=='PASS_FINAL_STATIC_C132_PROOF_SOURCE_GATE','preflight PASS')
for n,h in prior['frozen_pins'].items():ck(sha(A/n)==h,'preflight source unchanged '+n)
# Hash frozen old source artifacts only; do not replay their native/metadata evaluators.
hashspec=rd(A/'input_hashes.json')
for group in ['frozen_sources','selected_original_artifacts']:
 for path,h in hashspec[group].items():ck(sha(path)==h,'original pinned bytes '+path)
jets={k:rd(m['jet_path'])for k,m in entries.items()}
ck(len(jets)==325,'325 true source jets')
ck(cfg['selected_box_count']==len(cfg['selected_boxes'])==93,'93 frozen target boxes')
ids=[x['node_id']for x in cfg['selected_boxes']]
ck(set(ids)=={x.stem for x in(A/'box_results').glob('*.json')} and len(set(ids))==93,'93 complete saved local boxes')
ck(cfg['target']=='1/63500'and cfg['reference_target']=='1/65000','both exact targets')
ck(len(cfg['windows'])==5 and all(0<Q(w['a'])<Q(w['T'])for w in cfg['windows'])and all(Q(cfg['windows'][j+1]['T'])==Q(cfg['windows'][j]['a'])for j in range(4)),'five disjoint positive bands')
summary=[];newpasses=refpasses=0
for ni,item in enumerate(cfg['selected_boxes']):
 out=rd(A/'box_results'/(item['node_id']+'.json'));ck(out['node_id']==item['node_id']and out['box']==item['box'],'original local box preserved')
 ck([b['window']for b in item['bands']]==cfg['windows']==[b['window']for b in out['band_envelopes']],'each five bands exactly once')
 lower=Q(0)
 for ib,ob in zip(item['bands'],out['band_envelopes']):
  refs=ib['candidate_keys'];rows=ob['center_bounds'];keys=[r['semantic_key']for r in refs]
  ck(len(keys)==len(set(keys))<=16 and [r['semantic_key']for r in rows]==keys,'same band unique chosen centers')
  lows=[]
  for ref,row in zip(refs,rows):
   k=ref['semantic_key'];jet=jets[k];m=entries[k]
   ck(row['true_center']==ref['true_center']==jet['center']==m['expected']['center']and Q(jet['center'][0])==1,'original true center preserved')
   ck(row['window']==ib['window']==jet['window']==m['expected']['window'],'exact same original band')
   ck(row['source_artifact_key']==m['artifact_key']and row['source_proof']==jet['source_proof']and row['precision_bits']==512,'saved source reference')
   center=list(map(Q,jet['center'][1:]));bounds=[tuple(map(Q,b))for b in item['box']];L,U=map(Q,jet['F_interval']);grad=[tuple(map(Q,b))for b in jet['face_gradient_intervals']];a=Q(ib['window']['a'])
   ck(0<L<=U and len(bounds)==len(grad)==5,'positive authentic source interval')
   ck(len(row['corners'])==32,'32 saved corners')
   rowlow=[];rowhigh=[]
   for bits,corner in zip(itertools.product((0,1),repeat=5),row['corners']):
    d=[bounds[j][bit]-center[j]for j,bit in enumerate(bits)]
    b=sum(min(g[0]*v,g[1]*v)for g,v in zip(grad,d));f=min(U,max(L,b));q=sum(v*v for v in d)/32;e=b/f-q/(a*a)
    ck(corner['corner_bits']==list(bits)and list(map(Q,corner['displacement']))==d,'corner and displacement')
    ck(Q(corner['gradient_dot_lower'])==b and Q(corner['clamped_true_value'])==f and Q(corner['frame32_quadratic'])==q and Q(corner['exponent'])==e,'independent exact robust corner formula')
    fresh=arb(fmpq(f.numerator,f.denominator))*arb(fmpq(e.numerator,e.denominator)).exp();fl,fu=ends(fresh);sl,su=Q(corner['enclosure_lower']),Q(corner['enclosure_upper'])
    ck(sl<=fl<=fu<=su,'fresh768 enclosed by saved512 corner')
    rowlow.append(sl);rowhigh.append(su);corners+=1
   ck(Q(row['complete_box_lower'])==min(rowlow)and Q(row['robust_tangent_min_enclosure_upper'])==min(rowhigh),'whole-box min before same-band max')
   lows.append(min(rowlow))
  best=max(lows,default=Q(0));lower+=best
  ck(Q(ob['band_complete_box_lower'])==best and ob['missing_band_zero']==(len(rows)==0),'same-band max and missing zero')
  ck(ob['maximizing_keys']==[k for k,v in zip(keys,lows)if v==best],'same-band max provenance keys')
 ck(Q(out['summed_disjoint_band_lower'])==lower,'disjoint time-band positive sum')
 new=lower>=Q(cfg['target']);ref=lower>=Q(cfg['reference_target']);newpasses+=new;refpasses+=ref
 ck(out['passes_new_local_target']==new and out['passes_reference_local_target']==ref,'two exact local targets')
 ck(out['global_g_certificate']is False and out['new_C6_bound']is False,'local scope only')
 summary.append({k:v for k,v in out.items()if k!='band_envelopes'})
 if(ni+1)%10==0 or ni==92:progress({'phase':'SAVED_CORNER_VERIFY','pid':os.getpid(),'complete_boxes':ni+1,'corners768':corners,'assertions':checks})
ck(corners==cfg['actual_planned_corner_exp']==saved['corner_exp_used']==95232,'one finite exact corner count')
ck(saved['boxes']==summary and saved['selected_boxes']==93 and saved['root_id']==cfg['root_id'],'complete consolidated result')
ck(saved['new_target_local_passes']==newpasses==1 and saved['reference_target_local_passes']==refpasses==9,'consolidated exact1 and9 local passes')
ck(saved['all_selected_pass_new_target']is False and saved['root_new_target_closed']is False and saved['global_g_certificate']is False and saved['new_C6_bound']is False and saved['native_calls']==0,'no root/global/C6 claim')
result={'status':'PASS_ONE_COMPLETE_INDEPENDENT_C132_SAVED_ARB768_REVIEW','assertions':checks,'boxes':93,'new768_corner_exp':corners,'precision_bits':768,'new_target_local_passes':newpasses,'reference_target_local_passes':refpasses,'native_producer_old_baseline_evaluator_calls':0,'root_globalg_C6_proved':False,'method':'independent exact rational reconstruction; fresh768 directed enclosures; wholebox min then sameband max then disjoint-band sum'}
(P/'review_results.json').write_text(json.dumps(result,indent=2)+'\n');progress({'phase':'COMPLETE',**result,'pid':os.getpid()});print(json.dumps(result),flush=True)

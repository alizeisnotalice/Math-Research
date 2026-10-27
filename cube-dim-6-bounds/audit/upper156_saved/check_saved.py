"""Independent whole C156 saved certificate. Never imports producer, replays no old functional."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,permutations,product
from math import prod
import json,gzip,hashlib,os,sys,time,traceback
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
TARGET=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/audit/upper156_saved')
CHECKS=0
def need(x,m):
 global CHECKS
 CHECKS+=1
 if not x:raise ArithmeticError(m)
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def interval(raw):
 a,b=map(F,raw);need(0<a<=b,'POSITIVE_ORDERED_INTERVAL');return a,b
def shown(z):return list(map(str,z))
def determinant(matrix):
 # Independent literal Leibniz expansion, dimension fixed six, no polynomial/fan engine.
 total=0
 for p in permutations(range(6)):
  parity=sum(p[i]>p[j]for i in range(6)for j in range(i+1,6))
  total+=(-1 if parity%2 else 1)*prod(matrix[i][p[i]]for i in range(6))
 return total
def release_requirements():
 r=read(P/'registration.json')
 return {'round':'156_SAVED_INDEPENDENT','authorized':True,'scope':'ONE_COMPLETE_SAVED156_AUDIT_10416_PAIRS_NO_PRODUCER_OR_OLD_BASELINE','registration_sha256':sha(P/'registration.json'),'inputs_sha256':sha(P/'inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'source_sha256':r['source_sha256'],'limits':r['limits'],'workers':1,'nice_minimum':15,'time_or_CPU_limit':None,'automatic_retries':0}
def run():
 need(P==TARGET,'TERRA_ONLY_ROOT');reg=read(P/'registration.json');cfg=read(P/'inputs.json');rel=read(P/'parent_release.json')
 for k,v in release_requirements().items():need(rel.get(k)==v,'RELEASE:'+k)
 need(rel.get('main_approved')is True and rel.get('independent_review_approved')is True,'DUAL_GATE')
 for n,h in reg['source_sha256'].items():need(sha(P/n)==h,'SOURCE_PIN:'+n)
 need(sha(P/'inputs.json')==reg['inputs_sha256']and sha(P/'source_manifest.json')==reg['source_manifest_sha256']and cfg['limits']==reg['limits'],'FROZEN_PACKET')
 frozen=read(P/'source_manifest.json')['sources']
 for row in frozen:need(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'FROZEN_SOURCE:'+row['path'])
 need(not(P/'started.json').exists(),'ONE_AUDIT_NO_RESTART');need(str(Path(sys.executable).resolve())==cfg['python_executable'],'ACTUAL_PINNED_PYTHON')
 need(all(os.environ.get(n)=='1'for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']),'ONE_THREAD_ENV')
 os.nice(max(0,15-os.getpriority(os.PRIO_PROCESS,0)));save('started.json',{'pid':os.getpid(),'release_sha256':sha(P/'parent_release.json'),'counts_after_process_loss':'UNKNOWN except durable completed checkpoints','time_or_CPU_limit':None})
 start=time.monotonic();cpu=time.process_time();current='metadata';counts={'raw_pairs':0,'caps_evaluated':0,'exact_column_determinants':0,'four_seed_pairing_products':0,'cauchy_sqrt_cache_misses':0,'cauchy_sqrt_pruned_by_exact_range':0,'new_I0_closed':0,'new_positive_k_anchor':0};scalar={k:{'attempted':0,'completed':0,'raised_exception':0}for k in ['exp','erfc','sqrt_h','sqrt_pi','pi_constant']};inflight=None;done=[]
 try:
  root=Path(cfg['producer_root']);saved=read(root/'result.json');inp=read(root/'inputs.json');paths=cfg['upstream_paths'];old130=read(paths['results.json']);old146=read(paths['candidate_certificate.json']);up128=read(paths['upstream.json']);table=read(paths['cone_table.json'])['cones']
  need(sha(root/'result.json')==cfg['producer_result_sha256']=='f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f','ACTUAL_RESULT');need(inp['cones']==table,'ORIGINAL_CONE_TABLE')
  need(inp['h']==old130['h']==saved['h']==up128['h_values'][5]=='1/128'and inp['p_interval']==up128['p_intervals'][5]==saved['p_interval']and inp['kappa_interval']==old130['pi_factor_interval']==saved['kappa_interval'],'INHERITED_EXACT_SCALARS')
  need(old130['status']=='COMPLETE_FIXED_FINITE_BOX_PILOT'and old130['pair_terms']==496 and old130['cone_count']==21 and old130['probability_renormalized']is False,'OLD130_COMPLETE_TRUE_MOMENT_CERTIFICATE')
  pL,pU=interval(inp['p_interval']);need(pU<1,'P_PROBABILITY');h=F(1,128);p2=(pL*pL,pU*pU)
  # Only these fresh scalar operations are new oracle work. No erf, I30, old kernel or sqrt(P).
  import flint
  from flint import arb,ctx
  need(flint.__version__=='0.9.0'and flint.__FLINT_VERSION__=='3.6.0'and str(Path(flint.__file__).resolve())==cfg['flint_init_path'],'PINNED_RUNTIME');ctx.prec=768;ctx.threads=1
  def ar(x):x=F(x);return arb(x.numerator)/x.denominator
  def ends(x):need(x.is_finite(),'FINITE_ARBITRARY_INTERVAL');return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))
  def call(label,fn):
   nonlocal inflight
   need(scalar[label]['attempted']<{'exp':3,'erfc':2,'sqrt_h':1,'sqrt_pi':1,'pi_constant':1}[label],'FINITE_NEW_SCALAR_ATTEMPT_CAP:'+label)
   inflight={'operation':label,'state':'INFLIGHT_UNKNOWN'};scalar[label]['attempted']+=1
   try:value=fn()
   except BaseException:scalar[label]['raised_exception']+=1;inflight=None;raise
   scalar[label]['completed']+=1;inflight=None;return value
  current='fresh scalar768';s=call('sqrt_h',lambda:ar(h).sqrt());pi=call('pi_constant',arb.pi);sqrtpi=call('sqrt_pi',pi.sqrt)
  z=s*call('exp',lambda:ar(-h).exp())/sqrtpi
  # Algebraically reduced I0 avoids the producer's two dependent p occurrences.
  fresh_I=ar(h)+ar(h-F(1,2))*ar(pL).union(ar(pU))+z;counts['new_I0_closed']=1
  q=call('exp',lambda:(s+s).exp())*call('erfc',lambda:(s+ar(F(1,2))).erfc())
  fresh_K=call('erfc',lambda:(ar(F(1,2))-s).erfc())/2-q/2
  need(fresh_K>0,'POSITIVE_NEW_K');fresh_tau=q/fresh_K;fresh_C=fresh_K*call('exp',lambda:(s*fresh_tau).exp());counts['new_positive_k_anchor']=1
  saved_I=interval(saved['I0_interval']);saved_K=interval(saved['anchor_k_interval']);tau=[interval(x)for x in saved['tangent_tau_intervals']];C=[interval(x)for x in saved['tangent_C_intervals']]
  fresh={'I0':ends(fresh_I),'K_at_s':ends(fresh_K),'tau_at_s':ends(fresh_tau),'C_at_s':ends(fresh_C),'tau_at_0':((1-pU)/pU,(1-pL)/pL)}
  for name,old in [('I0',saved_I),('K_at_s',saved_K),('tau_at_s',tau[1]),('C_at_s',C[1]),('tau_at_0',tau[0])]:
   new=fresh[name];need(old[0]<=new[0]<=new[1]<=old[1],'NEW768_ENCLOSURE_WITHIN_SAVED512:'+name)
  need(C[0]==(pL,pU),'C0_EXACT_P');need(all(v['attempted']==v['completed']and v['raised_exception']==0 for v in scalar.values()),'SCALAR_COMPLETE')
  need({k:v['completed']for k,v in scalar.items()}=={'exp':3,'erfc':2,'sqrt_h':1,'sqrt_pi':1,'pi_constant':1},'FINITE_SCALAR_BUDGET')
  save('scalar_certificate.json',{'fresh_Arb_bits':768,'fresh_ranges':{k:shown(v)for k,v in fresh.items()},'all_fresh_ranges_within_saved512':True,'scalar_calls':scalar,'erf_calls':0,'new_Cauchy_root_calls':0})
  # Every stored Cauchy root is verified by exact squaring, never recomputed.
  roots={}
  for header in saved['cauchy_root_certificate_files']:
   path=root/header['file'];need(sha(path)==header['sha256'],'ROOT_FILE_HASH');data=json.loads(gzip.decompress(path.read_bytes()));need(len(data)==header['entries'],'ROOT_HEADER_COUNT')
   for key,bounds in data.items():
    P4=F(key);need(str(P4)==key and P4>0 and P4 not in roots,'CANONICAL_UNIQUE_ROOT_KEY');lo,hi=interval(bounds);need(lo*lo<=P4<=hi*hi,'EXACT_ROOT_SQUARE_CERTIFICATE');roots[P4]=(lo,hi)
  need(len(roots)==147,'FIXED_SAVED_ROOT_COUNT');used=set();i2=(saved_I[0]**2,saved_I[1]**2)
  def mul(ranges):return prod(x[0]for x in ranges),prod(x[1]for x in ranges)
  def flat(D):return p2 if not D else(min(p2[0],i2[0]/D),min(p2[1],i2[1]/D))
  def four(P4):
   if saved_I[0]**4>=16*P4*pU**4:counts['cauchy_sqrt_pruned_by_exact_range']+=1;return p2
   need(P4 in roots,'REQUIRED_PAID_ROOT');used.add(P4);counts['cauchy_sqrt_cache_misses']=len(used);lo,hi=roots[P4];return min(p2[0],i2[0]/(4*hi)),min(p2[1],i2[1]/(4*lo))
  def tangent(columns,r,t,normalized=False):
   # Independent whole-block formula, reciprocal product first, no producer helper.
   dlow=prod(1+tau[r][0]*b+tau[t][0]*c for b,c in columns);dhigh=prod(1+tau[r][1]*b+tau[t][1]*c for b,c in columns)
   v=C[r][0]*C[t][0]/dhigh,C[r][1]*C[t][1]/dlow
   return(v[0]/p2[1],v[1]/p2[0])if normalized else v
  body=[(1,)+v for v in product((-1,1),repeat=5)];headers={x['cone_id']:x for x in saved['pair_files']};need(len(saved['pair_files'])==len(headers)==21 and set(headers)==set(range(21)),'COMPLETE21_PAIR_FILES');need(len(table)==21 and [c['id']for c in table]==list(range(21)),'TRUE21_CONE_ORDER')
  totals={'flat':[F(0),F(0)],'finite_best':[F(0),F(0)]};alpha_sum=F(0)
  for cone in table:
   cid=cone['id'];current={'cone_id':cid};header=headers[cid];path=root/header['file'];need(sha(path)==header['sha256'],'PAIR_FILE_HASH');records=json.loads(gzip.decompress(path.read_bytes()));need(len(records)==header['raw_pairs']==496,'COMPLETE496_RAW_RECORDS')
   rays=cone['rays'];rates=cone['rates'];need(len(rays)==len(rates)==6 and all(len(v)==6 for v in rays),'SIX_DIMENSIONAL_FIXED_RAYS')
   absolute=[[abs(sum(v[k]*ray[k]for k in range(6)))for ray in rays]for v in body]
   need([sum(row[j]for row in absolute)for j in range(6)]==rates and all(type(r)is int and r>0 for r in rates),'ORIGINAL_RAY_L1_RATES')
   alpha=F(abs(determinant(rays)),prod(rates));need(alpha>0 and alpha==F(cone['alpha'])==F(header['alpha']),'INDEPENDENT_DETERMINANT_ALPHA');alpha_sum+=alpha
   nodes=[[F(absolute[i][j],rates[j])for j in range(6)]for i in range(32)];need([[str(x)for x in row]for row in nodes]==cone['nodes_by_body'],'ALL32_BODY_IDENTITIES')
   conetotals={'flat':[F(0),F(0)],'finite_best':[F(0),F(0)]}
   for savedrow,(i,j)in zip(records,combinations(range(32),2)):
    need(savedrow['body_pair_zero_based']==[i,j],'LEX_COMPLETE_PAIR_IDENTITY');columns=tuple(zip(nodes[i],nodes[j]));dets={key:abs(columns[key[0]][0]*columns[key[1]][1]-columns[key[1]][0]*columns[key[0]][1])for key in combinations(range(6),2)};counts['exact_column_determinants']+=15
    maxD=max(dets.values());f=flat(maxD);need(F(savedrow['best_seed_determinant'])==maxD and interval(savedrow['flat_budget_interval'])==f,'FLAT_SAVED_EXPRESSION')
    singles=[]
    for b,c in columns:
     candidates=[(F(1),F(1))]if b==c==0 else[(F(1),F(1))]+[tangent([(b,c)],r,t,True)for r,t in product(range(2),repeat=2)]
     singles.append((min(x[0]for x in candidates),min(x[1]for x in candidates)))
    candidates=[]
    def add(name,priority,value):need(0<value[0]<=value[1],'POSITIVE_COMPLETE_CAP');candidates.append((priority,name,value))
    z=mul(singles);add('all_single_columns',(0,), (p2[0]*z[0],p2[1]*z[1]))
    for subset,D in dets.items():
     if not D:continue
     seed=flat(D);rest=mul([singles[k]for k in range(6)if k not in subset]);add('seed2:'+','.join(map(str,subset)),(1,)+subset,(seed[0]*rest[0],seed[1]*rest[1]))
    for subset in combinations(range(6),4):
     a,b,c,d=subset;P4=max(dets[a,b]*dets[c,d],dets[a,c]*dets[b,d],dets[a,d]*dets[b,c]);counts['four_seed_pairing_products']+=3
     if not P4:continue
     seed=four(P4);rest=mul([singles[k]for k in range(6)if k not in subset]);add('seed4:'+','.join(map(str,subset)),(2,)+subset,(seed[0]*rest[0],seed[1]*rest[1]))
    for r,t in product(range(2),repeat=2):add('global_tangent:'+str(r)+str(t),(3,r,t),tangent(columns,r,t))
    need(len(candidates)==savedrow['caps_evaluated']<=35,'ALL_FINITE_CAPS_NO_OMISSION');low=min(candidates,key=lambda x:(x[2][0],x[0]));high=min(candidates,key=lambda x:(x[2][1],x[0]));best=low[2][0],high[2][1]
    need(best==interval(savedrow['finite_best_budget_interval'])and low[1]==savedrow['selected_lower_cap']and high[1]==savedrow['selected_upper_cap'],'MINIMUM_AND_STABLE_SELECTOR');need(best[1]<=f[1],'FINITE_NOT_WORSE_FLAT')
    counts['raw_pairs']+=1;counts['caps_evaluated']+=len(candidates);need(counts['raw_pairs']<=10416 and counts['caps_evaluated']<=364560,'FIXED_PAIR_CAP_BUDGET')
    for name,value in [('flat',f),('finite_best',best)]:
     for end in (0,1):conetotals[name][end]+=value[end]
   for name in totals:
    need(shown(conetotals[name])==header['unweighted_budget_sums'][name],'COMPLETE_CONE_SUM')
    for end in (0,1):totals[name][end]+=alpha*conetotals[name][end]
   done.append(cid);save('cone_'+str(cid).zfill(2)+'_checkpoint.json',{'cone_id':cid,'pairs_complete':496,'all_pair_intervals_selectors_match':True,'counts_so_far':dict(counts),'assertions':CHECKS,'source_gzip_sha256':header['sha256']})
  counts['cauchy_sqrt_cache_misses']=len(used);need(used==set(roots),'ALL_PAID_ROOTS_USED_NO_GAP');need(counts==saved['counts']and counts['raw_pairs']==10416 and counts['caps_evaluated']==229314,'ALL_ACTUAL_COUNTS')
  need(46080*alpha_sum==F(inp['R0'])and str(alpha_sum)==saved['source_alpha_sum'],'COMPLETE_ALPHA_AND_GROUP_SCALE')
  kappa=interval(inp['kappa_interval']);goals={'goal7':F(7,10**8),'goal8':F(8,10**8),'saved146_legal_cap':F(old146['current_valid_beta2_upper'])};need(goals['saved146_legal_cap']==F(inp['saved146_legal_beta2_upper']),'PRIOR146_VALID_UPPER')
  need(goals['goal7']<F(old146['true_beta2_upper_threshold_for_strict_g0']['sufficient_exact'])and goals['goal8']<F(old146['true_beta2_upper_threshold_for_strict_saved128_improvement']['sufficient_exact']),'SOURCE146_EXACT_SUFFICIENT_TARGET_INTERFACES')
  oldtrue=interval(old130['moment_interval']);budgets={}
  for name,value in totals.items():
   value=kappa[0]*value[0],kappa[1]*value[1];old=saved['budgets'][name];need(value==tuple(map(F,old['budget_expression_interval']))and value[1]==F(old['valid_true_beta2_upper']),'COMPLETE_NEW_BUDGET_ENDPOINTS');need(value[1]-value[0]<=F(inp['global_budget_width_gate']),'COMPLETE_WIDTH')
   upper={k:value[1]<v for k,v in goals.items()};excluded={k:value[0]>=v for k,v in goals.items()};need(upper==old['upper_comparisons']and excluded==old['this_fixed_budget_family_cannot_reach']and old['budget_lower_is_actual_mu2_lower']is False,'EXACT_ALL_TARGET_COMPARISONS')
   need(oldtrue[0]<=value[1],'CONSISTENCY_WITH130_TRUE_COMPLETE_MU2_LOWER')
   budgets[name]={'budget_expression_interval':shown(value),'certified_true_beta2_upper':str(value[1]),'upper_comparisons':upper,'fixed_family_only_lower_comparisons':excluded,'BL_is_actual_mu2_lower':False}
  need(budgets['finite_best']['upper_comparisons']=={'goal7':False,'goal8':False,'saved146_legal_cap':True},'EXPECTED_NEW_CANDIDATE_IMPROVEMENT_ONLY')
  for row in frozen:need(sha(row['path'])==row['sha256'],'FROZEN_SOURCE_AT_END')
  for n,hv in reg['source_sha256'].items():need(sha(P/n)==hv,'CHECKER_SOURCE_AT_END')
  certificate={'status':'PASS_COMPLETE_SAVED156_INDEPENDENT_SCALAR768_AND_RATIONAL_CAPS','producer_result_sha256':cfg['producer_result_sha256'],'certified_true_beta2_upper':budgets['finite_best']['certified_true_beta2_upper'],'budgets':budgets,'old130_true_complete_moment_interval':shown(oldtrue),'old130_lower_le_new_BU':True,'BL_is_actual_mu2_lower':False,'counts':counts,'scalar_calls':scalar,'new_erf_calls':0,'new_Cauchy_root_calls':0,'old_functional_or_producer_calls':0,'cones_complete':done,'checks':CHECKS,'wall_before_receipt':time.monotonic()-start,'CPU_before_receipt':time.process_time()-cpu,'release_sha256':sha(P/'parent_release.json'),'new_M_upper':False,'new_C6_bound':False,'new_global_g':False}
  save('certificate.json',certificate);save('exit_receipt.json',{'status':certificate['status'],'certificate_sha256':sha(P/'certificate.json'),'counts':counts,'scalar_calls':scalar,'checks':CHECKS,'automatic_retries':0});print(json.dumps({'status':certificate['status'],'raw_pairs':10416,'caps':counts['caps_evaluated'],'new_Cauchy_root_calls':0,'checks':CHECKS}))
 except BaseException:
  exact=inflight is None and all(v['attempted']==v['completed']+v['raised_exception']for v in scalar.values())
  save('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'current':current,'completed_cones':done,'observed_counts':counts,'scalar_observed':scalar,'scalar_counts_exact':exact,'scalar_final_counts':scalar if exact else None,'scalar_inflight_status':inflight,'process_crash_policy':'uncheckpointed progress unknown, never0','checks':CHECKS,'automatic_retries':0,'new_C6_bound':False});raise
if __name__=='__main__':run()

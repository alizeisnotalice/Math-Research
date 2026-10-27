import json,hashlib,time,resource,os,sys,collections
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0);st=time.monotonic();cpu=time.process_time();root=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper');p=root/'round044/compute';o=root/'round038/compute';prior=Path(__file__).with_name('lp_batch.json');pb=prior.read_bytes();pr=json.loads(pb);prev={r['case']:r for r in pr['rows']}; assert len(prev)==76
sha=lambda b:hashlib.sha256(b).hexdigest();counts=collections.Counter();rows=[];errors=[]
for f in sorted(p.glob('C_case_phasea*.json')):
 try:
  raw=f.read_bytes();v=json.loads(raw);cb=(p/v['source']).read_bytes();c=json.loads(cb);rb=(o/c['source_raw']).read_bytes();d=json.loads(rb);r=prev[c['source_raw']];assert sha(cb)==v['source_sha256']==r['case_sha256'] and sha(rb)==c['source_raw_sha256']==r['raw_sha256']
  has=v['exact_candidate'] is not None;counts[(str(v['solver_status']),has)]+=1
  if not has:continue
  assert v['solver_status']==0 and r['global_improves'] and r['global_proposal_feasible'];h=list(map(F,d['witness_budgets']));q=list(map(F,d['q_after']));qg=list(map(F,c['global_proposed_q']));C=d['maximal_cliques'];N=len(h);U=F(d['U']);assert h==list(map(F,c['h'])) and q==list(map(F,c['q_initial'])) and C==c['maximal_cliques']
  q1=list(map(F,v['exact_candidate']['q_feasible']));assert len(q)==len(qg)==len(q1)==N and min(q1)>=0
  assert sum((h[i]*q[i] for i in range(N)),F())==U==sum((h[i]*qg[i] for i in range(N)),F())==sum((h[i]*q1[i] for i in range(N)),F())
  load=lambda x:max(sum((x[i] for i in clique),F()) for clique in C)
  L0=load(q);Lg=load(qg);L1=load(q1);mid=(L0+Lg)/2;e=v['exact_candidate'];support=[i for i in range(N) if q1[i]!=q[i]]
  assert L0==F(c['L_initial']) and Lg==F(c['global_proposed_load']) and Lg<L0 and mid==F(e['exact_target']) and L1==F(e['exact_load']) and L1<=mid
  assert support==e['exact_changed_indices'] and len(support)==e['exact_support_count'];assert e['exact_payment_conserved'] and e['exact_load_below_target']==(L1<mid)
  rows.append({'case':c['source_raw'],'input_sha256':sha(rb),'candidate_sha256':sha(raw),'compact_sha256':sha(cb),'support':support,'support_count':len(support),'initial_load':str(L0),'global_feasible_load':str(Lg),'midpoint':str(mid),'candidate_load':str(L1),'target_slack':str(mid-L1),'checks':'pass'})
 except Exception as ex:errors.append({'file':f.name,'error':repr(ex)})
counts={str(k):v for k,v in counts.items()};out={'exact_candidates_pass':len(rows),'errors':errors,'record_categories':counts,'support_histogram':dict(collections.Counter(x['support_count'] for x in rows)),'prior_lp_batch_sha256':sha(pb),'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'results':rows}
Path(__file__).with_name('audit_c.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='results'},indent=2))

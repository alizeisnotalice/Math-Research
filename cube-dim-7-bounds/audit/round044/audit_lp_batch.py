import json,hashlib,time,resource,os,sys,collections
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0);start=time.monotonic();cpu=time.process_time();base=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper');p=base/'round044/compute';old=base/'round038/compute';prev=Path('/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round038/phasea/results')
def sha(b):return hashlib.sha256(b).hexdigest()
regraw=(p/'registration.json').read_bytes();reg=json.loads(regraw);regs={v['raw_filename']:v for v in reg['cases']};assert len(regs)==76
rows=[];errors=[]
for f in sorted(p.glob('case_phasea*.json')):
 try:
  blob=f.read_bytes();c=json.loads(blob);name=c['source_raw'];raw=(old/name).read_bytes();d=json.loads(raw);pa=(prev/(Path(name).stem+'.audit.json')).read_bytes();v=json.loads(pa);assert sha(raw)==c['source_raw_sha256']==regs[name]['raw_sha256']==v['input_sha256'] and v['status']=='pass'
  C=d['maximal_cliques'];assert c['maximal_cliques']==C and c['geometry']==d['geometry'];h=list(map(F,d['witness_budgets']));q=list(map(F,d['q_after']));assert h==list(map(F,c['h'])) and q==list(map(F,c['q_initial']));N=len(h);U=F(d['U']);assert U==F(c['U']) and min(q)>=0 and sum((h[i]*q[i] for i in range(N)),F())==U
  L=max(sum((q[i] for i in Cj),F()) for Cj in C);assert L==F(c['L_initial'])==F(d['L_after'])
  m=list(map(F,c['mass_primal']));z=list(map(F,c['mass_dual']));assert len(m)==len(C) and len(z)==N and min(m)>=0 and min(z)>=0
  assert all(sum((m[j] for j,Cj in enumerate(C) if i in Cj),F())>=h[i] for i in range(N));assert all(sum((z[i] for i in Cj),F())<=1 for Cj in C)
  ml=sum((h[i]*z[i] for i in range(N)),F());mu=sum(m,F());assert 0<ml<=mu and ml==F(c['M_lower']) and mu==F(c['M_upper'])
  qg=list(map(F,c['global_proposed_q']));neg=sum(x<0 for x in qg);pg=sum((h[i]*qg[i] for i in range(N)),F());Lg=max(sum((qg[i] for i in Cj),F()) for Cj in C);assert Lg==F(c['global_proposed_load']);feasible=neg==0 and pg==U
  if c['global_improves']:assert feasible and Lg<L
  else:assert not (feasible and Lg<L)
  rows.append({'case':name,'case_sha256':sha(blob),'raw_sha256':sha(raw),'prior_sha256':sha(pa),'primal_dual_feasible':True,'mass_zero_gap':ml==mu,'relative_mass_gap':float((mu-ml)/mu),'M_lower':str(ml),'M_upper':str(mu),'global_proposal_feasible':feasible,'global_negative_coordinates':neg,'global_payment_equal':pg==U,'global_improves':bool(c['global_improves']),'relative_load_improvement':float((L-Lg)/L) if feasible else None})
 except Exception as e:errors.append({'file':f.name,'error':repr(e)})
summary={'cases':len(rows),'errors':errors,'zero_gap':sum(x['mass_zero_gap'] for x in rows),'nonzero_gap':sum(not x['mass_zero_gap'] for x in rows),'global_feasible':sum(x['global_proposal_feasible'] for x in rows),'global_improves':sum(x['global_improves'] for x in rows),'cases_with_negative_coordinates':sum(bool(x['global_negative_coordinates']) for x in rows),'registration_sha256':sha(regraw),'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'rows':rows}
Path(__file__).with_name('lp_batch.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2))

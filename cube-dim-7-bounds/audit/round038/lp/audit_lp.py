#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,os,time,resource
sys.set_int_max_str_digits(0)
st=time.monotonic();cpu=time.process_time()
src=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round038/compute');prior=Path(__file__).resolve().parent.parent/'phasea/results'
files=['lp_cert_phasea_C_2_05_full_disjoint.json','lp_cert_phasea_C_1_06_full_disjoint.json','lp_cert_phasea_A_3_15_partial.json']
def sha(raw):return hashlib.sha256(raw).hexdigest()
out=[]
for name in files:
 blob=(src/name).read_bytes();cert=json.loads(blob);case=cert['source']['file'];caseblob=(src/case).read_bytes();d=json.loads(caseblob);auditblob=(prior/(Path(case).stem+'.audit.json')).read_bytes();audit=json.loads(auditblob)
 assert audit['status']=='pass' and sha(caseblob)==audit['input_sha256']==cert['source']['sha256'];assert len(d['maximal_cliques'])==audit['cliques']
 C=[tuple(c) for c in d['maximal_cliques']];h=list(map(F,d['witness_budgets']));m=list(map(F,cert['primal_m']));z=list(map(F,cert['dual_y']));N=len(h);assert len(m)==len(C) and len(z)==N and min(m)>=0 and min(z)>=0
 assert all(sum((m[k] for k,c in enumerate(C) if i in c),F())>=h[i] for i in range(N))
 assert all(sum((z[i] for i in c),F())<=1 for c in C)
 ml=sum((h[i]*z[i] for i in range(N)),F());mu=sum(m,F());assert 0<ml<=mu and ml==F(cert['M_lower']) and mu==F(cert['M_upper'])
 U=F(d['U']);L0=F(d['L_before']);Lf=F(d['L_after']);kl=U/mu;ku=U/ml
 assert kl==F(cert['Lopt_lower']) and ku==F(cert['Lopt_upper']) and L0==F(cert['L0']) and Lf==F(cert['Lfinal'])
 gap=L0-ku-L0/F(10000);failure=2*Lf-L0-ku
 assert len(C)>1 and gap>=0 and failure>0
 out.append({'case':case,'stop':d['stop'],'N':N,'clique_count':len(C),'certificate_sha256':sha(blob),'audited_case_sha256':sha(caseblob),'phasea_independent_result_sha256':sha(auditblob),'M_lower':str(ml),'M_upper':str(mu),'Lopt_lower':str(kl),'Lopt_upper':str(ku),'eligible_margin':str(gap),'H_failure_margin':str(failure),'eligible_margin_rel_L0':float(gap/L0),'failure_margin_rel_L0':float(failure/L0),'checks':{'primal_nonnegative':True,'dual_nonnegative':True,'primal_all_vertices':True,'dual_all_cliques':True,'objective_values':True,'bounds_direction':True,'eligible_gap':True,'strict_H_fail':True,'nonclique_graph':True}})
result={'scope':'three LP certificates, exact feasibility only; no LP solving or old step replay','cases':out,'audit_script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('audit_lp.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps([{'case':x['case'],'stop':x['stop'],'eligible_margin_rel_L0':x['eligible_margin_rel_L0'],'failure_margin_rel_L0':x['failure_margin_rel_L0']} for x in out],indent=2))

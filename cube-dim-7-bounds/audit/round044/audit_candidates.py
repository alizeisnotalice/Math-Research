import json,hashlib,time,resource,os,sys
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0);start=time.monotonic();cpu=time.process_time();base=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper');p=base/'round044/compute';old=base/'round038/compute';prev=Path('/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round038/phasea/results')
def sha(b):return hashlib.sha256(b).hexdigest()
out=[]
for suffix,triple in [('C_1_13_full_disjoint',[1,3,9]),('C_4_07_full_disjoint',[1,13,16])]:
 bf=p/('B_case_phasea_'+suffix+'.json');bb=bf.read_bytes();b=json.loads(bb);cf=p/b['source_file'];cb=cf.read_bytes();c=json.loads(cb);rr=(p/('retry_'+b['source_file'])).read_bytes();raw=(old/c['source_raw']).read_bytes();d=json.loads(raw);pa=(prev/(Path(c['source_raw']).stem+'.audit.json')).read_bytes();v=json.loads(pa)
 assert sha(cb)==b['source_sha256'] and sha(rr)==b['source_retry_sha256'] and sha(raw)==c['source_raw_sha256']==v['input_sha256'] and v['status']=='pass'
 assert c['maximal_cliques']==d['maximal_cliques'] and list(map(F,c['h']))==list(map(F,d['witness_budgets'])) and list(map(F,c['q_initial']))==list(map(F,d['q_after'])) and c['geometry']==d['geometry']
 h=list(map(F,d['witness_budgets']));q=list(map(F,d['q_after']));C=d['maximal_cliques'];U=F(d['U']);q1=list(map(F,b['result']['q_new']));assert b['result']['triple']==triple and len(q1)==len(q) and min(q1)>=0
 assert all(q1[i]==q[i] for i in range(len(q)) if i not in triple);assert sum((h[i]*q[i] for i in range(len(q))),F())==U==sum((h[i]*q1[i] for i in range(len(q))),F())
 L=max(sum((q[i] for i in c),F()) for c in C);L1=max(sum((q1[i] for i in c),F()) for c in C);assert L==F(d['L_after'])==F(c['L_initial']) and L1==F(b['result']['L_new']) and L1<L
 out.append({'case':c['source_raw'],'triple':triple,'old_load':str(L),'new_load':str(L1),'strict_gap':str(L-L1),'relative_gap':float((L-L1)/L),'old_triple_q':[str(q[i]) for i in triple],'new_triple_q':[str(q1[i]) for i in triple],'hashes':{'B':sha(bb),'compact':sha(cb),'retry':sha(rr),'raw':sha(raw),'independent_graph':sha(pa)},'checks':'pass'})
res={'scope':'two exact three-coordinate proposals, not search completeness','cases':out,'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};Path(__file__).with_name('candidates.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

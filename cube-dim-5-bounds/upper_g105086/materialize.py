#!/usr/bin/env python3
"""Materialize the independent partial g105 state; preserve old types by map."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse,hashlib,json,os,time,resource
from scan import SOURCE,SOURCE_SHA,LOCKS,ENC,sha,digest,write
AUDIT081=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-audit-cont/cube-dim-5-bounds/audit_kernel081')
AUDIT_REPORT='76276ff907887dd68d984aaabc7a91677f6c92f607d9e51c3e4d2d749fbd15ce'
AUDIT_MANIFEST='ee9ee4b87aa2b4cae4adfad438265ecc918c45a1fea6dfc91ecb9808e3c0800b'
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True);args=ap.parse_args();out=args.work_dir.resolve()
 assert sha(SOURCE)==SOURCE_SHA and sha(AUDIT081/'report.md')==AUDIT_REPORT and sha(AUDIT081/'manifest.json')==AUDIT_MANIFEST
 summary=json.load(open(out/'summary.json'));val=json.load(open(out/'author_validation.json'));assert val['status']=='author_new_target_zero_H_verification_passed'
 reg={'task':'upper_g105086_materialization','pid':os.getpid(),'source':str(SOURCE),'source_sha256':SOURCE_SHA,'target':'21/200000','epsilon':'1/64',
  'code':str(Path(__file__).resolve()),'code_sha256':sha(__file__),'summary_sha256':sha(out/'summary.json'),'author_validation_sha256':sha(out/'author_validation.json'),
  'source_value_audit081_report_sha256':AUDIT_REPORT,'source_value_audit081_manifest_sha256':AUDIT_MANIFEST,
  'new_state':'96104 new-target accepted certificates,1856 pending original closed boxes,all source roots/splits/180040 cache;source070 leaf type remains in provenance and immutable source',
  'new_H_physical_gradient_kernel_calls':0,'new_log_calls':0}
 write(out/'materialization_registration.json',reg);source=json.load(open(SOURCE));allrows={}
 for b in sorted(out.glob('batch[0-9][0-9][0-9][0-9]')):
  for row in map(json.loads,open(b/'certificates.jsonl')):
   i=row['source_leaf_index'];assert i not in allrows;allrows[i]=row
 assert set(allrows)==set(range(97960));passed=[];failed=[];locations={}
 for i in range(97960):
  row=allrows[i];leaf=source['passed'][i];old=leaf['certificate']
  assert row['node']==leaf['node'] and row['source_certificate_sha256']==digest(old) and row['source_leaf_type']==old.get('mechanism',old.get('quadratic'))
  if row['certificate']['passed']:
   locations[i]={'kind':'passed','index':len(passed)}
   passed.append({'node':row['node'],'certificate':row['certificate'],'source070_leaf_index':i,'source070_leaf_type':row['source_leaf_type'],'source070_certificate_sha256':row['source_certificate_sha256']})
  else:failed.append(row)
 failed.sort(key=lambda r:(F(r['H_interval'][0]),F(r['certificate']['minimum_log_margin']),*r['node']))
 queue=[];pending=[];order=[]
 for j,row in enumerate(failed):
  i=row['source_leaf_index'];locations[i]={'kind':'queue','index':j};queue.append(row['node'])
  pending.append({'node':row['node'],'certificate':row['certificate'],'source070_leaf_index':i,'source070_leaf_type':row['source_leaf_type'],'source070_certificate_sha256':row['source_certificate_sha256']})
  record=source['cache'][row['certificate']['anchor']]
  order.append({'node':row['node'],'source_leaf_index':i,'source_leaf_type':row['source_leaf_type'],'root':source['roots'][row['node'][0]],
   'box_lower':row['box_lower'],'box_upper':row['box_upper'],'anchor':row['certificate']['anchor'],'H_interval':row['H_interval'],
   'parent_cached_H_lower':row['H_interval'][0],'parent_full_frame_minimum_log_margin':row['certificate']['minimum_log_margin'],
   'point':record['point'],'new_target':'21/200000','epsilon':'1/64'})
 assert len(passed)==96104 and len(queue)==1856 and len(locations)==97960
 mapper=out/'state_leaf_mapper.jsonl'
 with open(mapper,'w') as f:
  for row in map(json.loads,open(out/'active_leaf_mapper.jsonl')):
   i=row['source_leaf_index'];f.write(ENC.encode({**row,'new_target':'21/200000','new_state_location':locations[i]})+'\n')
  f.flush();os.fsync(f.fileno())
 write(out/'pending_order089.json',{'source070_sha256':SOURCE_SHA,'target':'21/200000','epsilon':'1/64','ordering':'cached H_lower,failed new-target fullQ margin,root,path','records':order})
 new={'task':'upper_g105086','model':'alpha1','epsilon':'1/64','target':'21/200000','cache_model_epsilon':'1/64',
  'roots':source['roots'],'splits':source['splits'],'split_history_source_target':'1/10000','cache':source['cache'],
  'passed':passed,'queue':queue,'pending_certificates086':pending,'inflight':{},'resume_jobs':{},'active':None,'counterexample':None,
  'status':'partial_g105_cover_waiting_main_review','source070':{'path':str(SOURCE),'sha256':SOURCE_SHA,'leaf_mapper':str(mapper),'leaf_mapper_sha256':sha(mapper)},
  'qualification086':{'new_target':'21/200000','epsilon':'1/64','math':'075 Q_full=I4/16,16corners','accepted':96104,'pending':1856,
   'source_value_audit081':{'report':str(AUDIT081/'report.md'),'report_sha256':AUDIT_REPORT,'manifest_sha256':AUDIT_MANIFEST},
   'last_checkpoint_sha256':summary['last_checkpoint_sha256'],'new_target_author_validation_sha256':sha(out/'author_validation.json'),
   'g105_closed_cover':False,'global_g_or_C5_claim':False}}
 new['fingerprint086']=digest({'requalification_registration_sha256':sha(out/'requalification_registration.json'),'last_checkpoint_sha256':summary['last_checkpoint_sha256'],'state_leaf_mapper_sha256':sha(mapper),'target':'21/200000','source070_sha256':SOURCE_SHA})
 final=out/'final';final.mkdir(exist_ok=True);write(final/'state.json',new)
 # Exact byte reconstruction hash also protects all reused source objects.
 expected_file=hashlib.sha256()
 for chunk in ENC.iterencode(new):expected_file.update(chunk.encode())
 expected_file.update(b'\n');assert expected_file.hexdigest()==sha(final/'state.json')
 assert new['roots']==source['roots'] and new['splits']==source['splits'] and new['cache']==source['cache']
 paths=defaultdict(list)
 for leaf in passed:paths[leaf['node'][0]].append(leaf['node'][1])
 for r,p in queue:paths[r].append(p)
 for r in range(52360):
  ps=sorted(paths[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
 assert sha(SOURCE)==SOURCE_SHA
 result={'status':'materialized_partial_g105_state','pid':os.getpid(),'source070_sha256':SOURCE_SHA,'target':'21/200000','epsilon':'1/64',
  'state_path':str(final/'state.json'),'state_sha256':sha(final/'state.json'),'state_bytes':(final/'state.json').stat().st_size,
  'accepted':96104,'pending':1856,'cache':180040,'splits':45600,'roots':52360,'exact_Kraft':52360,'inflight':0,
  'state_leaf_mapper_sha256':sha(mapper),'pending_order089_sha256':sha(out/'pending_order089.json'),'source_leaf_map_bijection':True,
  'source_cache_roots_split_records_preserved':True,'source_old_leaf_types_preserved_in_map':True,
  'new_H_physical_gradient_kernel_log_calls':0,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
  'g105_closed_cover':False,'global_g_or_C5_claim':False}
 write(final/'summary.json',result);print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()

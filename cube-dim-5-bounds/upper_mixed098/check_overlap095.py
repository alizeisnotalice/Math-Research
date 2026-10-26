#!/usr/bin/env python3
"""Exact (cutoff,point) registry intersection only; never evaluates H."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,os,time
W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098')
D=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-alt-cont/cube-dim-5-bounds/upper_adaptive095')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 began=time.monotonic();reg=json.load(open(D/'registration.json'));assert sha(D/'jobs.json.gz')==reg['jobs_sha256']
 tasks=json.load(open(W/'execution_tasks_194.json'));jobs=json.loads(gzip.decompress((D/'jobs.json.gz').read_bytes()));cutoff=F(1,128)
 ours={tuple(map(F,r['point'])) for r in tasks};assert len(ours)==194 and all(F(r['epsilon'])==cutoff for r in tasks)
 endpoint_times={F(t) for t in reg['time_endpoints']};other={}
 if cutoff in endpoint_times:
  for j in jobs:
   reused={F(t) for t in j['reuse_H']}
   if cutoff not in reused:other[tuple(map(F,j['point']))]=j['key']
 shared=sorted(ours & set(other));result={'status':'exact_registry_intersection_empty' if not shared else 'intersection_requires_main_reuse_decision','pid':os.getpid(),
  'epsilon':'1/128','requests098':len(ours),'095_registered_children':len(jobs),'095_new_H128_point_count':len(other),
  '095_H128_time_is_registered':cutoff in endpoint_times,'095_registration_path':str(D/'registration.json'),'095_registration_sha256':sha(D/'registration.json'),
  '095_jobs_sha256':reg['jobs_sha256'],'095_time_endpoints':reg['time_endpoints'],'098_execution_tasks_sha256':sha(W/'execution_tasks_194.json'),
  'intersection_count':len(shared),'intersection':[{'point':list(map(str,p)),'095_job_key':other[p]} for p in shared],
  'comparison':'exact Fraction cutoff and all5 point coordinates;not display decimals','new_H_gradient_kernel_calls':0,'seconds':time.monotonic()-began}
 (W/'overlap095_verification.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='095_time_endpoints'},indent=2))
if __name__=='__main__':main()

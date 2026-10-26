#!/usr/bin/env python3
"""Freeze a completed, author-verified finite batch; no numerical calls."""
from pathlib import Path
import hashlib,json,subprocess
W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089')
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 c=json.load(open(W/'CURRENT.json'));n=c['batch'];b=W/f'batch{n:04d}';v=json.load(open(W/f'author_validation_batch{n:04d}.json'));s=json.load(open(b/'summary.json'))
 assert s['status']=='finite_prefix_complete' and v['status']=='author_no_H_prefix_verification_passed' and v['state_sha256']==sha(b/'state.json')==sha(b/'recovered_overlay.json')
 pids=v['all_production_pids']+[v['pid']];ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid,ppid,command'],capture_output=True,text=True);assert len(ps.stdout.strip().splitlines())==1
 pm=b/'point_manifest.json'
 if not pm.exists():
  points=[]
  for p in sorted(b.glob('point[0-9][0-9][0-9][0-9].json')):
   x=json.load(open(p));assert x['stage']=='complete';points.append({'path':str(p),'sha256':sha(p),'point':x['point'],'pid':x['pid'],'bits':x['record']['bits'],'record':x['record']})
  assert len(points)==s['new_starts'];pm.write_text(json.dumps({'target':'21/200000','epsilon':'1/64','new_points':len(points),'points':points},sort_keys=True,separators=(',',':'))+'\n')
 data=json.load(open(pm));assert data['new_points']==s['new_starts'] and all(sha(Path(p['path']))==p['sha256'] for p in data['points'])
 files=[W/'registration.json',W/'pre_split_registration.json',W/'pre_split_mapper.json',W/'pre_split_validation.json',W/f'author_validation_batch{n:04d}.json']
 files += [b/name for name in ['registration.json','initial_overlay.json','state.json','journal.jsonl','summary.json','validation.json','recovery_validation.json','point_manifest.json','pool_registration.json']]
 m={'task':f'upper_g105089_batch{n:04d}','source086_sha256':v['source086_sha256'],'target':'21/200000','epsilon':'1/64',
  'production_code_sha256':sha(R/'run.py'),'verification_code_sha256':sha(R/'verify.py'),'files':{str(p):sha(p) for p in files},
  'pre_split':1856,'selected':s['selected'],'new_H_records':s['new_starts'],'cache_hits':s['actual_cache_hits'],'new_passed':s['new_passed'],'new_split':s['new_split'],
  'remaining_queue':s['queue'],'all_cumulative_production_and_current_verifier_pids_os_absent':True,'ps_output':ps.stdout,'g105_closed_cover':not s['queue'],'global_g_or_C5_claim':False}
 path=b/'manifest.json';assert not path.exists();path.write_text(json.dumps(m,sort_keys=True,indent=2)+'\n');(R/f'batch{n:04d}_manifest.json').write_bytes(path.read_bytes())
 print(json.dumps({'batch':n,'manifest_sha256':sha(path),'author_validation_sha256':sha(W/f'author_validation_batch{n:04d}.json'),'state_sha256':sha(b/'state.json'),'new_H_records':s['new_starts'],'queue':s['queue']},indent=2))
if __name__=='__main__':main()

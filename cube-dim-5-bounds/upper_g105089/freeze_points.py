#!/usr/bin/env python3
"""Publish only the completed batch's new raw point records, no calculation."""
from pathlib import Path
import hashlib,json
W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 c=json.load(open(W/'CURRENT.json'));b=W/f"batch{c['batch']:04d}";s=json.load(open(b/'summary.json'));assert s['status']=='finite_prefix_complete'
 pm=b/'point_manifest.json'
 if not pm.exists():
  rows=[]
  for p in sorted(b.glob('point[0-9][0-9][0-9][0-9].json')):
   x=json.load(open(p));assert x['stage']=='complete';rows.append({'path':str(p),'sha256':sha(p),'point':x['point'],'pid':x['pid'],'bits':x['record']['bits'],'record':x['record']})
  assert len(rows)==s['new_starts'];pm.write_text(json.dumps({'target':'21/200000','epsilon':'1/64','new_points':len(rows),'points':rows},sort_keys=True,separators=(',',':'))+'\n')
 print(json.dumps({'batch':c['batch'],'new_records':s['new_starts'],'point_manifest_path':str(pm),'point_manifest_sha256':sha(pm),'registration_sha256':sha(b/'registration.json'),'state_sha256':sha(b/'state.json'),'journal_sha256':sha(b/'journal.jsonl')},indent=2))
if __name__=='__main__':main()

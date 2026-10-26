#!/usr/bin/env python3
"""Frozen zero-kernel signed-basis requalification and tree transformation.

No previous implementation is imported. All numerical work uses Fraction.
The immutable source contains all retired records, retained by source index.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import resource
import sys
import time

SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover052/batch0043/state.json')
SOURCE_HASH='065382625da1645e267025ef7c26d292a29a62b65cb7ab97f49dfa9c06995dbf'
LOCKS={
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/upper_multibasis056/manifest.json':'0f3261b09dabc9557903547c0670125392e2bf106a8f4e7c3d952005b3978ba4',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/upper_multibasis056/report.md':'7aee0b3a6f685f5e722f2abd95cbf303357ef7ef4d4dfe8659c0b34e487c5303',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-gprobe/cube-dim-5-bounds/audit_multibasis060/report.md':'b2d7797bb487472383797399a75c4ee9856d0a0141f2564dcea28898a41dc081',
}
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  while chunk:=f.read(8*1024*1024):h.update(chunk)
 return h.hexdigest()
def digest(obj):
 h=hashlib.sha256()
 for chunk in ENC.iterencode(obj):h.update(chunk.encode())
 return h.hexdigest()
def write(path,obj):
 temp=path.with_suffix(path.suffix+'.tmp')
 with open(temp,'w') as f:
  for chunk in ENC.iterencode(obj):f.write(chunk)
  f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(temp,path)
def key(node):return (node[0],node[1])
def physical(node,roots):
 lo=[F(x,32) for x in roots[node[0]]];hi=[x+F(1,32) for x in lo]
 for b in node[1]:
  widths=[y-x for x,y in zip(lo,hi)];j=max(range(4),key=lambda j:widths[j]);m=(lo[j]+hi[j])/2
  assert b in '01'
  if b=='0':hi[j]=m
  else:lo[j]=m
 return lo,hi
def partition(roots,passed,queue,splits):
 leaves=set();byroot=defaultdict(list)
 for n in [r['node'] for r in passed]+queue:
  k=key(n);assert k not in leaves;leaves.add(k);byroot[k[0]].append(k[1])
 assert len(byroot)==len(roots)
 for r in range(len(roots)):
  paths=sorted(byroot[r]);assert all(not b.startswith(a) for a,b in zip(paths,paths[1:]))
  depth=max(map(len,paths));assert sum(1<<(depth-len(p)) for p in paths)==1<<depth
 # Every retained split is an internal tree node, and its recorded split axis
 # equals deterministic longest-side splitting (ties select first axis).
 splitkeys=set()
 for s in splits:
  r,p=s['node'];assert (r,p) not in splitkeys;splitkeys.add((r,p))
  assert s['children']==[[r,p+'0'],[r,p+'1']]
  assert s['axis']==len(p)%4
 for r,p in leaves|splitkeys:
  if p:assert (r,p[:-1]) in splitkeys
 assert len(leaves)==len(roots)+len(splits)
 return {'antichain':True,'roots_with_Kraft_sum_one':len(roots),'leaf_count':len(leaves),'internal_count':len(splitkeys)}
def qualify(split,state,index):
 c=split['failed_certificate'];assert c and not c['passed']
 lo,hi=physical(split['node'],state['roots']);rec=state['cache'][c['anchor']]
 point=list(map(F,rec['point']));assert point[0]==1 and ','.join(map(str,point))==c['anchor']
 assert all(a<=x<=b for a,x,b in zip(lo,point[1:],hi))
 values=[tuple(map(F,v)) for v in rec['values']];assert len(values)==5
 Hlo,Hhi=values[0];assert F(1,10000)<Hlo<=Hhi
 gd=[tuple(map(F,v)) for v in c['loggradient_directed']];assert len(gd)==4
 for (a,b),(gl,gu) in zip(values[1:],gd):
  assert a<=b
  ratios=[a/Hlo,a/Hhi,b/Hlo,b/Hhi];assert gl<=min(ratios)<=max(ratios)<=gu
 Delta=F(c['Delta_lower']);assert Delta>0
 lower=[a-x for a,x in zip(lo,point[1:])];upper=[b-x for b,x in zip(hi,point[1:])]
 linear=lambda v:sum(min(gl*x,gu*x) for x,(gl,gu) in zip(v,gd))
 old=[]
 for bits in product((0,1),repeat=4):
  v=tuple((upper[j] if bits[j] else lower[j]) for j in range(4))
  q=sum(x*x for x in v)/4-sum(v)**2/36
  old.append((Delta+linear(v)-4096*q,bits))
 om,ob=min(old);assert om==F(c['minimum_log_margin'])<=0 and list(ob)==c['worst_corner_bits']
 sets=[sorted(set([a,F(0),b])) for a,b in zip(lower,upper)];rows=[]
 for v in product(*sets):
  norm=sum(x*x for x in v);q=norm/4-sum(abs(x) for x in v)**2/36
  assert F(5,36)*norm<=q<=F(2,9)*norm
  rows.append((Delta+linear(v)-4096*q,v))
 nm,nv=min(rows);assert nm>=om
 cert={'mechanism':'signed_basis_81','passed':nm>0,'minimum_log_margin':str(nm),
  'worst_offset':list(map(str,nv)),'checked_points':len(rows),
  'quadratic':'psi(v)=sum(v^2)/4-sum(abs(v))^2/36',
  'Delta_lower':c['Delta_lower'],'loggradient_directed':c['loggradient_directed'],'anchor':c['anchor']}
 return {'source_split_index':index,'node':split['node'],'source_failed_certificate_sha256':digest(c),
  'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),'anchor_inside':True,
  'old_minimum_log_margin':str(om),'certificate':cert}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--resume',action='store_true');args=parser.parse_args()
 out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);began=time.monotonic()
 assert sha(SOURCE)==SOURCE_HASH and SOURCE.stat().st_size==743364738
 for p,h in LOCKS.items():assert sha(p)==h
 codehash=sha(__file__);regpath=out/'registration.json'
 if not args.resume:
  assert not regpath.exists()
  reg={'task':'upper_cover061','batch':1,'registered_utc':datetime.now(timezone.utc).isoformat(),
   'parent_state':str(SOURCE),'parent_state_sha256':SOURCE_HASH,'parent_state_bytes':SOURCE.stat().st_size,
   'code_path':str(Path(__file__).resolve()),'code_sha256':codehash,'locks':LOCKS,
   'epsilon':'1/64','target':'1/10000','max_tested_split_parents':2048,'seconds':120,
   'deadline_scope':'qualification loop only; immutable source loading, hashing, transformation and audits separately bounded by finite source size',
   'selection':'sort split parents by (path length,root index,path), skip descendants of accepted ancestors',
   'checkpoint_every_tested_parents':16,'resources':{'workers':1,'threads':1,'kernel_calls':0,'MATLAB_calls':0},
   'zero_new_calls':['H','gradient','kernel','log','exp','Gamma'],'new_kernel_chain_authorized':False}
  write(regpath,reg)
 reg=json.load(open(regpath));assert reg['code_sha256']==codehash and reg['parent_state_sha256']==SOURCE_HASH
 fp=digest({'registration':sha(regpath),'parent':SOURCE_HASH,'code':codehash})
 state=json.load(open(SOURCE));assert state['epsilon']=='1/64' and state['target']=='1/10000'
 assert not state['inflight'] and not state['resume_jobs'] and not state.get('active')
 assert [len(state[k]) for k in ('passed','splits','queue')]==[84802,93920,61478]
 source_audit=partition(state['roots'],state['passed'],state['queue'],state['splits'])
 cache_before=digest(state['cache']);root_before=digest(state['roots'])
 source_field_hashes={k:digest(v) for k,v in state.items() if k not in ('cache','roots','passed','splits','queue')}
 ordered=sorted(range(len(state['splits'])),key=lambda i:(len(state['splits'][i]['node'][1]),*state['splits'][i]['node']))
 accepted=defaultdict(list);rows=[];cursor=0
 checkpoint=out/'checkpoint.json'
 if args.resume:
  cp=json.load(open(checkpoint));assert cp['fingerprint']==fp and cp['stage']!='complete'
  rows=cp['rows'];cursor=cp['cursor']
  for row in rows:
   if row['certificate']['passed']:accepted[row['node'][0]].append(row['node'][1])
 def ancestor(n,strict=False):
  r,p=n;matches=[q for q in accepted[r] if p.startswith(q) and (not strict or len(p)>len(q))]
  assert len(matches)<=1
  return None if not matches else [r,matches[0]]
 def save(stage):write(checkpoint,{'fingerprint':fp,'pid':os.getpid(),'stage':stage,'cursor':cursor,'tested':len(rows),'rows':rows,'elapsed_seconds':time.monotonic()-began})
 save('ready_for_qualification');loopstart=time.monotonic()
 while cursor<len(ordered) and len(rows)<2048 and time.monotonic()-loopstart<120:
  index=ordered[cursor];cursor+=1;s=state['splits'][index]
  if ancestor(s['node']) is not None:continue
  row=qualify(s,state,index);rows.append(row)
  if row['certificate']['passed']:accepted[row['node'][0]].append(row['node'][1])
  if len(rows)%16==0:save('qualification_prefix')
 qualification_seconds=time.monotonic()-loopstart;save('qualification_complete')
 oldpassed=state['passed'];oldqueue=state['queue'];oldsplits=state['splits']
 retired={'passed':[],'queue':[],'splits':[]};newpassed=[];newqueue=[];newsplits=[]
 preserved={}
 for name,items,dest in [('passed',oldpassed,newpassed),('queue',oldqueue,newqueue),('splits',oldsplits,newsplits)]:
  kept=[]
  for i,item in enumerate(items):
   n=item if name=='queue' else item['node'];a=ancestor(n,strict=name!='splits')
   if a is None:dest.append(item);kept.append([i,digest(item)])
   else:retired[name].append({'source_index':i,'node':n,'unique_new_parent':a})
  preserved[name]=digest(kept)
 for row in rows:
  if row['certificate']['passed']:newpassed.append({'node':row['node'],'certificate':row['certificate']})
 state['passed']=newpassed;state['queue']=newqueue;state['splits']=newsplits
 new_audit=partition(state['roots'],newpassed,newqueue,newsplits)
 # Re-read every unchanged object in live output against the retained source objects.
 for name,olditems,newitems in [('passed',oldpassed,newpassed[:len(newpassed)-sum(r['certificate']['passed'] for r in rows)]),('queue',oldqueue,newqueue),('splits',oldsplits,newsplits)]:
  source_kept=[(i,x) for i,x in enumerate(olditems) if ancestor(x if name=='queue' else x['node'],strict=name!='splits') is None]
  assert len(source_kept)==len(newitems)
  assert all(a==b for (_,a),b in zip(source_kept,newitems))
  assert digest([[i,digest(x)] for i,x in source_kept])==preserved[name]
 assert digest(state['cache'])==cache_before and digest(state['roots'])==root_before
 assert {k:digest(state[k]) for k in source_field_hashes}==source_field_hashes
 write(out/'certificates.json',rows);write(out/'retired.json',retired)
 state.update({'task061':'upper_cover061','fingerprint061':fp,'status061':'finite_zero_kernel_conversion_complete',
  'conversion061':{'source_state':str(SOURCE),'source_state_sha256':SOURCE_HASH,'source_fingerprint':state['fingerprint'],
   'registration_sha256':sha(regpath),'certificates_sha256':sha(out/'certificates.json'),'retired_sha256':sha(out/'retired.json'),
   'tested_parents':len(rows),'accepted_parents':sum(r['certificate']['passed'] for r in rows),'cursor':cursor,
   'cache_sha256':cache_before,'new_kernel_chain_started':False}})
 write(out/'state.json',state)
 assert sha(SOURCE)==SOURCE_HASH
 for p,h in LOCKS.items():assert sha(p)==h
 validation={'status':'passed','source':source_audit,'new':new_audit,'cache_entries':len(state['cache']),
  'cache_canonical_sha256_before':cache_before,'cache_canonical_sha256_after':digest(state['cache']),
  'roots_canonical_sha256':root_before,'unchanged_node_index_digest':preserved,'source_top_fields_preserved':True,
  'retired_counts':{k:len(v) for k,v in retired.items()},'every_retired_queue_unique_accepted_ancestor':True,
  'zero_new_H_gradient_kernel_log_exp_Gamma_calls':True,'new_kernel_chain_started':False,
  'state_sha256':sha(out/'state.json'),'source_sha256_after':sha(SOURCE),'global_g_or_C5_claim':False}
 write(out/'validation.json',validation)
 summary={'task':'upper_cover061','status':'complete_waiting_main_review','pid':os.getpid(),'workers':[],
  'tested_parents':len(rows),'accepted_parents':sum(r['certificate']['passed'] for r in rows),
  'tested_points':sum(r['certificate']['checked_points'] for r in rows),'unvisited_split_order_entries':len(ordered)-cursor,
  'passed':len(newpassed),'split':len(newsplits),'queue':len(newqueue),'retired':validation['retired_counts'],
  'qualification_seconds':qualification_seconds,'wall_seconds':time.monotonic()-began,
  'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'kernel_calls':0,
  'minimum_accepted_log_margin':str(min(F(r['certificate']['minimum_log_margin']) for r in rows if r['certificate']['passed'])),
  'state_sha256':validation['state_sha256'],'state_path':str(out/'state.json'),'exit_code':0}
 write(out/'summary.json',summary);save('complete');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()

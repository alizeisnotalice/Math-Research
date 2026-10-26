#!/usr/bin/env python3
"""Exact H-ball scan and closed-leaf mapper; never evaluates H or gradients."""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from fractions import Fraction as F
import argparse,hashlib,json,os,resource,time
from pathlib import Path
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover070/final/state.json')
SOURCE_SHA='8727af764ea088b1c30522f0204d6630ad6e0a4d804c7698c5ed2d646cedf343'
TARGET=F(21,200000)
LOCKS={
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover070/manifest.json':'5f8caa2fe34775f1dfec518f484870828054966d10fc13575b28a211a8c5fdd2',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/report.md':'3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/manifest.json':'1d5e270eb1664d714ec80f0b371420d40028b24c3943196b9ab41a0b3de1f9d1'}
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(obj):
 h=hashlib.sha256()
 for s in ENC.iterencode(obj):h.update(s.encode())
 return h.hexdigest()
def write(path,obj):
 temp=path.with_suffix(path.suffix+'.tmp')
 with open(temp,'w') as f:
  for s in ENC.iterencode(obj):f.write(s)
  f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(temp,path)
def physical(node,roots):
 r,p=node;lo=[];hi=[]
 for j,k in enumerate(roots[r]):
  bits=p[j::4];n=len(bits);v=int(bits,2) if bits else 0
  x=F(k*(1<<n)+v,32*(1<<n));lo.append(x);hi.append(x+F(1,32*(1<<n)))
 return lo,hi
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 assert not (out/'scan_registration.json').exists();assert sha(SOURCE)==SOURCE_SHA
 for p,h in LOCKS.items():assert sha(p)==h
 reg={'task':'upper_g105086_scan','registered_utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),
  'source':str(SOURCE),'source_sha256':SOURCE_SHA,'code':str(Path(__file__).resolve()),'code_sha256':sha(__file__),'locks':LOCKS,
  'epsilon':'1/64','new_target':'21/200000','source_cache_count':180040,'source_active_leaf_count':97960,
  'finite_scan_counts':{'cache_H_balls':180040,'active_leaf_H_balls':97960,'root_Kraft':52360},
  'ordering':'active leaf cached H_lower ascending, root/path ascending; exact Fraction comparisons',
  'stop':'any active leaf anchor H_upper<new_target is a strict fixed-H target obstruction; no new H or gradients authorized',
  'resources':{'processes':1,'workers':0,'wall_alarm':False,'new_H_gradient_kernel_calls':0,'new_log_exp_Gamma_calls':0}}
 write(out/'scan_registration.json',reg);state=json.load(open(SOURCE));assert len(state['cache'])==180040 and len(state['passed'])==97960
 assert not state['queue'] and not state['inflight'] and not state['resume_jobs'];roots=state['roots'];assert len(roots)==52360
 cached={};cacheminlo=None;cacheminhi=None;cache_obstacles=[]
 for anchor,record in state['cache'].items():
  lo,hi=map(F,record['values'][0]);assert 0<lo<=hi;assert ','.join(record['point'])==anchor
  point=list(map(F,record['point']));assert point[0]==1 and all(0<=x<=1 for x in point[1:])
  cached[anchor]=(lo,hi)
  if cacheminlo is None or (lo,anchor)<cacheminlo:cacheminlo=(lo,anchor)
  if cacheminhi is None or (hi,anchor)<cacheminhi:cacheminhi=(hi,anchor)
  if hi<TARGET:cache_obstacles.append({'anchor':anchor,'H_interval':record['values'][0],'point':record['point'],'bits':record['bits']})
 leaves=defaultdict(list);mapping=[];active_obstacles=[];types=Counter();seen=set()
 for index,leaf in enumerate(state['passed']):
  node=leaf['node'];assert tuple(node) not in seen;seen.add(tuple(node));r,p=node;leaves[r].append(p)
  cert=leaf['certificate'];anchor=cert['anchor'];record=state['cache'][anchor];hl,hu=cached[anchor];lo,hi=physical(node,roots)
  point=list(map(F,record['point']));assert all(a<=x<=b for a,x,b in zip(lo,point[1:],hi))
  typ=cert.get('mechanism',cert.get('quadratic'));types[typ]+=1
  row={'source_leaf_index':index,'node':node,'root':roots[r],'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),
   'anchor':anchor,'point':record['point'],'H_interval':record['values'][0],'cached_bits':record['bits'],
   'source_leaf_type':typ,'source_certificate_sha256':digest(cert)}
  mapping.append(row)
  if hu<TARGET:active_obstacles.append(row)
 for r in range(52360):
  ps=sorted(leaves[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
 splitkeys={tuple(s['node']) for s in state['splits']};assert len(splitkeys)==45600
 for node in seen:
  for d in range(len(node[1])):assert (node[0],node[1][:d]) in splitkeys
 assert len(state['passed'])==52360+len(state['splits'])
 mapper=out/'active_leaf_mapper.jsonl'
 with open(mapper,'w') as f:
  for row in mapping:f.write(ENC.encode(row)+'\n')
  f.flush();os.fsync(f.fileno())
 order=sorted(range(len(mapping)),key=lambda i:(cached[mapping[i]['anchor']][0],*mapping[i]['node']))
 write(out/'risk_order.json',{'source_sha256':SOURCE_SHA,'mapper_sha256':sha(mapper),'source_leaf_indices':order})
 minlo=min(mapping,key=lambda r:(cached[r['anchor']][0],*r['node']));minhi=min(mapping,key=lambda r:(cached[r['anchor']][1],*r['node']))
 def cache_witness(pair):
  value,anchor=pair;record=state['cache'][anchor]
  return {'value':str(value),'anchor':anchor,'point':record['point'],'H_interval':record['values'][0],'bits':record['bits']}
 result={'status':'strict_fixed_H_target_obstruction' if active_obstacles else 'no_active_anchor_upper_obstruction','pid':os.getpid(),
  'target':'21/200000','cache_scanned':len(cached),'active_leaves_scanned':len(mapping),'distinct_active_anchors':len({r['anchor'] for r in mapping}),
  'cache_minimum_H_lower':cache_witness(cacheminlo),'cache_minimum_H_upper':cache_witness(cacheminhi),
  'active_minimum_H_lower':minlo,'active_minimum_H_upper':minhi,
  'cache_strict_upper_obstacles':len(cache_obstacles),'active_strict_upper_obstacles':len(active_obstacles),
  'active_H_lower_not_above_target':sum(cached[r['anchor']][0]<=TARGET for r in mapping),'source_leaf_types':dict(types),
  'roots_exact_Kraft':52360,'mapper_sha256':sha(mapper),'risk_order_sha256':sha(out/'risk_order.json'),
  'new_H_gradient_kernel_log_exp_Gamma_calls':0,'source_sha256_after':sha(SOURCE),'wall_seconds':time.monotonic()-began,
  'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'global_g_or_C5_claim':False}
 write(out/'active_upper_obstacles.json',active_obstacles);write(out/'cache_upper_obstacles.json',cache_obstacles);write(out/'diagnostics.json',result)
 for p,h in LOCKS.items():assert sha(p)==h
 print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()

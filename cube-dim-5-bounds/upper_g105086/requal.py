#!/usr/bin/env python3
"""Finite 2048-case zero-H blocks for the new g=21/200000 target only."""
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
from pathlib import Path
import argparse,hashlib,json,os,resource,time
from flint import arb,fmpq,ctx
from scan import SOURCE,SOURCE_SHA,TARGET,LOCKS,ENC,sha,digest,write,physical

def ab(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def lower(x):return F(str(x.lower().fmpq()))
def upper(x):return F(str(x.upper().fmpq()))
def ball(pair):
 l,u=map(F,pair);mid=(l+u)/2;rad=(u-l)/2
 return arb(fmpq(mid.numerator,mid.denominator),fmpq(rad.numerator,rad.denominator))
def certify(index,state):
 leaf=state['passed'][index];old=leaf['certificate'];record=state['cache'][old['anchor']]
 values=[tuple(map(F,p)) for p in record['values']];hl,hu=values[0];assert TARGET<hl<=hu
 ctx.prec=512;hball=ball(record['values'][0]);gradient=[]
 for pair,(a,b) in zip(record['values'][1:],values[1:]):
  ratio=ball(pair)/hball;gl,gu=lower(ratio),upper(ratio)
  assert gl<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=gu
  gradient.append((gl,gu))
 delta=lower((ab(hl)/ab(TARGET)).log());assert delta>0
 lo,hi=physical(leaf['node'],state['roots']);point=list(map(F,record['point']));assert point[0]==1
 terms=[]
 for j,(gl,gu) in enumerate(gradient):
  v=(lo[j]-point[j+1],hi[j]-point[j+1]);assert v[0]<=0<=v[1]
  terms.append([min(gl*x,gu*x)-256*x*x for x in v])
 rows=[(delta+sum(terms[j][bits[j]] for j in range(4)),bits) for bits in product((0,1),repeat=4)]
 margin,bits=min(rows);assert margin==delta+sum(min(v) for v in terms)
 cert={'mechanism':'full_frame_16_new_target_g105','target':'21/200000','epsilon':'1/64','quadratic':'Q_full=I4/16',
  'anchor':old['anchor'],'Delta_lower':str(delta),'loggradient_directed':[[str(a),str(b)] for a,b in gradient],
  'minimum_log_margin':str(margin),'worst_corner_bits':list(bits),'checked_points':16,'passed':margin>0}
 return {'source_leaf_index':index,'node':leaf['node'],'source_leaf_type':old.get('mechanism',old.get('quadratic')),
  'source_certificate_sha256':digest(old),'H_interval':record['values'][0],'cached_bits':record['bits'],
  'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),'certificate':cert}

def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True);ap.add_argument('--blocks',type=int,default=1);args=ap.parse_args();out=args.work_dir.resolve()
 assert 1<=args.blocks<=48 and sha(SOURCE)==SOURCE_SHA
 for p,h in LOCKS.items():assert sha(p)==h
 scanreg=json.load(open(out/'scan_registration.json'));diag=json.load(open(out/'diagnostics.json'));risk=json.load(open(out/'risk_order.json'))
 assert scanreg['source_sha256']==SOURCE_SHA and diag['active_strict_upper_obstacles']==0 and diag['cache_strict_upper_obstacles']==0
 assert diag['active_H_lower_not_above_target']==0 and sha(out/'active_leaf_mapper.jsonl')==diag['mapper_sha256']==risk['mapper_sha256']
 order=risk['source_leaf_indices'];assert len(order)==97960 and set(order)==set(range(97960))
 codehash=sha(__file__);helper=sha(Path(__file__).with_name('scan.py'));runreg=out/'requalification_registration.json'
 if not runreg.exists():
  write(runreg,{'task':'upper_g105086_requalification','registered_utc':datetime.now(timezone.utc).isoformat(),
   'source':str(SOURCE),'source_sha256':SOURCE_SHA,'code_path':str(Path(__file__).resolve()),'code_sha256':codehash,'scan_py_sha256':helper,
   'locks':LOCKS,'diagnostics_sha256':sha(out/'diagnostics.json'),'mapper_sha256':diag['mapper_sha256'],'risk_order_sha256':sha(out/'risk_order.json'),
   'target':'21/200000','epsilon':'1/64','curvature':'075 Q_full=I4/16','log_gradient_precision_bits':512,
   'finite_source_active_leaves':97960,'finite_block_size':2048,'maximum_blocks':48,
   'algorithm':'fresh outward log(H_lower/g105) lower endpoint and outward physical-gradient/H ball ratios; all16 exact Fraction corner margins; no g100 replay',
   'preservation':'old leaf types/certificates/cache/source remain unchanged; each new certificate separate and source-indexed',
   'resources':{'processes':1,'workers':0,'wall_alarm':False,'new_H_gradient_kernel_calls':0},
   'checkpoint':'every64 rows fsync certificate-prefix and record exact byteoffset+SHA; every complete2048 block chains parent checkpoint and row-file SHA',
   'result_scope':'uncertified boxes are certificate failures, not H counterexamples; no new kernels permitted'})
 rr=json.load(open(runreg));assert rr['code_sha256']==codehash and rr['scan_py_sha256']==helper
 state=json.load(open(SOURCE));assert len(state['passed'])==97960 and len(state['cache'])==180040 and not state['queue']
 n=0;prev=None;cursor=0
 while (out/f'batch{n+1:04d}'/'checkpoint.json').exists():
  b=out/f'batch{n+1:04d}';cp=json.load(open(b/'checkpoint.json'));br=json.load(open(b/'registration.json'))
  assert cp['parent_checkpoint_sha256']==prev and cp['registration_sha256']==sha(b/'registration.json') and cp['certificates_sha256']==sha(b/'certificates.jsonl')
  assert br['selected_source_leaf_indices']==order[cursor:cursor+2048] and br['start_rank']==cursor
  cursor=cp['next_rank'];prev=sha(b/'checkpoint.json');n+=1
 attempt=out/f'program_{os.getpid()}_registration.json'
 write(attempt,{'pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),'run_registration_sha256':sha(runreg),
  'starting_complete_blocks':n,'starting_rank':cursor,'parent_checkpoint_sha256':prev,'finite_max_blocks_this_program':args.blocks,
  'finite_max_cases_this_program':min(97960-cursor,args.blocks*2048),'new_H_gradient_kernel_calls':0,'wall_alarm':False})
 for step in range(args.blocks):
  if cursor==len(order):break
  b=out/f'batch{n+1:04d}';b.mkdir(exist_ok=True);selected=order[cursor:cursor+2048];brpath=b/'registration.json'
  br={'task':'upper_g105086','block':n+1,'source_sha256':SOURCE_SHA,'run_registration_sha256':sha(runreg),
   'code_sha256':codehash,'parent_checkpoint_sha256':prev,'start_rank':cursor,'selected_source_leaf_indices':selected,'max_cases':2048,'target':'21/200000'}
  if brpath.exists():assert json.load(open(brpath))==br
  else:write(brpath,br)
  prefix=b/'prefix_checkpoint.json';temp=b/'certificates.jsonl.tmp';completed=0;passed=0;failed=0;minimum=None
  hasher=hashlib.sha256()
  if prefix.exists():
   cp=json.load(open(prefix));assert cp['registration_sha256']==sha(brpath);completed=cp['completed'];passed=cp['passed'];failed=cp['failed'];minimum=F(cp['minimum_margin']) if cp['minimum_margin'] else None
   with open(temp,'r+b') as f:
    data=f.read(cp['byteoffset']);assert hashlib.sha256(data).hexdigest()==cp['prefix_sha256'];hasher.update(data);f.truncate(cp['byteoffset'])
  blockstart=time.monotonic()
  with open(temp,'ab') as f:
   def saveprefix():
    f.flush();os.fsync(f.fileno());write(prefix,{'pid':os.getpid(),'registration_sha256':sha(brpath),'completed':completed,'passed':passed,'failed':failed,
     'minimum_margin':str(minimum) if minimum is not None else None,'byteoffset':f.tell(),'prefix_sha256':hasher.hexdigest()})
   saveprefix()
   try:
    for pos in range(completed,len(selected)):
     row=certify(selected[pos],state);data=(ENC.encode(row)+'\n').encode();f.write(data);hasher.update(data);completed+=1
     c=row['certificate'];passed+=c['passed'];failed+=not c['passed'];m=F(c['minimum_log_margin']);minimum=m if minimum is None else min(minimum,m)
     if completed%64==0:saveprefix()
    saveprefix()
   except BaseException:
    saveprefix();raise
  os.replace(temp,b/'certificates.jsonl');assert sha(b/'certificates.jsonl')==hasher.hexdigest();cursor+=len(selected)
  summary={'block':n+1,'pid':os.getpid(),'tested':completed,'passed':passed,'uncertified':failed,'corners':completed*16,
   'minimum_margin':str(minimum),'start_rank':br['start_rank'],'next_rank':cursor,'seconds':time.monotonic()-blockstart,'new_H_gradient_kernel_calls':0}
  write(b/'summary.json',summary);cp={'block':n+1,'parent_checkpoint_sha256':prev,'registration_sha256':sha(brpath),
   'certificates_sha256':sha(b/'certificates.jsonl'),'summary_sha256':sha(b/'summary.json'),'next_rank':cursor,'tested':completed,'passed':passed,'uncertified':failed}
  write(b/'checkpoint.json',cp);prev=sha(b/'checkpoint.json');n+=1;write(out/'CURRENT.json',{'block':n,'next_rank':cursor,'checkpoint_sha256':prev})
  print(json.dumps({'event':'block_complete','directory':str(b),'checkpoint_sha256':prev,**summary}),flush=True)
 if cursor==len(order):
  failedpath=out/'uncertified_boxes.jsonl';statuspath=out/'all_leaf_status.jsonl';counts={'passed':0,'uncertified':0};seen=set();minpassed=None
  with open(failedpath,'w') as ff,open(statuspath,'w') as sf:
   for k in range(1,n+1):
    b=out/f'batch{k:04d}'
    for line in open(b/'certificates.jsonl'):
     row=json.loads(line);index=row['source_leaf_index'];assert index not in seen;seen.add(index);c=row['certificate'];m=F(c['minimum_log_margin'])
     status={'source_leaf_index':index,'node':row['node'],'source_leaf_type':row['source_leaf_type'],'passed':c['passed'],'minimum_log_margin':c['minimum_log_margin'],'certificate_block':k}
     sf.write(ENC.encode(status)+'\n')
     if c['passed']:counts['passed']+=1;minpassed=m if minpassed is None else min(minpassed,m)
     else:counts['uncertified']+=1;ff.write(ENC.encode(row)+'\n')
  assert seen==set(range(97960)) and sum(counts.values())==97960
  final={'status':'all_source_leaves_requalified_for_g105','pid':os.getpid(),'source_sha256':SOURCE_SHA,'target':'21/200000','epsilon':'1/64','blocks':n,
   'tested_source_leaves':97960,'corners':97960*16,**counts,'minimum_passed_margin':str(minpassed),'last_checkpoint_sha256':prev,
   'uncertified_boxes_path':str(failedpath),'uncertified_boxes_sha256':sha(failedpath),'all_leaf_status_sha256':sha(statuspath),
   'mapper_sha256':diag['mapper_sha256'],'root_mapping_complete':True,'source_leaf_types_preserved':True,
   'closed_root_cover_for_new_target':counts['uncertified']==0,'new_H_gradient_kernel_calls':0,'source_sha256_after':sha(SOURCE),'global_g_or_C5_claim':False}
  write(out/'summary.json',final);print(json.dumps({'event':'all_complete',**final}),flush=True)
 assert sha(SOURCE)==SOURCE_SHA
 for p,h in LOCKS.items():assert sha(p)==h
 write(out/f'program_{os.getpid()}_receipt.json',{'pid':os.getpid(),'exit_code':0,'completed_blocks':n,'next_rank':cursor,'wall_seconds':time.monotonic()-began,
  'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'new_H_gradient_kernel_calls':0})
if __name__=='__main__':main()

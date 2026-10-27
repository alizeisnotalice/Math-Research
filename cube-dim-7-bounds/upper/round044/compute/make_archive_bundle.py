#!/usr/bin/env python3
"""Compact exact R044 result strings; no solver calls or new mathematics."""
import os,sys,time,json,hashlib,gzip,collections,resource
from pathlib import Path
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent;BASE=D.parent.parent/'round038/compute';start=time.monotonic();cpu=time.process_time();sha=lambda b:hashlib.sha256(b).hexdigest();canon=lambda x:json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()
sum0=json.loads((D/'summary.json').read_bytes());cases=[]; rawmap={}
for rec in sum0['results']:
 p=D/rec['file'];raw=p.read_bytes();assert sha(raw)==rec['sha256'];r=json.loads(raw);rrp=D/('retry_'+p.name);rrraw=rrp.read_bytes();rr=json.loads(rrraw);assert rr['source_sha256']==sha(raw)
 bp=D/('B_'+p.name);cp=D/('C_'+p.name);B=json.loads(bp.read_bytes()) if bp.exists() else None;C=json.loads(cp.read_bytes()) if cp.exists() else None
 if B:assert B['source_sha256']==sha(raw)
 if C:assert C['source_sha256']==sha(raw)
 v={k:r[k] for k in ('source_raw','source_raw_sha256','source_final_q_sha256','n','a','U','q_initial','L_initial','mass_primal','mass_dual','M_lower','M_upper','global_proposed_q','global_proposed_load','global_improves')}
 v['triple_retry']={'status_counts':dict(collections.Counter(x['status'] for x in rr['triple_logs'])),'best_triple':rr['best_triple'],'best_q':rr['best_q'],'L_best':rr['L_best'],'strict_improvement':rr['strict_improvement']}
 v['B']={'status':B['stop'],'attempted':B['attempted'],'candidate':B['result']} if B else {'status':'not_in_B','attempted':0,'candidate':None}
 v['C']={k:C[k] for k in ('solver_status','status','node_count_numerical','numerical_dual_bound','numerical_mip_gap','numerical_support','exact_candidate')} if C else {'status':'missing_C','solver_status':None,'exact_candidate':None}
 v['raw_hashes']={'initial':sha(raw),'retry':sha(rrraw),'B':sha(bp.read_bytes()) if B else None,'C':sha(cp.read_bytes()) if C else None};cases.append(v)
cases.sort(key=lambda x:x['source_raw']);assert len(cases)==76
payload={'format':'R044 compact v1','R038_compact_sha256':sha((BASE/'phasea_compact.json.gz').read_bytes()),'cases':cases};body=canon(payload)
with (D/'r044_compact.json.gz').open('wb') as f:
 with gzip.GzipFile(filename='',fileobj=f,mode='wb',compresslevel=9,mtime=0) as z:z.write(body)
registrations=['registration.json','triple_retry_registration.json','B_registration.json','C_registration.json'];summaries=['summary.json','triple_retry_summary.json','B_summary.json','C_summary.json'];monitors=['resource_monitor.jsonl','retry_resource_monitor.jsonl','B_monitor.jsonl','C_monitor.jsonl'];producers=['r044_batch.py','r044_triple_retry.py','r044_b.py','r044_c.py']
out={'bundle_sha256':sha((D/'r044_compact.json.gz').read_bytes()),'content_sha256':sha(body),'source_R038_compact_sha256':payload['R038_compact_sha256'],'case_count':76,'files_sha256':{n:sha((D/n).read_bytes()) for n in registrations+summaries+monitors+producers},'scope':'compact extracted exact strings; omitted candidate grids and numeric proposals; raw preserved separately; no new computation','source_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(D/'archive_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 76 cases uncompressed',len(body),'gzip',(D/'r044_compact.json.gz').stat().st_size)

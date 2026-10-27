#!/usr/bin/env python3
"""Extract numerical extreme observations; no quadrature rerun."""
import os,json,hashlib,time,gzip,resource,decimal
from pathlib import Path
D=Path(__file__).resolve().parent;decimal.getcontext().prec=100;Dec=decimal.Decimal;sha=lambda b:hashlib.sha256(b).hexdigest();start=time.monotonic();cpu=time.process_time();cases=[]
for p in sorted(D.glob('case_m*_d*.json')):
 raw=p.read_bytes();r=json.loads(raw);rows=r['results'];assert len(rows)==4096 and r['status']=='completed'
 cases.append({'m':r['m'],'d':r['d'],'raw_filename':p.name,'raw_sha256':sha(raw),'status':r['status'],'triples_completed':r['triples_completed'],'cross_integral_calls':r['cross_integral_calls'],'quadrature_parts':r['quadrature_parts'],'pid':r['pid'],'wall_seconds':r['wall_seconds'],'cpu_seconds':r['cpu_seconds'],'peak_rss_bytes':r['peak_rss_bytes'],'min_observation':min(rows,key=lambda z:Dec(z['S_cross_over_n2'])),'max_observation':max(rows,key=lambda z:Dec(z['S_cross_over_n2']))})
assert len(cases)==15
payload={'format':'R047 sampled-case-summary v1','cases':cases,'scope':'extreme rows selected from frozen numerical raw; no independent proof of maximality over omitted rows'};body=json.dumps(payload,sort_keys=True,separators=(',',':')).encode()
with (D/'archive_summary.json.gz').open('wb') as f:
 with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0,compresslevel=9) as z:z.write(body)
files=['registration.json','summary.json','monitor.jsonl','cross_face.py','micro.json','check_moments.py','check_moments.json','r047_primary_matlab_aggregation.m','r047_primary_matlab_aggregation.json'];out={'bundle_sha256':sha((D/'archive_summary.json.gz').read_bytes()),'content_sha256':sha(body),'raw_case_hashes':{c['raw_filename']:c['raw_sha256'] for c in cases},'source_files_sha256':{x:sha((D/x).read_bytes()) for x in files},'source_sha256':sha(Path(__file__).read_bytes()),'case_count':15,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'matlab_provenance':'MATLAB script and receipt copied unchanged from coordinator; no MATLAB call made here','scope':'derived summaries from frozen raw, no quadrature'};(D/'archive_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 15 summaries gzip bytes',(D/'archive_summary.json.gz').stat().st_size)

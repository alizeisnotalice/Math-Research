#!/usr/bin/env python3
"""Extract one-case worst observations, never rerun integration."""
import os,json,hashlib,time,gzip,resource,decimal
from fractions import Fraction as Q
from pathlib import Path
D=Path(__file__).resolve().parent;decimal.getcontext().prec=100;Dec=decimal.Decimal;sha=lambda b:hashlib.sha256(b).hexdigest();start=time.monotonic();cpu=time.process_time()
def qd(s):
 x=Q(s);return Dec(x.numerator)/Dec(x.denominator)
cases=[]
for p in sorted(D.glob('case_g*_rho*.json')):
 raw=p.read_bytes();r=json.loads(raw);g=qd(r['g']);rho=qd(r['rho'])
 def deficit(z):return -Dec(z['mean_minus_one'])*rho/(g*g*qd(z['h_over_g']))
 def moment(z):return Dec(z['second_moment'])*rho**2/g**3
 rows=r['results'];assert len(rows)==256 and r['status']=='completed'
 chosen={'worst_mean_deficit':max(rows,key=deficit),'worst_second_moment':max(rows,key=moment),'worst_direct_error':max(rows,key=lambda z:Dec(z['direct_abs_error'])),'largest_mean_minus_one':max(rows,key=lambda z:Dec(z['mean_minus_one']))}
 cases.append({'g':r['g'],'rho':r['rho'],'raw_filename':p.name,'raw_sha256':sha(raw),'status':r['status'],'combos':r['combos_completed'],'distinct_r':r['distinct_r'],'integral_calls':r['integral_calls'],'pid':r['pid'],'wall_seconds':r['wall_seconds'],'cpu_seconds':r['cpu_seconds'],'peak_rss_bytes':r['peak_rss_bytes'],'observations':chosen,'worst_deficit_normalized':str(deficit(chosen['worst_mean_deficit'])),'worst_second_normalized':str(moment(chosen['worst_second_moment']))})
assert len(cases)==15
payload={'format':'R046 sampled-case-summary v1','cases':cases,'scope':'saved numerical worst rows as reported in producer data; not an independent proof of maxima over omitted rows'};body=json.dumps(payload,sort_keys=True,separators=(',',':')).encode()
with (D/'archive_summary.json.gz').open('wb') as f:
 with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0,compresslevel=9) as z:z.write(body)
files=['registration.json','summary.json','monitor.jsonl','r046_row.py','symbolic_check.py','symbolic_check.json'];out={'bundle_sha256':sha((D/'archive_summary.json.gz').read_bytes()),'content_sha256':sha(body),'raw_case_hashes':{c['raw_filename']:c['raw_sha256'] for c in cases},'source_files_sha256':{x:sha((D/x).read_bytes()) for x in files},'source_sha256':sha(Path(__file__).read_bytes()),'case_count':15,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'derived summaries from frozen raw, no integration'};(D/'archive_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 15 summaries gzip bytes',(D/'archive_summary.json.gz').stat().st_size)

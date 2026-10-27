#!/usr/bin/env python3
"""Portable arithmetic check of archived worst observations; no integration rerun."""
import os,json,hashlib,gzip,time,resource,decimal
from fractions import Fraction as Q
from pathlib import Path
D=Path(__file__).resolve().parent;decimal.getcontext().prec=100;Dec=decimal.Decimal;sha=lambda b:hashlib.sha256(b).hexdigest();start=time.monotonic();cpu=time.process_time()
def qd(s):
 x=Q(s);return Dec(x.numerator)/Dec(x.denominator)
manifest=json.loads((D/'archive_manifest.json').read_bytes());packed=(D/'archive_summary.json.gz').read_bytes();assert sha(packed)==manifest['bundle_sha256'];raw=gzip.decompress(packed);assert sha(raw)==manifest['content_sha256'];data=json.loads(raw);reg=json.loads((D/'registration.json').read_bytes());assert len(data['cases'])==len(reg['cases'])==15
assert all(sha((D/x).read_bytes())==v for x,v in manifest['source_files_sha256'].items());assert reg['code_sha256']==sha((D/'r046_row.py').read_bytes());assert json.loads((D/'symbolic_check.json').read_bytes())['status']=='PASS'
valid={(x['g'],x['rho']) for x in reg['cases']}
for c in data['cases']:
 assert (c['g'],c['rho']) in valid and c['raw_sha256']==manifest['raw_case_hashes'][c['raw_filename']] and c['combos']==256 and c['status']=='completed';g=qd(c['g']);rho=qd(c['rho']);obs=c['observations']
 for row in obs.values():
  h=g*qd(row['h_over_g']);r=qd(row['r']);assert h>0 and Dec(1)<=r<=Dec(1)+2*h and all(row['checks'].values())
  assert Dec(row['mean_minus_one'])<=Dec('1e-65') and Dec(row['mean_minus_one'])>=-6*h*g/rho-Dec('1e-65')
  assert Dec(row['second_moment'])<=31*g**3/rho**2+Dec('1e-65')
 a=-Dec(obs['worst_mean_deficit']['mean_minus_one'])*rho/(g*g*qd(obs['worst_mean_deficit']['h_over_g']));b=Dec(obs['worst_second_moment']['second_moment'])*rho**2/g**3
 assert a==Dec(c['worst_deficit_normalized']) and b==Dec(c['worst_second_normalized']) and Dec(obs['worst_direct_error']['direct_abs_error'])<Dec('1e-65')
out={'status':'PASS','case_count':15,'bundle_sha256':manifest['bundle_sha256'],'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'checks stored worst-row arithmetic/registered parameters only; does not independently prove maximality over omitted raw rows or rerun numerical integration'};(D/'archive_verify_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 15 stored case summaries; no raw-row or integral replay')

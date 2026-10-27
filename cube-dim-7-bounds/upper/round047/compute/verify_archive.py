#!/usr/bin/env python3
"""Check archived extremes/identities, no moment integration or new sampling."""
import os,json,gzip,time,hashlib,resource,decimal
from pathlib import Path
from fractions import Fraction as Q
D=Path(__file__).resolve().parent;Dec=decimal.Decimal;sha=lambda b:hashlib.sha256(b).hexdigest();start=time.monotonic();cpu=time.process_time();manifest=json.loads((D/'archive_manifest.json').read_bytes());packed=(D/'archive_summary.json.gz').read_bytes();assert sha(packed)==manifest['bundle_sha256'];raw=gzip.decompress(packed);assert sha(raw)==manifest['content_sha256'];data=json.loads(raw);reg=json.loads((D/'registration.json').read_bytes());assert len(data['cases'])==len(reg['cases'])==15
assert all(sha((D/x).read_bytes())==h for x,h in manifest['source_files_sha256'].items());assert reg['code_sha256']==sha((D/'cross_face.py').read_bytes());assert json.loads((D/'micro.json').read_bytes())['status']=='PASS' and json.loads((D/'check_moments.json').read_bytes())['status']=='PASS'
mat=json.loads((D/'r047_primary_matlab_aggregation.json').read_bytes());assert all(Q(mat['n'+str(i)])>0 for i in range(2,8))
valid={(x['m'],x['d']) for x in reg['cases']};grid=reg['grid']
for c in data['cases']:
 assert (c['m'],c['d']) in valid and c['raw_sha256']==manifest['raw_case_hashes'][c['raw_filename']] and c['triples_completed']==4096 and c['cross_integral_calls']==2176 and c['status']=='completed'
 for z in (c['min_observation'],c['max_observation']):
  assert z['u_over_d'] in grid['u_over_d'] and z['theta1'] in grid['theta1'] and z['theta2'] in grid['theta2'] and Dec(0)<=Dec(z['S_cross_over_n2'])<Dec(4096)
 assert Dec(c['min_observation']['S_cross_over_n2'])<=Dec(c['max_observation']['S_cross_over_n2'])
out={'status':'PASS','case_count':15,'bundle_sha256':manifest['bundle_sha256'],'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'checks saved numerical extreme rows and provenance only; no proof of maximality over omitted raw rows, no integral rerun, no interval bound'};(D/'archive_verify_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 15 stored extrema summaries; no raw-row or integral replay')

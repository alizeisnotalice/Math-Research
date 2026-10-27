#!/usr/bin/env python3
"""Package frozen Phase A outputs without rerunning any computation."""
import json,gzip,hashlib,sys,time,resource,os
from pathlib import Path
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
canonical=lambda obj:json.dumps(obj,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()
start=time.monotonic();cpu=time.process_time();summary=json.loads((D/'phasea_summary.json').read_bytes());assert len(summary['results'])==79
cases=[]
for entry in summary['results']:
 path=D/entry['output'];raw=path.read_bytes();assert sha(raw)==entry['sha256'];r=json.loads(raw)
 v={k:r[k] for k in ('historical_identity','n','a','N','geometry_fingerprint','geometry','generation_retry','witness_budgets','maximal_cliques','graph_nodes','U','q_before','L_before','L_after','stop','pair_evaluations','pid','wall_seconds','cpu_seconds','peak_rss_bytes')}
 v.update({'raw_filename':path.name,'raw_sha256':sha(raw),'accepted_pairs':[[s['i'],s['j']] for s in r['steps']],'q_after_canonical_sha256':sha(canonical(r['q_after'])),'accepted_step_count':len(r['steps'])})
 cases.append(v)
cases.sort(key=lambda v:v['raw_filename']);payload={'format':'R038 compact replay v1','case_count':len(cases),'cases':cases};body=canonical(payload)
with (D/'phasea_compact.json.gz').open('wb') as f:
 with gzip.GzipFile(filename='',fileobj=f,mode='wb',compresslevel=9,mtime=0) as z:z.write(body)
certs=sorted(D.glob('lp_cert_phasea_*.json'));assert len(certs)==3
manifest={'bundle_filename':'phasea_compact.json.gz','compressed_sha256':sha((D/'phasea_compact.json.gz').read_bytes()),'canonical_content_sha256':sha(body),'case_count':len(cases),'step_count':sum(c['accepted_step_count'] for c in cases),'raw_summary_sha256':sha((D/'phasea_summary.json').read_bytes()),'registration_sha256':sha((D/'phasea_registration.json').read_bytes()),'historical_originals':[{'filename':c['raw_filename'],'sha256':c['raw_sha256'],'q_after_canonical_sha256':c['q_after_canonical_sha256']} for c in cases],'three_lp_certificates':[{'filename':c.name,'sha256':sha(c.read_bytes()),'source_case':json.loads(c.read_bytes())['source']['file']} for c in certs],'generator_source_sha256':sha(Path(__file__).read_bytes()),'scope':'compact cache and source mapping only; not a new run or new mathematical search','pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(D/'compact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('PASS compact cases',len(cases),'steps',manifest['step_count'],'raw canonical bytes',len(body),'gzip bytes',(D/'phasea_compact.json.gz').stat().st_size)

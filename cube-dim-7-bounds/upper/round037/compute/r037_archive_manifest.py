#!/usr/bin/env python3
"""Compact historical run ledger; no reexecution of production tasks."""
import json,hashlib,glob,time,os,resource
from pathlib import Path
p=Path(__file__).resolve();root=p.parent;start=time.monotonic();cpu=time.process_time();summary=root/'r037_summary.json';d=json.loads(summary.read_text());files=sorted(root.glob('r037_[ABC]_[0-4].jsonl'));assert len(files)==15
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
cases=[]
for file in files:
 for ln in file.read_text().splitlines():
  r=json.loads(ln);cases.append({k:r[k] for k in ('template','stratum','caseindex','seed_hex','status','pid','case_wall_seconds')});cases[-1]['patterns']=[{'pattern':v['pattern'],'status':v['status']} for v in r.get('patterns',[])]
cases.sort(key=lambda r:(r['template'],r['stratum'],r['caseindex']));assert len(cases)==240 and len({(c['template'],c['stratum'],c['caseindex']) for c in cases})==240
selected=['r037_archive_manifest.sha256','R037_REPORT.md','r037_summary.json','r037_archive_manifest.py','r037_archive_manifest.json','r037_joint.py','joint_geometry.py','joint_atoms.py','joint_core.py','joint_compress.py','r037_failure_witness.json','r037_failure_geometry.json','r037_no_move.json','r037_check_export.json','r037_failure_witness.py','r037_failure_geometry.py','r037_no_move.py','r037_check_export.py']
inputs=['r037_joint.py','joint_geometry.py','joint_atoms.py','joint_core.py','joint_compress.py','r037_failure_witness.json','r037_failure_geometry.json','r037_no_move.json','r037_check_export.json','r037_failure_witness.py','r037_failure_geometry.py','r037_no_move.py','r037_check_export.py']
audit=Path('/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round037/joint_counterexample');auditfiles=['README.md','audit.py','audit.json']
out={'description':'Compact case status ledger and historical source/JSONL digests; no production cases recomputed','historical_run_parent_pid':67010,'historical_run_session':81622,'historical_run_summary_sha256':sha(summary),'historical_source_sha256':d['hashes'],'raw_15_jsonl_sha256':{f.name:sha(f) for f in files},'artifact_sha256':{name:sha(root/name) for name in inputs},'independent_audit_sha256':{name:sha(audit/name) for name in auditfiles},'archival_whitelist_compute':selected,'archival_whitelist_audit':auditfiles,'excluded_from_archive':'15 full per-worker JSONL remain in private clone; compact per-case statuses and immutable hashes retained here; optional preserve externally if future exact q data needed','case_statuses':cases,'manifest_script_sha256':sha(p),'manifest_pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(root/'r037_archive_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 240 status rows 15 original JSONL hashes; no rerun')

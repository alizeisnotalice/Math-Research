from pathlib import Path
import gzip,json,hashlib,time,os,resource,collections
st=time.monotonic();cpu=time.process_time();src=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round070/compute');ar=src/'archive';here=Path(__file__).resolve().parent;sha=lambda b:hashlib.sha256(b).hexdigest();mb=(ar/'archive_manifest.json').read_bytes();m=json.loads(mb)
for f,h in m['sha256'].items():assert sha((ar/f).read_bytes())==h
package=json.loads(gzip.decompress((ar/'original_chunks.json.gz').read_bytes()));assert set(package)=={'pilot','production','width','width_e'}
counts={};
for mode,sn in [('pilot','pilot_summary.json'),('production','production_summary.json'),('width','width_summary.json'),('width_e','width_e_summary.json')]:
 s=json.loads((ar/sn).read_bytes());rs={x['file']:x for x in s['receipts']};rows=package[mode];assert len(rows)==len(rs) and {x['file'] for x in rows}==set(rs)
 for x in rows:
  data=x['original_utf8'].encode();assert sha(data)==x['sha256']==rs[x['file']]['sha256'];assert sha((src/x['file']).read_bytes())==sha(data)
 counts[mode]=len(rows)
out={'scope':'read-only independent archive bytes and receipt completeness, no numerical rerun','manifest_sha256':sha(mb),'archive_members':len(m['sha256']),'embedded_receipts':counts,'all_hashes_match':True,'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'archive_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

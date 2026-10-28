from pathlib import Path
import gzip,json,hashlib,time,resource,os
st=time.monotonic();cpu=time.process_time();src=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round070/compute');ar=src/'frontier4_archive';here=Path(__file__).resolve().parent;sha=lambda b:hashlib.sha256(b).hexdigest();mb=(ar/'archive_manifest.json').read_bytes();manifest=json.loads(mb)
for name,h in manifest['sha256'].items():assert sha((ar/name).read_bytes())==h
summary=json.loads((ar/'frontier4_summary.json').read_bytes());reg=(ar/'frontier4_registration.json').read_bytes();assert summary['registration_sha256']==sha(reg);rs={x['file']:x for x in summary['receipts']};chunks=json.loads(gzip.decompress((ar/'original_chunks.json.gz').read_bytes()));assert len(chunks)==len(rs)==1280 and {x['file'] for x in chunks}==set(rs);top=[];counts={'domain':0,'payment':0,'eligible':0}
for entry in chunks:
 data=entry['original_utf8'].encode();f=entry['file'];assert sha(data)==entry['sha256']==rs[f]['sha256']==sha((src/f).read_bytes());c=json.loads(data);assert c['status']==rs[f]['status']=='completed' and c['counts']==rs[f]['counts'];top.extend(c['top32'])
 for k in counts:counts[k]+=c['counts'][k]
top.sort(key=lambda x:x['cost_float'],reverse=True);assert top[:64]==summary['top64_float_proposals'] and sum(counts.values())==10485760
rc=json.loads((ar/'frontier4_recheck.json').read_bytes());assert rc['screen_summary_sha256']==sha((ar/'frontier4_summary.json').read_bytes()) and rc['source_sha256']==sha((ar/'frontier4_recheck.py').read_bytes())
out={'scope':'archive bytes, recorded counts and float merge only; no strict numerical certificate','manifest_sha256':sha(mb),'manifest_entries':len(manifest['sha256']),'receipt_count':len(chunks),'counts':counts,'top64_float_merge_match':True,'recheck_rows':len(rc['rows']),'checks':'PASS','pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'script_sha256':sha(Path(__file__).read_bytes())}
(here/'frontier_archive_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

"""Independent exact point/time cache reuse and zero overlap check; noH."""
from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,time,os
H=Path(__file__).resolve().parent;B=H.parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def point(p):return tuple(str(F(x)) for x in p)
def main():
 t=time.monotonic();reg=read(H/'registration.json');assert sha(H/'jobs.json.gz')==reg['jobs_sha256'];jobs=read(H/'jobs.json.gz');pts={point(j['point']) for j in jobs};assert len(pts)==len(jobs);old={point(r['point']) for r in read(B/'upper_global080/roots.json.gz')}
 for p in sorted((B/'upper_adapt082').glob('batch*_jobs.json.gz')):old.update(point(r['point']) for r in read(p))
 old.update(point(r['point']) for r in read(B/'upper_adaptive095/jobs.json.gz'));assert len(old)==69545 and not old&pts;gp=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098/final_state.json');assert sha(gp)==reg['source098_final_state_sha256'];cache=read(gp)['cache_by_epsilon'];lookup={ts:{point(r['point']):r for r in records.values()} for ts,records in cache.items()};assert len(lookup['1/64'])==184796 and len(lookup['1/128'])==194;hits={ts:len(pts&set(tab)) for ts,tab in lookup.items()};assert hits=={'1/128':0,'1/64':230};reuse=0
 for j in jobs:
  pk=point(j['point']);expected={ts:lookup[ts][pk] for ts in lookup if pk in lookup[ts]};assert set(expected)==set(j['reuse_H'])
  for ts,r in expected.items():assert j['reuse_H'][ts]=={'lower':r['values'][0][0],'upper':r['values'][0][1]};reuse+=1
 out={'status':'pass','all28427_new_point_keys_unique':True,'old_H73_points':69545,'old_H73_intersection':0,'H64_cache_points':184796,'H128_cache_points':194,'exact_point_time_hits':hits,'all230_reused_intervals_exact_source_copy':True,'same_time_keys_kept_separate':True,'reused_values':reuse,'new_H_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-t,'source_sha256':sha(__file__)};(H/'registry_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()

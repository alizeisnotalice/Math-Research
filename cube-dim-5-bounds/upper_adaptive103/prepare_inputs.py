from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import gzip,json,hashlib,time
import moment085 as M
import ratio091 as R
H=Path(__file__).resolve().parent;OLD=H.parent/'upper_adaptive095'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def gz(p,x):p.write_bytes(gzip.compress(serial(x).encode(),mtime=0))
@lru_cache(maxsize=32)
def geom(offsets):
 lo=[a for a,b in offsets];hi=[b for a,b in offsets];mu=[(a+b)/2 for a,b in offsets];m2=M.max_abs_square_moment(offsets);vol=F(1)
 for a,b in offsets:vol*=b-a
 return {'offsets':[[str(a),str(b)] for a,b in offsets],'M2':str(m2),'mu_infinity':str(max(map(abs,mu))),'h_infinity':str(max(map(abs,lo+hi))),'volume':str(vol)}
def main():
 began=time.monotonic();assert not (H/'jobs.json.gz').exists();rows=json.loads(gzip.decompress((H/'children_plan.json.gz').read_bytes()));h64=json.loads(gzip.decompress((H/'reuse_values.json.gz').read_bytes()));ts=json.loads((OLD/'registration.json').read_text())['time_endpoints'];xs=set()
 for r in rows:
  point=list(map(F,r['point']));L=r['level'];offs=tuple((F(i,2**L)-x,F(i+1,2**L)-x) for i,x in zip(r['indices'],point[1:]));r['geometry']=geom(offs);pk=','.join(r['point']);r['reuse_H']={};r['reuse_sources']=[]
  if r['reuse_H73']:
   p=Path(r['H73_packet_path']);assert sha(p)==r['H73_packet_sha256'];packet=json.loads(gzip.decompress(p.read_bytes()));assert packet['point']==r['point'];r['reuse_H']={e['t']:e['H'] for e in packet['endpoints']};r['reuse_sources'].append({'kind':'completed_H73','path':str(p),'sha256':sha(p)})
  if pk in h64:
   for cache_time,entry in h64[pk].items():
    val=entry['H'];old=r['reuse_H'].get(cache_time)
    if old is not None:
     lo=max(F(old['lower']),F(val['lower']));up=min(F(old['upper']),F(val['upper']));assert lo<=up;val={'lower':str(lo),'upper':str(up)}
    r['reuse_H'][cache_time]=val;r['reuse_sources'].append({'kind':entry['source'],'time':cache_time,'source_reuse_file_sha256':sha(H/'reuse_values.json.gz')})
  xs.update(F(r['geometry']['h_infinity'])/F(a) for a in ts)
 old=json.loads((OLD/'ratio_catalog.json').read_text());missing=sorted(x for x in xs if str(x) not in old);(H/'ratio_registration.json').write_text(json.dumps({'status':'new103_scalar_entries_registered_before_production','used_x':list(map(str,sorted(xs))),'new_x':list(map(str,missing)),'old_catalog_sha256':sha(OLD/'ratio_catalog.json'),'ratio_source_sha256':sha(H/'ratio091.py'),'bits':192},indent=2)+'\n');cat={k:v for k,v in old.items()}
 for x in missing:
  lo,up=R.ratio_interval(x,192);a,b=R._sqrt_interval(1+x*x/16,192);c,d=R._sqrt_interval(1+x*x/256,192);cat[str(x)]={'lower':str(lo),'upper':str(up),'sqrt1_lower':str(a),'sqrt1_upper':str(b),'sqrt2_lower':str(c),'sqrt2_upper':str(d)}
 (H/'ratio_catalog.json').write_text(serial(cat)+'\n');gz(H/'jobs.json.gz',rows);out={'children':len(rows),'unique_offset_geometries':geom.cache_info().currsize,'used_ratio_x':len(xs),'retained_old_catalog_entries':len(old),'new_ratio_entries':len(missing),'catalog_count':len(cat),'H_values_reused':sum(len(r['reuse_H']) for r in rows),'new_H_registered':sum(73-len(r['reuse_H']) for r in rows),'new_H_calls_so_far':0,'jobs_sha256':sha(H/'jobs.json.gz'),'catalog_sha256':sha(H/'ratio_catalog.json'),'seconds':time.monotonic()-began};(H/'input_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()

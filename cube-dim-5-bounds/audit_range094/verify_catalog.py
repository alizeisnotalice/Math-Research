#!/usr/bin/env python3
"""094: integer-only algebraic ratio audit; no producer import, no H."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,gzip,os,time
P=Path(__file__).resolve().parent
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_shift091')
R=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-alt-cont/cube-dim-5-bounds/upper_range092')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 t=Path(str(p)+'.tmp');t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)
def root(q,bits):
 d=1<<bits;n=isqrt((q.numerator*d*d)//q.denominator)
 assert n*n*q.denominator<=q.numerator*d*d<(n+1)**2*q.denominator
 return F(n,d),F(n+1,d)
def main():
 t=time.monotonic();reg=json.loads((P/'catalog_registration.json').read_text())
 assert not (P/'catalog_results.json').exists(),'completed finite catalog refuses replay'
 for r in reg['source_files']:assert sha(r['path'])==r['sha256'],r['path']
 for r in json.loads((S/'manifest.json').read_text())['files']:assert sha(S/r['path'])==r['sha256']
 cat=json.loads((R/'ratio_catalog.json').read_text());nodes=json.loads(gzip.decompress((R/'nodes.json.gz').read_bytes()));prod=json.loads((R/'registration.json').read_text());rr=json.loads((R/'ratio_registration.json').read_text())
 assert len(cat)==89 and len(nodes)==56139
 hs={F(n['h_infinity']) for n in nodes};ts=list(map(F,prod['time_endpoints']))
 assert len(ts)==73 and ts[0]==F(1,512) and ts[-1]==1 and all(a<b for a,b in zip(ts,ts[1:]))
 needed={h/a for h in hs for a in ts};assert needed==set(map(F,cat))
 old={F(1,8),F(1),F(2),F(4),F(8),F(16)};overlap=needed&old
 assert not overlap,'existing cases require certificate consumption, not replay'
 result=[]
 for x in sorted(needed):
  c=cat[str(x)];lo,hi=F(c['lower']),F(c['upper']);s,t1=root(16+x*x,224);v,w=root(256+x*x,224)
  low=(s+x)**2*(v+x)**8/16**9;up=(t1+x)**2*(w+x)**8/16**9
  assert 1<lo<=low<up<=hi,(x,'224 independent enclosure')
  for j,d in [(1,16),(2,256)]:
   a,b=F(c[f'sqrt{j}_lower']),F(c[f'sqrt{j}_upper']);q=1+x*x/d
   assert a*a<=q<b*b and b-a==F(1,1<<192)
  al,au=F(c['sqrt1_lower']),F(c['sqrt1_upper']);bl,bu=F(c['sqrt2_lower']),F(c['sqrt2_upper'])
  directlo=(al+x/4)**2*(bl+x/16)**8;directup=(au+x/4)**2*(bu+x/16)**8
  D=1<<192
  assert lo==F((directlo*D).__floor__(),D)
  assert hi==F((directup*D).__floor__()+1,D)
  result.append({'x':str(x),'independent_lower':str(low),'independent_upper':str(up),'published_upper':str(hi),'pass':True})
  save(P/'catalog_checkpoint.json',{'completed':len(result),'count':89,'x':str(x)})
 out={'status':'pass','pid':os.getpid(),'seconds':time.monotonic()-t,'root_bits':224,'ratios':89,'required_nodes':len(nodes),'distinct_h':list(map(str,sorted(hs))),'time_endpoints':len(ts),'old_six_overlap':[],'new_H_calls':0,'old_exp_log_asinh_calls':0,'Arb_calls':0,'source_registration_sha256':sha(P/'catalog_registration.json'),'code_sha256':sha(__file__),'certificates':result}
 save(P/'catalog_results.json',out);print(json.dumps({k:v for k,v in out.items() if k!='certificates'}))
if __name__=='__main__':main()

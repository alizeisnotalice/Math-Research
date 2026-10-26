#!/usr/bin/env python3
"""Finite, read-only source extraction and exact complete-root ordering."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial,prod,comb
import gzip,json,hashlib,time,os,resource,sys,signal
from flint import arb,fmpq,ctx
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_global080');SRC=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover066/final/state.json');P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_shorttime076');LOCK='42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b'
def sha(p):
 v=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):v.update(b)
 return v.hexdigest()
def rq(x):
 a=F(str(x));return fmpq(a.numerator,a.denominator)
def ep(x,up=False):return str((x.upper() if up else x.lower()).fmpq())
def alarm(sig,frame):raise TimeoutError('prepare590s alarm')
def main():
 began=time.monotonic();signal.signal(signal.SIGALRM,alarm);signal.alarm(590);ctx.prec=384;assert not (H/'roots.json.gz').exists();assert sha(SRC)==LOCK
 assert sha(P/'manifest.json')=='c9221c9b10e8a9c931cb81b015fd19f79572cb52c3af9e9a2238f954f84f1fb5'
 short=json.loads((P/'results.json').read_text());C3=F(short['derived_constants']['C3']['upper']);C1=F(short['derived_constants']['C1']['upper']);piL=F(short['derived_constants']['pi_outward_384_dyadic']['lower']);eps=F(1,512);T=F(1,16)
 P5=2/(arb(rq(piL))**3*(8*arb(2).sqrt())**5);cap=arb(rq(C3/2*(eps**(-2)-T**(-2))))+arb(rq(C1))*arb(rq(T/eps)).log()+P5/(4*arb(rq(T))**4);assert cap<3
 E=next(v['E_norm'] for v in short['cases'] if v['epsilon']=='1/512');scalar={'formula':'C3/2*(eps^-2-T^-2)+C1*log(T/eps)+P5/(4*T^4)','epsilon':str(eps),'T':str(T),'Hcap':'3','computed_cap_upper':ep(cap,True),'computed_cap_display':float(cap),'P5_upper':ep(P5,True),'C3_upper':str(C3),'C1_upper':str(C1),'pi_lower':str(piL),'E076_norm':E,'source076_results_sha256':sha(P/'results.json'),'source076_manifest_sha256':sha(P/'manifest.json'),'strict_cap_verified':True}
 (H/'cap_scalar.json').write_text(json.dumps(scalar,indent=2)+'\n')
 state=json.loads(SRC.read_text());assert state['epsilon']=='1/64';assert len(state['roots'])==52360;roots=state['roots'];assert roots==[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
 rootanchor={}
 for e in state['passed']:
  rid,path=e['node']
  if path=='':assert rid not in rootanchor;rootanchor[rid]=(e['certificate']['anchor'],'passed.root')
 for e in state['splits']:
  rid,path=e['node']
  if path=='':assert rid not in rootanchor;rootanchor[rid]=(e['failed_certificate']['anchor'],'split.root')
 assert len(rootanchor)==52360
 records=[]
 for rid,indices in enumerate(roots):
  key,origin=rootanchor[rid];rec=state['cache'][key];assert key==','.join(rec['point']);point=list(map(F,rec['point']));assert point[0]==1 and all(F(i,32)<x<F(i+1,32) for i,x in zip(indices,point[1:]));hl,hu=map(F,rec['values'][0]);assert 0<hl<=hu
  weight=factorial(4)//prod(factorial(v) for v in Counter(indices).values());records.append({'rootid':rid,'indices':indices,'canonical_nonincreasing':indices,'point':rec['point'],'cached_H1over64':{'lower':str(hl),'upper':str(hu)},'weight':weight,'source_anchor_key':key,'origin':origin})
  if rid%4096==0 and time.monotonic()-began>540:raise TimeoutError('finite extraction soft deadline')
 records.sort(key=lambda r:(-(F(r['weight'])**4*F(r['cached_H1over64']['upper'])**5),r['rootid']))
 for rank,r in enumerate(records):r['rank']=rank
 assert sum(r['weight'] for r in records)==32**4
 payload=json.dumps(records,sort_keys=True,separators=(',',':')).encode();(H/'roots.json.gz').write_bytes(gzip.compress(payload,mtime=0));(S/'roots.json').write_bytes(payload)
 result={'status':'complete','count':len(records),'sum_weights':32**4,'weighted_volume':'1','ordering':'descending exact weight^4*cached_H1over64_upper^5, then ascending rootid','source066':str(SRC),'source066_sha256':LOCK,'root_input_sha256':sha(H/'roots.json.gz'),'uncompressed_root_input_sha256':sha(S/'roots.json'),'root_input_bytes':(H/'roots.json.gz').stat().st_size,'prepare_code_sha256':sha(__file__),'source076_results_sha256':sha(P/'results.json'),'seconds':time.monotonic()-began,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'pid':os.getpid(),'first32_rootids':[r['rootid'] for r in records[:32]],'source_anchor_origins':dict(Counter(r['origin'] for r in records))}
 (H/'prepare_receipt.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps(result,indent=2));print('Hcap',ep(cap,True),float(cap))
if __name__=='__main__':main()

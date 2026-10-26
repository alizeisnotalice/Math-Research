#!/usr/bin/env python3
"""Standalone user checker. NOT RUN during151 packaging.

Default prints help. --decode-only checks bytes. --verify-certificates checks
saved interval-based certificates (including required g105 mass cells).
Optional --numeric-count explicitly reproduces only the requested jet range
using the embedded original046 Arb kernel; needs python-flint and is a
production-kernel reproduction, not an independent PF audit.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
import argparse,base64,gzip,hashlib,json,sys
H=Path(__file__).resolve().parent
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dataset(name,m):
 for part in m['datasets'][name]:
  with gzip.open(H/part['path'],'rt',encoding='utf-8') as f:
   for line in f:yield json.loads(line)
def verify_bytes(m):
 for name,h in m['files'].items():assert sha(H/name)==h
 for parts in m['datasets'].values():
  for p in parts:
   assert sha(H/p['path'])==p['sha256'] and (H/p['path']).stat().st_size<100000000;h=hashlib.sha256();n=rows=0
   with gzip.open(H/p['path'],'rb') as f:
    for b in f:h.update(b);n+=len(b);rows+=1
   assert h.hexdigest()==p['uncompressed_sha256'] and n==p['uncompressed_bytes'] and rows==p['records']
def pointkey(p):return ','.join(map(str,map(F,p)))
ROOTS=[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
def box(n):
 lo=[F(x,32) for x in ROOTS[n[0]]];hi=[x+F(1,32) for x in lo]
 for j,b in enumerate(n[1]):
  k=j%4;mid=(lo[k]+hi[k])/2
  if b=='0':hi[k]=mid
  else:assert b=='1';lo[k]=mid
 return lo,hi
def ceildiv(a,b):return -((-a)//b)
LOG2={}
def basic_log_interval(a,b,bits):
 # 1<=a/b<=2. Positive atanh series, all rounding directed in integer units.
 assert b<=a<=2*b;S=1<<bits;n=a-b;d=a+b;n2=n*n;d2=d*d;pl=S*n//d;pu=ceildiv(S*n,d);lower=upper=0;N=bits//3+24
 for k in range(N):
  lower+=(2*pl)//(2*k+1);upper+=ceildiv(2*pu,2*k+1);pl=pl*n2//d2;pu=ceildiv(pu*n2,d2)
 upper+=ceildiv(2*pu*d2,(2*N+1)*(d2-n2));return lower,upper
def log_interval(x,bits):
 assert x>0;a,b=x.numerator,x.denominator;k=a.bit_length()-b.bit_length()
 if k>=0:b<<=k
 else:a<<=-k
 while a<b:a*=2;k-=1
 while a>=2*b:b*=2;k+=1
 lo,hi=basic_log_interval(a,b,bits)
 if bits not in LOG2:LOG2[bits]=basic_log_interval(2,1,bits)
 l2,u2=LOG2[bits]
 return (F(lo+k*(l2 if k>=0 else u2),1<<bits),F(hi+k*(u2 if k>=0 else l2),1<<bits))
def check_delta(delta,Hlower,target):
 for bits in (512,768,1024,2048,4096):
  lo,hi=log_interval(Hlower/target,bits)
  if delta<=lo:return
  if delta>hi:raise AssertionError('saved Delta is above certified log upper')
 raise AssertionError('log comparison unresolved at4096bits;not a numerical counterexample')
def margins(point,delta,gamma,q,lo,hi):
 rows=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)];rows.append(delta+sum(min(a*x,b*x) for x,(a,b) in zip(v,gamma))-q*sum(x*x for x in v))
 return min(rows)
def numeric_core(core,jet,target,delta_name,gradient_name):
 assert core['epsilon']==jet['epsilon'] and pointkey(core['point'])==pointkey(jet['point']);values=[tuple(map(F,p)) for p in jet['values']];hl,hu=values[0];assert 0<hl<=hu;gamma=[tuple(map(F,p)) for p in core[gradient_name]];assert len(gamma)==4;delta=F(core[delta_name]);check_delta(delta,hl,target)
 for (a,b),(gl,gu) in zip(values[1:],gamma):assert a<=b and gl<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=gu
 return list(map(F,core['point'])),delta,gamma,F(core['penalty_coefficient'])
def check_certificates(m,meta,begin,count):
 jets={}
 for r in dataset('jets',m):
  h=r.pop('minimal_projection_sha256');assert digest(r)==h and r['jet_id'] not in jets and r['model']==meta['model'] and r['epsilon'] in ('1/64','1/128') and len(r['values'])==5;jets[r['jet_id']]=r
 coefficients={r['canonical_id']:r for r in dataset('coefficients',m)};assert len(coefficients)==107530
 for oid,r in coefficients.items():assert digest(r['core'])==oid and r['raw_record_sha256']==jets[r['jet_id']]['raw_record_sha256']
 leaves=list(dataset('leaves',m));assert len(leaves)==111201;chosen=leaves[begin:begin+count] if count else leaves[begin:];needed={r['canonical_id'] for r in chosen};provenance={r['source_alias_index']:r for r in dataset('provenance',m)};bycore=defaultdict(list)
 for row in provenance.values():bycore[row['canonical_id']].append(row)
 supports={tuple(r['node']):r for r in dataset('supports',m)};support_done=set();core_done={}
 for oid in needed:
  obj=coefficients[oid];c=obj['core'];jet=jets[obj['jet_id']];args=numeric_core(c,jet,F(11,100000),'Delta110_lower','loggradient_directed');eps=F(c['epsilon']);q=args[3];assert q>0
  if c['curvature_rule_class']=='075_full_frame':assert q==1/(16*eps*eps)
  else:
   assert c['curvature_rule_class']=='129_anchor_fixed_Jensen' and eps==F(1,128) and q==1024-768*F(21,200000)/F(jet['values'][0][1]);aliases=[r for r in bycore[oid] if r['curvature_proof'] is not None];assert aliases;node=tuple(aliases[0]['curvature_proof']['node']);support=supports[node]
   if node not in support_done:
    cert=support['certificate'];sj=jets[support['jet_id']];sc={'epsilon':'1/64','point':sj['point'],'Delta_lower':cert['Delta_lower'],'loggradient_directed':cert['loggradient_directed'],'penalty_coefficient':cert['penalty_coefficient']};a=numeric_core(sc,sj,F(21,200000),'Delta_lower','loggradient_directed');assert a[3]==256;lo,hi=box(support['node']);assert margins(*a,lo,hi)>0;support_done.add(node)
   lo,hi=box(node);assert all(a<=x<=b for a,x,b in zip(lo,args[0][1:],hi))
  core_done[oid]=args
 for leaf in chosen:
  oid=leaf['canonical_id'];assert leaf['epsilon']==coefficients[oid]['core']['epsilon'] and leaf['closed_box'];lo,hi=box(leaf['node']);assert list(map(str,lo))==leaf['box_lower'] and list(map(str,hi))==leaf['box_upper'];margin=margins(*core_done[oid],lo,hi);assert margin>0 and margin==F(leaf['claimed_minimum_log_margin']);p=provenance[leaf['source_alias_index']];assert p['canonical_id']==oid and p['node']==leaf['node']
 if begin==0 and len(chosen)==111201:
  by=defaultdict(list);seen=set()
  for r in leaves:n=tuple(r['node']);assert n not in seen;seen.add(n);by[n[0]].append(n[1])
  assert len(by)==52360
  for root in range(52360):
   paths=sorted(by[root]);assert all(not b.startswith(a) for a,b in zip(paths,paths[1:]));depth=max(map(len,paths));assert sum(1<<(depth-len(p)) for p in paths)==1<<depth
 return {'certificate_leaves_checked':len(chosen),'unique_coefficients_checked':len(needed),'g105_support_cells_checked':len(support_done),'full_root_cover_checked':begin==0 and len(chosen)==111201,'scope':'saved qualifiedinterval certificate reproduction;formal global/C5 acceptance remains external'}
def recheck_numeric(m,meta,begin,count):
 docs={r['document_id']:r for r in dataset('documents',m)}
 for r in docs.values():
  b=r['text'].encode('utf-8') if r['encoding']=='utf-8' else base64.b64decode(r['data']);assert hashlib.sha256(b).hexdigest()==r['original_sha256']
 namespace={'__name__':'portable046_user_recheck','__file__':'<embedded046>'};exec(compile(docs[meta['kernel_document_id']]['text'],'<embedded046>','exec'),namespace);table=json.loads(gzip.decompress(base64.b64decode(docs[meta['bodytable_document_id']]['data'])));kernels={};checked=0
 for i,r in enumerate(dataset('jets',m)):
  if not begin<=i<begin+count:continue
  e=F(r['epsilon'])
  if e not in kernels:kernels[e]=namespace['FirstJetKernel'](table,e)
  got=kernels[e].evaluate(r['point'],512,None,None,128);assert got['complete'] and got['all_cones']==3504
  for v,pair in zip(got['values'],r['values']):assert F(pair[0])<=namespace['endpoint'](v)<=namespace['endpoint'](v,True)<=F(pair[1])
  checked+=1
 return {'numeric_jets_reproduced':checked,'scope':'explicit selectedrange,originalArbkernel512bits;not independentPF audit'}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--decode-only',action='store_true');ap.add_argument('--verify-certificates',action='store_true');ap.add_argument('--begin',type=int,default=0);ap.add_argument('--count',type=int,default=0);ap.add_argument('--numeric-begin',type=int,default=0);ap.add_argument('--numeric-count',type=int,default=0);a=ap.parse_args()
 if not(a.decode_only or a.verify_certificates or a.numeric_count):ap.print_help();return
 m=json.load(open(H/'manifest.json'));meta=json.load(open(H/'metadata.json'));verify_bytes(m);result={'compressed_and_decoded_SHA_checks':'PASS'}
 if a.verify_certificates:result.update(check_certificates(m,meta,a.begin,a.count))
 if a.numeric_count:result.update(recheck_numeric(m,meta,a.numeric_begin,a.numeric_count))
 print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':main()

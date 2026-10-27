#!/usr/bin/env python3
"""Frozen C certificate audit: integer polynomial endpoint arithmetic and
incremental disjoint orthotope subtraction. No producer/audit math imports,
root isolation, parameter search, optimization or new lower bound.
"""
from pathlib import Path
from functools import lru_cache,cmp_to_key
from itertools import product
from math import prod,gcd
import argparse,hashlib,json,os,time
H=Path(__file__).resolve().parent
ZERO=(0,)*25

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def save(p,x):
 p=Path(p);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n');os.replace(q,p)
def ratio(s):
 p=s.split('/');assert len(p) in (1,2);n=int(p[0]);d=int(p[1]) if len(p)==2 else 1;assert d>0;return n,d
def text_ratio(n,d):
 assert d>0;g=gcd(n,d);n//=g;d//=g;return str(n) if d==1 else str(n)+'/'+str(d)
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def minus(a,b):return tuple(x-y for x,y in zip(a,b))
def scale(a,k):return tuple(x*k for x in a)
@lru_cache(None)
def mul(a,b):
 out=[0]*25
 for k,x in enumerate(a):
  if not x:continue
  i,j=divmod(k,5)
  for ell,y in enumerate(b):
   if not y:continue
   p,q=divmod(ell,5);s,t=i+p,j+q;out[(s%5)*5+t%5]+=x*y*(2**(s//5))*(3**(t//5))
 return tuple(out)
def mon(i,j=0,n=1):
 a=[0]*25;a[(i%5)*5+j%5]=n*2**(i//5)*3**(j//5);return tuple(a)
def constant(n):return mon(0,0,n)

class Order:
 def __init__(self,cert):
  self.D=int(cert['root_isolation_denominator']);p,q=cert['root_lower_numerators'];self.p=p;self.q=q;D=self.D
  assert 0<p and 0<q and p**5<=2*D**5<(p+1)**5 and q**5<=3*D**5<(q+1)**5
  assert 5*(p+1)<6*D
  self.lowmon=tuple(p**i*q**j*D**(8-i-j) for i in range(5) for j in range(5));self.highmon=tuple((p+1)**i*(q+1)**j*D**(8-i-j) for i in range(5) for j in range(5));self.comparisons=0
 @lru_cache(None)
 def bounds(self,a):
  lo=hi=0
  for c,l,u in zip(a,self.lowmon,self.highmon):
   lo+=c*(l if c>=0 else u);hi+=c*(u if c>=0 else l)
  return lo,hi
 @lru_cache(None)
 def compare(self,a,b):
  self.comparisons+=1;d=minus(a,b)
  if d==ZERO:return 0
  lo,hi=self.bounds(d)
  if lo>0:return 1
  if hi<0:return -1
  raise ArithmeticError('saved root intervals do not decide endpoint order')
 def low(self,*xs):return min(xs,key=cmp_to_key(self.compare))
 def high(self,*xs):return max(xs,key=cmp_to_key(self.compare))

def radius_scaled4(lengths):
 n=prod(lengths);m=n;e=f=0
 while m%2==0:m//=2;e+=1
 while m%3==0:m//=3;f+=1
 assert m==1
 r=mon(e,f,2);fifth=constant(1)
 for j in range(5):fifth=mul(fifth,r)
 assert fifth==constant(32*n) # (2r_physical)^5=n, r_physical=r/4
 return r

def atom_box(lengths,spacings,order):
 """Intersect the half cell with all literal atom-distance inequalities."""
 radius=radius_scaled4(lengths);out=[];atoms=[]
 for a,l in zip(spacings,lengths):
  m=l//2;indices=tuple(range(-m,m+1)) if l%2 else tuple(range(1-m,m+1));assert len(indices)==l and len(set(indices))==l and all(-2<=i<=2 for i in indices)
  assert all(c%2==0 for c in a);half=tuple(c//2 for c in a)
  lo=order.high(ZERO,*[minus(scale(a,i),radius) for i in indices]);hi=order.low(half,*[plus(scale(a,i),radius) for i in indices])
  assert order.compare(lo,hi)<0
  assert order.compare(ZERO,lo)<=0 and order.compare(hi,half)<=0
  for i in indices:
   z=scale(a,i)
   for x in (lo,hi):assert order.compare(minus(x,z),radius)<=0 and order.compare(minus(z,x),radius)<=0
  out.extend((lo,hi));atoms.append(indices)
 all_atoms=set(product(*atoms));assert len(all_atoms)==prod(lengths)
 return tuple(out),radius,len(all_atoms)

def intersect(a,b):
 out=[]
 for j in range(5):
  lo=max(a[2*j],b[2*j]);hi=min(a[2*j+1],b[2*j+1])
  if lo>=hi:return None
  out.extend((lo,hi))
 return tuple(out)
def subtract(a,b):
 """Return disjoint part of A outside B; at most ten strips, using only endpoint ranks."""
 cut=intersect(a,b)
 if cut is None:return [a]
 work=list(a);out=[]
 for j in range(5):
  if work[2*j]<cut[2*j]:
   c=work.copy();c[2*j+1]=cut[2*j];out.append(tuple(c));work[2*j]=cut[2*j]
  if cut[2*j+1]<work[2*j+1]:
   c=work.copy();c[2*j]=cut[2*j+1];out.append(tuple(c));work[2*j+1]=cut[2*j+1]
 assert tuple(work)==cut and len(out)<=10
 assert all(all(c[2*j]<c[2*j+1] for j in range(5)) and intersect(c,b) is None for c in out)
 for i,c in enumerate(out):assert all(intersect(c,d) is None for d in out[:i])
 return out

def main(args):
 plan=read(H/'plan.json');assert sha(H/'plan.json')==args.plan_sha and not (H/'result.json').exists();assert sha(__file__)==plan['code_sha256']
 for p,r in plan['sources'].items():assert Path(p).stat().st_size==r['bytes'] and sha(p)==r['sha256']
 start=time.monotonic();save(H/'execution.json',{'status':'running','pid':os.getpid(),'plan_sha256':args.plan_sha,'before_any_new_math':True})
 cert=read(plan['certificate']);audit=read(plan['old_audit_result']);lock=read(plan['old_lock']);order=Order(cert);D=order.D
 assert cert['spacing']==lock['spacing']==['1+u^3/2']+['(1+u)/2']*4
 spacings=[plus(constant(4),mon(3,0,2))]+[plus(constant(2),mon(1,0,2))]*4 # physical spacing is vector/4
 assert all(order.compare(a,ZERO)>0 and order.compare(a,constant(8))<0 for a in spacings)
 labels=[tuple(ls) for ls in audit['retained_witness_counts']];assert len(labels)==len(set(labels))==42 and all(all(1<=l<=4 for l in ls) for ls in labels)
 boxes=[];witness=[]
 for ls in labels:
  b,r,n=atom_box(ls,spacings,order);boxes.append(b);witness.append({'lengths':ls,'mass':n,'radius_scaled4':r,'box_scaled4':b})
 assert len(set(boxes))==42
 endpoints=[sorted({b[2*j+s] for b in boxes for s in (0,1)},key=cmp_to_key(order.compare)) for j in range(5)];rank=[{c:i for i,c in enumerate(xs)} for xs in endpoints]
 ranks=[tuple(rank[j][b[2*j+s]] for j in range(5) for s in (0,1)) for b in boxes]
 segment_counts=[sum(any(b[2*j]<=a and a+1<=b[2*j+1] for b in ranks) for a in range(len(endpoints[j])-1)) for j in range(5)]
 assert segment_counts==cert['coordinate_segment_counts']==[12,14,14,14,14]
 save(H/'witnesses.json',{'endpoint_scale':4,'spacings_scaled4':spacings,'witnesses':witness,'rank_endpoints':endpoints,'rank_boxes':ranks})
 pieces=[];subtractions=0
 for index,b in enumerate(ranks):
  extra=[b]
  for old in pieces:
   next_extra=[]
   for c in extra:next_extra.extend(subtract(c,old));subtractions+=1
   extra=next_extra
   if not extra:break
  pieces.extend(extra);assert len(pieces)<=prod(segment_counts)
  save(H/'checkpoint.json',{'verified_witnesses':index+1,'disjoint_pieces':pieces,'subtractions':subtractions,'code_sha256':sha(__file__),'plan_sha256':args.plan_sha})
 # Independent geometric consistency checks on the disjoint partition.
 for i,b in enumerate(pieces):
  assert all(intersect(b,c) is None for c in pieces[:i])
  assert any(all(a[2*j]<=b[2*j] and b[2*j+1]<=a[2*j+1] for j in range(5)) for a in ranks)
 @lru_cache(None)
 def width(j,a,b):return minus(endpoints[j][b],endpoints[j][a])
 union=ZERO
 for b in pieces:
  volume=constant(1)
  for j in range(5):volume=mul(volume,width(j,b[2*j],b[2*j+1]))
  union=plus(union,volume)
 # Coordinate widths have denominator4; the ratio32 cancels4^5 to denominator32.
 expected=[ratio(cert['coefficients'].get(str(i//5)+','+str(i%5),'0')) for i in range(25)]
 for c,(n,d) in zip(union,expected):assert c*d==32*n
 for c,i in zip(union,range(25)):
  n,d=ratio(audit['ratio_coefficients'].get(str(i//5)+','+str(i%5),'0'));assert c*d==32*n
 lo,hi=order.bounds(union);den=32*D**8;assert lo<=hi
 savedlo,savedhi=map(ratio,cert['ratio_interval']);assert lo*savedlo[1]==savedlo[0]*den and hi*savedhi[1]==savedhi[0]*den
 qn,qd=19943,10000;assert lo*qd>qn*den
 b_up=ratio(read(plan['previous_B_lock'])['exact_interval'][1]);assert lo*b_up[1]>b_up[0]*den
 # Current lower ledger is a frozen semantic projection, not a future upper ledger.
 ledger=read(plan['lower_ledger_projection']);assert ledger['certificate']=='lower_round003/a2_equals_a3_results.json' and ledger['simple_rational_strict_lower']=='19943/10000'
 assert ledger['exact_interval']==cert['ratio_interval']
 for p,r in plan['sources'].items():assert sha(p)==r['sha256']
 save(H/'disjoint_certificate.json',{'endpoint_scale':4,'rank_endpoints':endpoints,'rank_boxes':ranks,'disjoint_pieces':pieces,'volume_numerator_coefficients_over1024':union,'ratio_denominator':32})
 result={'status':'frozen_existing_C_cert_independent_integer_disjoint_union_pass','pid':os.getpid(),'seconds':time.monotonic()-start,'plan_sha256':args.plan_sha,'retained_witnesses':42,'literal_atom_occurrences':sum(w['mass'] for w in witness),'disjoint_piece_count':len(pieces),'subtraction_calls':subtractions,'endpoint_comparisons':order.comparisons,'coordinate_segment_counts':segment_counts,'ratio_coefficients':{str(i//5)+','+str(i%5):text_ratio(c,32) for i,c in enumerate(union)},'ratio_interval':[text_ratio(lo,den),text_ratio(hi,den)],'strict_margin_above19943_over10000':text_ratio(lo*qd-qn*den,den*qd),'same_all25_coefficients':True,'same_original_ratio_interval':True,'all_sources_before_after_unchanged':True,'producer_or_old_audit_math_imports':False,'new_root_isolations':0,'new_parameter_search_or_lower_advancement':False,'erosion_scalar_replayed':False,'shared_arithmetic_limit':'Python arbitrary precision int plus gcd; quotient relations u^5=2,v^5=3 and original root brackets checked by integer fifth powers','disjoint_certificate_sha256':sha(H/'disjoint_certificate.json'),'witnesses_sha256':sha(H/'witnesses.json')}
 save(H/'result.json',result);print(json.dumps({'status':result['status'],'pid':result['pid'],'seconds':result['seconds'],'pieces':len(pieces),'result_sha256':sha(H/'result.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--plan-sha',required=True);main(p.parse_args())

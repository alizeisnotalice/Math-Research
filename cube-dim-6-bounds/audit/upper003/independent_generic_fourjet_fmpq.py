#!/usr/bin/env python3
"""Exact FLINT-fmpq four-jet reconstruction of one fixed round003 point."""
from __future__ import annotations
import hashlib, itertools, json, math, platform, time
from pathlib import Path
from flint import arb, ctx, fmpq
import flint

HERE=Path(__file__).resolve().parent
SRC=Path(__file__).resolve().parents[2] / 'upper/round003/direct_kernel'
FAN,POINTS,CPP,AUTHOR=(SRC/'fan_integer_input.txt',SRC/'face_points_generic.txt',SRC/'direct_kernel.cpp',SRC/'direct_kernel_generic_pilot_results.json')
Z=(fmpq(0),)*4
def hs(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def va(a,b): return tuple(x+y for x,y in zip(a,b,strict=True))
def vs(a,b): return tuple(x-y for x,y in zip(a,b,strict=True))
def vm(a,c): return tuple(x*c for x in a)
def cones():
 t=iter(FAN.read_text().split()); n=int(next(t));out=[]
 for _ in range(n):
  det=int(next(t));rs=[];prod=1
  for _ in range(6):
   q=int(next(t));prod*=q;rs.append((q,tuple(int(next(t)) for _ in range(6))))
  out.append((rs,fmpq(det,prod)))
 assert n==21 and not list(t);return out
def firstpoint():
 a=list(map(int,POINTS.read_text().split()));assert a[1]==32 and len(a)==194;return a[0],tuple(a[2:8])
def dd4(v):
 v=sorted(v);t=[[Z for _ in range(6)]for _ in range(6)]
 for i,x in enumerate(v):
  if not x:
   raise AssertionError('zero node')
  t[i][0]=(fmpq(0),fmpq(0),fmpq(0),fmpq(1 if x>0 else -1))
 for d in range(1,6):
  for i in range(d,6):
   if v[i]==v[i-d]:raise AssertionError('repeated node')
   t[i][d]=vm(vs(t[i][d-1],t[i-1][d-1]),fmpq(1)/(v[i]-v[i-d]))
 return t[5][5]
def main():
 begun=time.perf_counter();cs=cones();den,x=firstpoint();T=Z;terms=zeros=repeats=0;maxn=maxd=maxcross=maxproduct=maxl1=0;mingap=None
 for perm in itertools.permutations(range(6)):
  y=tuple(x[i] for i in perm)
  for signs in itertools.product((-1,1),repeat=6):
   tx=tuple(s*z for s,z in zip(signs,y,strict=True))
   for rs,alpha in cs:
    raw=[];v=[]
    for q,r in rs:
     pieces=[a*b for a,b in zip(tx,r,strict=True)];n=sum(pieces);d=den*q;raw.append((n,d));v.append(fmpq(n,d));maxn=max(maxn,abs(n));maxd=max(maxd,d);maxproduct=max(maxproduct,*(abs(a) for a in pieces));maxl1=max(maxl1,sum(abs(a) for a in pieces))
    v.sort();terms+=1
    z=any(not b for b in v);rep=any(v[i]==v[i-1] for i in range(1,6));zeros+=z;repeats+=sum(v[i]==v[i-1] for i in range(1,6))
    if z or rep:raise AssertionError('fixed high-denominator point unexpectedly degenerate')
    for a,b in itertools.product(raw,repeat=2):maxcross=max(maxcross,abs(a[0])*b[1])
    for i in range(1,6):
     g=v[i]-v[i-1];mingap=g if mingap is None or g<mingap else mingap
    T=va(T,vm(dd4(v),alpha))
 assert terms==967680 and T[:3]==Z[:3]
 q=T[3]/128;ctx.prec=512;gb=(arb(int(q.numer()))/int(q.denom()))/arb.pi()**5
 au=json.loads(AUTHOR.read_text())['points'][0]
 out={'status':'PASS_exact_fmpq_generic_fourjet','scope':'one fixed fully generic point only; exact point diagnostic, no global/cost/C6 claim','inputs':{'fan_integer_input_sha256':hs(FAN),'points_sha256':hs(POINTS),'cpp_sha256':hs(CPP),'author_results_sha256':hs(AUTHOR)},'point_numerators':list(x),'denominator':den,'terms':terms,'zero_node_terms':zeros,'repeated_node_pairs':repeats,'fourjet':{'epsilon^-5':str(T[0]),'epsilon^-3':str(T[1]),'epsilon^-1':str(T[2]),'sign':str(T[3])},'G_exact':'('+str(q)+')/pi^5','G_arb_512':gb.str(80,more=True),'integer_bounds':{'max_abs_raw_n':maxn,'max_raw_d':maxd,'max_abs_multiplication_summand':maxproduct,'max_l1_dot_bound_before_signed_accumulation':maxl1,'max_raw_cross':maxcross,'int64_limit':2**63-1,'int128_limit':2**127-1,'min_exact_adjacent_gap':str(mingap)},'producer_long_double_G':au['G_direct_long_double'],'producer_cancellation_ratio':au['cancellation_ratio'],'longdouble_scope':'diagnostic only: this host reports LDBL_MANT_DIG=53 and sizeof(long double)=8; Kahan does not bound DD or total roundoff.','elapsed_seconds':time.perf_counter()-begun,'software':{'python':platform.python_version(),'python_flint':flint.__version__},'source_sha256':hs(Path(__file__))}
 (HERE/'independent_generic_fourjet_fmpq_results.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()

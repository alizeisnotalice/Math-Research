#!/usr/bin/env python3
"""G30A: bounded interval sign/monotonicity attempt, never infer from samples."""
import os,sys,json,time,hashlib,resource,concurrent.futures
from fractions import Fraction as Q
from collections import deque
from pathlib import Path
from math import comb
from flint import arb, arb_series, fmpq, ctx

def ab(q): return arb(fmpq(q.numerator,q.denominator))
def ball(a,b): return arb(ab((a+b)/2),ab((b-a)/2))
def calc(n,a,b):
 t=arb_series([ball(a,b),arb(1)],prec=3); q=1+1/(2*t); h=(1+1/t)/2
 S=arb_series([0],prec=3)
 for j in range(1,n): S+=Q(comb(n-1,j),j).numerator * (q**(n-1-j)) * ((-1)**j) * (1-t**(-j)) / (Q(comb(n-1,j),j).denominator * 2**j)
 B=1-h**n+n*t.log()-n*S/2-n*q**(n-1)*t.log()/2
 A=(1+1/t)**n
 tt=arb_series([t[0],t[1]],prec=2); AA=arb_series([A[0],A[1]],prec=2); BB=arb_series([B[0],B[1]],prec=2); bp=arb_series([B[1],2*B[2]],prec=2)
 H=(tt*(tt+1)*bp/n)**2-AA*BB
 return H[0],H[1]
def sign(x):
 if x>0:return 1
 if x<0:return -1
 return 0
def run(n):
 ctx.prec=512; start=time.monotonic(); cpu=time.process_time(); d=Q(1,4*n*2**n); points=[1+d]; x=d
 while x<1:
  x=min(2*x,Q(1));points.append(1+x)
 while points[-1]<4*n: points.append(min(2*points[-1],Q(4*n)))
 queue=deque(zip(points,points[1:])); count=0; signed=0; monotone=0; roots=[]; reasons={}; failures=[]
 while queue and count<4096 and time.monotonic()-start<45:
  a,b=queue.popleft();count+=1
  try:
   H,Hp=calc(n,a,b); s=sign(H)
   if s: signed+=1;continue
   sp=sign(Hp)
   if sp:
    sa=sign(calc(n,a,a)[0]); sb=sign(calc(n,b,b)[0])
    if sa and sb:
     if sa!=sb: roots.append([str(a),str(b),sa,sb,sp])
     else: monotone+=1
     continue
   mid=(a+b)/2
   queue.append((a,mid));queue.append((mid,b))
  except Exception as exc: failures.append([str(a),str(b),repr(exc)])
 stop='complete' if not queue else ('box limit' if count>=4096 else '45-second per-dimension wall limit')
 return {'n':n,'precision_bits':ctx.prec,'pid':os.getpid(),'boxes_evaluated':count,'strict_sign_boxes':signed,'same_sign_monotone_boxes':monotone,'isolated_roots':roots,'unresolved_boxes':[[str(a),str(b)] for a,b in queue]+failures,'stop':stop,'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'[1+1/(4 n 2^n),4n]; tail and endpoint excluded subject to separate proofs'}
if __name__=='__main__':
 with concurrent.futures.ProcessPoolExecutor(max_workers=3) as ex: rows=list(ex.map(run,(4,7,64)))
 out={'batch':'G30A','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limits':'4096 boxes and 45 seconds per dimension; 3 workers; 512 bits','note':'Results are a bounded attempt; isolated boxes require strict interval derivative and opposite endpoint signs. Not a global conclusion absent independent proof for outside scope and zero unresolved boxes.','results':rows}
 Path(__file__).with_name('g30a_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([{k:v for k,v in r.items() if k not in ('unresolved_boxes',)} | {'unresolved_count':len(r['unresolved_boxes'])} for r in rows],indent=2))

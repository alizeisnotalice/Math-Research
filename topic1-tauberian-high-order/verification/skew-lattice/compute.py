#!/usr/bin/env python3
"""Exact rational interval certificate for the stated R10 formula, not its geometry."""
from fractions import Fraction as F
from math import comb
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse,json,os,time
STEPS=200

def bound(x,mode):
 with localcontext() as c:
  c.prec=35;c.rounding=mode
  return str(Decimal(x.numerator)/Decimal(x.denominator))

def cases():
 first=[(n,k,d) for n in (2,3,4) for k in (F(n),F(2*n),F(4*n),F(16*n,3),F(7*n)) for d in (8,16,32,64,128,256,512,1024)]
 more=[(n,F(16*n,3),d) for n in range(5,17) for d in (16,32,64,128,256,512,1024)]
 assert len(first+more)==204 and len(set(first+more))==204
 return first+more

def evaluate(spec,root=None):
 n,k,den=spec;u=F(1,den);v=2*u;a=v**n;delta=u/(1+2*u);r=1/(1+2*u);e=k*u**4;m=(2*delta)**n;det=1-e**2;K=4*u**2*r**(2-n);width=delta-e
 poly=lambda h:K*sum((F(comb(n,j))*r**(n-j)*h**j for j in range(1,n+1)),F(0))
 guards={'e_lt_delta_over_4':e<delta/4,'r0_gt_2delta':r>2*delta,'positive_separation':1-2*e-delta>0,'positive_determinant_gap':det>m,'outer_gap_positive':width>0,'double_budget_at_outer':width>0 and poly(width)>2*e**2}
 if not all(guards.values()):return {'n':n,'kappa':str(k),'u':str(u),'status':'guard_failed','guards':guards},None
 if root is None:
  lo=F(0);hi=width
  for _ in range(STEPS):
   mid=(lo+hi)/2
   if poly(mid)<e**2:lo=mid
   else:hi=mid
  q=lo*2**STEPS/width;assert q.denominator==1;idx=q.numerator
 else:
  assert root['scale']=='delta-e' and root['denominator_power2']==STEPS
  idx=int(root['lower_numerator_hex'],16);assert int(root['upper_numerator_hex'],16)==idx+1 and 0<=idx<2**STEPS
  lo=width*idx/2**STEPS;hi=width*(idx+1)/2**STEPS
 def integral(x):return sum((F(comb(n,j))*r**(n-j)*x**(j+1)/F(j+1) for j in range(1,n+1)),F(0))
 Al=2*(e*lo-K/e*integral(hi));Au=2*(e*hi-K/e*integral(lo));f=(2*delta)**(n-2);bl=2*Al*f;bu=2*Au*f
 if n>=3:
  bl+=4*(n-2)*lo*(e+lo)**2*(2*delta)**(n-3);bu+=4*(n-2)*hi*(e+hi)**2*(2*delta)**(n-3)
 common=(det-m)*(1-m)*a*v**8;gl=(m*e**2-bu*(1-m))/common;gu=(m*e**2-bl*(1-m))/common
 guards.update(root_enclosed=poly(lo)<e**2<=poly(hi),root_inside=hi<width,A_lower_positive=Al>0,A_interval_order=Al<=Au,b_interval_order=bl<=bu,G_interval_order=gl<=gu)
 sign='positive' if gl>0 else ('negative' if gu<0 else 'unresolved')
 result={'n':n,'kappa':str(k),'u':str(u),'status':'certified' if all(guards.values()) else 'arithmetic_guard_failed','guards':guards,'root':{'scale':'delta-e','denominator_power2':STEPS,'lower_numerator_hex':hex(idx),'upper_numerator_hex':hex(idx+1)},'G':[bound(gl,ROUND_FLOOR),bound(gu,ROUND_CEILING)],'G_sign':sign,'predicted_limit':str((k*k-k**3/F(8*n))/256)}
 return result,(gl,gu)

def worker(spec):
 r,_=evaluate(spec);return os.getpid(),r

def verify(doc):
 expected=cases();records=doc['cases'];assert len(records)==len(expected)
 for record,spec in zip(records,expected):
  assert (record['n'],F(record['kappa']),F(record['u']))==(spec[0],spec[1],F(1,spec[2]))
  assert record.get('root') is not None or record['status']=='guard_failed'
  fresh,exact=evaluate(spec,record.get('root'));assert fresh['status']==record['status'] and fresh['guards']==record['guards']
  if exact is not None:
   assert fresh['root']==record['root'] and fresh['G_sign']==record['G_sign']
   assert F(Decimal(record['G'][0]))<=exact[0]<=exact[1]<=F(Decimal(record['G'][1]))
   assert fresh['predicted_limit']==record['predicted_limit']
 signs={s:sum(r.get('G_sign')==s for r in records) for s in ('positive','negative','unresolved')}
 assert signs==doc['summary']['signs'] and sum(r['status']=='guard_failed' for r in records)==doc['summary']['guard_failed']
 assert sum(r['status']=='arithmetic_guard_failed' for r in records)==doc['summary']['arithmetic_guard_failed']
 return signs

def main():
 p=argparse.ArgumentParser();sub=p.add_subparsers(dest='mode',required=True);g=sub.add_parser('generate');g.add_argument('--output-dir',type=Path,required=True);v=sub.add_parser('verify');v.add_argument('--input',type=Path,required=True);args=p.parse_args()
 if args.mode=='verify':
  start=time.perf_counter();d=json.loads(args.input.read_text());signs=verify(d);print(json.dumps({'verify':'PASS','cases':len(d['cases']),'signs':signs,'seconds':time.perf_counter()-start}));return
 start=time.perf_counter();specs=cases()
 with ProcessPoolExecutor(max_workers=15) as pool:results=list(pool.map(worker,specs))
 pids=sorted({pid for pid,_ in results});rs=[r for _,r in results];signs={s:sum(r.get('G_sign')==s for r in rs) for s in ('positive','negative','unresolved')}
 data={'summary':{'cases':len(rs),'requested_workers':15,'actual_worker_pids':pids,'actual_worker_count':len(pids),'elapsed_seconds':time.perf_counter()-start,'signs':signs,'guard_failed':sum(r['status']=='guard_failed' for r in rs),'arithmetic_guard_failed':sum(r['status']=='arithmetic_guard_failed' for r in rs),'bisection_steps':STEPS},'cases':rs}
 args.output_dir.mkdir(parents=True,exist_ok=True);target=args.output_dir/'certificates.json';target.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'output':str(target),'summary':data['summary']}))
if __name__=='__main__':main()

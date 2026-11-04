#!/usr/bin/env python3
from fractions import Fraction as F
from math import comb
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse,json,os,time
STEPS=400

def bound(x,mode):
 with localcontext() as c:
  c.prec=35;c.rounding=mode
  return str(Decimal(x.numerator)/Decimal(x.denominator))

def cases():return [(n,F(16*n,3),2**j) for n in range(2,17) for j in (8,10,12,14,16,18,20)]

def evaluate(spec,roots=None):
 n,k,den=spec;u=F(1,den);v=2*u;a=v**n;delta=u/(1+2*u);r=1/(1+2*u);e=k*u**4;m=(2*delta)**n;det=(1-e)**(n-1)*(1+(n-1)*e);K=4*u**2*r**(2-n);width=delta-e;M=3**n
 poly=lambda h:K*sum((F(comb(n,j))*r**(n-j)*h**j for j in range(1,n+1)),F(0))
 guards={'u_le_quarter':u<=F(1,4),'strain':(n-1)*k*u**3<F(1,3),'budget_guard':M*k**2*u**5<F(4*n,9),'positive_width':width>0,'e_positive_and_lt_delta_over_2':0<e<delta/2,'strain_e_lt_1_over_12':(n-1)*e<F(1,12),'r0_at_least_2_over_3':r>=F(2,3),'r0_gt_2delta':r>2*delta,'det_positive_and_gt_m':det>m>0,'det_minus_m_gt_29_over_36':det-m>F(29,36),'outer_double_M_budget':width>0 and poly(width)>2*M*e**2}
 rec={'n':n,'kappa':str(k),'u':str(u),'M':M,'guards':guards}
 if not all(guards.values()):rec['status']='guard_failed';return rec,None
 def bracket(mult,name):
  if roots is None:
   lo=F(0);hi=width
   for _ in range(STEPS):
    mid=(lo+hi)/2
    if poly(mid)<mult*e**2:lo=mid
    else:hi=mid
   idx=lo*2**STEPS/width;assert idx.denominator==1;idx=idx.numerator
  else:
   data=roots[name];assert data['scale']=='delta-e' and data['denominator_power2']==STEPS
   idx=int(data['lower_numerator_hex'],16);assert int(data['upper_numerator_hex'],16)==idx+1 and 0<=idx<2**STEPS
   lo=width*idx/2**STEPS;hi=width*(idx+1)/2**STEPS
  return lo,hi,{'scale':'delta-e','denominator_power2':STEPS,'lower_numerator_hex':hex(idx),'upper_numerator_hex':hex(idx+1)}
 l,h,sroot=bracket(1,'hstar');lm,hm,mroot=bracket(M,'hM')
 integ=lambda x:sum((F(comb(n,j))*r**(n-j)*x**(j+1)/F(j+1) for j in range(1,n+1)),F(0))
 Al=2*(e*l-K/e*integ(h));Au=2*(e*h-K/e*integ(l));pairl=2*Al*(2*delta)**(n-2);pairu=2*Au*(2*delta)**(n-2)
 bl=comb(n,2)*pairl;bu=comb(n,2)*pairu
 if n>=3:
  bl+=8*comb(n,3)*(e+lm)**3*(2*delta)**(n-3);bu+=8*comb(n,3)*(e+hm)**3*(2*delta)**(n-3)
 common=(det-m)*(1-m)*a*v**8;gl=(m*(1-det)-bu*(1-m))/common;gu=(m*(1-det)-bl*(1-m))/common
 guards.update(hstar_enclosed=poly(l)<e**2<=poly(h),hM_enclosed=poly(lm)<M*e**2<=poly(hm),both_inside=h<width and hm<width,root_order=h<lm,A_lower_positive=Al>0,A_interval_order=Al<=Au,B_interval_order=bl<=bu,G_interval_order=gl<=gu)
 sign='positive' if gl>0 else ('negative' if gu<0 else 'unresolved')
 rec.update(status='certified' if all(guards.values()) else 'arithmetic_guard_failed',roots={'hstar':sroot,'hM':mroot},G=[bound(gl,ROUND_FLOOR),bound(gu,ROUND_CEILING)],G_sign=sign,predicted_limit=str(F(comb(n,2)*n*n,27)))
 return rec,(gl,gu)

def worker(spec):return os.getpid(),evaluate(spec)[0]
def verify(doc):
 expected=cases();records=doc['cases'];assert len(records)==len(expected)
 for record,spec in zip(records,expected):
  assert (record['n'],F(record['kappa']),F(record['u']))==(spec[0],spec[1],F(1,spec[2]));assert record.get('roots') is not None or record['status']=='guard_failed'
  fresh,exact=evaluate(spec,record.get('roots'));assert fresh['status']==record['status'] and fresh['guards']==record['guards']
  if exact is not None:
   assert fresh['roots']==record['roots'] and fresh['G_sign']==record['G_sign'] and fresh['predicted_limit']==record['predicted_limit']
   assert F(Decimal(record['G'][0]))<=exact[0]<=exact[1]<=F(Decimal(record['G'][1]))
 signs={s:sum(r.get('G_sign')==s for r in records) for s in ('positive','negative','unresolved')};assert signs==doc['summary']['signs'];assert sum(r['status']=='guard_failed' for r in records)==doc['summary']['guard_failed'];assert sum(r['status']=='arithmetic_guard_failed' for r in records)==doc['summary']['arithmetic_guard_failed'];return signs

def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='mode',required=True);g=sp.add_parser('generate');g.add_argument('--output-dir',type=Path,required=True);v=sp.add_parser('verify');v.add_argument('--input',type=Path,required=True);args=p.parse_args();start=time.perf_counter()
 if args.mode=='verify':
  doc=json.loads(args.input.read_text());signs=verify(doc);print(json.dumps({'verify':'PASS','cases':len(doc['cases']),'signs':signs,'seconds':time.perf_counter()-start}));return
 with ProcessPoolExecutor(max_workers=15) as pool:res=list(pool.map(worker,cases()))
 pids=sorted({pid for pid,_ in res});records=[r for _,r in res];signs={s:sum(r.get('G_sign')==s for r in records) for s in ('positive','negative','unresolved')}
 out={'summary':{'cases':len(records),'requested_workers':15,'actual_worker_pids':pids,'actual_worker_count':len(pids),'elapsed_seconds':time.perf_counter()-start,'signs':signs,'guard_failed':sum(r['status']=='guard_failed' for r in records),'arithmetic_guard_failed':sum(r['status']=='arithmetic_guard_failed' for r in records),'bisection_steps_per_root':STEPS},'cases':records};args.output_dir.mkdir(parents=True,exist_ok=True);target=args.output_dir/'certificates.json';target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'output':str(target),'summary':out['summary']}))
if __name__=='__main__':main()

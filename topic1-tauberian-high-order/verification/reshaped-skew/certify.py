from fractions import Fraction as F
from math import comb
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json,os,time
BASE=Path(__file__).resolve().parent;OUT=BASE/'fraction';STEPS=320

def bound(x,mode):
 with localcontext() as c:c.prec=30;c.rounding=mode;return str(Decimal(x.numerator)/Decimal(x.denominator))

def case(n,k,den,root=None,return_exact=False):
 u=F(1,den);lam=F(11,10)*k**3/(8*n*n);v=2*u;a=v**n;opt=u/(1+2*u);w=opt+lam*u**6;rho=1-2*w;e=k*u**3;det=1-e**2;outer=(2*w)**n;M=a*rho**n;D=outer-M;K=a/(2*w)**(n-2);width=w-e
 poly=lambda h:K*sum((F(comb(n,j))*rho**(n-j)*h**j for j in range(1,n+1)),F(0))
 g={'positive_width':width>0,'e_lt_w_over_4':0<e<w/4,'rho_gt_2w':rho>2*w,'separation_positive':1-2*e-w>0,'det_gt_outer':det>outer,'det_gt_M':det>M>0,'D_nonnegative':D>=0,'double_budget':width>0 and rho>0 and poly(width)>2*e**2}
 out={'n':n,'kappa':str(k),'lambda':str(lam),'u':str(u),'guards':g,'M':str(M),'bisection_steps':STEPS}
 if not all(g.values()):
  out.update(status='guard_failed',failed=[x for x in g if not g[x]]);return (out,None) if return_exact else out
 if root is None:
  lo=F(0);hi=width
  for _ in range(STEPS):
   mid=(lo+hi)/2
   if poly(mid)<e**2:lo=mid
   else:hi=mid
 else:
  assert root['scale']=='w-e' and root['pow2']==STEPS
  index=int(root['lower_hex'],16);assert int(root['upper_hex'],16)==index+1 and 0<=index<2**STEPS
  lo=width*index/2**STEPS;hi=width*(index+1)/2**STEPS
 integ=lambda x:sum((F(comb(n,j))*rho**(n-j)*x**(j+1)/F(j+1) for j in range(1,n+1)),F(0))
 Al=2*(e*lo-K/e*integ(hi));Au=2*(e*hi-K/e*integ(lo));bl=2*Al*(2*w)**(n-2);bu=2*Au*(2*w)**(n-2)
 if n>=3:
  bl+=4*(n-2)*lo*(e+lo)**2*(2*w)**(n-3);bu+=4*(n-2)*hi*(e+hi)**2*(2*w)**(n-3)
 low=D-bu;high=D-bl;lower='positive' if low>0 else ('negative' if high<0 else 'unresolved')
 m0=(2*opt)**n;G=(M-det*m0)/((det-M)*(1-m0)*a*v**6);idx=lo*2**STEPS/width;assert idx.denominator==1
 g.update(root_enclosed=poly(lo)<e**2<=poly(hi),root_inside=hi<width,A_lower_positive=Al>0,A_order=Al<=Au,B_order=bl<=bu)
 out.update(status='certified' if all(g.values()) else 'arithmetic_guard_failed',budget_eligible=(lower=='positive'),admissible_G_sign=(('positive' if G>0 else ('negative' if G<0 else 'zero')) if lower=='positive' else 'not_applicable'),root={'scale':'w-e','pow2':STEPS,'lower_hex':hex(idx.numerator),'upper_hex':hex(idx.numerator+1)},D_minus_B=[bound(low,ROUND_FLOOR),bound(high,ROUND_CEILING)],budget_sign=lower,G=[bound(G,ROUND_FLOOR),bound(G,ROUND_CEILING)],G_sign='positive' if G>0 else ('negative' if G<0 else 'zero'),predicted_limit=str((k*k-2*n*lam)/64))
 return (out,(low,high,G)) if return_exact else out

def worker(sp):
 n,k=sp;t=time.perf_counter();r=[case(n,k,2**j) for j in (3,4,5,6,7,8,9,10)]
 path=OUT/f'n{n}_k{str(k).replace("/","_")}.json';path.write_text(json.dumps({'pid':os.getpid(),'n':n,'kappa':str(k),'elapsed':time.perf_counter()-t,'cases':r},indent=2)+'\n')
 return {'n':n,'kappa':str(k),'pid':os.getpid(),'elapsed':time.perf_counter()-t,'guard_failed':sum(x['status']=='guard_failed' for x in r),'budget_positive':sum(x.get('budget_sign')=='positive' for x in r),'budget_negative':sum(x.get('budget_sign')=='negative' for x in r),'G_positive':sum(x.get('G_sign')=='positive' for x in r),'G_negative':sum(x.get('G_sign')=='negative' for x in r)}
if __name__=='__main__':
 OUT.mkdir(exist_ok=True);t=time.perf_counter();spec=[(n,k) for n in (2,3,4) for k in (F(n),F(2*n),F(8*n,3),F(3*n),F(7*n,2))];print(json.dumps({'master':os.getpid(),'workers_requested':15,'branches':len(spec),'lambda':'11/10*kappa^3/(8*n^2)','u_denominators':[2**j for j in range(3,11)]}),flush=True)
 with ProcessPoolExecutor(max_workers=15) as pool:results=list(pool.map(worker,spec))
 z={'master':os.getpid(),'elapsed':time.perf_counter()-t,'worker_pids':sorted(set(x['pid'] for x in results)),'branches':results};(BASE/'fraction_summary.json').write_text(json.dumps(z,indent=2)+'\n');rows=[r for p in sorted(OUT.glob('n*_k*.json')) for r in json.loads(p.read_text())['cases']];assert len(rows)==120;(BASE/'cases.json').write_text(json.dumps({'cases':rows},indent=2)+'\n');print(json.dumps(z),flush=True)

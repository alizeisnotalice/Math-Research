from fractions import Fraction as F
from math import comb
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path
import argparse,json,os,time,multiprocessing as mp
BASE=Path(__file__).resolve().parent;STEPS=400

def dec(x,rounding):
 with localcontext() as ctx:ctx.prec=35;ctx.rounding=rounding;return str(Decimal(x.numerator)/Decimal(x.denominator))
def outward(lo,hi):return [dec(lo,ROUND_FLOOR),dec(hi,ROUND_CEILING)]
def laurent_neighbor(n,w,e):
 dp={0:F(1)}
 for _ in range(n):
  nxt={}
  for j,val in dp.items():
   for step,mul in ((0,2*w),(-1,e),(1,e)):nxt[j+step]=nxt.get(j+step,F(0))+val*mul
  dp=nxt
 return dp[0]-(2*w)**n

def case(n,j=20,roots=None):
 u=F(1,2**j);k=F(8*n,3);pairs=comb(n,2);lam=F(11,10)*pairs*k**3/(8*n*n);v=2*u;a=v**n;d0=u/(1+2*u);w=d0+lam*u**6;rho=1-2*w;e=k*u**3;det=(1-e)**(n-1)*(1+(n-1)*e);M=a*rho**n;mout=(2*w)**n;m0=(2*d0)**n;K=a/(2*w)**(n-2);width=w-e;mult=3**n
 poly=lambda h:K*sum((F(comb(n,i))*rho**(n-i)*h**i for i in range(1,n+1)),F(0))
 guards={'u_le_quarter':u<=F(1,4),'w_lt_quarter':0<w<F(1,4),'strain':(n-1)*e<F(1,12),'e_lt_w_over_2':0<e<w/2,'rho_gt_2w':rho>2*w,'det_gt_outer':det>mout,'det_gt_M':det>M>0,'M_gt_half_outer':2*M>mout,'double_high_budget':width>0 and poly(width)>2*mult*e**2}
 row={'n':n,'u':str(u),'kappa':str(k),'lambda':str(lam),'bisection_steps_per_root':STEPS,'guards':guards}
 if not all(guards.values()):row.update(status='guard_failed',failed=[x for x in guards if not guards[x]]);return row,None
 def bracket(mul,name):
  if roots is None:
   lo=F(0);hi=width
   for _ in range(STEPS):
    mid=(lo+hi)/2
    if poly(mid)<mul*e**2:lo=mid
    else:hi=mid
   idx=lo*2**STEPS/width;assert idx.denominator==1;idx=idx.numerator
  else:
   t=roots[name];assert t['scale']=='w-e' and t['pow2']==STEPS
   idx=int(t['lower_hex'],16);assert int(t['upper_hex'],16)==idx+1 and 0<=idx<2**STEPS
   lo=width*idx/2**STEPS;hi=width*(idx+1)/2**STEPS
  return lo,hi,{'scale':'w-e','pow2':STEPS,'lower_hex':hex(idx),'upper_hex':hex(idx+1)}
 l,h,sroot=bracket(1,'hstar');lh,hh,hroot=bracket(mult,'hH')
 integ=lambda x:sum((F(comb(n,i))*rho**(n-i)*x**(i+1)/F(i+1) for i in range(1,n+1)),F(0))
 Al=2*(e*l-K/e*integ(h));Au=2*(e*h-K/e*integ(l));bl=pairs*2*Al*(2*w)**(n-2);bu=pairs*2*Au*(2*w)**(n-2)
 if n>=3:bl+=8*comb(n,3)*(e+lh)**3*(2*w)**(n-3);bu+=8*comb(n,3)*(e+hh)**3*(2*w)**(n-3)
 D=mout-M;dl=D-bu;dh=D-bl
 neighbor=sum((F(comb(n,q)*comb(q,q//2))*e**q*(2*w)**(n-q) for q in range(2,n+1,2)),F(0));dp=laurent_neighbor(n,w,e)
 radius=rho+w;defect_upper=(M+neighbor)/(2*radius)**n;strict_margin=u**n-defect_upper;N=2**(7*j);inverse_norm=(1+(2*n-3)*e)/((1-e)*(1+(n-1)*e));fd=F(2*N-5,2*N)**n*M/(det-M)-m0/(1-m0)
 guards.update(hstar_enclosed=poly(l)<e**2<=poly(h),hH_enclosed=poly(lh)<mult*e**2<=poly(hh),roots_order=h<lh,hH_lt_width=hh<width,A_positive=Al>0,A_order=Al<=Au,B_order=bl<=bu,budget_sufficient=dl>=0,neighbor_DP_identity=neighbor==dp,center_strict=strict_margin>0,center_witness_in_finite_box=inverse_norm*radius<N,finite_ratio_strict=fd>0)
 row.update(status='certified' if all(guards.values()) else 'arithmetic_guard_failed',exact_parameters={'w':str(w),'rho':str(rho),'e':str(e),'det':str(det),'M':str(M),'mout':str(mout)},roots={'hstar':sroot,'hH':hroot},A_bounds=outward(Al,Au),B_bounds=outward(bl,bu),D_minus_B=outward(dl,dh),N=str(N),center={'coordinates':['0']*n,'radius':str(radius),'inverse_operator_norm':str(inverse_norm),'neighbor_outer_bound':str(neighbor),'defect_upper':str(defect_upper),'strict_margin':str(strict_margin),'interpretation':'actual defect <= upper; actual mean >= 1-upper, not equality'},finite_halo_minus_Ln={'exact_fraction':str(fd),'strictly_positive':fd>0},definitions={'T':'(1-e)I+eJ','G':'outer cube minus bad set','A':'G intersect {x1<c}, where |A|=M; c>0 by outer half-volume bound','E_N':'T([-N,N]^n) minus union over z in Z^n of (Tz+A)','finite_bound':'1+((2N-5)/(2N))^n M/(det-M)'})
 return row,(dl,dh,strict_margin,fd,Al,Au,bl,bu)

def dispatch(n,barrier,queue):
 start=time.perf_counter();barrier.wait();attempts=[]
 for j in (20,22):
  r,_=case(n,j);attempts.append(r)
  if r['status']=='certified':break
 queue.put({'n':n,'pid':os.getpid(),'elapsed':time.perf_counter()-start,'attempts':attempts})

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--output-dir',type=Path,default=BASE/'results');args=ap.parse_args();start=time.perf_counter()
 if args.pilot:
  r,_=case(2);print(json.dumps({'pilot_n':2,'status':r['status'],'guards':r['guards']}));return
 ctx=mp.get_context('spawn');barrier=ctx.Barrier(15);queue=ctx.Queue();ps=[ctx.Process(target=dispatch,args=(n,barrier,queue)) for n in range(2,17)]
 for p in ps:p.start()
 res=sorted((queue.get() for _ in ps),key=lambda x:x['n'])
 for p in ps:p.join();assert p.exitcode==0
 summary={'master':os.getpid(),'cases':len(res),'worker_pids':sorted(x['pid'] for x in res),'actual_worker_count':len(set(x['pid'] for x in res)),'elapsed_seconds':time.perf_counter()-start,'certified':sum(x['attempts'][-1]['status']=='certified' for x in res),'failed':sum(x['attempts'][-1]['status']!='certified' for x in res)}
 args.output_dir.mkdir(parents=True,exist_ok=True);(args.output_dir/'certificates.json').write_text(json.dumps({'summary':summary,'cases':res},indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()

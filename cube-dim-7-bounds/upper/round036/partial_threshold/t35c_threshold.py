#!/usr/bin/env python3
"""Arb covers sufficient scalar envelopes; NOT a standalone proof of R035."""
import os,sys,json,time,hashlib,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,fmpq,ctx
ctx.prec=512

def ab(q):return arb(fmpq(q.numerator,q.denominator))
def ball(a,b):return arb(ab((a+b)/2),ab((b-a)/2))
def margins(z):
 e=(-z).exp();sqrtz=z.sqrt();d=sqrtz*(-z/2).exp();h=(-z/2).exp()/(100*sqrtz);q=1-3*d/4
 mmean=z*(3*z*e+3*d/4)/(4*q)
 row=(ab(Q(13,100))+48*d*e**4+4*d*e)/(2*q)
 col=ab(Q(17,200))+d+64*d*e**4
 vmax=(d+h+8*e**5)/2
 perr=mmean+row+col+vmax+vmax/(1-vmax)*(z/4+mmean+row+col)
 bad=ab(Q(2000000,3))*z**2*sqrtz*(-z/4).exp()+ab(Q(400,121))*sqrtz*(-z/4).exp()+51200*(-5*z/4).exp()
 tau=7200*sqrtz*(-z/4).exp()
 rem=131072*((3*z+1)*(-z/2).exp()/sqrtz+(9*z*z+1)*(-3*z/2).exp()/sqrtz)
 return {'d_lt_quarter':ab(Q(1,4))-d,'x_ge_3L_over4':ab(Q(1,4))-3*z*e,'q_gt_half':q-ab(Q(1,2)),'h_zeta_lt_d':d-h-8*e**5,'target_rho_gt_4':q/(2*e)-4,'corner_target_face_gap':1/(1+3*d/4)-2*e-ab(Q(1,2)),'target_cheb_below_quarter':ab(Q(1,4))-ab(Q(20000,3))*d**3,'source_Ld_gt_one':d/e-1,'target_nd_gt_two':d/(2*e**2)-2,'max_loss_lt_one':1-vmax,'uniform_product_error_lt_two':2-perr,'source_bad_below_one36':ab(Q(1,36))-bad,'tau_below_half_row':1-tau,'W_remainder_below_delta_half':1-rem}
def run(th):
 start=time.monotonic();cpu=time.process_time();fail=[];queue=[(Q(th),Q(1024))] if th<1024 else [];boxes=0;passed=0;names=list(margins(arb(th)))
 point=margins(arb(th));badpoints=[k for k,v in point.items() if not v>0]
 if badpoints:return {'threshold':th,'status':'sufficient envelope fails at endpoint; not a theorem counterexample','failed_at_endpoint':badpoints,'pid':os.getpid(),'boxes':0,'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 while queue and boxes<8192:
  a,b=queue.pop();boxes+=1;v=margins(ball(a,b));bad=[k for k,z in v.items() if not z>0]
  if not bad:passed+=1;continue
  if b-a<Q(1,1024):fail.append({'a':str(a),'b':str(b),'unresolved':bad});continue
  mid=(a+b)/2;queue.extend([(a,mid),(mid,b)])
 return {'threshold':th,'status':'finite interval covered; tail requires analytic monotonic proof' if not queue and not fail else 'unresolved','names':names,'boxes':boxes,'passed_boxes':passed,'unresolved':fail+[{ 'a':str(a),'b':str(b),'reason':'box limit'} for a,b in queue],'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
if __name__=='__main__':
 ctx.prec=512;path=Path(__file__);start=time.monotonic();print('START',os.getpid(),'workers<=14 thresholds 64 128 256 512 1024',flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=5) as ex:res=list(ex.map(run,(64,128,256,512,1024)))
 out={'batch':'T35C','precision_bits':512,'boxes_per_threshold_cap':8192,'script_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'parent_pid':os.getpid(),'wall_seconds':time.monotonic()-start,'scope':'sufficient envelopes only; neither tail monotonicity nor full geometric proof independently audited by this script','results':res};path.with_name('t35c_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(res,indent=2))

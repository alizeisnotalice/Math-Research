"""PREPARED ONLY: no scientific run without a separately bound parent release.

Single process. Default SIGALRM177 stays armed through receipt, print and exit.
No old program, native kernel, MATLAB, peak, time grid or child is launched.
"""
import signal
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 177.0)

import sys
sys.dont_write_bytecode = True
import datetime,hashlib,json,os,resource,time,traceback
from pathlib import Path
from fractions import Fraction as Q
resource.setrlimit(resource.RLIMIT_CPU,(60,60))
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
release_path=P/'parent_release.json'
if not release_path.exists():raise SystemExit('104 PREPARED_NOT_RELEASED: refusing mathematics')
release=json.loads(release_path.read_text())
if not(release.get('authorized')is True and release.get('round')=='104' and release.get('script_sha256')==sha(Path(__file__)) and release.get('registration_sha256')==sha(P/'registration.json')):
 raise SystemExit('104 release missing or mismatched; refusing mathematics')
if any((P/n).exists()for n in ['results.json','run_receipt.json','started.json','checkpoints.json']):
 raise SystemExit('104 one attempt only; existing run artifact')
started=time.monotonic();started_cpu=time.process_time()
start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
(P/'started.json').write_text(json.dumps({'round':'104','status':'RELEASED_STARTED','pid':os.getpid(),'started_utc':start_utc,'SIGALRM_handler':'SIG_DFL','fuse177_active':True,'CPU_seconds':60,'children':0,'release_sha256':sha(release_path)},indent=2)+'\n')
completed_count=0

def run_math():
 import sympy as sp
 from flint import arb,fmpq,ctx
 ctx.prec=384
 reg=json.loads((P/'registration.json').read_text())
 assert sha(P/'upstream_constants.json')==reg['upstream_constants_sha256']
 assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
 old=json.loads((P/'upstream_constants.json').read_text())
 def A(x):
  x=Q(x);return arb(fmpq(x.numerator,x.denominator))
 def rat(x):return Q(int(x.p),int(x.q))
 def ends(x):
  l,e=x.lower().man_exp();u,f=x.upper().man_exp();e=int(e);f=int(f);g=min(e,f)
  return {'lower_mantissa':str(int(l)*2**(e-g)),'upper_mantissa':str(int(u)*2**(f-g)),'common_binary_exponent':str(g)}
 def saved(x,side):return Q(int(x[side+'_mantissa']))*Q(2)**int(x['common_binary_exponent'])
 def exact_integral(f,u,l,r):
  # Fixed positive rational pieces; elementary real log/rational primitives.
  rr=Q(0);logs={}
  for term in sp.Add.make_args(sp.apart(f,u)):
   num,den=sp.fraction(term);coef,factors=sp.factor_list(den,u)
   if not factors:
    assert not num.has(u);rr+=rat(sp.Rational(num,coef))*(r-l);continue
   assert len(factors)==1 and not num.has(u)
   pol,k=factors[0];pol=sp.Poly(pol,u);assert pol.degree()==1
   aa=rat(pol.nth(1));bb=rat(pol.nth(0));nn=rat(sp.Rational(num,coef))
   vl=aa*l+bb;vr=aa*r+bb;assert vl*vr>0
   if k==1:
    ratio=abs(vr/vl);logs[ratio]=logs.get(ratio,Q(0))+nn/aa
   else:rr+=nn/(aa*(1-k))*(vr**(1-k)-vl**(1-k))
  val=A(rr)
  for ratio,coef in logs.items():val+=A(coef)*A(ratio).log()
  return rr,logs,val
 def merge_integrals(items):
  rr=Q(0);logs={}
  for rpart,ls,_ in items:
   rr+=rpart
   for ratio,coef in ls.items():logs[ratio]=logs.get(ratio,Q(0))+coef
  val=A(rr)
  for ratio,coef in logs.items():val+=A(coef)*A(ratio).log()
  return rr,logs,val
 def serial_int(data):
  rr,logs,val=data
  return {'rational_part':str(rr),'logs':[{'coefficient':str(c),'positive_rational_argument':str(r)}for r,c in logs.items()if c],'interval':ends(val)}
 def radial(lam,delta):
  def F(w):return -w**-5/Q(5)+Q(3,4)*w**-4-w**-3+Q(1,2)*w**-2
  return (F(Q(1))-F(1-lam*delta))/lam**4
 cfg=[];pi=arb.pi();q=Q(6,5);delta=Q(1,64);v4=pi*pi/2;c=A(32)/A(3).sqrt()
 def done(row):
  global completed_count
  cfg.append(row)
  completed_count=len(cfg)
  (P/'checkpoints.json').write_text(json.dumps({'round':'104','completed_configurations':cfg,'pid':os.getpid(),'remaining_fuse_seconds':signal.getitimer(signal.ITIMER_REAL)[0],'SIGALRM_handler_is_default':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL},indent=2)+'\n')

 # 1. Exact ordered-pair/coarea structure. No quotient peak is computed.
 full=[192,480,320];removed=[12,30,20];counts=[a-b for a,b in zip(full,removed)]
 assert counts==[180,450,300] and sum(counts)==31*30==930
 dets=[20,32,36]
 Sigma=sum((n/A(d).sqrt()for n,d in zip(counts,dets)),A(0))
 D=512*A(2).sqrt()*Sigma/(pi**4*c**3)
 assert (512*A(2).sqrt()/(pi**4*c**3)-3*A(6).sqrt()/(64*pi**4)).contains(0)
 gaps=[]
 for k in range(1,6):
  g=1/(4*A(2).sqrt())-A(k).sqrt()*A(delta)
  assert g.lower()>A(delta).upper()
  gaps.append({'k':k,'all930_pair_distance_lower':ends(g)})
 done({'configuration':1,'kind':'restored four-dimensional positive source/coarea interface','ordered_pair_counts':counts,'J2_squared':dets,'pair_count':930,'Sigma_inv_J2':ends(Sigma),'D_equals_3sqrt6Sigma_over64pi4':ends(D),'all_pair_gaps':gaps,'Crem_missing_peak_loss_and_third_deletion_used':False,'actual_full_symbol_coercivity_reused_from098':'32/sqrt3'})

 # 2. One fixed core, original normal radius, true mean/rho joint budget.
 u=sp.symbols('u');left=Q(3,8);mid=Q(1,2);right=Q(5,8)
 def Zcore(mu):return sp.cancel(8/(9*mu**2*(1-u))+256/(1233*(sp.Rational(7,8)-u)**2*(1-u)))
 Zl=Zcore(u);Zr=Zcore(1-u)
 core_int=merge_integrals([exact_integral(Zl,u,left,mid),exact_integral(Zr,u,mid,right)])
 core_M=max(rat(Zl.subs(u,sp.Rational(3,8))),rat(Zl.subs(u,sp.Rational(1,2))),rat(Zr.subs(u,sp.Rational(5,8))))
 assert core_int[2].lower()>0
 Jr7=Q(old['radial7_exact'])
 Pcore=A(3).sqrt()/(128*pi**4)*(A(Q(5,14))*v4*A(delta)**A(Q(14,5)))**A(Q(5,6))
 Bcore=D*(4*v4*A(Jr7)*A(core_M)**A(Q(1,5))*core_int[2])**A(Q(5,6))
 Ncore=Pcore+Bcore;Icore=Ncore**A(q);Icore10=10*Icore
 g0=Q(1,65000);target=64*A(g0/6)**A(q);paid=A(Q(old['paid_CEH_I_upper']))
 S_core_plan=target-paid-Icore10
 old100norm=A(saved(old['norm100_core'],'upper'))
 done({'configuration':2,'kind':'restored two-deletion core all-time joint spatial budget','u_interval':['3/8','5/8'],'rho_max':'1/64','Zcore_max_exact':str(core_M),'Zcore_integral':serial_int(core_int),'radial7_exact_reused':str(Jr7),'large_branch_norm':ends(Pcore),'complement_norm':ends(Bcore),'full_G_core_norm':ends(Ncore),'I_one_core':ends(Icore),'I_ten_cores':ends(Icore10),'ratio_to_saved100_same_core_norm':ends(Ncore/old100norm),'paid_CEH_plus_core_plan_I':ends(paid+Icore10),'S_paid_CEH_plus_core_other_unpaid_zero':ends(S_core_plan),'global_I_status':'UNPAID: other middle, far, and global g0 remain unproved','operator_breakthrough':False})

 # 3. Conservative all31 full-middle, full-time analytic joint budget.
 l=Q(1,4);r=Q(3,4)
 def Zall(mu):return sp.cancel(24/mu**3+512/(75*mu*(sp.Rational(7,8)-u)**2))
 Al=Zall(u);Ar=Zall(1-u)
 all_int=merge_integrals([exact_integral(Al,u,l,mid),exact_integral(Ar,u,mid,r)])
 all_M=max(rat(Al.subs(u,sp.Rational(1,4))),rat(Al.subs(u,sp.Rational(1,2))),rat(Ar.subs(u,sp.Rational(3,4))))
 Jr20=radial(Q(20),delta);assert Jr20>0
 rows=[];I_D=A(0)
 for k,n in zip(range(1,6),[1,5,10,10,5]):
  m=6-k;shape=A(Q(m,k)).sqrt()
  Pk=256/(A(k*m).sqrt()*pi**4*c**3)*(A(Q(5,7))*v4*shape*A(delta)**A(Q(14,5)))**A(Q(5,6))
  Bk=D/(k*k*m)*(4*v4*shape*A(Jr20)*A(all_M)**A(Q(1,5))*all_int[2])**A(Q(5,6))
  Nk=Pk+Bk;Ik=n*Nk**A(q);I_D+=Ik
  rows.append({'k':k,'multiplicity':n,'Pnorm':ends(Pk),'Bnorm':ends(Bk),'full_norm_one':ends(Nk),'weighted_I':ends(Ik)})
 oldD=(A(Q(old['short051_rational_upper']))+A(Q(old['long074_rational_upper'])))**A(q)
 plan=paid+I_D;oldplan=paid+oldD
 done({'configuration':3,'kind':'all31 complete middle all-time cost diagnostic','u_interval':['1/4','3/4'],'rho_max':'1/64','common_rho_factor':'(1-20rho)^-5','Zall_integral':serial_int(all_int),'Zall_max_exact':str(all_M),'radial20_exact':str(Jr20),'cases':rows,'new_middle_I_upper':ends(I_D),'old051short_plus074long_same_middle_I_upper_expression':ends(oldD),'new_over_old_middle_expression':ends(I_D/oldD),'paid_CEH_plus_new_middle_I_plan':ends(plan),'S_new_plan_far_optimistically_zero':ends(target-plan),'S_old_plan_far_optimistically_zero':ends(target-oldplan),'U_new_plan_far_optimistically_zero':ends(plan/A(g0/6)**A(q)),'g0_status':'unproved globally','far_status':'UNPAID, zero only for optimistic certificate-expression diagnostic','actual_I_or_C6_lower_bound':False})

 # 4. Exact source witness refutes an antiparallel-only large-normal cover.
 v=sp.Matrix([1]*6);w=sp.Matrix([1,1,1,-1,-1,-1]);U=sp.Matrix.hstack(v,w)
 Proj=sp.eye(6)-U*(U.T*U).inv()*U.T
 dirs=[sp.Matrix([1,-1,-1,-1,-1,-1]),sp.Matrix([1,-1,1,1,1,1]),sp.Matrix([1,1,-1,1,1,1])]
 coeffs=[-sp.Rational(3,8),sp.Rational(3,8),sp.Rational(3,8)]
 assert all(abs(a)>sp.Rational(1,4) for a in coeffs)
 source=sp.Rational(5,8)*w+sum((a*d for a,d in zip(coeffs,dirs)),sp.zeros(6,1))
 assert source==sp.Matrix([1,1,1,sp.Rational(1,2),sp.Rational(1,2),sp.Rational(1,2)])
 normal=[a*Proj*d for a,d in zip(coeffs,dirs)]
 assert sum(normal,sp.zeros(6,1))==sp.zeros(6,1)
 for i in range(3):
  assert (normal[i].T*normal[i])[0]>0
  for j in range(i):
   assert normal[i]+normal[j]!=sp.zeros(6,1)
   assert 2*(normal[i].T*normal[j])[0]+(normal[i].T*normal[i])[0]==0
 done({'configuration':4,'kind':'strict actual-source counterexample to antiparallel-pair-only cover','nonzero_normal_large_parameter_threshold':'1/4','a_v':'0','a_w':'5/8','other_directions':[list(map(int,d))for d in dirs],'other_parameters':['-3/8','3/8','3/8'],'physical_point':['1','1','1','1/2','1/2','1/2'],'three_nonparallel_normal_contributions_sum_zero':True,'no_two_normal_contributions_cancel':True,'scope':'refutes an antiparallel-only large-label pair cover, not a valid complete fractional cover or unique-largest-label partition'})
 assert len(cfg)==reg['configuration_cap']==4
 return {'status':'PASS_FOUR_CONFIGURATIONS_SOURCE_PROOF_REVIEW_SEPARATE','round':'104','configurations':cfg,'precision_bits':384,'source_code_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json'),'upstream_constants_sha256':sha(P/'upstream_constants.json'),'new_native_kernel_peak_grid_replays':0,'children':0,'Cn_improved':False}

status='FIRST_FAILURE_STOP';failure=None
try:
 result=run_math()
 (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
 status='COMPLETED'
except Exception:
 failure=traceback.format_exc()
finally:
 # The177-second default fuse is intentionally still armed here and at exit.
 receipt={'status':status,'round':'104','pid':os.getpid(),'completed_configurations':completed_count,'started_utc':start_utc,'receipt_constructed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_before_final_receipt_write':time.monotonic()-started,'CPU_since_start_marker_before_final_receipt_write':time.process_time()-started_cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'SIGALRM_handler_is_default':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL,'fuse_not_cancelled':True,'fuse_remaining_before_final_receipt_write':signal.getitimer(signal.ITIMER_REAL)[0],'fuse_scope':'continues through receipt.write_text, print and normal process exit','children_created':0,'cleanup_required':'none; single process','retries':0,'failure':failure,'release_sha256':sha(release_path)}
 (P/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':status,'completed_configurations':completed_count,'pid':os.getpid(),'default_fuse_active_at_receipt':receipt['SIGALRM_handler_is_default']}))
if status!='COMPLETED':raise SystemExit(1)

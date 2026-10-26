"""Four new finite configurations, no old kernel/function/grid program replay."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,os,time
import sympy as sp
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def A(x):
 x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def rat(x):return Q(int(x.p),int(x.q))
def ends(x):
 l,e=x.lower().man_exp();u,f=x.upper().man_exp();e=int(e);f=int(f);g=min(e,f)
 return {'lower_mantissa':str(int(l)*2**(e-g)),'upper_mantissa':str(int(u)*2**(f-g)),'common_binary_exponent':str(g)}
def frozen(x,side):return Q(int(x[side+'_mantissa']))*Q(2)**int(x['common_binary_exponent'])
def J(a,b,c,pi):return pi*(a+b+c)/(2*a*b*c*(a+b)*(a+c)*(b+c))
def integrate_rational(f,u,l,r):
 """Exact real definite integral of this fixed positive rational function."""
 rpart=Q(0);logs={}
 for term in sp.Add.make_args(sp.apart(f,u)):
  num,den=sp.fraction(term);coef,facs=sp.factor_list(den,u)
  if not facs:
   assert not num.has(u);rpart+=rat(sp.Rational(num,coef))*(r-l);continue
  assert len(facs)==1 and not num.has(u)
  pol,k=facs[0];pol=sp.Poly(pol,u);assert pol.degree()==1
  aa=rat(pol.nth(1));bb=rat(pol.nth(0));nn=rat(sp.Rational(num,coef))
  vl=aa*l+bb;vr=aa*r+bb;assert vl*vr>0
  if k==1:
   ratio=abs(vr/vl);logs[ratio]=logs.get(ratio,Q(0))+nn/aa
  else:rpart+=nn/(aa*(1-k))*(vr**(1-k)-vl**(1-k))
 return rpart,logs
def eval_integral(rpart,logs):
 z=A(rpart)
 for ratio,coef in logs.items():z+=A(coef)*A(ratio).log()
 return z
def main():
 if (P/'results.json').exists():raise RuntimeError('refuse repeat')
 st=time.monotonic();cpu=time.process_time();ctx.prec=384
 reg=json.loads((P/'registration.json').read_text());assert sha(P/'upstream_constants.json')==reg['upstream_constants_sha256']
 old=json.loads((P/'upstream_constants.json').read_text());pi=arb.pi();cfg=[]
 eta=A(3).sqrt()/64;hv0=A(3).sqrt()/2-eta;hv1=A(Q(3,32)).sqrt()-eta
 hw0=15*A(3).sqrt()/16-eta;hw1=A(Q(3,128)).sqrt()-eta;g=1/(4*A(2).sqrt())-eta
 Crem=A(Q(44743079,2970))/pi**6;Krem=A(Q(1000,729))*Crem*(1+A(Q(6,928)).sqrt())
 prod=928*960*992
 # 1. New middle application of the existing062 three-resolvent formula.
 Bfull=Krem/prod*(J(g/A(928).sqrt(),hv0/A(992).sqrt(),hv1/A(960).sqrt(),pi)+J(g/A(928).sqrt(),hw0/A(992).sqrt(),hw1/A(960).sqrt(),pi))
 B098lo=frozen(old['B098_coarse'],'lower');B098hi=frozen(old['B098_coarse'],'upper')
 assert Bfull.upper()<A(B098lo).lower()
 cfg.append({'configuration':1,'kind':'new closed-form full-time middle complement','B_middle_full':ends(Bfull),'ratio_to_saved098_B0':ends(Bfull/A(B098lo)),'old_B0_source_dyadic':old['B098_coarse'],'formula_origin':'062 theorem section5; reused, not new theorem or replay','new_time_grid':False})

 # 2. One new joint u/rho rational majorant and its exact spatial integrals.
 u=sp.symbols('u');a=sp.Rational(7,310)*(1-u)
 bvleft=sp.Rational(6,155)*u;bvright=sp.Rational(6,155)*(1-u)
 bw=sp.Rational(6,155)*(sp.Rational(7,8)-u);cv=sp.Rational(3,112);cw=sp.Rational(7,128)
 def RR(a,b,c):return sp.cancel((a+b+c)/(a*b*c*(a+b)*(a+c)*(b+c)))
 RVl=RR(a,bvleft,cv);RVr=RR(a,bvright,cv);RW=RR(a,bw,cw)
 aa,bb,cc=sp.symbols('aa bb cc',positive=True)
 assert sp.factor(RR(aa,bb,cc)-(1/(bb*cc)+1/(aa*cc)+1/(aa*bb))/((aa+bb)*(aa+cc)*(bb+cc)))==0
 left=Q(3,8);mid=Q(1,2);right=Q(5,8);delta=Q(1,64)
 Rpoints=[rat((RVl+RW).subs(u,sp.Rational(left.numerator,left.denominator))),rat((RVl+RW).subs(u,sp.Rational(1,2))),rat((RVr+RW).subs(u,sp.Rational(right.numerator,right.denominator)))]
 Rmax=max(Rpoints)
 rational=Q(0);logs={}
 for f,l,r in [(RVl,left,mid),(RVr,mid,right),(RW,left,right)]:
  rr,ll=integrate_rational(f,u,l,r);rational+=rr
  for ratio,coef in ll.items():logs[ratio]=logs.get(ratio,Q(0))+coef
 Uintegral=eval_integral(rational,logs)
 assert Uintegral.lower()>0 and Uintegral.upper()<A(Rmax*(right-left)).lower()
 def F(w):return -w**-5/Q(5)+Q(3,4)*w**-4-w**-3+Q(1,2)*w**-2
 radial=(F(Q(1))-F(1-7*delta))/7**4
 assert radial>0
 cfg.append({'configuration':2,'kind':'new analytic joint-u-rho majorant','u_interval':['3/8','5/8'],'rho_max':'1/64','normalized_gap_lower_functions':{'a':'7(1-u)/310','bv_left':'6u/155','bv_right':'6(1-u)/155','bw':'6(7/8-u)/155','cv':'3/112','cw':'7/128'},'common_rho_factor':'(1-7rho)^-5','convex_piece_endpoint_values_exact':[str(z)for z in Rpoints],'R_max_exact':str(Rmax),'u_integral_rational_part':str(rational),'u_integral_log_terms':[{'coefficient':str(c),'log_positive_rational':str(r)}for r,c in logs.items()if c],'u_integral':ends(Uintegral),'radial_integral_exact':str(radial),'u_or_rho_quadrature_grid':False})

 # 3. Full-time G spatial budgets on the same non-thin new domain.
 q=Q(6,5);v4=pi*pi/2;volume=v4/4*A(delta)**4
 sing=v4*A(Q(5,14))*A(delta)**A(Q(14,5))
 K=A(3).sqrt()/(128*pi**4);Pnorm=K*sing**A(Q(5,6))
 C=pi*Krem/(2*prod)
 Bnorm=C*(4*v4*A(radial)*A(Rmax)**A(Q(1,5))*Uintegral)**A(Q(5,6))
 Bnew=Pnorm+Bnorm;Ione=Bnew**A(q);Iten=10*Ione
 B098uniform=Pnorm+512*volume**A(Q(5,6))
 B098actual=Pnorm+A(B098hi)*volume**A(Q(5,6))
 Bcloseduniform=Pnorm+Bfull*volume**A(Q(5,6))
 gate=Q(4,48750);g0=Q(1,65000)
 assert 64**5*(g0/6)**6==gate**6
 Itarget=64*A(g0/6)**A(q);Soptimistic=Itarget-Iten
 Uoptimistic=Iten/A(g0/6)**A(q)
 gamma_needed=(A(gate)/A(10)**A(Q(5,6))-Pnorm)/Bnorm
 assert Bnew.upper()<B098actual.lower()
 tinyI=A(frozen(old['I098_tiny_Q'],'upper'))
 cfg.append({'configuration':3,'kind':'new complete local G budget / same-domain comparison / ten-domain S diagnostic','volume_one_upper':ends(volume),'Pnorm':ends(Pnorm),'joint_Bnorm':ends(Bnorm),'new_full_G_norm':ends(Bnew),'new_I_one':ends(Ione),'new_I_ten':ends(Iten),'old098_512_same_domain_norm':ends(B098uniform),'old098_saved_B0_same_domain_norm':ends(B098actual),'closed_middle_uniform_same_domain_norm':ends(Bcloseduniform),'new_over_098_saved_B0':ends(Bnew/B098actual),'background_fraction_of_norm':ends(Bnorm/Bnew),'conditional_whole_B_gate':'4/48750','conditional_whole_I_target':ends(Itarget),'S_ten_optimistic':ends(Soptimistic),'U_ten_optimistic':ends(Uoptimistic),'constructed_ten_domain_budget_fails':bool(Soptimistic.upper()<0),'Crem_factor_needed_even_if_other_regions_zero':ends(gamma_needed),'tiny098_I_upper_over_full_target':ends(tinyI/Itarget),'actual_I_or_C6_lower_bound':False,'g0_global_status':'UNPAID','complementary_face_budget':'UNPAID'})

 # 4. Lawful actual omitted-column moment coercivity vs saved exact peaks.
 x=sp.symbols('x',real=True);f=(32-x)**3/(96-x*x)
 derivative=sp.factor(sp.diff(f,x))
 assert sp.factor(derivative+(32-x)**2*(288-64*x-x*x)/(96-x*x)**2)==0
 xstar=-32+4*A(82).sqrt();c31=((32-xstar)**3/(96-xstar**2)).sqrt()
 counts=[435,1740,435,580,1305];dets=[64,96,128,160,192]
 Rvals=[Q(old['projection_R'][z])for z in ['4110','3210','2220','3111','2211']]
 assert 6*sum(counts)==26970
 assert 24576*sum((n*r for n,r in zip(counts,Rvals)),Q(0))==Q(44743079,2970)
 coarea_sum=sum((n/A(d).sqrt()for n,d in zip(counts,dets)),A(0))
 Cmoment=196608*coarea_sum/(pi**5*c31**3)
 Cold=A(Q(1000,729))*Crem
 cfg.append({'configuration':4,'kind':'actual31-column moment coercivity and legal aggregate peak replacement test','moment_second':'32|xi|^2-|b.xi|^2','moment_fourth_upper':'96|xi|^4-|b.xi|^4','scalar_x_range':'0<=x<=6','minimizer_exact':'-32+4sqrt(82)','c31':ends(c31),'ordered_count':26970,'saved_R_and_coarea_sum_reused':True,'new_lawful_Crem_after_omission':ends(Cmoment),'old_exact_peak_Crem_after_omission':ends(Cold),'new_over_old':ends(Cmoment/Cold),'automatic_replacement_improves':bool(Cmoment.upper()<Cold.lower()),'old_model_peak_not_mixed':True})
 assert len(cfg)==reg['configuration_cap']==4
 result={'status':'PASS_FOUR_FINITE_CONFIGURATIONS_NO_GLOBAL_BOUND','configurations':cfg,'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'pid':os.getpid(),'precision_bits':384,'code_sha256':sha(Path(__file__)),'upstream_constants_sha256':sha(P/'upstream_constants.json'),'native_G_H_p_calls':0,'old_program_replays':0,'time_grid':0,'matlab_calls':0,'Cn_improved':False}
 (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'configurations':4,'seconds':result['wall_seconds'],'new_middle_Bfull':str(Bfull),'new_local_norm':str(Bnew),'ten_domain_S':str(Soptimistic),'Crem_replacement_ratio':str(Cmoment/Cold)}))
if __name__=='__main__':main()

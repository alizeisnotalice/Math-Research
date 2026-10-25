"""Six configurations, gated on pinned A041 freeze and F041 PASS."""
import json,hashlib,time,os
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def B(x):
 x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def ends(v):
 l,e=v.lower().man_exp();u,f=v.upper().man_exp();e,f=int(e),int(f);g=min(e,f)
 return dict(lower_mantissa=str(int(l)*2**(e-g)),upper_mantissa=str(int(u)*2**(f-g)),common_binary_exponent=str(g))
def formal():
 import sympy as S
 a,b,c,e,k,V=S.symbols('a b c e k V',positive=True)
 I2=lambda x,y:S.pi/(2*x*y*(x+y))
 J=lambda x,y,z:S.pi*(x+y+z)/(2*x*y*z*(x+y)*(x+z)*(y+z))
 r=S.factor((I2(a,c)-I2(b,c))/(b*b-a*a)-J(a,b,c));assert r==0
 c1=dict(id=1,name='formal_three_quadratic_integral',partial_fraction_residual=str(r))
 r2=S.factor(-S.diff(I2(a,c),a)/(2*a)-J(a,a,c));r3=S.factor(J(a,a,a)-3*S.pi/(16*a**5));assert r2==r3==0
 c2=dict(id=2,name='formal_double_and_triple_root_values',double_root_residual=str(r2),triple_root_residual=str(r3),continuity='dominated by(t^2+m^2)^-3 on compact positive parameter sets')
 R=4/(3*e**3)-16/e+32*S.atan(1/(2*e));dr=S.factor(S.diff(R,e)+1/(e**4*(S.Rational(1,4)+e*e)));lim=S.limit(R,e,S.oo);assert dr==lim==0
 loose=S.Rational(6,7)*k**S.Rational(1,6)*S.pi**S.Rational(-5,6)*e**S.Rational(7,6)+V**S.Rational(5,6)*k/(3*e**3)
 opt=S.simplify(S.diff(loose,e).subs(V,e**5/(k*S.pi)));assert opt==0
 c3=dict(id=3,name='formal_P_tail_and_single_loose_optimizer',tail_derivative_residual=str(dr),tail_at_infinity=str(lim),optimizer_residual=str(opt))
 return [c1,c2,c3],S.__version__
def main():
 if (P/'results.json').exists():raise SystemExit('Refuse repeat')
 release=P/'dependency_release.json'
 if not release.exists():raise SystemExit('WAIT: A041 frozen and F041 PASS not released; no calculations')
 dep=json.loads(release.read_text());assert dep['status']=='A041_FROZEN_F041_PASS'
 for f,h in dep['sha256'].items():assert sha(Path(f))==h,('dependency changed',f)
 src=json.loads(Path(dep['A041_results']).read_text());old=src['configurations']
 assert next(v for v in old if v['id']==2)['Crem_times_pi6_exact']=='44743079/2970'
 assert next(v for v in old if v['id']==4)['R1_exact']=='67320683/3251298631680'
 tic=time.monotonic();rows,sympy_version=formal();ctx.prec=640
 rho=Q(1,64);pi=arb.pi();R1=Q(67320683,3251298631680);D=B(R1)/(32*pi**6);Cbig=4*D
 Crem=B(Q(44743079,2970))/pi**6;V=4*B(6).sqrt()*pi*pi/45*B(rho)**5;eta=B(6).sqrt()*B(rho)
 cs=[B(6).sqrt()/2-eta,B(Q(10,3)).sqrt()/2-eta,1/B(2).sqrt()-eta];As=[992,960,928]
 assert all(c.lower()>0 for c in cs)
 eps=(Cbig*pi*V)**B(Q(1,5));vpow=V**B(Q(5,6))
 rows.append(dict(id=4,name='one_fixed_cap_and_unique_analytic_epsilon',rho=str(rho),D=ends(D),Cbig=ends(Cbig),Crem=ends(Crem),volume=ends(V),c=[ends(v) for v in cs],A=As,epsilon_star=ends(eps)))
 short=B(Q(6,7))*Cbig**B(Q(1,6))*pi**B(Q(-5,6))*eps**B(Q(7,6))
 R=4/(3*eps**3)-16/eps+32*(1/(2*eps)).atan();tail=vpow*D*R;loose_tail=vpow*Cbig/(3*eps**3)
 assert R.lower()>0 and (loose_tail-tail).lower()>0
 rows.append(dict(id=5,name='P_short_plus_complete_tail',short_norm_upper=ends(short),tail_integral=ends(R),tail_norm_upper=ends(tail),loose_tail_norm_upper=ends(loose_tail),whole_P_norm_upper=ends(short+tail)))
 aa=[cs[i]/B(As[i]).sqrt() for i in range(3)];a,b,c=aa
 J=pi*(a+b+c)/(2*a*b*c*(a+b)*(a+c)*(b+c));JB=J/(As[0]*As[1]*As[2]);bad=vpow*Crem*JB;total=short+tail+bad
 assert JB.lower()>0 and bad.lower()>0
 rows.append(dict(id=6,name='B_all_time_and_fully_paid_cap',scaled_poles=[ends(x) for x in aa],B_density_integral_factor=ends(JB),whole_B_norm_upper=ends(bad),complete_cap_norm_upper=ends(total),cap_integral_Gq_upper=ends(total**B(Q(6,5))),fraction_of_conditional_whole_face_4_over48750=ends(total/B(Q(4,48750))),whole_face_certified=False))
 assert len(rows)==6
 out=dict(status='COMPLETE_FULL_TIME_CAP_ONLY',precision_bits=640,configurations=rows,sympy=sympy_version,pid=os.getpid(),elapsed_seconds=time.monotonic()-tic,source_sha256=sha(Path(__file__)),kernel_calls=0,rho_values=1,epsilon_values=1,whole_face_I_certified=False,g_floor_certified=False,Hnorm_goal_certified=False,new_C6_bound=False)
 (P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],configurations=6,seconds=out['elapsed_seconds'],epsilon_display=str(eps),Pshort_display=str(short),Ptail_display=str(tail),Bfull_display=str(bad),total_display=str(total))))
if __name__=='__main__':main()

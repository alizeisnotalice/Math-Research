"""Four registered small configurations; no G/H/p kernel call or old replay."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib, itertools, json, os, time
import sympy as s
from flint import arb, fmpq, ctx

P=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def A(x):
    x=Q(x)
    return arb(fmpq(x.numerator,x.denominator))
def ends(x):
    l,e=x.lower().man_exp();u,f=x.upper().man_exp();e=int(e);f=int(f);g=min(e,f)
    return {'lower_mantissa':str(int(l)*2**(e-g)), 'upper_mantissa':str(int(u)*2**(f-g)), 'common_binary_exponent':str(g)}
def dyadic(x,which):
    return Q(int(x[which+'_mantissa']))*Q(2)**int(x['common_binary_exponent'])

def main():
    if (P/'results.json').exists(): raise RuntimeError('one attempt only: result already exists')
    st=time.monotonic();cpu=time.process_time();ctx.prec=384
    reg=json.loads((P/'registration.json').read_text())
    assert sha(P/'upstream_constants.json')==reg['upstream_constants_sha256']
    frozen=json.loads((P/'upstream_constants.json').read_text())
    configs=[];pi=arb.pi()

    # 1. One actual quotient, exact coarea, and its Rademacher moments.
    V=[s.Matrix([1,*z]) for z in itertools.product([-1,1],repeat=5)]
    v=s.Matrix([1]*6);w=s.Matrix([1,1,1,-1,-1,-1]);a=s.Matrix([1,1,-1,1,1,-1])
    N=s.Matrix.hstack(v,w);Proj=s.eye(6)-N*(N.T*N).inv()*N.T
    b=Proj*a;Gram=sum((z*z.T for z in V),s.zeros(6))
    assert Gram==32*s.eye(6)
    assert (N.T*N).det()==36 and s.Matrix.hstack(v,w,a).T*s.Matrix.hstack(v,w,a)==s.Matrix([[6,0,2],[0,6,0],[2,0,6]])
    assert s.simplify((b.T*b)[0]-s.Rational(16,3))==0
    assert Proj*Gram*Proj==32*Proj
    x=s.symbols('x0:6');m2=0;m4=0
    for z in V:
        lin=sum(z[j]*x[j] for j in range(6));m2+=lin**2;m4+=lin**4
    n2=sum(z*z for z in x)
    assert s.expand(m2/32-n2)==0
    assert s.expand(m4/32-(3*n2**2-2*sum(z**4 for z in x)))==0
    t,r,u=s.symbols('t r u',positive=True)
    assert s.simplify((t**2*t**-4).subs(t,r/u)*r/u**2-1/r)==0
    configs.append({'configuration':1,'kind':'actual k3 quotient / coarea / moments / Mellin substitution','J2_squared':'36','projected_a_squared':'16/3','triple_Gram_determinant':'192','frame_Gram':'32I6','quotient_Gram':'32I4','nonzero_quotient_columns':30,'second_moment':'sum x_i^2','fourth_moment':'3(sum x_i^2)^2-2sum x_i^4','coercivity':'32/sqrt(3)','time_to_line_prefactor':'1/r'})

    # 2. New time-first coefficient; stored projection peak is only read.
    R=Q(frozen['R3_2211']);q3=s.sqrt(192)*s.Rational(R.numerator,R.denominator)/(2*s.pi)**3
    assert s.simplify(q3-s.sqrt(3)*s.Rational(R.numerator,R.denominator)/s.pi**3)==0
    K=A(3).sqrt()/(128*pi**4)
    assert K.upper()<A(Q(1,7000)).lower()
    Q3e=A(192).sqrt()*A(R)/(8*pi**3)
    cone_factor=A(Q(20,19))**3
    Kcone=1024*cone_factor*Q3e/(12*pi**2)
    M3=frozen['M3_saved_dyadic']
    M3lo=dyadic(M3,'lower');M3hi=dyadic(M3,'upper')
    assert M3lo/2>Q(1,5000)
    configs.append({'configuration':2,'kind':'new direct time integral coefficients / saved peak conversion','H_identity':'Q3(e)/(2*r)','uniform_H_coefficient':ends(1/(2*pi**2*(A(32)/A(3).sqrt())**3)),'uniform_P_G_coefficient':ends(K),'rational_K_upper':'1/7000','direction_e0_Q3':ends(Q3e),'cone_radius_chordal':'1/64','cone_symbol_loss':'1/20','cone_P_G_coefficient':ends(Kcone),'saved_M3_half_lower_exact':str(M3lo/2),'saved_M3_half_above':'1/5000','new_no_kernel_value':True})

    # 3. All-time complementary branch, and one fixed inner tube.
    eta=A(3).sqrt()/64
    hv0=A(3).sqrt()/2-eta;hv1=A(Q(3,32)).sqrt()-eta
    hw0=15*A(3).sqrt()/16-eta;hw1=A(Q(3,128)).sqrt()-eta
    g=1/(4*A(2).sqrt())-eta
    assert min(hv0.lower(),hv1.lower(),hw0.lower(),hw1.lower(),g.lower())>0
    Crem=A(Q(44743079,2970))/pi**6
    Krem=A(Q(1000,729))*Crem*(1+A(Q(6,928)).sqrt())
    B0=pi*Krem/(2*A(928).sqrt()*g)*(1/(hv0**2*hv1**2)+1/(hw0**2*hw1**2))
    assert B0.upper()<512
    delta=Q(1,2**24);delta_outer=Q(1,64)
    # u is in [31/64,33/64], with |x-x(u)| <= sqrt(3)/64.
    Rmax=(A(3*(1+Q(33,64)**2)).sqrt()+eta)**2
    dvmax=(A(Q(3,2)*(1-Q(31,64))**2).sqrt()+eta)**2
    assert Rmax.upper()<4 and dvmax.upper()<A(Q(1,2)).lower()
    ratio=5000*(Q(1,7000)+512*delta)
    assert ratio<1
    configs.append({'configuration':3,'kind':'full complement and strict local estimate comparison','B0':ends(B0),'rational_B0_upper':'512','outer_rho':'1/64','inner_rho':str(delta),'mean_B_interval':['31/64','33/64'],'r_squared_upper_check':ends(Rmax),'dv_squared_upper_check':ends(dvmax),'old_Phi_lower':'1/(5000*rho)','new_G_upper':'1/(7000*rho)+512','new_over_old_uniform_upper_exact':str(ratio),'strictly_better_than_old_Phi':True,'all_positive_times_paid_for_local_G':True})

    # 4. Spatial q budget only after full time integration. No full-face I.
    q=Q(6,5);v4=pi*pi/2
    vol=v4/32*A(delta)**4
    singular=v4*A(Q(5,112))*A(delta)**A(Q(14,5))
    B=A(Q(1,7000))*singular**A(Q(5,6))+512*vol**A(Q(5,6))
    I=B**A(q)
    gate=A(Q(4,48750));g0=Q(1,65000)
    Itarget=64*A(g0/6)**A(q)
    assert I.upper()<A(Q(1,10**23)).lower()
    configs.append({'configuration':4,'kind':'one new local full-G spatial budget and S interface','volume_upper':ends(vol),'rho_power_integral_upper':ends(singular),'local_B_upper':ends(B),'local_I_upper':ends(I),'local_I_rational_upper':'1/10^23','conditional_full_B_gate_exact':'4/48750','local_B_over_gate':ends(B/gate),'conditional_I_target':ends(Itarget),'S_status':'UNPAID: g_L and complete complementary I absent','same_model':'original equal-weight body32 alpha=1'})
    assert len(configs)==reg['configuration_cap']==4
    result={'status':'PASS_FOUR_CONFIGURATIONS_THEOREM_REQUIRES_REVIEW','configurations':configs,'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'precision_bits':384,'pid':os.getpid(),'new_G_H_p_calls':0,'old_native_or_baseline_replays':0,'new_time_grid':0,'matlab_calls':0,'source_code_sha256':sha(Path(__file__)),'upstream_constants_sha256':sha(P/'upstream_constants.json'),'C6_improved':False}
    (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'configurations':4,'seconds':result['wall_seconds'],'B0':str(B0),'K':str(K),'local_B':str(B),'local_I':str(I),'local_ratio_exact':str(ratio)}))

if __name__=='__main__':main()

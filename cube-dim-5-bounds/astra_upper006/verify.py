"""Exact scalar pilot; never promotes a corner value to a global g bound."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, floor
import hashlib
import json
import sys
import time
from flint import arb, ctx

ROOT=Path(__file__).resolve().parent


def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def main():
    began=time.monotonic();ctx.prec=256
    source=ROOT/'corner_moments.json'
    data=json.loads(source.read_text())
    assert data['completed_W_power']==7
    moments={int(k):list(map(F,v)) for k,v in data['moments'].items()}
    old_path=ROOT.parent/'astra_upper004'/'corner_checkpoint.json'
    old=json.loads(old_path.read_text())
    assert moments[1]==list(map(F,old['moments']['E_d_W_power'][:8]))
    assert moments[2]==list(map(F,old['moments']['E_d_squared_W_power'][:8]))
    b=list(map(F,old['moments']['E_b_corner_W_power']))
    for k in range(1,5):
        assert moments[k][1]==0
        assert all(moments[k][n]>=0 for n in range(8))
        assert 0<moments[k][0]<=1
        for n in (0,2,4,6):
            if k>1:assert moments[k][n]<=moments[k-1][n]
    for k in (2,3):
        for n in (0,2,4,6):
            assert moments[k][n]**2<=moments[k-1][n]*moments[k+1][n]
    m1,m2=moments[1][0],moments[2][0]
    tau=m2/m1;s=F(5,11);rho=1-s
    assert tau==F(15800608,26277615) and F(1,2)<tau<1
    a=[F(1)]
    for k in range(1,4):a.append(a[-1]*(rho+k-1)/k)
    assert a==[F(1),F(6,11),F(51,121),F(476,1331)]
    # T^11=tau^-6: the two rational inequalities below independently
    # certify the rounded endpoint without relying on floating evaluation.
    T=ar(tau)**(-ar(rho))
    Tlo=F(floor(float(T)*10**12)-1,10**12)
    Thi=Tlo+F(3,10**12)
    assert Tlo>0 and Tlo**11*tau**6<1<Thi**11*tau**6
    assert ar(Tlo)<T<ar(Thi)
    normalized=[F(0)]*5
    for k in range(4):
        for j in range(k+1):normalized[j+1]+=a[k]*comb(k,j)*(-1/tau)**j
    # Positive normalized polynomial, proved using u in (-1,1):
    # ell_1/(d T)>=s, ell_3-ell_1=d T a2 u²(1+(28/33)u)>=0.
    assert 1-rho==s and a[3]/a[2]==F(28,33)<1
    eta=F(1,8)
    Ubase=sum((F(comb(7,n))*eta**n*sum((normalized[k]*moments[k][n] for k in range(1,5)),F(0)) for n in range(8)),F(0))
    U=Tlo*Ubase
    V=sum((F(comb(11,n))*eta**n*b[n] for n in range(12)),F(0))
    assert Ubase>0 and U>0 and V>0
    # T-dependent expectation gains are also exact rational after cancellation.
    first_normalized=[F(0),1+rho,-rho/tau]
    U1base=sum((F(comb(7,n))*eta**n*sum((first_normalized[k]*moments[k][n] for k in (1,2)),F(0)) for n in range(8)),F(0))
    assert Ubase>U1base>0
    optimal_a=(ar(U)/ar(V))**(arb(1)/4)
    scale=F(floor(float(optimal_a)*1000),1000)
    assert scale**4<=U/V<(scale+F(1,1000))**4
    value=F(11,4)*U*scale**7-F(7,4)*V*scale**11
    constants_path=ROOT.parent/'astra_upper003'/'verification.json'
    prior=json.loads(constants_path.read_text())
    cL=F(prior['constants']['c_lower'])
    c=3/(2**(arb(15)/2)*arb.pi()**(arb(5)/2))
    assert ar(cL)<c
    exact_corner_lower=cL*value
    target=F(4583,10**8)
    # All-a corner ceiling is for this fixed eta and rounded ell3 only.
    ceiling=c*ar(U)**(arb(11)/4)/ar(V)**(arb(7)/4)
    output={'status':'one_necessary_corner_pilot_only','new_global_g':False,'new_C5_upper_bound':False,
            'tau':str(tau),'s':str(s),'Taylor_u_coefficients':[str(x) for x in a],
            'normalized_d_coefficients':[str(x) for x in normalized],
            'T_lower':str(Tlo),'T_upper':str(Thi),'T_Arb':str(T),
            'T_endpoints_pure_rational_power_check':True,
            'moments':{str(k):[str(x) for x in v] for k,v in moments.items()},
            'E_ell3_normalized_minus_ell1_normalized_at_f7':str(Ubase-U1base),
            'eta':str(eta),'a_selected':str(scale),'a_star_Arb':str(optimal_a),
            'exact_U':str(U),'exact_V':str(V),'corner_polynomial_value':str(value),
            'c_lower':str(cL),'exact_corner_g_lower':str(exact_corner_lower),
            'corner_g_lower_Arb':str(ar(exact_corner_lower)),
            'all_a_corner_ceiling_Arb':str(ceiling),
            'necessary_screen_target':str(target),'passes_necessary_screen':exact_corner_lower>target,
            'crosschecks':{'all_16_old_mixed_moments_equal':True,'odd_first_moment_zero':True,
                           'nonnegative_moments':True,'even_weighted_d_monotone_and_logconvex':True},
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,old_path,constants_path,ROOT/'registration.md',ROOT/'corner_moments.cpp',Path(__file__)]},
            'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ['eta','a_selected','corner_g_lower_Arb','all_a_corner_ceiling_Arb','passes_necessary_screen','runtime']},indent=2))


if __name__=='__main__':main()

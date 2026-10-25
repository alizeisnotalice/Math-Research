#!/usr/bin/env python3
"""Certify two fixed Young-shape budgets by a determinant-Laplace Cauchy bound."""
import json,os,time
from fractions import Fraction as F
import sympy as sp
from flint import arb,ctx
from reciprocal_linear_screen import P,RHOS,ar,dump,digest

START=time.monotonic()

def main():
    reg=P/'cauchy_registration.json'
    assert json.loads(reg.read_text())['status']=='registered_before_computation'
    ctx.prec=256
    # det(A) Laplace factor 1/(2 sqrt(t))^5, then the H^-5 Mellin integral.
    integral_constant=sp.simplify(2*sp.gamma(5)/(32*sp.gamma(5)))
    assert integral_constant==sp.Rational(1,16)
    inverse_first=sp.gamma(3)/sp.gamma(2)
    assert inverse_first==2
    zinv=16*sp.gamma(sp.Rational(3,2))
    cauchy_constant=sp.simplify(zinv**2*integral_constant*inverse_first)
    assert cauchy_constant==8*sp.pi
    result=[]
    c1=1/(32*arb.pi()**3)
    supplied_target=ar(F(15198363,100000000000))
    weaker_target=ar(F(15198,100000000))
    for rho in RHOS:
        C=ar(1-rho).sqrt()+5*ar(1-rho/3).sqrt()+10
        transverse=12+4/ar(1-rho/3).sqrt()
        parallel=2+9/ar(1-rho/3).sqrt()+5/ar(1-rho).sqrt()
        detD=transverse**4*parallel
        Uupper=(8*arb.pi()*detD/C**7).sqrt()
        V=ar(F(128,7))*transverse**4/C**5
        gupper=c1*Uupper**3/V**2
        Ucrit=(supplied_target*V**2/c1)**(arb(1)/3)
        assert Uupper<Ucrit
        assert gupper<weaker_target<supplied_target
        result.append({'rho':str(rho),'D_parallel':str(parallel),'det_D':str(detD),
                       'U_upper':str(Uupper),'V':str(V),'Ucrit':str(Ucrit),
                       'all_a_corner_Young_upper':str(gupper),
                       'U_upper_below_Ucrit':bool(Uupper<Ucrit),
                       'strictly_below_weaker_target_15198e_minus8':bool(gupper<weaker_target),
                       'g_upper_float':float(gupper)})
        dump('cauchy_checkpoint.json',{'pid':os.getpid(),'completed_rhos':len(result),'rows':result})
    out={'status':'both_fixed_rho_Young_families_strictly_excluded_at_corner',
         'proof_constants':{'E_detA_H_minus5':'detD/(16 C^5)','E_H_minus1':'2/C^2','U_squared_upper':'8 pi detD/C^7'},
         'rows':result,'source_script':digest(__file__),'registration':digest(reg),
         'shared_source':digest(P/'reciprocal_linear_screen.py'),
         'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-START,'Arb_bits':ctx.prec},
         'new_global_g':False,'new_C5_upper_bound':False,
         'not_an_upper_bound_on_Green':True}
    dump('cauchy_upper.json',out);print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()

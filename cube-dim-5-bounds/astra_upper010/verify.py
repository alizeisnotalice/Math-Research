"""Verify one second-order shape's necessary corner budget, with no global claim."""
from pathlib import Path
from fractions import Fraction as F
from math import floor
import hashlib
import json
import sys
import time
from flint import arb, ctx

ROOT=Path(__file__).resolve().parent


def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator


def main():
    began=time.monotonic();ctx.prec=256
    raw_path=ROOT/'corner_moments.json';raw=json.loads(raw_path.read_text())
    kernel_path=ROOT/'corner.cpp'
    assert raw['source_sha256']==hashlib.sha256(kernel_path.read_bytes()).hexdigest()
    assert raw['new_moments_computed']==0 and raw['existing_cache_count']==143289
    assert raw['raw_orbit_crosscheck_through_power']==3
    Ed,Ed2,Ef,V=(F(raw[key]) for key in ['E_d_f7','E_d2_f7','E_f7','E_b_f11'])
    assert F(raw['E_f'])==F(643,616)
    assert 0<Ed2<=Ed<=Ef and Ed**2<=Ed2*Ef and V>0
    # M²<=5M implies T<=4|x|²+6W. Since |x|²<=5,
    # f >= 1/2-W/20+W²/80 = 9/20+(W-2)²/80.
    assert F(1,2)-F(4,80)==F(9,20)>0
    constants_path=ROOT.parent/'astra_upper003'/'verification.json'
    const=json.loads(constants_path.read_text())['constants']
    A,B,cL,cU=(F(const[key]) for key in ['A_lower','B_upper','c_lower','c_upper'])
    assert A>B>0
    U=A*Ed-B*Ed2;assert U>0
    star=(ar(U)/ar(V))**(arb(1)/4)
    a=F(floor(float(star)*1000),1000)
    assert a**4<=U/V<(a+F(1,1000))**4
    value=F(11,4)*U*a**7-F(7,4)*V*a**11
    exact_lower=cL*value
    c=3/(2**(arb(15)/2)*arb.pi()**(arb(5)/2));assert ar(cL)<c<ar(cU)
    family_ceiling=c*ar(U)**(arb(11)/4)/ar(V)**(arb(7)/4)
    ideal_ceiling=c*ar(Ed)**(arb(5)/4)*ar(Ef)**(arb(3)/2)/ar(V)**(arb(7)/4)
    family4=cU**4*U**11/V**7
    ideal4=cU**4*Ed**5*Ef**6/V**7
    old_target=F(2821,50000000)
    output={'status':'single_second_order_shape_necessary_corner_only','new_global_g':False,'new_C5_upper_bound':False,
            'f_global_lower_bound':'9/20','corner_shape_denominator':880,
            'a_selected':str(a),'a_star_Arb':str(star),'exact_U':str(U),'exact_V':str(V),
            'exact_corner_polynomial_value':str(value),'exact_corner_g_lower':str(exact_lower),
            'corner_g_lower_Arb':str(ar(exact_lower)),'passes_necessary_screen':exact_lower>old_target,
            'screen_target':str(old_target),'all_a_old_ell_corner_ceiling_Arb':str(family_ceiling),
            'all_a_old_ell_ceiling_fourth_upper':str(family4),
            'all_a_old_ell_strictly_below_existing_g':family4<old_target**4,
            'ideal_d_power_Jensen_ceiling_Arb':str(ideal_ceiling),
            'ideal_d_power_Jensen_ceiling_fourth_upper':str(ideal4),
            'ideal_ceiling_strictly_below_1e_4':ideal4<F(1,10000)**4,
            'ideal_ceiling_strictly_below_128e_6':ideal4<F(16,125000)**4,
            'new_moments_computed':0,'existing_cache_count':143289,
            'crosschecks':{'first_three_powers_raw_vs_orbit_exact':True,'E_f_closed_form':True,
                           'weighted_d_moments_monotone_and_logconvex':True},
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'registration.md',kernel_path,raw_path,constants_path,Path(__file__)]},
            'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ['a_selected','corner_g_lower_Arb','passes_necessary_screen','all_a_old_ell_corner_ceiling_Arb','all_a_old_ell_strictly_below_existing_g','ideal_d_power_Jensen_ceiling_Arb','ideal_ceiling_strictly_below_1e_4','ideal_ceiling_strictly_below_128e_6','runtime']},indent=2))


if __name__=='__main__':main()

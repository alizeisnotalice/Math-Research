"""Weighted Jensen ceiling for ideal determinant power, not a Green upper bound."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
import hashlib
import importlib.util
import json
import sys
import time
from flint import arb, ctx

ROOT=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def ar(x):
    x=F(x);return arb(x.numerator)/x.denominator


def choose(x,n):
    q=F(1)
    for j in range(n):q*=F(x-j,j+1)
    return q


def rising(x,n):
    q=F(1)
    for j in range(n):q*=x+j
    return q


def raw_W(order):
    beta=F(3,4)
    denominator=[sum((F(count,16)*choose(beta,n)*(-weight)**n for count,weight in [(1,10),(5,2),(10,-2)]),F(0)) for n in range(order+1)]
    reciprocal=[F(1)]
    for n in range(1,order+1):reciprocal.append(-sum((denominator[j]*reciprocal[n-j] for j in range(1,n+1)),F(0)))
    return [reciprocal[n]*factorial(n)/rising(beta,n) for n in range(order+1)]


def main():
    began=time.monotonic();ctx.prec=256
    helper_path=ROOT.parent/'astra_upper005'/'checkpoint.py';helper=load('locked005_ceiling_helper',helper_path)
    old_path=ROOT.parent/'astra_upper005'/'checkpoint.json';old=json.loads(old_path.read_text())
    d=list(map(F,old['moments_to_order_22']['E_d_W_power']))
    b=list(map(F,old['moments_to_order_22']['E_b_corner_W_power']))
    w=raw_W(22);assert w[:4]==[F(1),F(0),F(10,7),F(300,77)]
    mixed_path=ROOT.parent/'astra_upper004'/'corner_checkpoint.py';mixed=load('locked004_mixed_helper',mixed_path);mixed.ORDER=22
    assert tuple(w)==mixed.mixed((0,0,0),(0,0,0))
    cU=F(old['constants']['c_upper']);cL=F(old['constants']['c_lower'])
    c=3/(2**(arb(15)/2)*arb.pi()**(arb(5)/2));assert c<ar(cU)
    A=[F(comb(7,n),2**n)*d[n] for n in range(8)]
    B=[F(comb(7,n),2**n)*w[n] for n in range(8)]
    V=[F(comb(11,n),2**n)*b[n] for n in range(12)]
    right=helper.pmul(helper.ppow(A,5),helper.ppow(B,6));left=helper.ppow(V,7)
    linear=[]
    for target in [F(1,10000),F(16,125000)]:
        residual=[target**4*x-cU**4*y for x,y in zip(left,right)]
        assert len(left)==len(right)==len(residual)==78
        result=helper.certify_positive(residual)
        result.update({'target':str(target),'degree':77,'scope':'all a>=0, 0<=eta<=1/2, ideal d^(5/11)'})
        # Include exact residual even if this particular Bernstein budget fails.
        result['power_coefficients']=[str(x) for x in residual]
        linear.append(result)
        if result['certified']:break
    # Diagnostics at the already selected eta=1/8 and the family endpoint:
    # a negative exact residual distinguishes a loose ceiling from failure
    # caused only by the finite subdivision budget. No new shape is searched.
    linear_diagnostics=[]
    for eta in [F(1,8),F(1,2)]:
        aa=sum((F(comb(7,n))*eta**n*d[n] for n in range(8)),F(0))
        bb=sum((F(comb(7,n))*eta**n*w[n] for n in range(8)),F(0))
        vv=sum((F(comb(11,n))*eta**n*b[n] for n in range(12)),F(0))
        lower4=cL**4*aa**5*bb**6/vv**7
        upper4=cU**4*aa**5*bb**6/vv**7
        linear_diagnostics.append({'eta':str(eta),'exact_A':str(aa),'exact_B':str(bb),'exact_V':str(vv),
                                   'ceiling_Arb':str(c*ar(aa)**(arb(5)/4)*ar(bb)**(arb(3)/2)/ar(vv)**(arb(7)/4)),
                                   'ceiling_fourth_lower':str(lower4),'ceiling_fourth_upper':str(upper4),
                                   'ceiling_strictly_above_128e_6':lower4>F(16,125000)**4,
                                   'ceiling_strictly_below_128e_6':upper4<F(16,125000)**4})
    quadratic=[]
    for kappa in [F(1,512),F(1,256),F(1,128)]:
        f=[F(1),F(1,8),-kappa];f7=helper.ppow(f,7);f11=helper.ppow(f,11)
        aa=sum((coef*d[n] for n,coef in enumerate(f7)),F(0))
        bb=sum((coef*w[n] for n,coef in enumerate(f7)),F(0))
        vv=sum((coef*b[n] for n,coef in enumerate(f11)),F(0))
        assert aa>0 and bb>0 and vv>0
        fourth=cU**4*aa**5*bb**6/vv**7
        ceiling=c*ar(aa)**(arb(5)/4)*ar(bb)**(arb(3)/2)/ar(vv)**(arb(7)/4)
        quadratic.append({'kappa':str(kappa),'exact_A':str(aa),'exact_B':str(bb),'exact_V':str(vv),
                          'upper_ceiling_fourth_power':str(fourth),'ceiling_Arb':str(ceiling),
                          'strictly_below_1e_4':fourth<F(1,10000)**4,
                          'strictly_below_128e_6':fourth<F(16,125000)**4})
    output={'status':'weighted_Jensen_ideal_d_power_necessary_corner_ceiling','new_g_bound':False,'new_C5_upper_bound':False,
            'not_an_upper_bound_for_G':True,'raw_W_moments':[str(x) for x in w],
            'raw_moment_crosschecks_passed':True,'c_upper':str(cU),'linear_family':linear,
            'linear_diagnostics':linear_diagnostics,'three_fixed_quadratic_shapes':quadratic,
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'registration.md',Path(__file__),helper_path,old_path,mixed_path]},
            'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'linear':[{k:v for k,v in x.items() if k not in ['power_coefficients','partition']}|{'leaf_count':len(x.get('partition',[]))} for x in linear],
                      'linear_diagnostics':[{k:v for k,v in x.items() if k in ['eta','ceiling_Arb','ceiling_strictly_above_128e_6','ceiling_strictly_below_128e_6']} for x in linear_diagnostics],
                      'quadratic':[{k:v for k,v in x.items() if k in ['kappa','ceiling_Arb','strictly_below_1e_4','strictly_below_128e_6']} for x in quadratic],
                      'runtime':output['runtime']},indent=2))


if __name__=='__main__':main()

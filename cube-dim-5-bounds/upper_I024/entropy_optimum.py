#!/usr/bin/env python3
"""One closed log-sum optimum batch; reads022 types, never reruns its geometry."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
import time

import flint
from flint import arb, ctx

HERE=Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def iroot(n,k):
    lo,hi=0,1<<((n.bit_length()+k-1)//k)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**k<=n:
            lo=mid
        else:
            hi=mid
    return hi if hi**k==n else lo


def root(x,k,digits=100):
    x=F(x)
    scale=10**digits
    j=iroot(x.numerator*scale**k//x.denominator,k)
    lo=F(j,scale)
    hi=lo if lo**k==x else F(j+1,scale)
    assert lo**k<=x<=hi**k
    return lo,hi


def power_interval(bounds,p,k):
    return root(bounds[0]**p,k)[0],root(bounds[1]**p,k)[1]


def decimal(x,up=False,digits=48):
    scale=10**digits
    n=(x*scale).__ceil__() if up else (x*scale).__floor__()
    sign='-' if n<0 else ''
    whole,rem=divmod(abs(n),scale)
    return f'{sign}{whole}.{rem:0{digits}d}'


def encode(bounds):
    return {'lower_fraction':str(bounds[0]),'upper_fraction':str(bounds[1]),
            'decimal_outward':[decimal(bounds[0]),decimal(bounds[1],True)]}


def arb_bounds(x):
    return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))


def parse(s):
    c,r=map(F,s[1:-1].split(' +/- '))
    return c-r,c+r


def main():
    began=time.monotonic()
    registration=json.loads((HERE/'registration.json').read_text())
    for source in registration['sources']:
        assert sha(Path(source['source']))==source['sha256']
    source_path=next(Path(x['source']) for x in registration['sources'] if x['source'].endswith('upper_I022/results.json'))
    old=json.loads(source_path.read_text())
    oldJ=parse(old['arb_384bit']['J_upper'])
    oldI=parse(old['arb_384bit']['complete_I_upper'])

    # Algebraic Z, independently enclosed with integer power/root comparisons.
    sqrt2,sqrt3=root(F(2),2),root(F(3),2)
    cubert34=root(F(3,4),3)
    A=(1+sqrt2[0],1+sqrt2[1])
    Apower=power_interval(A,4,3)
    B=(7+8*sqrt2[0]+3*sqrt3[0],7+8*sqrt2[1]+3*sqrt3[1])
    Bpower=power_interval(B,4,3)
    coeff=((5+3*cubert34[0])/16,(5+3*cubert34[1])/16)
    Z=(F(3,2)+coeff[0]*Apower[0]+Bpower[0]/256,
       F(3,2)+coeff[1]*Apower[1]+Bpower[1]/256)
    assert 0<Z[0]<Z[1]
    merged_thetas=[F(1,8),F(3,16),F(5,16),F(3,8)]
    Rpower=[(Bpower[0]/32,Bpower[1]/32),
            (cubert34[0]*Apower[0],cubert34[1]*Apower[1]),
            Apower,(F(4),F(4))]
    Zcheck=(sum(w*r[0] for w,r in zip(merged_thetas,Rpower)),
            sum(w*r[1] for w,r in zip(merged_thetas,Rpower)))
    assert Zcheck[0]==Z[0] and Zcheck[1]==Z[1]
    q_bounds=[(theta*value[0]/Z[1],theta*value[1]/Z[0])
              for theta,value in zip(merged_thetas,Rpower)]
    assert sum(x[0] for x in q_bounds)<=1<=sum(x[1] for x in q_bounds)

    ctx.prec=384
    # Direct sum over five022 types is only a cross-check of the new algebraic Z.
    Zdirect=arb(0)
    assert sum(F(row['theta']) for row in old['geometry']['probability_types'])==1
    for row in old['geometry']['probability_types']:
        R=arb(1)
        for histogram in row['probability_histograms']:
            R*=sum((int(count)*ar(F(prob))**ar(F(3,4)) for prob,count in histogram.items()),arb(0))
        Zdirect+=ar(F(row['theta']))*R**ar(F(4,3))
    direct_bounds=arb_bounds(Zdirect)
    assert Z[0]<=direct_bounds[0]<=direct_bounds[1]<=Z[1]

    gamma=ar(F(7,4)).gamma()
    moment=ar(F(6,5)).gamma()/ar(F(23,20)).gamma()
    # Jopt^20=2^20 Gamma^20 m^100/Z^3.
    NJ=arb_bounds(ar(2)**20*gamma**20*moment**100)
    J20=NJ[0]/Z[1]**3,NJ[1]/Z[0]**3
    Jopt=root(J20[0],20)[0],root(J20[1],20)[1]
    # Iopt^14=3^20 Gamma^20 m^100/[2^130 60^14 pi^22 Z^3].
    NI=arb_bounds(ar(3)**20*gamma**20*moment**100/(ar(2)**130*ar(60)**14*arb.pi()**22))
    I14=NI[0]/Z[1]**3,NI[1]/Z[0]**3
    Iopt=root(I14[0],14)[0],root(I14[1],14)[1]
    cap=F((Iopt[1]*10**13).__ceil__(),10**13)
    assert Iopt[1]<cap
    checks={'Jopt_budget_improves022':Jopt[1]<oldJ[0],
            'Iopt_budget_improves022':Iopt[1]<oldI[0],
            'family_minimum_I_budget_exceeds_11e_minus6':Iopt[0]>F(11,10**6),
            'family_minimum_I_budget_exceeds_118e_minus7':Iopt[0]>F(118,10**7)}
    out={'status':'all_flag_weight_AMGM_ordinary_Holder_budget_optimum_certified',
         'actual_J_or_I_lower_bound':False,'new_C5_upper_bound':False,
         'script_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),
         'source_inputs':registration['sources'],
         'family':{'base_flags':4560,'all_weights_allowed':'q>=0, sumq=1; not limited to symmetric type weights',
                   'theta_base':'Cflag/(60W)','Z_formula':'3/2+(5+3*(3/4)^(1/3))/16*(1+sqrt(2))^(4/3)+(7+8sqrt(2)+3sqrt(3))^(4/3)/256',
                   'qopt':'theta_f*R_f^(4/3)/Z','optimal_J_budget':'2Gamma(7/4)*m(3/20)^5*Z^(-3/20)',
                   'scope':'1; ordinary Holder; fixed022 normalized factorization'},
         'Z_exact_interval':encode(Z),'Z_direct_Arb_crosscheck':Zdirect.str(90),
         'merged_type_qopt_intervals':[{'base_theta':str(theta),'R_power_4_over_3':encode(value),'qopt_mass':encode(q)}
                                      for theta,value,q in zip(merged_thetas,Rpower,q_bounds)],
         'J_optimal_budget_interval':encode(Jopt),'I_optimal_budget_interval':encode(Iopt),
         'strict_rational_I_upper':str(cap),'checks':checks,
         'proof_identity':'log budget(q)=log[2Gamma*m^5]-b log Z+b KL(q||qopt)',
         'runtime':{'python':sys.executable,'flint':flint.__version__,'Arb_bits':ctx.prec,
                    'Fraction_decimal_root_digits':100,'seconds':time.monotonic()-began,
                    'geometry_recomputed':False,'parameter_optimization_calls':0,'spatial_or_probability_samples':0}}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'Z':encode(Z)['decimal_outward'],
                      'Jopt':encode(Jopt)['decimal_outward'],'Iopt':encode(Iopt)['decimal_outward'],
                      'rational_I_upper':str(cap),'checks':checks,'seconds':out['runtime']['seconds']},indent=2))


if __name__=='__main__':
    main()

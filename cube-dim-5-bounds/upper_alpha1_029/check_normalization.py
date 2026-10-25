#!/usr/bin/env python3
"""Necessary alpha1 qualification checks only; no Kanter or beta3/4 replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys,time

import flint
from flint import arb,ctx
import sympy as sp

HERE=Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def iroot(n,k):
    lo,hi=0,1<<((n.bit_length()+k-1)//k)
    while lo+1<hi:
        m=(lo+hi)//2
        if m**k<=n:
            lo=m
        else:
            hi=m
    return hi if hi**k==n else lo


def root(x,k,digits=150):
    x=F(x)
    scale=10**digits
    j=iroot(x.numerator*scale**k//x.denominator,k)
    lo=F(j,scale)
    hi=lo if lo**k==x else F(j+1,scale)
    assert lo**k<=x<=hi**k
    return lo,hi


def arb_ends(x):
    return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))


def parse(s):
    c,r=map(F,s[1:-1].split(' +/- '))
    return c-r,c+r


def decimal(x,up=False,digits=60):
    scale=10**digits
    j=(x*scale).__ceil__() if up else (x*scale).__floor__()
    sign='-' if j<0 else ''
    a,b=divmod(abs(j),scale)
    return f'{sign}{a}.{b:0{digits}d}'


def encode(bounds):
    return {'lower_fraction':str(bounds[0]),'upper_fraction':str(bounds[1]),
            'decimal_outward':[decimal(bounds[0]),decimal(bounds[1],True)]}


def main():
    began=time.monotonic()
    reg=json.loads((HERE/'registration.json').read_text())
    resolved={}
    for row in reg['sources']:
        current=HERE.parents[1]/row['path']
        original=Path(row['source'])
        selected=next((p for p in (current,original) if p.is_file() and sha(p)==row['sha256']),None)
        assert selected is not None,row['path']
        resolved[row['path']]=selected
    prefix='cube-dim-5-bounds/astra_upper014/'
    directory=resolved[prefix+'manifest.json'].parent
    manifest=json.loads(resolved[prefix+'manifest.json'].read_text())
    assert sha(resolved[prefix+'manifest.json'])=='7a65aa85088c1beaaaa8a6ac60974449d9175f5f0d37ef9c0c6327207dd8bea3'
    for row in manifest['white_list']:
        assert sha(directory/row['path'])==row['sha256']
    flags=json.loads(resolved[prefix+'flag_types.json'].read_text())
    geometry_candidates=[HERE.parent/'upper_I022/results.json',Path(flags['source_path'])]
    assert any(p.is_file() and sha(p)==flags['source_sha256'] for p in geometry_candidates)
    budget=json.loads(resolved[prefix+'budget.json'].read_text())
    source=next(row for row in budget['cases'] if row['alpha']=='1')

    alpha,beta,n=F(1),F(1,2),F(5)
    q=n/(n-alpha)
    r=(n-alpha)/2
    b=beta/n
    assert (q,r,b)==(F(5,4),F(2),F(1,10))
    assert (q-1)/(2*q)==b and 5*b==beta
    expected=[F(81,8),F(27),F(36),F(36),F(64)]
    type_rows=[]
    Z=F(0)
    for row,wanted in zip(flags['types'],expected):
        R=sp.Integer(1)
        for histogram in row['probability_histograms']:
            assert sum(F(p)*int(c) for p,c in histogram.items())==1
            R*=sum(sp.Rational(F(p).numerator,F(p).denominator)**sp.Rational(1,2)*int(c)
                   for p,c in histogram.items())
        assert 2+sum(sum(map(int,h.values())) for h in row['probability_histograms'])==16
        square=sp.simplify(R*R)
        assert square==sp.Rational(wanted.numerator,wanted.denominator)
        theta=F(row['theta'])
        Z+=theta*wanted
        type_rows.append({'theta':str(theta),'R_squared':str(wanted)})
    assert Z==F(2661,64)
    optimized=[F(row['theta'])*value/Z for row,value in zip(flags['types'],expected)]
    assert optimized==[F(27,887),F(108,887),F(144,887),F(96,887),F(512,887)]
    assert sum(optimized)==1

    # Symbolic alpha1 identities, independent of the general-alpha source formula.
    t,a=sp.symbols('t a',positive=True)
    time_integral=sp.integrate(t**(-5)*sp.exp(-a/(4*t**2)),(t,0,sp.oo))
    assert sp.simplify(time_integral-8/a**2)==0
    kappa=sp.simplify((4*sp.pi)**(-sp.Rational(5,2))*time_integral*a**2)
    Ztilt=sp.simplify(1/(16*sp.gamma(sp.Rational(3,2))))
    c=sp.simplify(kappa*Ztilt)
    assert sp.simplify(c-1/(32*sp.pi**3))==0
    z=sp.symbols('z',positive=True)
    generator_integral=sp.integrate((1-sp.cos(z))/z**2,(z,0,sp.oo))
    assert sp.simplify(generator_integral-sp.pi/2)==0
    h,g,I=sp.symbols('h g I',positive=True)
    mass_cost=4*I*h**(-sp.Rational(1,4))/(g-h)
    log_derivative=sp.simplify(sp.diff(sp.log(mass_cost),h))
    assert sp.simplify(log_derivative-(5*h-g)/(4*h*(g-h)))==0
    assert sp.simplify(mass_cost.subs(h,g/5)**4-3125*I**4/g**5)==0
    vp=(1,1,1,1,1)
    vm=(1,1,1,1,-1)
    point=tuple((F(x)+F(y))/2 for x,y in zip(vp,vm))
    assert point==(1,1,1,1,0)
    powers=[2*k-5 for k in range(6)]
    assert powers==[-5,-3,-1,1,3,5]

    ctx.prec=512
    pi=arb.pi()
    # Only beta1/2: unshifted Mellin ratio m=2Gamma(1/5)/Gamma(1/10).
    m=2*ar(F(1,5)).gamma()/ar(F(1,10)).gamma()
    m50=m**50
    I8=arb_ends(m50/(ar(2)**34*ar(2661)*ar(60)**8*pi**9))
    Ibound=root(I8[0],8)[0],root(I8[1],8)[1]
    J10=arb_ends(ar(2)**10*pi**5*m50/ar(Z))
    Jbound=root(J10[0],10)[0],root(J10[1],10)[1]
    cball=arb_ends(1/(32*pi**3))
    needed10=F(5)**10*I8[0]/32**8,F(5)**10*I8[1]/32**8
    needed=root(needed10[0],10)[0],root(needed10[1],10)[1]
    checks={}
    for name,own,saved in [('m',arb_ends(m),source['Gamma_ratio']),
                           ('J',Jbound,source['J_upper_expression']['ball']),
                           ('c1',cball,source['c']['ball']),
                           ('I',Ibound,source['I_upper_expression']['ball']),
                           ('g_required_exact_I_formula',needed,source['g_required_for_cost_32_using_exact_expression']['ball'])]:
        original=parse(saved)
        assert original[0]<=own[0]<=own[1]<=original[1],name
        checks[name]={'independent_interval':encode(own),'contained_in_locked_source_ball':True}
    Icap=F(36112162239,500000000000000)
    assert Ibound[1]<Icap
    safe_g=root(F(3125)*Icap**4/32**4,5)

    out={'status':'alpha1_normalization_and_source_qualification_checks_passed_g_joint_deferred',
         'script_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),
         'sources':reg['sources'],'astra014_manifest_SHA_verified':True,'astra014_whitelist_hashes_checked_only':7,
         'Kanter_or_diagnostic_execution':False,'beta3_over4_formula_execution':False,
         'exact_parameters':{'alpha':'1','beta':'1/2','q':'5/4','r':'2','negative_det_b':'1/10',
                             'R_squared':list(map(str,expected)),'Zflag':str(Z),
                             'optimized_type_weights':list(map(str,optimized)),
                             'L2_generator_domain':'H1','form_domain':'H1/2','large_cutoff_power':'-1',
                             'small_cutoff_L1_power':'4','kernel_scale_denominator_power':'1',
                             'Green_homogeneity':'-4','rank_heat_powers_k0_to5':powers},
         'symbolic_checks':{'Gaussian_time_integral':str(time_integral),'kappa1':str(kappa),
                             'ET_minus_beta':str(Ztilt),'c1':str(c),
                             'Cauchy_generator_integral_0_inf':str(generator_integral),
                             'obstacle_time_integral_factor':'5/4','cost_log_derivative':str(log_derivative),
                             'optimal_positive_part_level':'g/5','cost_fourth_power':'3125*I^4/g^5',
                             '11110_pair_representation':{'v':vp,'w':vm,'coefficients':['1/2','1/2']}},
         'independent_512bit_checks':checks,'usable_strict_rational_I_upper':str(Icap),
         'required_g_using_rational_I_cap':encode(safe_g),
         'new_g_certificate_read_or_joined':False,'new_C5_upper_bound':False,
         'current_C5_upper':32,
         'runtime':{'python':sys.executable,'flint':flint.__version__,'sympy':sp.__version__,
                    'Arb_bits':ctx.prec,'Fraction_root_decimal_digits':150,
                    'seconds':time.monotonic()-began,'spatial_or_probability_samples':0}}
    (HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'R_squared':out['exact_parameters']['R_squared'],
                      'Zflag':str(Z),'I_interval':encode(Ibound)['decimal_outward'],
                      'safe_required_g':encode(safe_g)['decimal_outward'],
                      'g_joint_deferred':True,'seconds':out['runtime']['seconds']},indent=2))


if __name__=='__main__':
    main()

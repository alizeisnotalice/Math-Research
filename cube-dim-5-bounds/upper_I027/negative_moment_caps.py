#!/usr/bin/env python3
"""Exactly three new R negative determinant moment caps; no old formula replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys,time

import flint
from flint import arb,ctx

HERE=Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def endpoints(x):
    return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))


def decimal(x,up=False,digits=50):
    scale=10**digits
    n=(x*scale).__ceil__() if up else (x*scale).__floor__()
    sign='-' if n<0 else ''
    w,r=divmod(abs(n),scale)
    return f'{sign}{w}.{r:0{digits}d}'


def encode(bounds):
    return {'lower_fraction':str(bounds[0]),'upper_fraction':str(bounds[1]),
            'decimal_outward':[decimal(bounds[0]),decimal(bounds[1],True)]}


def parse(s):
    c,r=map(F,s[1:-1].split(' +/- '))
    return c-r,c+r


def main():
    began=time.monotonic()
    reg=json.loads((HERE/'registration.json').read_text())
    for row in reg['sources']:
        assert sha(Path(row['source']))==row['sha256']
    zsource=next(Path(row['source']) for row in reg['sources'] if row['source'].endswith('upper_I024/results.json'))
    zdata=json.loads(zsource.read_text())['Z_exact_interval']
    Zlo,Zhi=F(zdata['lower_fraction']),F(zdata['upper_fraction'])
    assert 0<Zlo<Zhi
    oldsource=next(Path(row['source']) for row in reg['sources'] if row['source'].endswith('round011_shape_integral_results.json'))
    old=json.loads(oldsource.read_text())
    oldQ=parse(old['arb_256bit']['negative_determinant_moment_upper_Q'])

    beta=F(3,4)
    rows=[]
    ctx.prec=384
    gamma_common=ar(F(7,4)).gamma()
    for gamma in [F(1,5),F(3,14),F(1,4)]:
        k=5*gamma-beta
        assert 0<k<beta
        scale_unsimplified=1+k/beta-5*gamma
        scale_simplified=5*gamma/3
        assert scale_unsimplified==scale_simplified
        # Unshifted Mellin formula, rather than reproducing the old shifted-Gamma path.
        moment=ar(gamma/beta).gamma()/(ar(beta)*ar(gamma).gamma())
        prefactor=(gamma_common*ar(16)**ar(scale_simplified)
                   *ar(1-k/beta).gamma()/ar(1-k).gamma()*moment**5)
        lower_expression=prefactor*ar(Zhi)**(-ar(gamma))
        upper_expression=prefactor*ar(Zlo)**(-ar(gamma))
        bounds=(endpoints(lower_expression)[0],endpoints(upper_expression)[1])
        assert bounds[0]>0 and bounds[0]<bounds[1]
        cap=F((bounds[1]*10**12).__ceil__(),10**12)
        assert bounds[1]<cap
        row={'gamma':str(gamma),'k_5gamma_minus_beta':str(k),'k_over_beta':str(k/beta),
             'net_16_power':str(scale_simplified),'negative_S_moment_Arb':moment.str(90),
             'cap_interval':encode(bounds),'strict_rational_R_negative_moment_cap':str(cap),
             'positive_T_moment_exists':True,'all_original_negative_moments_exist':True}
        if gamma==F(3,14):
            assert bounds[1]<oldQ[0]
            ratio=(bounds[0]/oldQ[1],bounds[1]/oldQ[0])
            row['old_locked_Q3_over14_printed_interval']=encode(oldQ)
            row['new_cap_strictly_below_old_lower_endpoint']=True
            row['new_over_old_cap_ratio_interval']=encode(ratio)
        rows.append(row)

    out={'status':'three_new_R_negative_determinant_moment_caps_derived_finite_batch_passed',
         'script_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),
         'source_inputs':reg['sources'],'beta':'3/4','Zflag_source_interval':encode((Zlo,Zhi)),
         'rows':rows,'geometry_recomputed':False,'old_Q_formula_recomputed':False,
         'majorant_optimization_calls':0,'new_I_bound':False,'new_C5_upper_bound':False,
         'scope':'The caps bound E_R d^(-gamma) for the same complete shape law. They change true moment input; no claim about a scalar-majorant family under the new caps.',
         'runtime':{'python':sys.executable,'flint':flint.__version__,'Arb_bits':ctx.prec,
                    'seconds':time.monotonic()-began,'parameter_points':3,'spatial_or_probability_samples':0}}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'rows':[{'gamma':r['gamma'],
                      'decimal_cap_interval':r['cap_interval']['decimal_outward'],
                      'strict_rational_cap':r['strict_rational_R_negative_moment_cap']} for r in rows],
                      'seconds':out['runtime']['seconds']},indent=2))


if __name__=='__main__':
    main()

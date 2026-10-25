"""Arb256 deterministic parameter-tail bounds; no density sampling or grid.

All geometric distances are computed with exact rational arithmetic first.
See report.md for domination arguments preceding the certified scalar algebra.
"""
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import platform
import time
from flint import arb, ctx
import sympy as sp

ctx.prec = 256
HERE = Path(__file__).resolve().parent
START = time.monotonic()
REG = json.loads((HERE/'registration.json').read_text())
INPUTS = REG['inputs']
X = tuple(map(F, INPUTS['point']))
V = [(1,)+s for s in itertools.product((-1,1),repeat=5)]


def aq(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def min_certified(values):
    # Arb comparisons prove the selected candidate does not exceed any other.
    index=min(range(len(values)),key=lambda i: float(values[i].mid()))
    chosen=values[index]
    assert all(i==index or chosen < values[i] for i in range(len(values)))
    return index,chosen


def geometry():
    rows=[]
    xx=dot(X,X)
    for vi,v in enumerate(V):
        plane=[]; quotient=[]
        z=tuple(X[j]-v[j]*X[0] for j in range(1,6))
        for wi,w in enumerate(V):
            if vi==wi: continue
            vw=dot(v,w);vx=dot(v,X);wx=dot(w,X)
            plane.append((wi,xx-(6*vx*vx-2*vw*vx*wx+6*wx*wx)/(36-vw*vw)))
            u=tuple(w[j]-v[j] for j in range(1,6))
            quotient.append((wi,dot(z,z)-dot(z,u)**2/dot(u,u)))
        wp,d2=min(plane,key=lambda r:r[1])
        wq,e2=min(quotient,key=lambda r:r[1])
        assert d2>0 and e2>0
        rows.append({'v_index':vi,'plane_minimizer':wp,'distance_to_pair_planes_squared':str(d2),
                     'quotient_minimizer':wq,'quotient_distance_to_lines_squared':str(e2)})
    return rows


def main():
    assert REG['input_sha256']==hashlib.sha256(json.dumps(INPUTS,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert len(INPUTS['small_cutoffs'])*len(INPUTS['large_cutoffs'])==16
    rows=geometry()
    a=sp.Rational(3,2);nu_exact=sp.Rational(9,4)
    k6=a*sp.gamma(nu_exact)/(2**(a+1)*sp.pi**3)
    k5=a*sp.gamma(nu_exact)/(2**sp.Rational(3,2)*sp.pi**sp.Rational(5,2))
    assert sp.simplify(k6/k5-1/(2*sp.sqrt(sp.pi)))==0
    base_ids=[0,3,5,14,22,25]
    base_in=[];base_out=[]
    for v in V:
        flip_in=[v[j]*V[0][j] for j in range(6)]
        bi=[tuple(w[j]*flip_in[j] for j in range(6)) for w in [V[k] for k in base_ids]]
        flip_out=[v[j]*V[1][j] for j in range(6)]
        bo=[tuple(w[j]*flip_out[j] for j in range(6)) for w in [V[k] for k in base_ids]]
        assert v in bi and v not in bo
        assert abs(sp.Matrix(bi).det())==160 and abs(sp.Matrix(bo).det())==160
        quot=[tuple(w[j]-v[j] for j in range(1,6)) for w in bi if w!=v]
        assert abs(sp.Matrix(quot).det())==160
        base_in.append([V.index(w) for w in bi]);base_out.append([V.index(w) for w in bo])
    M=(arb(2)/3).gamma()/((arb(3)/2)*arb.pi())
    A=(arb(5)/2).gamma()/(arb.pi()*(arb.pi()/3).sin()**(arb(5)/2))
    full_long=M**6/160
    quotient_long=M**5/160
    sum_H=arb(0);sum_B=arb(0)
    for r in rows:
        d=aq(r['distance_to_pair_planes_squared']).sqrt()
        e=aq(r['quotient_distance_to_lines_squared']).sqrt()
        # Full stable heat kernel on the affine line x+R*v, after splitting
        # one direction into two independent half-time stable increments.
        # A det160 basis excluding v exists; active bases avoid half-weight v
        # unless it is a restricted (large) parameter. Thus no inverse half
        # weight factor is needed in the unrestricted density sup bounds.
        short_H=31*A**2*M**4*(31*arb(6).sqrt()/d)**5
        H=arb(10)/3*short_H**(arb(9)/10)*full_long**(arb(1)/10)
        # 5-dimensional quotient has31 nonzero directions; every pair is
        # independent; max direction norm<=2*sqrt5; active basis det>=32.
        short_B=arb(31*30)/32*A**2*M**3*(60*arb(5).sqrt()/e)**5
        B=arb(10)/3*short_B**(arb(9)/10)*quotient_long**(arb(1)/10)
        sum_H+=H;sum_B+=B
        r.update({'H_full_affine_line_budget':H.str(65),'B_quotient_fractional_potential_budget':B.str(65)})
    big_constant=sum_B/(2*arb.pi().sqrt())
    small=[]
    for ls in INPUTS['small_cutoffs']:
        ell=aq(ls)
        moments=[]
        for order in INPUTS['negative_moment_orders']:
            # U is positive beta-stable at time1/2; S=U+W.
            r=arb(order)
            moment=arb(2)**(4*r/3)*(1+4*r/3).gamma()/(1+r).gamma()
            moments.append(ell**order*moment)
        mi,pneg=min_certified(moments)
        chernoff=(-arb(27)/(4096*ell**3)).exp()
        ci,pchosen=min_certified([arb(1),pneg,chernoff])
        small.append({'ell':ls,'selected_negative_moment_order':INPUTS['negative_moment_orders'][mi],
                      'negative_moment_probability_bound':pneg.str(65),'chernoff_probability_bound':chernoff.str(65),
                      'selected_probability_method':['1','negative_moment','chernoff'][ci],
                      'small_tail_upper':(sum_H*pchosen).str(65),
                      'small_tail_upper_negative_moment_only':(sum_H*pneg).str(65),
                      '_small_bound':sum_H*pchosen,'_moment_small':sum_H*pneg})
    big=[]
    for rs in INPUTS['large_cutoffs']:
        R=aq(rs)
        moment_bounds=[]
        for os in INPUTS['positive_moment_orders']:
            r=aq(os)
            moment=(1-4*r/3).gamma()/(1-r).gamma()
            moment_bounds.append(moment*R**(-r-arb(1)/2))
        mi,chosen=min_certified(moment_bounds)
        big.append({'R':rs,'selected_positive_moment_order':INPUTS['positive_moment_orders'][mi],
                    'truncated_negative_half_moment_bound':chosen.str(65),
                    'large_tail_upper':(big_constant*chosen).str(65),'_big_bound':big_constant*chosen})
    grid=[]
    for sm,bg in itertools.product(small,big):
        total=sm['_small_bound']+bg['_big_bound']
        grid.append({'ell':sm['ell'],'R':bg['R'],'total_removed_parameter_contribution_upper':total.str(65),
                     'certified_below_1e_minus6':bool(total<aq(INPUTS['target'])),
                     'moment_only_removed_contribution_upper':(sm['_moment_small']+bg['_big_bound']).str(65),
                     'moment_only_certified_below_1e_minus6':bool(sm['_moment_small']+bg['_big_bound']<aq(INPUTS['target']))})
    for r in small:r.pop('_small_bound');r.pop('_moment_small')
    for r in big:r.pop('_big_bound')
    ell=aq('1/32');R=arb(10)**18;norm2=aq(dot(X,X));nu=arb(9)/4
    K6=(arb(3)/2)*nu.gamma()/(arb(2)**(arb(5)/2)*arb.pi()**3)
    coarse_lower=K6*(32*R)**(-3)*(32*ell/norm2)**nu
    coarse_upper=K6*(32*ell)**(-3)*(32*R/norm2)**nu
    assert coarse_upper-coarse_lower > 10**6
    output={'id':REG['id'],'status':'DETERMINISTIC_TRUNCATION_ONLY; INTERIOR_POINT_VALUE_NOT_EVALUATED',
            'no_new_C6_bound':True,'kernel_samples':0,'grid_cells':0,'matlab_calls':0,
            'input_sha256':REG['input_sha256'],'registration_sha256':hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),
            'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':platform.python_version(),'sympy':sp.__version__,'arb_precision_bits':ctx.prec,
            'gaussian_potential_constant_ratio_exact':'1/(2*sqrt(pi))',
            'exact_basis_certificate':{'seed_indices':base_ids,'bases_containing_each_v':base_in,
                                      'bases_excluding_each_v':base_out,'all_absolute_determinants':160},
            'geometric_distances_exact':rows,
            'constants_arb256':{'M':M.str(65),'A':A.str(65),'sum_H':sum_H.str(65),
                                'sum_B_over_2sqrtpi':big_constant.str(65)},
            'small_cutoff_budgets':small,'large_cutoff_budgets':big,'sixteen_pairs':grid,
            'successful_pair_count':sum(r['certified_below_1e_minus6'] for r in grid),
            'coarse_single_box_integrand_envelope_not_integral_certificate':{
                'ell':'1/32','R':'1000000000000000000','norm_squared_exact':str(dot(X,X)),
                'lower':coarse_lower.str(65),'upper':coarse_upper.str(65),
                'limitation':'Uniform eigenvalue enclosure alone leaves enormous width; not an impossibility proof for better interior quadrature.'},
            'remaining_gap':'No certified quadrature of the32-dimensional finite box; a tail budget does not certify a Green point interval.',
            'wall_seconds':time.monotonic()-START}
    assert output['wall_seconds']<1200
    (HERE/'truncation_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'successes':output['successful_pair_count'],'constants':output['constants_arb256'],
                      'pairs':grid,'wall_seconds':output['wall_seconds']},indent=2))


if __name__=='__main__':main()

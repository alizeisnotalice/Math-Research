#!/usr/bin/env python3
"""One new exact line certificate for saved C4 data; no derivative recomputation."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial, prod
import json
from pathlib import Path
import sys


def multiply(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def product(polys):
    out = [F(1)]
    for p in polys:
        out = multiply(out, p)
    return out


def add(p, q, scale=F(1)):
    out = p[:]+[F(0)]*max(0,len(q)-len(p))
    for i, c in enumerate(q):
        out[i] += scale*c
    return out


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    root = Path(__file__).resolve().parent
    reg = json.loads((root/'registration.json').read_text())
    require(len(reg['configurations']) == 1, 'One fixed configuration only')
    cfg = reg['configurations'][0]
    author = json.loads(Path(reg['author_results_path']).read_text())
    data = next(c for c in author['configurations'] if c['id']=='Q003-amplitude-C4')
    center = list(map(F,cfg['center']))
    direction = list(map(F,cfg['direction']))
    affine = [list(map(F,row)) for row in cfg['scalar_center_affine_coefficients']]
    line = [(sum(a*y for a,y in zip(row,center)),
             sum(a*v for a,v in zip(row,direction))) for row in affine]
    # Exact pi^3*a(y) = (1/8) sum_j product_(l!=j)(d_j-d_l)^(-2).
    # A common denominator consists of the six squared scalar-center gaps.
    # Each numerator term contains precisely the three edges not incident to j.
    edges = list(combinations(range(4),2))
    squares = {}
    margins = []
    smax = F(cfg['line_parameter_halfwidth'])
    for j,l in edges:
        d = [line[j][0]-line[l][0],line[j][1]-line[l][1]]
        margin = abs(d[0])-smax*abs(d[1])
        require(margin>0,'Line remains in distinct-center chamber')
        margins.append({'scalar_center_pair':[j,l],'absolute_gap_lower':str(margin)})
        squares[(j,l)] = multiply(d,d)
    denominator = product(squares.values())
    numerator = [F(0)]
    for j in range(4):
        complement = [squares[e] for e in edges if j not in e]
        require(len(complement)==3,'Complement has three edges')
        numerator = add(numerator,product(complement),F(1,8))
    require(denominator[0]>0,'Nonzero positive denominator at line center')

    jets = data['point_partials']
    require(len(jets)==35,'Saved C4 table has 35 entries')
    seen = set()
    polynomial = [F(0)]*5
    for row in jets:
        alpha = tuple(row['multiindex'])
        require(len(alpha)==3 and min(alpha)>=0 and sum(alpha)<=4,'Valid saved multiindex')
        require(alpha not in seen,'Unique saved multiindex')
        seen.add(alpha)
        polynomial[sum(alpha)] += (F(row['value'])*prod(v**k for v,k in zip(direction,alpha))
                                  /prod(factorial(k) for k in alpha))
    difference = add(numerator,multiply(denominator,polynomial),F(-1))
    require(all(c==0 for c in difference[:5]),'Saved directional Taylor polynomial matches rational identity through order4')
    require(any(c!=0 for c in difference[5:]),'Nontrivial remainder certificate')
    require(all(abs(v)*smax <= F(cfg['tangential_box_halfwidth']) for v in direction),
            'Whole new line segment lies in registered author box')
    output = {'schema':'B085-fixed-direction-certificate-results-v1','status':'PASS',
              'actual_configurations':1,'direction':cfg['direction'],
              'line_parameter_interval':['-'+str(smax),str(smax)],
              'scope':'One new line contraction of the saved 35 C4 jets. Exact polynomial identity; no recomputation of the 35 derivatives and no density/Green evaluation.',
              'pi_cubed_amplitude_rational_numerator':list(map(str,numerator)),
              'pi_cubed_amplitude_rational_denominator':list(map(str,denominator)),
              'saved_jets_directional_Taylor_coefficients_order0to4':list(map(str,polynomial)),
              'numerator_minus_denominator_times_Taylor_coefficients':list(map(str,difference)),
              'first_five_residual_coefficients_zero':True,'line_gap_margins':margins,
              'saved_mixed_partials_read':35,'new_directional_coefficients_checked':5,
              'body_pair_or_triple_spans_enumerated':0,'density_evaluations':0,
              'full6D_kernel_calls':0,'body32_full_kernel_calls':0,
              'dependencies':['Python standard library Fraction'],'python':sys.version}
    (root/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','configurations':1,'new_directional_coefficients_checked':5,'kernel_calls':0}))


if __name__=='__main__':
    main()

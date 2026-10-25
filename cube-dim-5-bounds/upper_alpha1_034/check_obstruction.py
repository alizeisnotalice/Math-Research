#!/usr/bin/env python3
"""Read saved full-Green balls; exact continuity-gate arithmetic, no point replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json

HERE=Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def parse(s):
    c,r=map(F,s[1:-1].split(' +/- '))
    return c-r,c+r


def root(x,k,digits=100):
    x=F(x);scale=10**digits
    n=x.numerator*scale**k//x.denominator
    lo,hi=0,1<<((n.bit_length()+k-1)//k)
    while lo+1<hi:
        m=(lo+hi)//2
        if m**k<=n:lo=m
        else:hi=m
    j=hi if hi**k==n else lo
    a=F(j,scale);b=a if a**k==x else F(j+1,scale)
    assert a**k<=x<=b**k
    return a,b


def dec(x,up=False,digits=50):
    scale=10**digits;n=(x*scale).__ceil__() if up else (x*scale).__floor__()
    a,b=divmod(n,scale)
    return f'{a}.{b:0{digits}d}'


def encode(x):
    return {'lower_fraction':str(x[0]),'upper_fraction':str(x[1]),
            'decimal_outward':[dec(x[0]),dec(x[1],True)]}


def main():
    reg=json.loads((HERE/'registration.json').read_text())
    for row in reg['sources']:assert sha(Path(row['source']))==row['sha256']
    sources={Path(r['source']).name:Path(r['source']) for r in reg['sources']}
    old=json.loads(sources['green_generic.json'].read_text())['points'][0]
    new=json.loads(sources['batch_result.json'].read_text())['selected_value']
    budget=json.loads(sources['budget.json'].read_text())
    alpha1=next(c for c in budget['cases'] if c['alpha']=='1')
    U=parse(alpha1['I_upper_expression']['ball'])
    Icap=F(alpha1['I_upper_expression']['rational_upper'])
    floor=json.loads(sources['jensen_floor.json'].read_text())
    floorI=parse(floor['I_budget_floor'])
    rows=[]
    for label,row,field in [('saved1009',old,'G'),('new104729',new,'G_ball')]:
        point=tuple(map(F,row['point']))
        assert point[0]==1 and all(0<x<1 for x in point[1:])
        assert len({abs(x) for x in point})==5
        ball=parse(row[field])
        if 'G_upper_dyadic' in row:
            lower,upper=F(row['G_lower_dyadic']),F(row['G_upper_dyadic'])
            assert ball[0]<=lower<=upper<=ball[1]
        else:
            lower,upper=ball
        # Necessary I_face for any same-G lower g<=G(point)<=upper to give cost<32.
        gate=root(F(32)**4*upper**5/F(3125),4)
        assert U[0]>gate[1]
        assert F(3125)*U[0]**4>F(32)**4*upper**5
        cost_min=root(F(3125)*U[0]**4/upper**5,4)
        assert cost_min[0]>32
        if label=='new104729':assert floorI[0]>gate[1]
        rows.append({'label':label,'point':list(map(str,point)),
                     'source_full_G_field':field,'H_epsilon_used':False,
                     'saved_full_G_interval':encode((lower,upper)),
                     'outside_S2_reason':'five distinct positive coordinate magnitudes; pair planes allow at most two',
                     'essential_lower_g_must_not_exceed_G_point':True,
                     'necessary_I_face_gate_using_G_upper':encode(gate),
                     '014_formula_budget_minimum_cost_using_G_upper':encode(cost_min),
                     '014_formula_cannot_certify_cost_below32':True})
    assert F(4,5)*10*16==128
    out={'status':'saved_full_Green_continuity_obstruction_and_face_normalization_passed',
         'script_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),
         'sources':reg['sources'],'point_rows':rows,
         'I_face_definition':'integral_[0,1]^4 G(1,y)^(5/4)dy',
         'positive_part_mass_at_level1':'128*I_face',
         'direct_unitface_QMC_weight':'1/N, no factor128 and no sorting factor',
         '033_default_budget_floor_exceeds_new_necessary_gate':True,
         'new_Green_or_old_fan_evaluations':0,'new_I_or_C5_bound':False,
         'scope':'Continuity plus read-only full-G point certificate gives an essential-g obstruction. It excludes current014 budget, not the actual I or every Cauchy mechanism.'}
    (HERE/'obstruction_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'new_gate':rows[-1]['necessary_I_face_gate_using_G_upper']['decimal_outward'],
                      '014_minimum_cost':rows[-1]['014_formula_budget_minimum_cost_using_G_upper']['decimal_outward'],
                      'point_evaluations':0,'mass_factor':128},indent=2))


if __name__=='__main__':main()

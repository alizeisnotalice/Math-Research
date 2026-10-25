"""Exact stationary-root isolation and at most two rational K3 candidates."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, gcd, lcm, floor, ceil
from functools import reduce
import hashlib
import importlib.util
import json
import sys
import time
import sympy as sp
from flint import arb, ctx

ROOT=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def ar(x):
    x=F(x);return arb(x.numerator)/x.denominator


def trim(p):
    while len(p)>1 and p[-1]==0:p.pop()
    return p


def primitive(p):
    p=trim(p);den=reduce(lcm,(q.denominator for q in p),1)
    integers=[q.numerator*(den//q.denominator) for q in p]
    common=reduce(gcd,(abs(q) for q in integers),0)
    return [F(q,common) for q in integers]


def derivative(p):return [n*p[n] for n in range(1,len(p))] or [F(0)]


def peval(p,x):
    y=F(0)
    for q in reversed(p):y=y*x+q
    return y


def remainder(a,b):
    a=list(a)
    while len(a)>=len(b) and a!=[F(0)]:
        scale=a[-1]/b[-1];offset=len(a)-len(b)
        for i,q in enumerate(b):a[offset+i]-=scale*q
        trim(a)
    return a


def sturm(p):
    seq=[primitive(p),primitive(derivative(p))]
    while True:
        r=remainder(seq[-2],seq[-1])
        if r==[F(0)]:break
        seq.append(primitive([-q for q in r]))
    return seq


def variations(seq,x):
    signs=[1 if y>0 else -1 for p in seq if (y:=peval(p,x))!=0]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def interval_eval(p,lo,hi):
    a=b=F(0)
    for q in reversed(p):
        products=[a*lo,a*hi,b*lo,b*hi];a=min(products)+q;b=max(products)+q
    return a,b


def score_interval(U,V,lo,hi):
    ul,uh=interval_eval(U,lo,hi);vl,vh=interval_eval(V,lo,hi)
    assert ul>0 and vl>0
    return ul**11/vh**7,uh**11/vl**7


def main():
    began=time.monotonic();ctx.prec=256
    helper_path=ROOT.parent/'astra_upper005'/'checkpoint.py';helper=load('locked005_poly_helper',helper_path)
    bern_path=ROOT.parent/'astra_upper004'/'verify_global.py';bern=load('locked004_bernstein_helper',bern_path)
    corner_path=ROOT.parent/'astra_upper006'/'verification.json';corner=json.loads(corner_path.read_text())
    old004_path=ROOT.parent/'astra_upper004'/'corner_checkpoint.json';old004=json.loads(old004_path.read_text())
    constants_path=ROOT.parent/'astra_upper003'/'verification.json';constants=json.loads(constants_path.read_text())['constants']
    normalized=list(map(F,corner['normalized_d_coefficients']))
    dm={k:list(map(F,corner['moments'][str(k)])) for k in range(1,5)}
    bm=list(map(F,old004['moments']['E_b_corner_W_power']))
    U=[F(comb(7,n))*sum((normalized[k]*dm[k][n] for k in range(1,5)),F(0)) for n in range(8)]
    V=[F(comb(11,n))*bm[n] for n in range(12)]
    left=helper.pmul(derivative(U),V);right=helper.pmul(U,derivative(V))
    D=primitive([11*a-7*b for a,b in zip(left,right)])
    assert len(D)-1<=17 and peval(D,F(0))!=0 and peval(D,F(1,2))!=0
    seq=sturm(D);count=variations(seq,F(0))-variations(seq,F(1,2))
    x=sp.Symbol('eta');poly=sp.Poly(sum(sp.Rational(q.numerator,q.denominator)*x**i for i,q in enumerate(D)),x)
    isolated=poly.intervals(eps=sp.Rational(1,10**14),inf=0,sup=sp.Rational(1,2))
    intervals=[]
    for (lo,hi),multiplicity in isolated:
        lo,hi=F(lo),F(hi)
        assert lo<hi and multiplicity==1
        assert variations(seq,lo)-variations(seq,hi)==1
        intervals.append((lo,hi))
    assert len(intervals)==count
    critical=[(F(0),F(0))]+intervals+[(F(1,2),F(1,2))]
    scores=[score_interval(U,V,lo,hi) for lo,hi in critical]
    winner=max(range(len(scores)),key=lambda i:scores[i][0])
    assert all(scores[winner][0]>scores[i][1] for i in range(len(scores)) if i!=winner)
    # Positivity is analytic from the polynomial minorant and f>=0.
    # Endpoint intervals and the isolated stationary intervals cover all maxima.
    cU,cL=F(constants['c_upper']),F(constants['c_lower'])
    TL,TU=F(corner['T_lower']),F(corner['T_upper'])
    upper4=cU**4*TU**11*scores[winner][1]
    display_upper=ar(upper4)**(arb(1)/4)
    tick=ceil(float(display_upper)*10**8)
    while F(tick,10**8)**4<=upper4:tick+=1
    rational_ceiling=F(tick,10**8)
    assert upper4<rational_ceiling**4
    lo,hi=critical[winner]
    if lo==hi:
        candidates_eta=[lo]
    else:
        assert floor(lo*1000)==floor(hi*1000)
        integer=floor(lo*1000);candidates_eta=[F(integer,1000),F(integer+1,1000)]
    space_path=ROOT.parent/'astra_upper007'/'space_moments.json';space=json.loads(space_path.read_text())['families']
    families={name:[{tuple(map(int,key.split(','))):F(value) for key,value in p.items()} for p in family] for name,family in space.items()}
    records=[]
    for eta in candidates_eta:
        u=TL*peval(U,eta);v=peval(V,eta)
        star=(ar(u)/ar(v))**(arb(1)/4);a=F(floor(float(star)*1000),1000)
        assert a**4<=u/v<(a+F(1,1000))**4
        corner_value=F(11,4)*u*a**7-F(7,4)*v*a**11
        row={'eta':str(eta),'a':str(a),'exact_corner_value':str(corner_value),'passes_necessary_screen':cL*corner_value>F(5563,10**8)}
        if row['passes_necessary_screen']:
            power={}
            for n in range(8):
                factor=F(11,4)*a**7*comb(7,n)*eta**n*TL
                for k in range(1,5):
                    for key,value in families['d'+str(k)][n].items():power[key]=power.get(key,F(0))+factor*normalized[k]*value
            for n in range(12):
                factor=-F(7,4)*a**11*comb(11,n)*eta**n
                for key,value in families['b'][n].items():power[key]=power.get(key,F(0))+factor*value
            power={k:v for k,v in power.items() if v}
            assert bern.corner_value(power)==corner_value
            face,degrees,tensor,index=bern.bernstein_face(power)
            assert degrees==(6,6,6,6) and len(tensor)==2401
            minimum=tensor[index];lower=cL*minimum;simple=F(lower.numerator*10**8//lower.denominator,10**8)
            row.update({'positive_global_certificate':minimum>0,'global_g_exact_rational':str(lower),'global_g_simple_rational':str(simple),
                        'float_display_only':float(lower),'minimum_Bernstein_index':index,'minimum_Bernstein_value':str(minimum),
                        'Bernstein_coefficient_count':len(tensor),'Bernstein_coefficients_sha256':hashlib.sha256(json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True).encode()).hexdigest(),
                        'symmetric_power_coefficients':{','.join(map(str,k)):str(v) for k,v in sorted(power.items())}})
        records.append(row)
    output={'status':'exact_K3_linear_family_optimization_and_registered_rational_candidates','new_C5_upper_bound':False,
            'normalized_U_power_coefficients':[str(q) for q in U],'V_power_coefficients':[str(q) for q in V],
            'primitive_derivative_coefficients':[str(q) for q in D],'derivative_degree':len(D)-1,
            'sturm_sequence_degrees':[len(p)-1 for p in seq],'root_count_in_zero_to_half':count,
            'stationary_root_intervals':[[str(lo),str(hi)] for lo,hi in intervals],
            'critical_points':[[str(lo),str(hi)] for lo,hi in critical],
            'critical_score_intervals':[[str(lo),str(hi)] for lo,hi in scores],
            'winning_critical_index':winner,'ideal_K3_all_eta_all_a_ceiling':str(rational_ceiling),
            'ideal_K3_ceiling_fourth_upper':str(upper4),'ceiling_display_Arb':str(display_upper),
            'strictly_below_1e_4':upper4<F(1,10000)**4,'not_an_upper_bound_for_G':True,
            'candidates':records,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'registration.md',Path(__file__),helper_path,bern_path,corner_path,old004_path,constants_path,space_path]},
            'runtime':{'python':sys.executable,'sympy':sp.__version__,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ['derivative_degree','root_count_in_zero_to_half','stationary_root_intervals','winning_critical_index','ideal_K3_all_eta_all_a_ceiling','ceiling_display_Arb','strictly_below_1e_4','runtime']}|{'candidates':[{k:v for k,v in row.items() if k!='symmetric_power_coefficients'} for row in records]},indent=2))


if __name__=='__main__':main()

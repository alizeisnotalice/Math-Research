"""Eight fixed Gamma peak cases and exact algebra checks. No MC or cubature."""
from pathlib import Path
import hashlib
import itertools
import json
import math
import platform
import time
import sympy as sp
from flint import arb,ctx
ctx.prec=256
HERE=Path(__file__).resolve().parent
START=time.monotonic()
REG=json.loads((HERE/'registration.json').read_text())
INPUT=REG['inputs']


def peak(dim,c):
    sphere=2*arb.pi()**(arb(dim)/2)/(arb(dim)/2).gamma()
    return sphere*(arb(2*dim)/3).gamma()/((2*arb.pi())**dim*(arb(3)/2)*c**(arb(2*dim)/3))


def exact_checks():
    z=sp.symbols('z0:6',real=True)
    e2=sum(sum(si*zi for si,zi in zip(s,z))**2 for s in itertools.product((-1,1),repeat=6))/64
    e4=sum(sum(si*zi for si,zi in zip(s,z))**4 for s in itertools.product((-1,1),repeat=6))/64
    assert sp.expand(e2-sum(t*t for t in z))==0
    assert sp.expand(e4-3*sum(t*t for t in z)**2+2*sum(t**4 for t in z))==0
    assert sp.Rational(3,2)*sp.Rational(4,5)+4*sp.Rational(1,5)==2
    q=sp.eye(5)+sp.ones(5,5)
    assert q.det()==6
    pair=[]
    for k in range(1,6):
        v=[1]*6;w=[1]+[-1]*k+[1]*(5-k)
        active=1
        rows=[]
        for j in range(1,6):
            if j==active:continue
            a=sp.Rational(v[j]+w[j],2)
            b=sp.Rational(v[j]-w[j],2)
            row=[sp.Integer(0)]*6;row[j]=1;row[0]=-a;row[active]=-b
            rows.append(row)
        T=sp.Matrix(rows); gramdet=(T*T.T).det()
        assert gramdet==k*(6-k)
        assert T*sp.Matrix(v)==sp.zeros(4,1) and T*sp.Matrix(w)==sp.zeros(4,1)
        pair.append({'hamming_k':k,'quotient_gram_determinant':int(gramdet),
                     'ordered_pair_multiplicity_per_v':math.comb(5,k),'active_2by2_absolute_determinant':2})
    # Exact moment at a fixed vector; c=32 is an actual false strengthening.
    probe=[1,1,0,0,0,0]
    moment=sum(abs(sum(si*xi for si,xi in zip(s,probe)))**sp.Rational(3,2)
               for s in itertools.product((-1,1),repeat=6))/64
    quotient_probe=[1,-1,0,0,0]
    quotient_moment=sum(abs(sum(si*xi for si,xi in zip(s,quotient_probe))-sum(quotient_probe))**sp.Rational(3,2)
                        for s in itertools.product((-1,1),repeat=5))/32
    assert sp.simplify(moment-sp.sqrt(2))==0
    assert sp.simplify(quotient_moment-sp.sqrt(2))==0
    norm_power=2**sp.Rational(3,4)
    assert sp.simplify(moment/norm_power-2**sp.Rational(-1,4))==0
    assert moment**4 < norm_power**4
    return {'second_moment_identity':str(e2),'fourth_moment_identity':'3*||z||^4-2*sum(z_i^4)',
            'holder_moment_weights':['4/5','1/5'],'free_quotient_gram_det':6,
            'pair_quotient_geometry':pair,
            'false_c32_counterexamples':{'full_vector':[1,1,0,0,0,0],
              'free_quotient_vector':[1,-1,0,0,0],
              'normalized_moment_exact':'2^(-1/4)<1',
              'proposed_lower_comparison':'2^(-1/4)>3^(-1/4)'}}


def main():
    assert len(INPUT['peak_cases'])==8
    assert REG['input_sha256']==hashlib.sha256(json.dumps(INPUT,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    exact=exact_checks()
    c=32/arb(3)**(arb(1)/4)
    cs={m:32/(3-arb(2)/m)**(arb(1)/4) for m in (4,5,6)}
    constants={};values={}
    for dim_s,label in INPUT['peak_cases']:
        dim=int(dim_s)
        used=c if label=='common' else cs[int(label.removeprefix('dimension'))]
        name=f'd{dim}_{label}'
        values[name]=peak(dim,used)
        constants[name]={'coercivity_c':used.str(65),'peak_upper':values[name].str(65)}
    assert values['d6_common'].overlaps(arb(3)/(2**24*arb.pi()**3))
    assert values['d6_dimension6'].overlaps(1/(6*2**20*arb.pi()**3))
    M=(arb(2)/3).gamma()/((arb(3)/2)*arb.pi())
    A=(arb(5)/2).gamma()/(arb.pi()*(arb.pi()/3).sin()**(arb(5)/2))
    old_one=A*M**5*32**(arb(5)/2)
    old_two=31*A**2*M**4*32**(arb(5)/2)*(31*arb(6).sqrt())**(arb(5)/2)
    common_factor=A**2*32**(arb(5)/2)*(31*arb(6).sqrt())**(arb(5)/2)
    geom_sum=sum(arb(math.comb(5,k))/arb(k*(6-k)).sqrt() for k in range(1,6))
    one_factor=32*A*32**(arb(5)/2)
    cap_input={
      'old_F':(old_one,old_two),
      'common_free':(one_factor*values['d5_common'],496*common_factor*values['d4_common']),
      'common_coarea':(one_factor*values['d5_common']/arb(6).sqrt(),16*geom_sum*common_factor*values['d4_common']),
      'dimension_free':(one_factor*values['d5_dimension5'],496*common_factor*values['d4_dimension4']),
      'dimension6_coarea':(one_factor*values['d5_dimension6']/arb(6).sqrt(),16*geom_sum*common_factor*values['d4_dimension6'])}
    oldD=arb(15)/4*old_one**(arb(1)/5)*old_two**(arb(4)/5)
    oldcap=(6*oldD)**(arb(4)/3)*arb.pi()**2/28
    cap={}
    for label,(one,two) in cap_input.items():
        D=arb(15)/4*one**(arb(1)/5)*two**(arb(4)/5)
        capcoef=(6*D)**(arb(4)/3)*arb.pi()**2/28
        radius=(arb(1)/10**6/capcoef)**(arb(3)/7)
        cap[label]={key:value.str(65) for key,value in {
          'C1':one,'C2':two,'D':D,'cap_coefficient':capcoef,
          'one_big_improvement':old_one/one,'two_big_improvement':old_two/two,
          'D_improvement':oldD/D,'cap_coefficient_improvement':oldcap/capcoef,
          'delta_for_I_budget_1e_minus6':radius,
          'I_cap_budget_at_delta_1e_minus3':capcoef*(arb(1)/1000)**(arb(7)/3)}.items()}
        cap[label]['cap_at_delta_1e_minus3_exceeds_1e_minus6']=bool(capcoef*(arb(1)/1000)**(arb(7)/3)>arb(1)/10**6)
        if label=='dimension6_coarea':
            assert radius<arb(1)/10**5
            assert capcoef>5*10**5
    product6=M**6/160;product5=M**5/160
    result={'id':REG['id'],'status':'ANALYTIC_COERCIVITY_AND_PEAK_IMPROVEMENT; NO_C6_BOUND',
            'input_sha256':REG['input_sha256'],'registration_sha256':hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),
            'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'exact_checks':exact,'eight_peak_cases':constants,
            'peak_improvements':{'six_product_over_radial':(product6/values['d6_common']).str(65),
              'five_product_over_free_radial':(product5/values['d5_common']).str(65),
              'five_product_over_free_coarea_radial':(product5/(values['d5_common']/arb(6).sqrt())).str(65)},
            'cap_comparisons':cap,'pair_coarea_sum':geom_sum.str(65),
            'runtime':{'python':platform.python_version(),'sympy':sp.__version__,'arb_bits':ctx.prec,
                       'kernel_evaluations':0,'MC_samples':0,'high_dim_quadrature_nodes':0,
                       'matlab_calls':0,'wall_seconds':time.monotonic()-START}}
    assert result['runtime']['wall_seconds']<1200
    (HERE/'coercivity_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'peaks':constants,'improvements':result['peak_improvements'],'cap':cap,
                      'wall_seconds':result['runtime']['wall_seconds']},indent=2))


if __name__=='__main__':main()

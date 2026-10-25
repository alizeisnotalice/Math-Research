#!/usr/bin/env python3
"""Eight formal scalar identities; saved-table validation, never a full replay."""
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import json
import math
import os
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from flint import fmpq as Q, fmpq_poly as P
import compress

HERE=Path(__file__).resolve().parent
CASES=[(1,2,3,4,5,6),(-3,-2,-1,1,2,3),(0,1,2,3,4,5),
       (0,0,-2,-1,1,3),(1,1,-1,-1,2,3),(-2,-2,-2,1,3,4),
       (1,1,1,1,-1,2),(0,0,0,0,0,0)]

def cmul(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]

def remaining_denominator(counts,q,order):
    out=P([1])
    for c,m in counts.items():
        multiplicity=m-(order if c==q else 0)
        assert multiplicity>=0
        out*=P([c,1])**multiplicity
    return out

def scalar():
    assert not (HERE/'scalar_checks.json').exists(),'eight scalar cases are a single registered run'
    tic=perf_counter();results=[]
    for index,nodes in enumerate(CASES):
        bs=tuple(map(Q,nodes));counts=Counter(bs)
        coefficients=compress.pole_coefficients(bs)
        real=P([0]);imag=P([0])
        for (b,j),a in coefficients.items():
            product=(P([1]),P([0]))
            for c,m in counts.items():
                for _ in range(m-(j if c==b else 0)):
                    product=cmul(product,(P([0,1]),P([-c])))
            phase=(j-6)%4
            cr,ci=((a,Q(0)),(Q(0),a),(-a,Q(0)),(Q(0),-a))[phase]
            rr,ii=cmul(product,(P([cr]),P([ci])))
            real+=rr;imag+=ii
        assert real==P([1]) and imag==P([0]),'complex rational identity failure'
        # Independent direct real numerator: Re prod(T+i*b) is z^3-e2*z^2+e4*z-e6.
        elementary=P([1])
        for b in bs:elementary*=P([1,b])
        target=P([-elementary[6],elementary[4],-elementary[2],1])
        squares=Counter(b*b for b in bs)
        actual=P([0])
        folded=compress.fold_positive(coefficients)
        for (b,j),a in folded.items():
            actual+=P(compress.numerator_coefficients(b,j))*a*remaining_denominator(squares,b*b,j)
        for j in (2,4,6):
            a=coefficients.get((Q(0),j),Q(0))
            if a:
                actual+=(1 if (j//2-3)%2==0 else -1)*a*remaining_denominator(squares,Q(0),j//2)
        assert actual==target,'real z-polynomial identity failure'
        results.append({'case':index+1,'nodes':list(nodes),'complex_cross_product_is_one':True,
                        'real_z_numerator_identity':True,'maximum_formal_degree':6})
    out={'status':'PASS_EIGHT_FORMAL_PARTIAL_IDENTITIES','pid':os.getpid(),'cases_count':8,
         'numeric_kernel_calls':0,'positive_time_values':0,'full_symbolic_passes':0,
         'elapsed_seconds':perf_counter()-tic,'cases':results,
         'compress_source_sha256':compress.sha(HERE/'compress.py'),'checker_sha256':compress.sha(HERE/'check.py')}
    compress.dump(HERE/'scalar_checks.json',out)
    print(json.dumps(out,indent=2))

def saved():
    tic=perf_counter()
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    for row in manifest['inputs']:assert compress.sha(row['path'])==row['sha256']
    record=json.loads((HERE/'run_receipt.json').read_text())
    result=json.loads((HERE/'results.json').read_text())
    data=json.loads((HERE/'partial_tables.json').read_text())
    assert record['status']=='SYMBOLIC_FULL_COMPLETE_STAGE1_STOPPED'
    assert record['tuples_processed']==40320 and record['logical_terms']==967680
    assert record['point_orbit_cone_terms']==80640 and record['symbolic_full_passes']==1
    assert record['numeric_kernel_calls']==record['positive_time_values']==0
    assert compress.sha(HERE/'partial_tables.json')==result['partial_tables_sha256']
    assert compress.sha(HERE/'scalar_checks.json')==record['scalar_checks_sha256']
    signed={(F(r['b']),r['order']):F(r['coefficient'])for r in data['signed_A']}
    assert len(signed)==len(data['signed_A']) and all(signed.values())
    folded={}
    for b in sorted({abs(b) for b,j in signed if b}):
        for j in range(1,7):
            c=signed.get((b,j),F(0))+(-1)**j*signed.get((-b,j),F(0))
            if c:folded[b,j]=c
    saved_fold={(F(r['positive_b']),r['order']):F(r['coefficient'])for r in data['positive_fold_C']}
    assert folded==saved_fold
    for row in data['zero_pole_channels']:
        j=row['order'];a=F(row['A'])
        assert a==signed.get((F(0),j),F(0))
        if j%2==0:assert a==0
    # Independently expand each tiny numerator about z=-b^2 using Fraction.
    real=defaultdict(F);powers={};highest_witnesses=[]
    for (b,j),c in folded.items():
        q=b*b;powers[q]=max(powers.get(q,0),j)
        at_pole=F(0)
        for r in range(j//2+1):
            coefficient=(1 if (j-r-3)%2==0 else -1)*math.comb(j,2*r)*b**(j-2*r)
            at_pole+=coefficient*(-q)**r
            for s in range(r+1):
                real[q,j-s]+=c*coefficient*math.comb(r,s)*(-q)**(r-s)
        assert at_pole==(-1)**(j+1)*2**(j-1)*b**j and at_pole!=0
    real={key:a for key,a in real.items()if a}
    saved_real={(F(r['b_squared']),r['order']):F(r['coefficient'])for r in data['real_z_principal']}
    assert real==saved_real
    for q,j in sorted(powers.items()):
        assert real[q,j]!=0
        assert all((qq!=q or k<=j)for qq,k in real)
        highest_witnesses.append({'b_squared':str(q),'order':j,'highest_coefficient':str(real[q,j])})
    saved_factors={F(r['b_squared']):r['exponent']for r in data['denominator_factors']}
    assert saved_factors==powers and sum(powers.values())==result['reduced_denominator_degree']
    gamma_signed=sum((a*(1 if b>0 else -1)for (b,j),a in signed.items()if b and j==1),F(0))
    gamma_fold=sum((c for (b,j),c in folded.items()if j==1),F(0))
    assert gamma_fold==gamma_signed==F(result['Gamma_G'])
    assert F(result['G_times_pi5_exact'])==gamma_fold/128
    # Independent integral basis: ∫(T²+b²)^(-j)dT=(pi/2)*binom(2j-2,j-1)/(4^(j-1)b^(2j-1)).
    gamma_real=F(0)
    for (q,j),a in real.items():
        bn,bd=math.isqrt(q.numerator),math.isqrt(q.denominator)
        assert bn*bn==q.numerator and bd*bd==q.denominator
        b=F(bn,bd)
        gamma_real+=a*math.comb(2*j-2,j-1)/(4**(j-1)*b**(2*j-1))
    assert gamma_real==gamma_fold
    out={'status':'PASS_SAVED_PARTIAL_TABLE_IDENTITIES_AND_EXACT_G','pid':os.getpid(),
         'elapsed_seconds':perf_counter()-tic,'additional_symbolic_full_passes':0,
         'additional_numeric_kernel_calls':0,'positive_time_values':0,
         'input_tuple_enumeration_repeated':False,'scalar_cases_repeated':False,
         'nonzero_signed_slots':len(signed),'positive_fold_slots':len(folded),
         'real_principal_slots':len(real),'true_denominator_degree':sum(powers.values()),
         'true_factor_count':len(powers),'zero_even_cancellation':True,
         'signed_fold_real_integral_Gamma_all_equal':True,
         'highest_pole_witnesses':highest_witnesses,
         'source_hashes':[{'path':str(HERE/f),'sha256':compress.sha(HERE/f)}for f in
                          ('input_manifest.json','scalar_checks.json','run_receipt.json','partial_tables.json','results.json','check.py')]}
    compress.dump(HERE/'checked_results.json',out)
    print(json.dumps({k:v for k,v in out.items()if k not in ('highest_pole_witnesses','source_hashes')},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=('scalar','saved'))
    (scalar if parser.parse_args().mode=='scalar'else saved)()

#!/usr/bin/env python3
"""One r=3/2 global budget, with no random or spatial computations."""
import hashlib,json,os,sys,time
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
from flint import arb,ctx

P=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_upper014')
START=time.monotonic()

def digest(path):
    p=Path(path);return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}

def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator

def main():
    reg=P/'registration.json'
    assert json.loads(reg.read_text())['status']=='registered_before_computation'
    lock=json.loads((SOURCE/'manifest.json').read_text())
    for entry in lock['white_list']:
        assert digest(SOURCE/entry['path'])['sha256']==entry['sha256'],entry['path']
    types=json.loads((SOURCE/'flag_types.json').read_text())['types']
    Z=sp.Integer(0);scales=[]
    for typ in types:
        R=sp.Integer(1)
        for block in typ['probability_histograms']:
            assert sum(sp.Rational(p)*count for p,count in block.items())==1
            R*=sum(sp.sqrt(sp.Rational(p))*count for p,count in block.items())
        R2=sp.simplify(R*R);scales.append(str(R2))
        Z+=sp.Rational(typ['theta'])*R2
    assert Z==sp.Rational(2661,64)
    assert sum(sp.Rational(t['theta']) for t in types)==1
    ctx.prec=256
    r=ar(F(3,2));beta=ar(F(1,2));gam=ar(F(1,6));k=ar(F(1,3))
    assert k<beta
    trace=arb(16)**(k/beta)*(1-k/beta).gamma()/(1-k).gamma()
    flag=arb(16)**(-5*gam)*((1+gam/beta).gamma()/(1+gam).gamma())**5*ar(F(2661,64))**(-gam)
    moment=16*(arb(3)/2).gamma()*trace*flag
    simplified_moment=64*(arb(2)/3).gamma()**3/(arb(3).sqrt()*arb.pi()*arb(2661)**(arb(1)/6))
    assert (moment-simplified_moment).contains(0)
    c1=1/(32*arb.pi()**3)
    K=c1**r*arb.pi()**2/60*arb(5).sqrt()*moment**r
    simplified_K=arb(10).sqrt()*(arb(2)/3).gamma()**(arb(9)/2)/(90*arb.pi()**4*arb(887)**(arb(1)/4))
    assert (K-simplified_K).contains(0)
    eps=ar(F(1,100000))
    v_exact_budget=eps**6/K**5
    checks={q:bool(K<ar(q)) for q in ('1/3840','1/3800','13/50000')}
    assert checks['1/3840']
    cap=F(1,3840)
    v_sufficient=F(1,100000)**6/cap**5
    assert ar(cap)**5*ar(v_sufficient)==eps**6 or (ar(cap)**5*ar(v_sufficient)-eps**6).contains(0)
    coefficient=K**(arb(5)/6)
    error_at_1e12=coefficient*ar(F(1,10**12))**(arb(1)/6)
    out={'status':'one_exponent_global_L_3over2_budget_certified',
         'parameters':{'alpha':'1','beta':'1/2','r':'3/2','gamma':'1/6','trace_k':'1/3'},
         'exact_beta_half_flag_R_squared':scales,'exact_flag_Z':'2661/64',
         'trace_positive_moment':str(trace),'detA_negative_moment_upper':str(flag),
         'R_det_negative_moment_upper':str(moment),
         'moment_simplified':'64 Gamma(2/3)^3/(sqrt(3) pi 2661^(1/6))',
         'K_integral_G_3over2_upper':str(K),
         'K_simplified':'sqrt(10) Gamma(2/3)^(9/2)/(90 pi^4 887^(1/4))',
         'rational_K_checks':checks,'strict_rational_K_upper':str(cap),
         'Holder_coefficient_K_5over6':str(coefficient),
         'target_complement_error':'1/100000','max_volume_using_exact_K_expression':str(v_exact_budget),
         'sufficient_volume_using_rational_K':str(v_sufficient),
         'Holder_error_bound_at_volume_1e_minus12':str(error_at_1e12),
         'source_lock_all_7_hashes_verified':True,
         'source_files':[digest(SOURCE/x) for x in ('manifest.json','flag_types.json')],
         'script':digest(__file__),'registration':digest(reg),
         'runtime':{'pid':os.getpid(),'python':sys.executable,'Arb_bits':ctx.prec,'seconds':time.monotonic()-START},
         'new_global_g':False,'new_C5_upper_bound':False,
         'claim_scope':'Entire positive unit face L^(3/2), including all singular-neighborhood contributions as a.e. integrals'}
    (P/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Registered one-point order-maximum conditional Jensen certificate."""
from pathlib import Path
from fractions import Fraction
from itertools import product
from collections import Counter
from flint import arb,ctx
import hashlib,json,sys,time

HERE=Path(__file__).resolve().parent

def stats(z,terms=64):
    """F, truncated mean m, truncated harmonic mean a; z>=1."""
    beta=arb(3)/4
    assert z>=1
    root_half=(arb(1)/2).sqrt()
    starts=[beta.gamma()*root_half/arb.pi()*z**(-beta),
            (2*beta).gamma()/(2*arb.pi())*z**(-2*beta),
            (3*beta).gamma()*root_half/(6*arb.pi())*z**(-3*beta)]
    tail=arb(0);moment=arb(0);inverse_tail=arb(0)
    for residue,term in enumerate(starts,1):
        tt=arb(0);mm=arb(0);nn=arb(0)
        for j in range(terms):
            k=residue+4*j
            tt+=term
            mm+=term*(beta*k)/(1-beta*k)*z
            nn+=term*(beta*k)/(1+beta*k)/z
            term=-term*(beta*k)*(beta*k+1)*(beta*k+2)*z**(-3)/((k+1)*(k+2)*(k+3)*(k+4))
        k=residue+4*terms
        assert z**(-3)/k<1
        tail+=tt+arb(0).union(term)
        moment+=mm+arb(0).union(term*(beta*k)/(1-beta*k)*z)
        inverse_tail+=nn+arb(0).union(term*(beta*k)/(1+beta*k)/z)
    F=1-tail;n0=(1/beta).gamma()/beta
    inverse_truncated=n0-inverse_tail
    m=moment/F;a=F/inverse_truncated
    assert F>0 and F<1 and m>a and m<z and a>0
    return F,m,a

def ar(q):
    q=Fraction(q)
    return arb(q.numerator)/q.denominator

def phi(p,q,m,a,z):
    return (16*m)**(-2)*(11*m+5*z)**(-arb(1)/2)*(q/(16*a)+p/(11*a+5*z))**(-arb(7)/4)

def main():
    started=time.monotonic();ctx.prec=256
    groups=Counter(Fraction(sum(v[:4])**2,5) for v in [(1,)+t for t in product((-1,1),repeat=4)])
    assert groups==Counter({Fraction(16,5):2,Fraction(4,5):8,Fraction(0):6})
    kappa=3*(arb(7)/4).gamma()/(arb(2)**(arb(7)/2)*arb.pi()**(arb(5)/2))
    beta=arb(3)/4;n0=(1/beta).gamma()/beta;a0=1/n0
    endpoints=[arb(1)]
    for _ in range(1024):endpoints.append(endpoints[-1]*65/64)
    moments=[]
    for z in endpoints:
        if time.monotonic()-started>120:raise TimeoutError('Registered120s checkpoint reached')
        moments.append(stats(z))
    outputs=[];total_lower=arb(0);total_upper_finite=arb(0);total_high=arb(0)
    F1=moments[0][0];FU,mU,aU=moments[-1];U=endpoints[-1]
    weights=[(r[0]**16-l[0]**16)/16 for l,r in zip(moments,moments[1:])]
    assert all(w>0 for w in weights)
    for p_exact,count in sorted(groups.items()):
        p=ar(p_exact);q=4-p;lower=arb(0);upper=arb(0)
        for i,weight in enumerate(weights):
            F_lo,m_lo,a_lo=moments[i];F_hi,m_hi,a_hi=moments[i+1]
            lo,hi=endpoints[i:i+2]
            low=(16*m_hi)**(-2)*(11*m_hi+5*hi)**(-arb(1)/2)*(q/(16*a_lo)+p/(11*a_lo+5*lo))**(-arb(7)/4)
            high=(16*m_lo)**(-2)*(11*m_lo+5*lo)**(-arb(1)/2)*(q/(16*a_hi)+p/(11*a_hi+5*hi))**(-arb(7)/4)
            assert high>low and low>0
            lower+=weight*low;upper+=weight*high
        high_tail=(1-FU**16)/16*(16*mU)**(-2)*(5*U)**(-arb(1)/2)*(q/(16*a0))**(-arb(7)/4)
        total_lower+=count*lower;total_upper_finite+=count*upper;total_high+=count*high_tail
        outputs.append({'p':str(p_exact),'q':str(4-p_exact),'multiplicity':count,'Green_lower_contribution':str(kappa*count*lower),'Green_finite_upper_contribution':str(kappa*count*upper),'Green_high_tail_upper':str(kappa*count*high_tail)})
    # B_a<=16zI, m>=a>=F/n0. At |x|²=4, phi<=C*n0²*z^(5/4)/F².
    C=arb(4)**(arb(7)/4)/(arb(256)*arb(5).sqrt())
    low_tail=16*C*n0*n0*F1**14/14
    Glo=kappa*total_lower;Ghi=kappa*(total_upper_finite+low_tail+total_high)
    result={'status':'ORDER_MAX_SINGLE_POINT_FORMULA_CERTIFIED_INTERVAL','point':[1,1,1,1,0],'entropy_penalty':0,'partition_is_sum_over16_disjoint_events':True,'precision_bits':ctx.prec,
      'projection_groups':outputs,'partition':{'bins':1024,'ratio':'65/64','left':'1','right':str(U),'mass_max_below1':str(F1**16),'mass_max_in_partition':str(FU**16-F1**16),'mass_max_above_U':str(1-FU**16),'all_weights_positive':True},
      'conditional_stats_at1':{'F':str(moments[0][0]),'m':str(moments[0][1]),'a':str(moments[0][2])},
      'Green_lower_Arb':str(Glo),'Green_formula_upper_Arb':str(Ghi),'Green_low_omission_upper':str(kappa*low_tail),'Green_high_omission_upper':str(kappa*total_high),
      'strictly_above_52e_minus6':Glo>arb(52)/1000000,'certifies_1_over10000':Glo>arb(1)/10000,'formula_strictly_below_1_over10000':Ghi<arb(1)/10000,
      'scope':'Only this maximum-order conditional Jensen formula at11110. Upper endpoint is not an upper bound onactualG. No spatial grid or top2 layer evaluated.',
      'seconds':time.monotonic()-started,'python':sys.executable,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()

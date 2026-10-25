"""Exact necessary-corner ceiling plus independent multivariate Stieltjes check."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import comb,factorial
import hashlib,json,time,sys
from flint import arb,ctx

P=Path(__file__).resolve().parent
EDGES=[(i,j) for i in range(5) for j in range(i+1,5)]

def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator

def coefficient_check(k):
    # Direct truncated multivariate geometric inverse of the average
    # sqrt(1-sum h_e v_i v_j), with ordinary Taylor coefficients.
    used=[i for i,n in enumerate(k) if n];cap=tuple(k[i] for i in used);zero=(0,)*len(used)
    states=list(product(*(range(n+1) for n in cap)))
    V=[(1,)+s for s in product((-1,1),repeat=4)]
    def multiply(a,b):
        c={}
        for u,cu in a.items():
            for v,cv in b.items():
                w=tuple(x+y for x,y in zip(u,v))
                if all(x<=y for x,y in zip(w,cap)):c[w]=c.get(w,F(0))+cu*cv
        return {k:v for k,v in c.items() if v}
    S={}
    for state in states:
        n=sum(state)
        if not n:continue
        char=sum(__import__('functools').reduce(lambda a,b:a*b,((v[EDGES[e][0]]*v[EDGES[e][1]])**power for e,power in zip(used,state)),1) for v in V)
        falling=F(1)
        for j in range(n):falling*=F(1,2)-j
        divisor=1
        for j in state:divisor*=factorial(j)
        value=falling*((-1)**n)*F(char,16)/divisor
        if value:S[state]=-value
    total={zero:F(1)};term={zero:F(1)}
    for _ in range(sum(cap)):
        term=multiply(term,S)
        for key,value in term.items():total[key]=total.get(key,F(0))+value
    rising=F(1)
    for j in range(sum(cap)):rising*=F(1,2)+j
    factor=1
    for n in cap:factor*=factorial(n)
    return total.get(cap,F(0))*factor/rising

def main():
    began=time.monotonic();ctx.prec=384
    raw=json.loads((P/'corner_screen.json').read_text());assert raw['beta']=='1/2'
    for test in raw['independent_check_targets']:
        assert coefficient_check(test['counts'])==F(test['moment'])
    mom={key:list(map(F,v)) for key,v in raw['moments'].items()}
    def weighted(name,n):return sum((F(comb(n,k))*2**(n-k)*mom[name][k] for k in range(n+1)),F(0))/12**n
    one=weighted('one',2);d=weighted('d',2);V=weighted('b',3);assert min(one,d,V)>0
    U2=one*d;c=1/(32*arb.pi()**3);ceiling=c*ar(U2)**(arb(3)/2)/ar(V)**2
    # Independent rational cap: pi>333/106 is enough. This comparison
    # never reconstructs a certificate from the displayed decimal.
    pi_lower=F(333,106);simple_cap=F(21,200000)
    atan5_lower=sum(((-1)**j*F(1,5**(2*j+1)*(2*j+1)) for j in range(6)),F(0))
    machin_lower=16*atan5_lower-4*F(1,239)
    assert machin_lower>pi_lower
    cap_squared_upper=U2**3/(32**2*pi_lower**6*V**4)
    assert cap_squared_upper<simple_cap**2
    old=json.loads((P.parent/'astra_upper014'/'budget.json').read_text())['cases'][0]
    target_lo=F(old['g_required_for_cost_32_using_exact_expression']['rational_lower'])
    target_hi=F(old['g_required_for_cost_32_using_exact_expression']['rational_upper'])
    assert simple_cap<target_lo
    exclusion=ceiling<ar(target_lo);nonexclusion=ceiling>ar(target_hi);assert exclusion or nonexclusion
    out={'status':'exact_fixed_family_necessary_corner_screen','new_global_g':False,'new_C5_bound':False,
         'shape':'u=(W+2)/12','E_u_squared':str(one),'E_d_u_squared':str(d),'E_B_u_cubed':str(V),
         'positive_U_squared_upper':str(U2),'optimized_corner_ceiling':str(ceiling),
         'pure_rational_corner_cap':str(simple_cap),'pi_rational_lower':str(pi_lower),
         'Machin_exact_lower':str(machin_lower),
         'corner_ceiling_squared_rational_upper':str(cap_squared_upper),
         'pure_rational_cap_below_014_threshold_lower':True,
         'excludes_required_014_threshold':bool(exclusion),'passes_only_necessary_test':bool(nonexclusion),
         'ceiling_below_1e_4':bool(ceiling<ar(F(1,10000))),
         'independent_multivariate_checks':len(raw['independent_check_targets']),
         'max_graph_degree':7,'same_model_c1':str(c),
         'scope':'An upper bound on the fixed Young certificate family at the corner; not on true G.',
         'source_sha256':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [P/'registration.md',P/'graph_moments_beta_half.hpp',P/'corner_screen.cpp',P/'corner_screen.json',Path(__file__)]},
         'runtime':{'python':sys.executable,'arb_bits':ctx.prec,'seconds':time.monotonic()-began}}
    (P/'screen_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()

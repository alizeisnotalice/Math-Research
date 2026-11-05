"""R064 independent exact algebra checks; no geometric certification."""
from fractions import Fraction as F
from math import comb
import json

N=9
coefficients={}
for h in range(1,N+1):
    for key,value in (((h,N-h),1),((h-1,N-h+1),-1)):
        coefficients[key]=coefficients.get(key,0)+value
coefficients={key:value for key,value in coefficients.items() if value}
assert coefficients=={(9,0):1,(0,9):-1}
alpha,beta=F(1,8),F(1,4)
ledger=[8+alpha,8-alpha+beta,9-alpha-3*beta,
        10-5*alpha,11-2*alpha,F(17,2)-3*alpha,9-beta]
guards=[1-alpha,1-alpha-3*beta,2-5*alpha]
assert min(ledger)==F(65,8) and min(guards)>0

def rational_envelope(n):
    size=n+9
    def const(c):return [F(c)]+[F(0)]*size
    def shift(p,k):return [F(0)]*k+p[:size+1-k]
    def add(*ps):return [sum(p[j] for p in ps) for j in range(size+1)]
    def scale(p,c):return [c*x for x in p]
    def mul(p,q):return [sum(p[k]*q[j-k] for k in range(j+1)) for j in range(size+1)]
    def reciprocal(p):
        out=const(1/p[0])
        for j in range(1,size+1):out[j]=-sum(p[k]*out[j-k] for k in range(1,j+1))/p[0]
        return out
    def binomial(k):
        out=const(1)
        for j in range(1,size+1):out[j]=out[j-1]*(k-j+1)/j
        return out
    K=F((n-1)*n**3,54)
    c8=-F(13,18)-F(n*(n-2),3)+F(400,243)
    Q=scale(add(shift(binomial(-4),6),scale(shift(binomial(-5),7),F(10,9)),
                scale(shift(binomial(-2),8),c8)),K)
    if n==2:Q=add(Q,shift(const(F(8,27)),8))
    A=binomial(n);a=shift(const(1),n)
    actual=mul(A,reciprocal(add(A,scale(a,-1),scale(mul(a,Q),-1))))
    baseline=mul(A,reciprocal(add(A,scale(a,-1))))
    Jhat=F((n-1)*n**3*(2609+1971*n+81*n*n),26244)+(F(16,27) if n==2 else 0)
    P8=add(baseline,shift(const(K),n+6),shift(const(-(n+F(26,9))*K),n+7),shift(const(Jhat),n+8))
    difference=add(actual,scale(P8,-1))
    assert all(x==0 for x in difference[:n+9])
    return {'n':n,'difference_zero_through_degree':n+8,'degree_n_plus_9':str(difference[n+9])}

output={'nine_layer_identity_nonzero_terms':[{ 'f_index':k[0],'tau_degree':k[1],'coefficient':c}
         for k,c in coefficients.items()],
        'ledger_exponents':list(map(str,ledger)), 'guard_exponents':list(map(str,guards)),
        'N9_prefix_exponents':['65/8','9','10','little-o(9)','10'],
        'N8_dangerous_prefix':'little-o(8), insufficient for O(v^(65/8))',
        'rational_envelope': [rational_envelope(n) for n in range(2,10)],
        'max_all_m_target_exponent':str(1+F(65,16)),
        'parameter_error_exponent':'6',
        'natural_rate_t_exponent':'1/(8n)',
        'geometry_certified':False,'all_assertions_passed':True}
print(json.dumps(output,ensure_ascii=False,indent=2))

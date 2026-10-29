#!/usr/bin/env python3
"""Stdlib-only exact certificate checks and independent small-case test."""
import json,math,subprocess,time,hashlib,pathlib,sys
HERE=pathlib.Path(__file__).resolve().parent

def rem(a,b):
    while a.bit_length()>=b.bit_length(): a^=b<<(a.bit_length()-b.bit_length())
    return a

def mul(a,b,p):
    z=0
    while b:
        if b&1:z^=a
        a<<=1;b>>=1
    return rem(z,p)

def power(a,n,p):
    z=1
    while n:
        if n&1:z=mul(z,a,p)
        a=mul(a,a,p);n>>=1
    return z

def gcd(a,b):
    while b:a,b=b,rem(a,b)
    return a

def factor_trial(n):
    out=[];t=2
    while t*t<=n:
        if n%t==0:
            e=0
            while n%t==0:n//=t;e+=1
            out.append([t,e])
        t=3 if t==2 else t+2
    if n>1:out.append([n,1])
    return out

def primality_trial(n):
    bound=math.isqrt(n)
    first=next((a for a in range(2,bound+1) if n%a==0),None)
    return {'integer':n,'is_prime':n>=2 and first is None,'trial_bound':bound,'first_divisor':first}

def field_certificate(d,p,g):
    assert p.bit_length()==d+1
    x=2;frob=[]
    for j in range(1,d+1):
        x=mul(x,x,p);frob.append(x)
    factors_d=factor_trial(d)
    gcds=[{'prime_divisor_of_d':r,'frobenius_index':d//r,'gcd':gcd(p,frob[d//r-1]^2)} for r,e in factors_d]
    irreducible=frob[-1]==2 and all(t['gcd']==1 for t in gcds)
    K=(1<<d)-1;factors_K=factor_trial(K)
    primes=[primality_trial(r) for r,e in factors_K]
    assert math.prod(r**e for r,e in factors_K)==K and all(t['is_prime'] for t in primes)
    powers=[{'prime_factor':r,'exponent':K//r,'residue':power(g,K//r,p)} for r,e in factors_K]
    assert irreducible and power(g,K,p)==1 and all(t['residue']!=1 for t in powers)
    return {'d':d,'modulus':p,'modulus_hex':hex(p),'representation':'bit i is coefficient of X^i; modulus includes leading X^d','degree_factorization':factors_d,'degree_primality':[primality_trial(r) for r,e in factors_d], 'rabin_frobenius_at_d':frob[-1],'rabin_gcd_tests':gcds,'irreducible':irreducible,'K':K,'K_factorization':factors_K,'K_prime_certificates':primes,'generator':g,'generator_power_K':power(g,K,p),'generator_quotient_powers':powers,'generator_is_primitive':True}

def forms(basis):
    W=[0]
    for x in basis:W+= [y^x for y in W]
    return W

def py_check(d,p,g,basis):
    W=forms(basis);V=W[1:];K=(1<<d)-1
    products=sorted(mul(x,y,p) for i,x in enumerate(V) for y in V[i:])
    ratios=sorted(mul(x,power(y,K-1,p),p) for x in V for y in V if x!=y)
    shift=sorted(mul(g,x,p) for x in ratios);gi=power(g,K-1,p)
    def h(a):
        z=14695981039346656037
        for x in a:
            for v in x.to_bytes(8,'little'):z=((z^v)*1099511628211)&((1<<64)-1)
        return f'{z:016x}'
    return {'product_count':len(products),'product_duplicate_adjacent_count':sum(x==y for x,y in zip(products,products[1:])), 'product_sorted_fnv1a64':h(products),'ratio_count':len(ratios),'ratio_duplicate_adjacent_count':sum(x==y for x,y in zip(ratios,ratios[1:])),'ratio_sorted_fnv1a64':h(ratios),'S_intersect_gS_count':len(set(ratios)&set(shift)),'short_forbidden_hits':sum(t in ratios for t in [g,gi,mul(g,g,p),mul(gi,gi,p)])}

def call(d,p,b,g,seed,attempts,basis=None):
    args=[str(HERE/'gf_check'),str(d),str(p),str(b),str(g),str(seed),str(attempts)]+[str(x) for x in (basis or [])]
    out=subprocess.run(args,text=True,capture_output=True,check=True)
    return [json.loads(line) for line in out.stdout.splitlines()]

def main():
    start=time.perf_counter()
    smallcert=field_certificate(8,0x11d,2)
    cases=[]
    # Include a passing case and a deliberate dependent basis failure.
    for basis in [[1,6],[1,2],[1,1],[17,129]]:
        c=call(8,0x11d,2,2,0,1,basis)[0]
        if c['rank']==2:
            py=py_check(8,0x11d,2,basis)
            assert all(c[k]==v for k,v in py.items()),(c,py)
        else:py={'rank_failure_independently_verified':len(set(forms(basis)))<4}
        cases.append({'basis':basis,'C':c,'independent_Python':py})
    cert=field_certificate(47,(1<<47)|(1<<5)|1,2)
    cost=call(47,cert['modulus'],8,2,123,1)[0]
    result={'small_field_certificate':smallcert,'small_double_implementation_checks':cases,'target_field_certificate':cert,'degree47_b8_cost_probe':cost,'seconds':time.perf_counter()-start,'python':sys.version}
    (HERE/(sys.argv[1] if len(sys.argv)>1 else 'PREFLIGHT.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()

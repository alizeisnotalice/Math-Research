#!/usr/bin/env python3
"""Registered exact Hadamard partition and 512-bit analytic moment check. No MC."""
from pathlib import Path
from itertools import product
from fractions import Fraction
import hashlib,json,platform,time
import sympy as sp
import flint
from flint import arb, fmpq, ctx
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE_RELS=['AGENTS.md','cube-dim-7-bounds/upper/round001/gaussian/proof.md',
 'cube-dim-7-bounds/upper/round001/theory/proof.md']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def snapshot(): return {s:sha(ROOT/s) for s in SOURCE_RELS}
def rat(a,b=1): return arb(fmpq(a,b))
def enclosure(x):
    return {'arb':str(x),'lower_rational':str(x.lower().fmpq()),'upper_rational':str(x.upper().fmpq())}
start=time.perf_counter(); before=snapshot()
V=tuple((1,)+s for s in product((-1,1),repeat=6))
# Rows0,...,6 of the8x8 character table, columns u in F2^3.
code=tuple(sum((((r&u).bit_count()%2)<< (r-1)) for r in range(1,7)) for u in range(8))
assert len(set(code))==8 and all(a^b in code for a in code for b in code)
unused=set(range(64)); cosets=[]
while unused:
    q=min(unused); coset=tuple(sorted(q^c for c in code))
    assert set(coset)<=unused
    cosets.append(coset); unused-=set(coset)
assert len(cosets)==8 and len(set(sum((list(c) for c in cosets),[])))==64
blocks=[]; determinants=[]
for coset in cosets:
    H=sp.Matrix.hstack(*(sp.Matrix((1,)+tuple((-1)**((c>>k)&1) for k in range(6))) for c in coset))
    assert H*H.T==8*sp.eye(7)
    assert all(tuple(H[:,j]) in V for j in range(8))
    minors=[int(H[:,[j for j in range(8) if j!=omit]].det()) for omit in range(8)]
    assert all(abs(d)==512 for d in minors)
    blocks.append([[int(H[i,j]) for j in range(8)] for i in range(7)])
    determinants.append(minors)
assert len(set(tuple(b[i][j] for i in range(7)) for b in blocks for j in range(8)))==64
# Cauchy-Binet prefactor and degree, plus exact equality case all S_i=1.
assert Fraction(512**2,4**7)==16
assert sum(512**2 for _ in range(8))==8**7
assert sp.Rational(1,2)/8/7==sp.Rational(1,112)
# Gamma moment: Y_i=T P_i gives T^(8a-a)=T^(7a).
a=sp.symbols('a',positive=True)
assert sp.expand(8*a-a)==7*a
# Evaluate only preregistered analytic endpoints.
ctx.prec=512
pi=arb.pi(); aa=rat(1,112)
moment=(rat(2)**(-4*aa))*(rat(4)+7*aa).gamma()/(rat(4)+8*aa).gamma()*((rat(1,2)+aa).gamma()/rat(1,2).gamma())**8
MG=(rat(8)**rat(-1,2))*moment**8
M0=rat(1,4)*(rat(65,128).gamma()/rat(1,2).gamma())**64
assert MG>0 and M0-MG>0
I_factor=rat(105)*(rat(2)**rat(25,6))*(pi**rat(13,12))
IG=MG**rat(7,6)/I_factor
I0=M0**rat(7,6)/I_factor
assert I0-IG>0
upper=Fraction(str(IG.upper().fmpq()))
den=10**18
coarse=Fraction((upper.numerator*den+upper.denominator-1)//upper.denominator,den)
assert rat(coarse.numerator,coarse.denominator)-IG>0
assert I0-rat(coarse.numerator,coarse.denominator)>0
assert snapshot()==before
out={'status':'PASS','claim':'improved complete analytic I upper budget; no new C7 bound',
 'environment':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'arb_bits':ctx.prec},
 'finite_counts':{'directions':64,'codewords':8,'cosets':8,'canonical_minors':8,'additional_minors':56,'total_minors':64},
 'code_masks':list(code),'coset_masks':[list(c) for c in cosets],
 'blocks':blocks,'signed_7_minors':determinants,
 'formulas':{'block_det':'16*(sum_8 Y_i)/(product_8 Y_i)',
 'block_negative_moment':'2^(-4a)*Gamma(4+7a)/Gamma(4+8a)*(Gamma(1/2+a)/Gamma(1/2))^8',
 'a':'1/112','M_group':'8^(-1/2)*(E det(Sigma_g)^(-1/112))^8',
 'M_Gamma':'(1/4)*(Gamma(65/128)/Gamma(1/2))^64',
 'I_from_M':'M^(7/6)/(105*2^(25/6)*pi^(13/12))'},
 'intervals':{'block_moment_a':enclosure(moment),'M_group':enclosure(MG),'M_Gamma':enclosure(M0),
 'M_baseline_minus_group':enclosure(M0-MG),'M_group_over_baseline':enclosure(MG/M0),
 'I_group':enclosure(IG),'I_Gamma':enclosure(I0),'I_baseline_minus_group':enclosure(I0-IG),
 'rational_I_upper_gap':enclosure(rat(coarse.numerator,coarse.denominator)-IG)},
 'rigorous_rational_I_upper':str(coarse),
 'rigorous_dyadic_I_upper':str(IG.upper().fmpq()),
 'sources':before,'registration_sha256':sha(HERE/'registration.md'),
 'check_sha256':sha(Path(__file__)),'elapsed_seconds':time.perf_counter()-start}
(HERE/'sources.json').write_text(json.dumps({'root':str(ROOT),'sha256':before},ensure_ascii=False,indent=2)+'\n')
(HERE/'results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':out['status'],'counts':out['finite_counts'],'M_group':str(MG),'M_Gamma':str(M0),'ratio':str(MG/M0),'I_group':str(IG),'I_Gamma':str(I0),'rational_I_upper':str(coarse),'seconds':out['elapsed_seconds']},indent=2))

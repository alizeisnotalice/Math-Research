#!/usr/bin/env python3
"""Fixed finite structural checks; no Green-node or g/I certification."""
from pathlib import Path
from itertools import product, combinations
from collections import Counter
from math import comb
import hashlib, json, platform, time
import sympy as sp
from flint import arb, ctx
import flint

HERE = Path(__file__).resolve().parent
EVIDENCE = Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/evidence-20261025/cube-finite-dim-bounds')
SOURCES = ['handoff/README_资料包说明.md', 'audit/second_review/theory/audit.md',
           'audit/second_review/theory/paper_finite_source.md', 'audit/second_review/n4/review.md',
           'audit/034_audit.md', 'audit/second_review/n3/report.md']
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def snapshots(): return {s: sha(EVIDENCE/s) for s in SOURCES}
start = time.perf_counter(); before = snapshots()
V = tuple((1,)+s for s in product((-1,1), repeat=6))
assert len(V)==64 and len(set(V))==64
M = sp.Matrix.hstack(*(sp.Matrix(v) for v in V))
assert M*M.T == 64*sp.eye(7)
# Delete first row/column of Sylvester H8; normalize column signs, then one column to all1.
H8 = sp.Matrix([[(-1)**((i&j).bit_count()) for j in range(8)] for i in range(8)])
assert H8*H8.T==8*sp.eye(8)
B = H8[1:,1:]
for j in range(7): B[:,j] *= B[0,j]
B = sp.diag(*B[:,0])*B
assert B[:,0] == sp.ones(7,1) and abs(B.det())==512
basis_checks=[]
for v in V:
    Bv=sp.diag(*v)*B
    cols=[tuple(Bv[:,j]) for j in range(7)]
    assert cols[0]==v and all(c in V for c in cols) and abs(Bv.det())==512
    basis_checks.append([V.index(c) for c in cols])
# Every triple is independent. Its V-closure canonically identifies its span.
seen_triples={}; closure_hist=Counter(); spans=0; membership_tests=0
for T in combinations(range(64),3):
    if T in seen_triples: continue
    a,b,c=(V[i] for i in T)
    d1=tuple((b[k]-a[k])//2 for k in range(7))
    d2=tuple((c[k]-a[k])//2 for k in range(7))
    pivot=None
    for i,j in combinations(range(1,7),2):
        det=d1[i]*d2[j]-d1[j]*d2[i]
        if det: pivot=(i,j,det); break
    assert pivot is not None, ('dependent triple',T)
    i,j,det=pivot; closure=[]
    for z,w in enumerate(V):
        e=tuple((w[k]-a[k])//2 for k in range(7))
        an=e[i]*d2[j]-e[j]*d2[i]
        bn=d1[i]*e[j]-d1[j]*e[i]
        membership_tests+=1
        if all(an*d1[k]+bn*d2[k]==det*e[k] for k in range(7)):
            closure.append(z)
    assert len(closure) in (3,4) and set(T)<=set(closure)
    key=tuple(closure); spans+=1; closure_hist[len(closure)]+=1
    for U in combinations(closure,3):
        assert U not in seen_triples or seen_triples[U]==key
        seen_triples[U]=key
assert len(seen_triples)==comb(64,3)==41664
assert closure_hist==Counter({3:22400,4:4816})
dependent4_formula=sum(comb(6,d)*2**(6-d)*comb(2**(d-1),2) for d in range(2,7))
assert dependent4_formula==closure_hist[4]==4816
assert spans==27216 and comb(64,4)-dependent4_formula==630560
# Exact n7 exponents, radial identity, threshold derivative and cost.
n=sp.Integer(7); q=n/(n-1); g,c,J,I=sp.symbols('g c J I',positive=True)
assert q==sp.Rational(7,6)
assert sp.simplify(((2-n)*(q-1)+1)/q)==sp.Rational(1,7)
assert sp.Rational(1,7)+1==sp.Rational(8,7)
assert -(n-1)/q==-sp.Rational(36,7)
r,Gy=sp.symbols('r Gy',positive=True)
R=(Gy/c)**(1/(n-1))
radial=sp.integrate((Gy*r**(1-n)-c)*r**(n-1),(r,0,R))
assert sp.simplify(radial-(n-1)/n*Gy**q*c**(-1/(n-1)))==0
mass=(n-1)/n*c**(-1/(n-1))*J/(2**n*(g-c))
assert sp.simplify(sp.diff(sp.log(mass),c)-(n*c-g)/((n-1)*c*(g-c)))==0
assert sp.simplify(mass.subs({c:g/7,J:7*128*I})-7**q*I/g**q)==0
Cf=sp.Integer(64)**3/sp.Integer(512)
assert Cf==512
# Algebraic factor after raising E/epsilon^(8/7) to7: (7/8)^7 * 2^9 / pi^13.
assert sp.Rational(512)*sp.Integer(64)**6/sp.Integer(2)**36==512
ctx.prec=256
pi=arb.pi(); E_coeff=arb(7)/8 * (arb(2)**(arb(9)/7)) / (pi**(arb(13)/7))
def enc(x): return {'enclosure':str(x),'lower':str(x.lower()),'upper':str(x.upper())}
checks={str(e):enc(E_coeff*(arb(e.numerator)/e.denominator)**(arb(8)/7)) for e in (sp.Rational(1,100),sp.Rational(1,1000))}
after=snapshots(); assert before==after
out={'status':'PASS','scope':'finite exact structure and analytic formula checks; no g/I/C7 certificate',
 'environment':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'arb_bits':ctx.prec},
 'counts':{'directions':64,'pairs':comb(64,2),'triples':41664,'distinct_rank3_spans':spans,'rank3_closure_histogram':dict(closure_hist),'dependent_quadruples':dependent4_formula,'independent_quadruples':630560,'membership_tests':membership_tests},
 'basis':{'matrix':[[int(B[i,j]) for j in range(7)] for i in range(7)],'abs_determinant':512,'all_64_target_completion_indices':basis_checks},
 'formulas':{'q':'7/6','k0123_density_lower_exponents':[-7,-5,-3,-1],'off_S3_density_upper_exponent':1,'short_time_norm_exponent':'8/7','orthant_norm_factor':'2^(-36/7)','face_sup':'512*pi^(-7)*t^(-5)','E_coefficient':'(7/8)*2^(9/7)*pi^(-13/7)','cost':'7^(7/6)*I/g^(7/6)'},
 'arb_expression_enclosures':{'E_coefficient':enc(E_coeff),'E_epsilon':checks},
 'sources':before,'check_sha256':sha(Path(__file__)),'elapsed_seconds':time.perf_counter()-start}
(HERE/'sources.json').write_text(json.dumps({'root':str(EVIDENCE),'sha256':before},ensure_ascii=False,indent=2)+'\n')
(HERE/'results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','counts','formulas','arb_expression_enclosures','elapsed_seconds']},ensure_ascii=False,indent=2))

#!/usr/bin/env python3
"""Independent finite audit; does not import/run any author implementation."""
from pathlib import Path
from itertools import product, combinations
from collections import Counter
from math import comb
import hashlib, json, platform, time
import sympy as s
import flint
from flint import arb, ctx
HERE=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-theory/cube-dim-7-bounds/upper/round001/theory')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def snapshot(): return {p.name:sha(p) for p in sorted(AUTHOR.iterdir()) if p.is_file()}
start=time.perf_counter(); before=snapshot()
manifest=json.loads((AUTHOR/'manifest.json').read_text())
for name,v in manifest['files'].items():
    assert before[name]==v['sha256'] and (AUTHOR/name).stat().st_size==v['bytes']
receipt=json.loads((AUTHOR/'run_receipt.json').read_text())
for name,key in [('check.py','check_sha256'),('results.json','results_sha256'),('run.stdout','stdout_sha256'),('run.stderr','stderr_sha256')]: assert before[name]==receipt[key]
provenance=json.loads((AUTHOR/'sources.json').read_text())
for name,h in provenance['sha256'].items(): assert sha(Path(provenance['root'])/name)==h
V=[(1,)+v for v in product((-1,1),repeat=6)]
M=s.Matrix(7,64,lambda i,j:V[j][i]); assert M*M.T==64*s.eye(7)
# Independent Fano-plane incidence construction: BB^T=8I-J, determinant^2=8^6.
raw=s.Matrix(7,7,lambda i,j: -1 if ((int(i+1)&int(j+1)).bit_count()%2)==0 else 1)
assert raw*raw.T==8*s.eye(7)-s.ones(7)
B=raw*s.diag(*list(raw[0,:])); B=s.diag(*list(B[:,0]))*B
assert B[0,:]==s.ones(1,7) and B[:,0]==s.ones(7,1) and abs(B.det())==512
for v in V:
    T=s.diag(*v)*B
    assert tuple(T[:,0])==v and all(tuple(T[:,j]) in V for j in range(7)) and abs(T.det())==512
# Derive exact Gram polynomial symbolically, then independently sweep ALL quadruples.
a,b,c,d,e,f=s.symbols('a b c d e f')
gram4=s.Matrix([[7,a,b,c],[a,7,d,e],[b,d,7,f],[c,e,f,7]])
formula=7**4-49*(a*a+b*b+c*c+d*d+e*e+f*f)+14*(a*b*d+a*c*e+b*c*f+d*e*f)+a*a*f*f+b*b*e*e+c*c*d*d-2*(a*b*e*f+a*c*d*f+b*c*d*e)
assert s.expand(gram4.det()-formula)==0
G=[[sum(x*y for x,y in zip(u,v)) for v in V] for u in V]
triples=0
for i,j,k in combinations(range(64),3):
    aa,bb,dd=G[i][j],G[i][k],G[j][k]
    assert 343+2*aa*bb*dd-7*(aa*aa+bb*bb+dd*dd)>0
    triples+=1
hist=Counter()
for i,j,k,l in combinations(range(64),4):
    aa,bb,cc,dd,ee,ff=G[i][j],G[i][k],G[i][l],G[j][k],G[j][l],G[k][l]
    det=2401-49*(aa*aa+bb*bb+cc*cc+dd*dd+ee*ee+ff*ff)+14*(aa*bb*dd+aa*cc*ee+bb*cc*ff+dd*ee*ff)+aa*aa*ff*ff+bb*bb*ee*ee+cc*cc*dd*dd-2*(aa*bb*ee*ff+aa*cc*dd*ff+bb*cc*dd*ee)
    assert det>=0
    hist[det]+=1
assert sum(hist.values())==635376 and hist[0]==4816 and sum(v for k,v in hist.items() if k)>0
# Distinct rank3 spans inferred using separately reviewed <=4 vertex geometric lemma.
assert triples==41664 and triples-3*hist[0]==27216 and triples-4*hist[0]==22400
# Alternative exact midpoint bucketing; pair count per midpoint determines parallelograms.
midpoints=Counter(tuple(x+y for x,y in zip(u[1:],v[1:])) for u,v in combinations(V,2))
dep_by_midpoints=sum(comb(n,2) for n in midpoints.values())
assert dep_by_midpoints==hist[0]
n=s.Integer(7); q=n/(n-1); r,A,cc,gg,II=s.symbols('r A c g I',positive=True)
radial=s.integrate((A*r**(-6)-cc)*r**6,(r,0,(A/cc)**s.Rational(1,6)))
assert s.simplify(radial-s.Rational(6,7)*A**q*cc**(-s.Rational(1,6)))==0
cost=s.Rational(6,7)*cc**(-s.Rational(1,6))*7*128*II/(128*(gg-cc))
assert s.simplify(s.diff(s.log(cost),cc)-(7*cc-gg)/(6*cc*(gg-cc)))==0
assert s.simplify(cost.subs(cc,gg/7)-7**q*II/gg**q)==0
interp=s.simplify((-5*(q-1)+1)/q); assert interp==s.Rational(1,7)
assert [2*k-7 for k in range(4)]==[-7,-5,-3,-1]
assert 2*4-7==1 and 1+interp==s.Rational(8,7)
# Raising coefficient to seventh avoids irrational simplification ambiguity.
assert s.Rational(7,8)**7*512*64**6/s.Integer(2)**36==s.Rational(7,8)**7*2**9
ctx.prec=256
E=arb(7)/8*arb(2)**(arb(9)/7)/arb.pi()**(arb(13)/7)
def enc(x):return {'enclosure':str(x),'lower':str(x.lower()),'upper':str(x.upper())}
out={'status':'PASS independent finite audit; analytical review is in review.md','environment':{'python':platform.python_version(),'sympy':s.__version__,'python_flint':flint.__version__,'arb_bits':ctx.prec},'source_hashes':before,'author_manifest_and_receipt_consistent':True,'historical_provenance_hashes_match':True,'directions':64,'gram':'64 I7','independent_Fano_basis':[[int(B[i,j]) for j in range(7)]for i in range(7)],'all_64_det512_completions':True,'triples_exact_Gram_sweep':triples,'quadruples_exact_Gram_sweep':sum(hist.values()),'quadruple_Gram_determinant_histogram':dict(sorted(hist.items())),'dependent_quadruples':hist[0],'independent_quadruples':635376-hist[0],'midpoint_check':dep_by_midpoints,'distinct_rank3_spans_from_geometry':27216,'three_vertex_rank3_spans':22400,'four_vertex_rank3_spans':4816,'formulas':{'q':'7/6','density_lower_exponents':[-7,-5,-3,-1],'off_S3_density_upper':'O(t)','short_time_norm':'epsilon^(8/7)','orthant_norm_factor':'2^(-36/7)','cost':'7^(7/6) I/g0^(7/6)'},'E_coefficient':enc(E),'E_epsilon':{str(d):enc(E*(arb(1)/d)**(arb(8)/7))for d in (100,1000)},'author_implementation_imported_or_executed':False,'g_I_C7_certification':False,'elapsed_seconds':time.perf_counter()-start}
after=snapshot(); assert before==after
out['source_hashes_before_after_equal']=True
(HERE/'results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ('status','quadruples_exact_Gram_sweep','dependent_quadruples','independent_quadruples','midpoint_check','E_coefficient','E_epsilon','elapsed_seconds')},indent=2))

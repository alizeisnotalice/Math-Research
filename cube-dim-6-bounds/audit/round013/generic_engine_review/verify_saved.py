#!/usr/bin/env python3
"""One saved-data review; no imports or calls of any producer/evaluator."""
import sys
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib, json, os, time
import sympy as sp
from flint import arb, ctx, fmpq, fmpz
ctx.prec=1536
start=time.monotonic()
HERE=Path(__file__).resolve().parent
C=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round013/generic_green_engine')
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
manifest=read(C/'manifest.json'); hashes={}
for name,row in manifest['files'].items():
    p=C/name; assert sha(p)==row['sha256'],name;hashes[str(p)]=sha(p)
for name,digest in manifest['readonly_sources'].items():
    p=Path(name); assert sha(p)==digest,name;hashes[str(p)]=digest
build=read(C/'build_record.json')
assert sha(Path(manifest['validated_binary']))==manifest['binary_sha256']==build['binary_sha256']
hashes[manifest['validated_binary']]=manifest['binary_sha256']
for name,digest in build['runtime_dylib_sha256'].items():
    assert sha(Path(name))==digest;hashes[name]=digest
assert sha(C/'generic_green.cpp')==build['source_sha256']=='516d69876325810868728ff450db7a8844dcd6f7258223e6495db819b9bda7e2'

def bounds(row):
    e=int(row['common_binary_exponent']);s=F(2)**e
    lo=F(int(row['lower_mantissa']))*s;hi=F(int(row['upper_mantissa']))*s
    assert lo<=hi
    return lo,hi
def ar(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))
def encloses_exact(row,value):
    lo,hi=bounds(row);assert lo<=value<=hi
def sq(x):return sp.Rational(x.numerator,x.denominator)
def fq(x):return F(int(sp.numer(x)),int(sp.denom(x)))

# Independent Hermite interpolation polynomial Q matches all sign jets
# needed by the original DD5 and every singly-copied DD6. Its monomial
# divided differences and derivatives are complete homogeneous polynomials.
def polynomial_reference(z):
    counts=Counter(z);conditions=[(v,r)for v,m in counts.items()for r in range(m+1)]
    degree=len(conditions)-1
    M=sp.Matrix([[0 if k<r else sp.factorial(k)/sp.factorial(k-r)*sq(v)**(k-r)
                  for k in range(degree+1)]for v,r in conditions])
    b=sp.Matrix([(1 if v>0 else -1)if r==0 else 0 for v,r in conditions])
    co=[fq(v)for v in M.inv()*b]
    h=[F(1)]+[F(0)]*max(0,degree-5)
    for v in z:
        for k in range(1,len(h)):h[k]+=v*h[k-1]
    value=sum((co[k]*h[k-5]for k in range(5,degree+1)),F(0))
    partial=[]
    for v in z:
        partial.append(sum((co[k]*sum((v**(s-1)*h[k-5-s]for s in range(1,k-4)),F(0))for k in range(6,degree+1)),F(0)))
    assert sum(partial)==0
    assert sum((a*b for a,b in zip(z,partial)),F(0))==-5*value
    if len(set(z))==6:
        weights=[]
        for i,v in enumerate(z):
            prod=F(1)
            for j,w in enumerate(z):
                if i!=j:prod*=v-w
            weights.append(F(1 if v>0 else -1)/prod)
        assert sum(weights)==value
        for i,v in enumerate(z):assert -sum(((weights[i]+weights[j])/(v-z[j])for j in range(6)if j!=i),F(0))==partial[i]
    return value,partial,degree

reg=read(C/'registration.json');an=read(C/'anchors_native.json')
assert len(reg['anchors'])==len(an['anchors'])==6
cases=[]
for i,(raw,row)in enumerate(zip(reg['anchors'],an['anchors'])):
    z=[F(n,d)for n,d in raw]
    if F(0)in z:
        assert row=={'unsupported_zero_node':True};cases.append({'index':i,'zero_rejected':True});continue
    val,grad,deg=polynomial_reference(z)
    encloses_exact(row['value'],val)
    for r,g in zip(row['partials'],grad):encloses_exact(r,g)
    cases.append({'index':i,'Hermite_polynomial_degree_bound':deg,'value_exact':str(val),'partials_exact':[str(x)for x in grad],'native_contains_exact':True})

n=read(C/'native_result.json');gref=read(C/'gradient_exact_reference.json');old=read(OLD);nv=read(C/'native_validation.json')
assert old['point_numerators']==gref['point_numerators'] and old['denominator']==gref['denominator']
assert list(map(int,nv['command'][3:]))==[old['denominator'],*old['point_numerators']]
assert gref['full_group_terms']==n['full_group_terms']==967680 and n['computed_sign_pairs']==483840
assert not gref['global_negative_pairing_used'] and n['zero_node_terms']==0 and n['confluent_terms']==0
G=arb(fmpq(fmpz(old['Q_over_pi5_numerator']),fmpz(old['Q_over_pi5_denominator'])))/arb.pi()**5
gx=arb(fmpq(fmpz(gref['gradient_Q_over_pi5_numerator']),fmpz(gref['gradient_Q_over_pi5_denominator'])))/arb.pi()**5
for row,ref in ((n['G'],G),(n['gradient'][1],gx)):
    lo,hi=bounds(row);assert ar(lo)<ref.lower() and ref.upper()<ar(hi)
x=[F(v,old['denominator'])for v in old['point_numerators']]
loG,hiG=bounds(n['G']);lo=5*loG;hi=5*hiG
for xi,row in zip(x,n['gradient']):
    a,b=bounds(row);lo+=min(xi*a,xi*b);hi+=max(xi*a,xi*b)
assert lo<=0<=hi
a,b=bounds(n['Euler_x_grad_plus5G']);assert a<=0<=b

# Contract-level integer bounds; no full-group node enumeration.
M64=2**63-1;M128=2**127-1
assert 6*M64*10**6<M128
assert M64*10**6<M128
assert M64*M64<2**126
out={'status':'PASS_FROZEN_GENERIC_NONZERO_ENGINE_SAVED_DATA_AND_FORMULAS',
 'kernel_evaluations':0,'full_reference_executions':0,'binary_launched':False,
 'saved_data_runs':1,'scalar_configurations':6,'cases':cases,
 'reference_precision_bits':1536,'native_G_strictly_contains_exact_reference_enclosure':True,
 'native_x2_gradient_strictly_contains_independent_exact_reference_enclosure':True,
 'other_five_full_gradient_components_independent_fmpq_reference':False,
 'Euler_rebuilt_using_exact_fraction_endpoints':True,'native_Euler_contains_zero':True,
 'group_terms':967680,'computed_negative_pairs':483840,
 'int128_six_term_dot_max_exact':str(6*M64*10**6),
 'int128_node_cross_product_max_exact':str(M64*M64),
 'source_hash_checks':hashes,
 'limitations':['only verified fan and canonical integer CLI inputs','zero nodes unsupported as a whole','unreadable or empty fan first token not guarded before int n use','stoll permits numeric prefixes; caller must validate syntax','specific checked FLINT3.6 ABI; no portability certification','no uniform point/box gradient bound or C6 improvement'],
 'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,'python':sys.version,'sympy':sp.__version__,
 'script_sha256':sha(Path(__file__)),'exit_code':0}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k not in ('source_hash_checks','cases')},indent=2))

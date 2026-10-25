#!/usr/bin/env python3
"""Independent six-scalar/saved-endpoint audit; never launches a kernel."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib,json,os,time
import sympy as sp
from flint import arb,ctx,fmpq,fmpz
ctx.prec=1536
start=time.monotonic()
HERE=Path(__file__).resolve().parent
C=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round014/generic_short_engine')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
m=read(C/'manifest.json');hashes={}
for name,row in m['files'].items():
 p=C/name;assert sha(p)==row['sha256'],name;hashes[str(p)]=sha(p)
for name,digest in m['readonly_sources'].items():
 assert sha(Path(name))==digest,name;hashes[name]=digest
build=read(C/'build_record.json');binary=Path(m['validated_binary'])
assert sha(binary)==build['binary_sha256']==m['binary_sha256'];hashes[str(binary)]=sha(binary)
assert sha(C/'generic_short.cpp')==build['source_sha256']=='247a1f6c6967ad52f34320b999818214488700651b46cb8e5feb3eac7064953b'
for name,digest in build['runtime_dylib_sha256'].items():
 assert sha(Path(name))==digest;hashes[name]=digest
hashes[str(C/'manifest.json')]=sha(C/'manifest.json')
def bounds(row):
 s=F(2)**int(row['common_binary_exponent']);a=F(int(row['lower_mantissa']))*s;b=F(int(row['upper_mantissa']))*s
 assert a<=b;return a,b
def A(q):
 q=F(q);return arb(fmpq(q.numerator,q.denominator))
def S(q):return sp.Rational(q.numerator,q.denominator)
def Q(q):return F(int(sp.numer(q)),int(sp.denom(q)))
def contained(row,ref):
 a,b=bounds(row);assert A(a)<ref.lower() and ref.upper()<A(b)
def contains_rational(row,q):
 a,b=bounds(row);assert a<=q<=b
def complexmul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])

reg=read(C/'registration.json');nativeanchors=read(C/'anchors_native.json');savedrefs=read(C/'anchor_validation.json')
t=F(reg['T']);assert t==F(4,125)==F(nativeanchors['T'])==F(savedrefs['T'])
zsym=sp.Symbol('z');ts=S(t)
phi_positive_jets=[None]+[sp.diff(-ts/(ts**2+zsym**2),zsym,k-1)/sp.factorial(k) for k in range(1,7)]
densityjets=[sp.diff(zsym/(ts**2+zsym**2),zsym,k)/sp.factorial(k) for k in range(6)]
def inverse_hermite(conditions):
 d=len(conditions)-1
 M=sp.Matrix([[0 if k<r else sp.binomial(k,r)*S(v)**(k-r) for k in range(d+1)]for v,r in conditions])
 return [[Q(a)for a in row]for row in M.inv().tolist()]
def scalar_refs(z):
 counts=Counter(z);conditions=[(v,k)for v,n in counts.items()for k in range(n+1)]
 inv=inverse_hermite(conditions);rhs=[]
 for v,k in conditions:
  rhs.append((arb.pi()/2 if v>0 else-arb.pi()/2)-(A(v)/A(t)).atan() if k==0 else A(Q(phi_positive_jets[k].subs(zsym,S(v)))))
 co=[sum((A(a)*b for a,b in zip(row,rhs)),arb(0))for row in inv]
 d=len(co)-1;h=[F(1)]+[F(0)]*(d-5)
 for v in z:
  for k in range(1,len(h)):h[k]+=v*h[k-1]
 val=sum((co[k]*A(h[k-5])for k in range(5,d+1)),arb(0));grad=[]
 for v in z:
  grad.append(sum((co[k]*A(sum((v**(s-1)*h[k-5-s]for s in range(1,k-4)),F(0)))for k in range(6,d+1)),arb(0)))
 # Rational density is independently a degree<=5 Hermite leading coefficient.
 dc=[(v,k)for v,n in counts.items()for k in range(n)];di=inverse_hermite(dc)
 dr=[Q(densityjets[k].subs(zsym,S(v)))for v,k in dc]
 density=sum((a*b for a,b in zip(di[5],dr)),F(0))
 product=(F(1),F(0))
 for v in z:product=complexmul(product,(t/(t*t+v*v),v/(t*t+v*v)))
 assert density==product[0]
 e=5*val+A(t*density)+sum((A(v)*g for v,g in zip(z,grad)),arb(0));assert e.contains(0)
 return val,grad,density

assert len(reg['anchors'])==len(nativeanchors['anchors'])==len(savedrefs['cases'])==6
cases=[]
for i,(raw,nr,sr)in enumerate(zip(reg['anchors'],nativeanchors['anchors'],savedrefs['cases'])):
 z=[F(a,b)for a,b in raw]
 if F(0)in z:
  assert nr=={'unsupported_zero_node':True} and sr['zero_rejected'];cases.append({'index':i,'zero_rejected':True});continue
 val,grad,p=scalar_refs(z);contained(nr['value'],val)
 for n,g in zip(nr['partials'],grad):contained(n,g)
 contains_rational(nr['p_cone'],p)
 for n,r in [(nr['value'],sr['phi_DD5_reference_dyadic']),*zip(nr['partials'],sr['partial_DD6_reference_dyadic'])]:
  a,b=bounds(n);c,d=bounds(r);assert a<=c<=d<=b
 pr=sr['p_DD5_equals_complexprod_exact'];assert p==F(int(pr['numerator']),int(pr['denominator']))
 # Exact Fraction interval Euler check for the saved native scalar outputs.
 lo,hi=bounds(nr['value']);lo*=5;hi*=5
 for factor,row in [(t,nr['p_cone']),*zip(z,nr['partials'])]:
  a,b=bounds(row);lo+=min(factor*a,factor*b);hi+=max(factor*a,factor*b)
 assert lo<=0<=hi
 cases.append({'index':i,'native_contains_independent_1536bit_phi_and_six_partial_references':True,'rational_density_exact':str(p),'density_interpolation_equals_complex_product':True,'native_and_reference_joint_Euler':True})

v=read(C/'validation_results.json');n=read(C/'native_result.json');paths=v['input_sha256']
for name,digest in paths.items():assert sha(Path(name))==digest
def source(suffix):return read(next(Path(p)for p in paths if p.endswith(suffix)))
g=source('exact_local_reference_results.json');hh=source('g_time_moments/results.json');dg=source('gradient_exact_reference.json');dh=source('short_minimum_diagnostic/gradient_results.json')
assert g['point_numerators']==hh['center_numerators']==dg['point_numerators']==dh['center_numerators']
assert g['denominator']==hh['denominator']==dg['denominator']==dh['denominator']
assert list(map(int,v['command'][3:]))==[4,125,g['denominator'],*g['point_numerators']]
hr=next(r for r in hh['time_nodes']if F(r['tau'])==t);dhr=next(r for r in dh['gradient_passes']if F(r['tau'])==t)
G=arb(fmpq(fmpz(g['Q_over_pi5_numerator']),fmpz(g['Q_over_pi5_denominator'])))/arb.pi()**5
DG=arb(fmpq(fmpz(dg['gradient_Q_over_pi5_numerator']),fmpz(dg['gradient_Q_over_pi5_denominator'])))/arb.pi()**5
Kref=G-arb(hr['H_arb384']);Dref=DG-arb(dhr['gradient_x2_to_x6'][0])
for row,ref in [(n['K'],Kref),(n['gradient'][1],Dref)]:
 a,b=bounds(row);assert ref.lower()<A(a) and A(b)<ref.upper()
assert hr['complete_group_terms']==dhr['complete_group_terms']==n['full_group_terms']==967680
assert n['computed_sign_pairs']==483840 and n['zero_node_terms']==0
x=[F(a,g['denominator'])for a in g['point_numerators']]
lo,hi=bounds(n['K']);lo*=5;hi*=5
for factor,row in [(t,n['p_T']),*zip(x,n['gradient'])]:
 a,b=bounds(row);lo+=min(factor*a,factor*b);hi+=max(factor*a,factor*b)
assert lo<=0<=hi
a,b=bounds(n['Euler_x_grad_plus5K_plusTp']);assert a<=0<=b
assert bounds(n['p_T'])[0]>0
out={'status':'PASS_TRUSTED_INPUT_GENERIC_SHORT_ENGINE_REVIEW','kernel_evaluations':0,'binary_launched':False,'saved_data_runs':1,'scalar_configurations':6,'reference_precision_bits':1536,'cases':cases,'native_K_inside_old_G_minus_H_certificate':True,'native_x2_gradient_inside_old_DG_minus_DH_certificate':True,'full_point_old_comparisons_are_consistency_checks_not_exact_value_enclosures_of_native':True,'full_point_p_independent_full_reference':False,'full_point_Euler_exact_Fraction_intervals':True,'p_full_point_strictly_positive':True,'source_hash_checks':hashes,'limitations':['verified fan hash and canonical integer CLI strings required','empty/unreadable fan parser and numeric-prefix parser limits inherited from013','all exact zero nodes reject whole generic evaluation','fixed checked FLINT3.6 ABI only','no fullcube or continuous box certificate'],'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,'python':sys.version,'sympy':sp.__version__,'script_sha256':sha(Path(__file__)),'exit_code':0}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:val for k,val in out.items()if k not in('cases','source_hash_checks')},indent=2))

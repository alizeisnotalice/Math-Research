#!/usr/bin/env python3
"""Portable read-only four-moment audit. Prints JSON; no panel regeneration."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,time,platform
import sympy as sp
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent;CUBE=HERE.parents[2]
p=argparse.ArgumentParser()
for key,rel in [('numeric','upper/round006/radau4_mellin'),('theory','upper/round006/theory'),('r004','upper/round004/pair_mellin'),('r005','upper/round005/radau_mellin')]:p.add_argument('--'+key+'-dir',type=Path,default=CUBE/rel)
a=p.parse_args();NUM=a.numeric_dir.resolve();SOL=a.theory_dir.resolve();OLD4=a.r004_dir.resolve();OLD5=a.r005_dir.resolve()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[x for x in sorted(NUM.iterdir())if x.is_file()]+[x for x in sorted(SOL.iterdir())if x.is_file()]+[x/n for x in (OLD4,OLD5)for n in ('certificate.json','manifest.json')]
def snap():return {str(p):sha(p)for p in files}
start=time.perf_counter();before=snap()
for n,h in json.loads((NUM/'manifest.json').read_text())['files'].items():assert sha(NUM/n)==h,n
for line in (SOL/'manifest.sha256').read_text().splitlines():
    h,n=line.split(maxsplit=1);assert sha(SOL/Path(n).name)==h
old=[]
for k,src in ((4,OLD4),(5,OLD5)):
    copy=NUM/f'r00{k}-certificate.json';raw=copy.read_bytes();manifest=json.loads((NUM/f'r00{k}-manifest.json').read_text())
    assert sha(copy)==manifest['files']['certificate.json']
    assert raw==(src/'certificate.json').read_bytes() and sha(src/'certificate.json')==json.loads((src/'manifest.json').read_text())['files']['certificate.json']
    old.append(json.loads(raw))
O4,O5=old;C=json.loads((NUM/'certificate.json').read_text())
assert C['reused_q_1_8']==O4['q'] and C['reused_q_R005']==O5['new_q']
def dy(d):return Q(int(d[0]))*Q(2)**int(d[1])
def box(d):return (dy(d['lower']),dy(d['upper']))
def contains(d,iv):assert dy(d['lower'])<=iv[0]<=iv[1]<=dy(d['upper'])
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    vals=[x*y for x in a for y in b];return(min(vals),max(vals))
def div(a,b):
    assert b[0]*b[1]>0
    return mul(a,(1/b[1],1/b[0]))
def const(x):return(Q(x),Q(x))
def scale(a,n):return mul(a,const(n))
def powi(a,n):
    if n==0:return const(1)
    if n%2==0 and a[0]<0<a[1]:return(Q(0),max(a[0]**n,a[1]**n))
    v=[a[0]**n,a[1]**n];return(min(v),max(v))
ctx.prec=768
def A(q):q=q if isinstance(q,Q) else Q(str(q));return arb(fmpq(q.numerator,q.denominator))
def ep(x):return(Q(str(x.lower().fmpq())),Q(str(x.upper().fmpq())))
def stringify(iv):return[str(x)for x in iv]
def sqrtiv(iv):assert iv[0]>0;return(ep(A(iv[0]).sqrt())[0],ep(A(iv[1]).sqrt())[1])
def compact(iv):return stringify((ep(A(iv[0]))[0],ep(A(iv[1]))[1]))
def m(v):return(-A(Q(4,7))*v*A(2).log()+(4+v).lgamma()-(4+A(Q(8,7))*v).lgamma()+8*((A(Q(1,2))+v/7).lgamma()-A(Q(1,2)).lgamma())).exp()
def corr(d,n):
    t=1+d;v=d
    for _ in range(n):v=2*v*(t+1);t=2*t*t-1
    root=t.sqrt().sqrt();return v/(root*(root**3+root**2+root+1))
# Exact moment/Chebyshev/old-variable identities, no floating roots.
x=sp.symbols('x');R=sp.symbols('R1:5');u1=R[0]-1;u2=R[1]/2-2*R[0]+sp.Rational(3,2);u3=R[2]/4-3*R[1]/2+15*R[0]/4-sp.Rational(5,2);u4=R[3]/8-R[2]+7*R[1]/2-7*R[0]+sp.Rational(35,8)
assert sp.simplify(4*u1+2*u2-(R[1]-1))==0
assert sp.simplify(16*u2+16*u3+4*u4-(R[3]/2-2*R[1]+sp.Rational(3,2)))==0
assert sp.expand(sp.chebyshevt(8,sp.chebyshevt(2,x))-sp.chebyshevt(16,x))==0
assert sp.Poly(sp.chebyshevt(16,x),x).LC()==2**15
P=[json.loads(s)for s in (NUM/'panels.jsonl').read_text().splitlines()];counts=[2048,1024,512,512]
assert len(P)==4096 and {(r['group'],r['panel'])for r in P}=={(g,j)for g,n in enumerate(counts)for j in range(n)} and all(r['status']=='completed'for r in P)
new=[]
for g,(sn,n)in enumerate(zip((1,3,5,7),counts)):
    L=Q(0);U=Q(0)
    for row in P:
        if row['group']==g:
            l,u=box(row['integral']);assert 0<l<=u;L+=l;U+=u
    srat=Q(sn,32);ar=srat/2;rh=3*ar/4;hr=Q(16,n);eta=hr/rh
    assert 2*hr*n==32 and 0<hr<rh<ar
    assert eta==[Q(2,3),Q(4,9),Q(8,15),Q(8,21)][g]
    d=C['new_q'][g];assert Q(d['s'])==srat
    contains(d['polynomial_lower'],(L,L));contains(d['polynomial_upper'],(U,U));contains(d['eta'],const(eta))
    s=A(srat);av=A(ar);rho=A(rh);norm=arb.pi()*s.gamma()
    star=((av-rho).lgamma()+(av+rho).lgamma()).exp()*m(av-rho)*m(av+rho)
    E=32*star*A(eta)**129/(1-A(eta));tail=m(av)**2*(2*(av+16).lgamma()).exp()/(norm*31*A(32)**31)
    qlo=((A(L)-E)/norm).lower();qhi=((A(U)+E)/norm+tail).upper();assert qlo>0 and qhi>qlo
    contains(d['q'],(Q(str(qlo.fmpq())),Q(str(qhi.fmpq()))))
    for k,v in [('Mstar',star),('quadrature_abs_error',E),('q_tail_bound',tail)]:contains(d[k],ep(v))
    new.append({'s':str(srat),'panels':n,'eta':str(eta),'exact_polynomial_sum':stringify((L,U)),'raw_error':stringify(ep(E)),'tail':stringify(ep(tail)),'q':stringify((Q(str(qlo.fmpq())),Q(str(qhi.fmpq()))))})
# Fully rational interval arithmetic, avoiding midpoint/rounded-ball substitution in roots.
AA=box(O4['q']);newq=[box(d['q'])for d in C['new_q']];prev={Q(d['s']):box(d['q'])for d in O5['new_q']}
RR=[div(mul(newq[1],newq[2]),powi(AA,2)),div(mul(prev[Q(1,16)],prev[Q(3,16)]),powi(AA,2)),div(mul(newq[0],newq[3]),powi(AA,2)),div(prev[Q(1,4)],powi(AA,2))]
for d,iv in zip(C['R'],RR):contains(d,iv)
r1,r2,r3,r4=RR
mm=[sub(r1,const(1)),add(sub(scale(r2,Q(1,2)),scale(r1,2)),const(Q(3,2))),sub(add(sub(scale(r3,Q(1,4)),scale(r2,Q(3,2))),scale(r1,Q(15,4))),const(Q(5,2))),add(sub(add(sub(scale(r4,Q(1,8)),r3),scale(r2,Q(7,2))),scale(r1,7)),const(Q(35,8)))]
for d,iv in zip(C['moments'],mm):contains(d,iv);assert iv[0]>0
m1,m2,m3,m4=mm
DD=sub(mul(m1,m3),powi(m2,2));assert DD[0]>0;contains(C['D'],DD)
alpha=div(sub(mul(m2,m3),mul(m1,m4)),DD);beta=div(sub(mul(m2,m4),powi(m3,2)),DD);disc=sub(powi(alpha,2),scale(beta,4))
assert alpha[1]<0 and beta[0]>0 and disc[0]>0
for key,iv in [('alpha',alpha),('beta',beta),('discriminant',disc)]:contains(C[key],iv)
root2=scale(add(neg(alpha),sqrtiv(disc)),Q(1,2));root1=div(beta,root2)
assert 0<root1[0] and root1[1]<root2[0]
for d,iv in zip(C['roots'],(root1,root2)):contains(d,iv)
w1=div(sub(mul(m1,root2),m2),mul(root1,sub(root2,root1)))
w2=div(sub(m2,mul(m1,root1)),mul(root2,sub(root2,root1)))
w0=sub(sub(const(1),w1),w2);weights=[w0,w1,w2]
for d,iv in zip(C['weights'],weights):contains(d,iv);assert iv[0]>0
mass=add(add(w0,w1),w2);assert mass[0]<=1<=mass[1]
for k,mom in enumerate([const(1)]+mm):
    val=add(mul(w1,powi(root1,k)),mul(w2,powi(root2,k)))
    if k==0:val=add(w0,val)
    res=sub(val,mom);assert res[0]<=0<=res[1]
    assert dy(C['moment_residuals'][k]['lower'])<=0<=dy(C['moment_residuals'][k]['upper'])
# Lower weights and lower roots are mathematically safe only after the preceding enclosures.
wl1=dy(C['weights'][1]['lower']);wl2=dy(C['weights'][2]['lower']);rl1=dy(C['roots'][0]['lower']);rl2=dy(C['roots'][1]['lower'])
Delta=A(wl1)*corr(A(rl1),4)+A(wl2)*corr(A(rl2),4)
dl=dy(C['delta_lower']['lower']);assert dl==dy(C['delta_lower']['upper']) and 0<dl<ep(Delta)[0]<1
F=1/(2*arb.pi()**A(Q(7,2)))**A(Q(7,6))*(2*arb.pi()**A(Q(7,2))/A(Q(7,2)).gamma())/896
B=1-A(dl);M=A(AA[1])**4*B**2/2;I=M**A(Q(7,6))*F
for k,v in [('B_upper',B),('M_upper',M),('I_upper',I)]:contains(C[k],ep(v))
# Independently valid lower bound of the *true-moment* R005 deterministic expression.
mu5=sub(r2,const(1));nu5=add(sub(scale(r4,Q(1,2)),scale(r2,2)),const(Q(3,2)))
assert mu5[0]>0 and nu5[0]>mu5[1]**2
Delta5=A(mu5[1])**2/A(nu5[0])*corr(A(nu5[0]/mu5[1]),3);d5upper=ep(Delta5)[1];assert 0<d5upper<1
M5safe=A(AA[0])**4*(1-A(d5upper))**2/2;I5safe=M5safe**A(Q(7,6))*F
old_declared_lower=dy(C['R005_exact_budget_expression_lower']['lower'])
assert old_declared_lower<=ep(I5safe)[0]
RI=Q(C['rational_I_upper_numerator'],C['rational_I_upper_denominator']);assert RI==Q(7723627498344,10**18)
assert ep(I)[1]<dy(C['I_upper']['upper'])<RI<old_declared_lower<=ep(I5safe)[0]
assert snap()==before
out={'status':'PASS independent four-moment Radau proof and scalar interface; no repeated panels','environment':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'bits':768},'new_panels_recomputed':0,'ledger_unique_counts':counts,'old_certificates_verbatim_reused':True,'four_new_q_reassemblies':new,'exact_rational_moment_intervals':[stringify(z)for z in mm],'positive_D_and_root_weight_gates':True,'independent_roots_outward_dyadic':[compact(root1),compact(root2)],'independent_weights_outward_dyadic':[compact(w)for w in weights],'declared_correction_lower':str(dl),'independent_correction_at_declared_lower_inputs':stringify(ep(Delta)),'I_upper_expression':stringify(ep(I)),'independent_R005_true_budget_safe_lower_expression':stringify(ep(I5safe)),'declared_R005_budget_lower':str(old_declared_lower),'rational_complete_I_upper':str(RI),'strict_comparison':'I768upper < declared dyadic upper <7723627498344/10^18 < declared R005 lower <= independent valid R005 lower','source_hashes':before,'source_before_after_equal':True,'author_checker_executed':False,'new_MC':0,'global_g_proved':False,'new_C7_bound':False,'elapsed_seconds':time.perf_counter()-start}
print(json.dumps(out,indent=2))

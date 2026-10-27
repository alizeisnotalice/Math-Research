#!/usr/bin/env python3
"""Portable stdout-only Radau proof-interface audit; never recomputes panels."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,time,platform
import sympy as sp
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent;CUBE=HERE.parents[2]
p=argparse.ArgumentParser();p.add_argument('--numeric-dir',type=Path,default=CUBE/'upper/round005/radau_mellin');p.add_argument('--theory-dir',type=Path,default=CUBE/'upper/round005/theory');p.add_argument('--r004-dir',type=Path,default=CUBE/'upper/round004/pair_mellin');args=p.parse_args()
NUM=args.numeric_dir.resolve();SOL=args.theory_dir.resolve();OLD=args.r004_dir.resolve()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[x for x in sorted(NUM.iterdir())if x.is_file()]+[x for x in sorted(SOL.iterdir())if x.is_file()]+[OLD/'certificate.json',OLD/'manifest.json']
def snap():return {str(x):sha(x)for x in files}
start=time.perf_counter();before=snap()
for name,h in json.loads((NUM/'manifest.json').read_text())['files'].items():assert sha(NUM/name)==h,name
for line in (SOL/'manifest.sha256').read_text().splitlines():
    h,name=line.split(maxsplit=1);assert sha(SOL/Path(name).name)==h
oldraw=(NUM/'r004-certificate.json').read_bytes();old=json.loads(oldraw)
assert sha(NUM/'r004-certificate.json')==json.loads((NUM/'r004-manifest.json').read_text())['files']['certificate.json']
assert oldraw==(OLD/'certificate.json').read_bytes()
assert sha(OLD/'certificate.json')==json.loads((OLD/'manifest.json').read_text())['files']['certificate.json']
cert=json.loads((NUM/'certificate.json').read_text())
assert cert['reused_q_1_8']==old['q']
# Exact Chebyshev recurrence and global proof algebra helpers; no roots approximated.
x=sp.symbols('x');t2=2*x*x-1;t4=2*t2*t2-1;t8=2*t4*t4-1
assert sp.expand(t8-sp.chebyshevt(8,x))==0 and t8.subs(x,1)==1 and sp.Poly(t8,x).LC()==128
CFG=[(Q(1,16),Q(1,64),Q(3,128),1024),(Q(3,16),Q(1,32),Q(3,64),512),(Q(1,4),Q(1,16),Q(3,32),256)]
P=[json.loads(s)for s in (NUM/'panels.jsonl').read_text().splitlines()]
assert len(P)==1792 and {(j['group'],j['panel'])for j in P}=={(g,j)for g,c in enumerate(CFG)for j in range(c[3])} and all(j['status']=='completed'for j in P)
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
def A(q):q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def ep(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def m(v):return (-A(Q(4,7))*v*A(2).log()+(4+v).lgamma()-(4+A(Q(8,7))*v).lgamma()+8*((A(Q(1,2))+v/7).lgamma()-A(Q(1,2)).lgamma())).exp()
def cert_contains(d,lo,hi):assert dy(d['lower'])<=Q(str(lo))<=Q(str(hi))<=dy(d['upper'])
ctx.prec=768
new=[]
for g,(ss,hh,rr,n)in enumerate(CFG):
    assert 2*hh*n==32 and hh/rr==Q(2,3) and 0<hh<rr<ss/2
    L=Q(0);U=Q(0)
    for row in P:
        if row['group']==g:
            l,u=dy(row['integral']['lower']),dy(row['integral']['upper']);assert 0<l<=u;L+=l;U+=u
    d=cert['new_q'][g];assert Q(d['s'])==ss
    cert_contains(d['polynomial_lower'],L,L);cert_contains(d['polynomial_upper'],U,U)
    s=A(ss);a=s/2;rho=A(rr);norm=arb.pi()*s.gamma()
    star=((a-rho).lgamma()+(a+rho).lgamma()).exp()*m(a-rho)*m(a+rho)
    E=32*star*A(Q(2,3))**129*3
    tail=m(a)**2*(2*(a+16).lgamma()).exp()/(norm*31*A(32)**31)
    qlo=((A(L)-E)/norm).lower();qhi=((A(U)+E)/norm+tail).upper();assert 0<qlo<qhi
    cert_contains(d['q'],qlo.fmpq(),qhi.fmpq())
    for key,v in [('Mstar',star),('quadrature_abs_error',E),('q_tail_bound',tail)]:cert_contains(d[key],*ep(v))
    new.append({'s':str(ss),'panels':n,'raw_exact_polynomial_sum':[str(L),str(U)],'Mstar':ep(star),'raw_Taylor_error':ep(E),'q_tail':ep(tail),'q_interval':[str(qlo.fmpq()),str(qhi.fmpq())]})
# Exact rational endpoint safety: delta increases in mu and decreases in nu.
AL,AU=dy(old['q']['lower']),dy(old['q']['upper'])
BL,BU=dy(cert['new_q'][0]['q']['lower']),dy(cert['new_q'][0]['q']['upper'])
CL,CU=dy(cert['new_q'][1]['q']['lower']),dy(cert['new_q'][1]['q']['upper'])
DL,DU=dy(cert['new_q'][2]['q']['lower']),dy(cert['new_q'][2]['q']['upper'])
R1L=BL*CL/AU**2;R1U=BU*CU/AL**2;R2L=DL/AU**2;R2U=DU/AL**2
cert_contains(cert['R1'],R1L,R1U);cert_contains(cert['R2'],R2L,R2U)
muL=R1L-1;muU=R1U-1;nuL=R2L/2-2*R1U+Q(3,2);nuU=R2U/2-2*R1L+Q(3,2)
cert_contains(cert['mu'],muL,muU);cert_contains(cert['nu'],nuL,nuU)
assert 0<muL and nuU>=muL**2 and nuL>muU**2
muused=dy(cert['mu']['lower']);nuused=dy(cert['nu']['upper'])
assert 0<muused<=muL and nuused>=nuU
# Independent stable positive recurrence for 1-T8^(-1/4), no subtractive exp/log.
mu=A(muused);nu=A(nuused);d=nu/mu;xx=1+d
t2minus=2*d*(xx+1);T2=1+t2minus;t4minus=2*t2minus*(T2+1);T4=1+t4minus;t8minus=2*t4minus*(T4+1);T8=1+t8minus
R=T8.sqrt().sqrt();delta=(mu**2/nu)*t8minus/(R*(R**3+R**2+R+1))
dl=dy(cert['delta_lower']['lower']);assert dl==dy(cert['delta_lower']['upper'])
assert 0<dl<Q(ep(delta)[0]) and dl<1
# Same exact declared inputs as Terra; M_old is a safe outward upper budget.
Mold=dy(old['M_pair']['upper']);assert Mold>=AU**4/2
M=A(Mold)*(1-A(dl))**2
C=1/(2*arb.pi()**A(Q(7,2)));H=2*arb.pi()**A(Q(7,2))/A(Q(7,2)).gamma()
I=M**A(Q(7,6))*C**A(Q(7,6))*H/896
cert_contains(cert['M_new_upper'],*ep(M));cert_contains(cert['I_new_upper'],*ep(I))
RI=Q(cert['rational_I_upper_numerator'],cert['rational_I_upper_denominator'])
assert RI==Q(7730422817512,10**18)
assert Q(ep(I)[1])<dy(cert['I_new_upper']['upper'])<RI<dy(old['I_pair']['lower'])
assert snap()==before
out={'status':'PASS independent Radau proof/scalar-interface audit; no repeated panels','environment':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'bits':768},'new_panel_counts':[1024,512,256],'new_panels_recomputed':0,'exact_R004_certificate_verbatim_reuse':True,'new_pair_moment_scalar_reassemblies':new,'exact_mu_interval':[str(muL),str(muU)],'exact_nu_interval':[str(nuL),str(nuU)],'safe_mu_lower_used':str(muused),'safe_nu_upper_used':str(nuused),'independent_delta_at_declared_inputs':ep(delta),'declared_delta_lower':str(dl),'declared_delta_below_independent_lower':True,'M_upper_expression':ep(M),'I_upper_expression':ep(I),'rational_complete_I_upper':str(RI),'strict_comparison':'I768upper < declared dyadic upper <7730422817512/10^18 < reused R004 budget lower','source_hashes':before,'source_before_after_equal':True,'author_checker_executed':False,'new_MC':0,'new_Green_nodes':0,'global_g_proved':False,'new_C7_bound':False,'elapsed_seconds':time.perf_counter()-start}
print(json.dumps(out,indent=2))

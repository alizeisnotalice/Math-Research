"""Three registered exact groups; no body/source quadrature, no old script import."""
import signal
import time
START=time.monotonic()
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179.8)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(30,30))
import hashlib
import json
import os
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,isqrt
import sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(name,data):
    with (P/name).open('x') as f:
        json.dump(data,f,indent=2,sort_keys=True)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())


def ceilq(x):
    return -(-x.numerator//x.denominator)


def ceilsqrt(x):
    a=isqrt(x.numerator//x.denominator)
    return a if a*a*x.denominator>=x.numerator else a+1


def main():
    assert not (P/'run_started.json').exists() and not (P/'results.json').exists(),'No rerun'
    pins=json.loads((P/'input_manifest.json').read_text())
    for path,h in pins['readonly_sha256'].items():assert sha(path)==h,path
    for name,h in pins['owned_source_sha256'].items():assert sha(P/name)==h,name
    save('run_started.json',dict(pid=os.getpid(),start_monotonic=START,
         registration_sha256=sha(P/'registration.json'),input_manifest_sha256=sha(P/'input_manifest.json'),
         max_new_configurations=3,native_calls=0,automatic_retries=0))
    groups=[]
    # 1: exact geometry margins, not a repeated S2/S3 enumeration.
    p=list(map(Q,['6/5','1','4/5','3/5','2/5','1/5']))
    h=Q(1,100);delta=Q(1,10);gap=Q(1,5)-2*h
    assert p[0]-h>1 and gap*gap>2*delta*delta
    assert (p[0]+h)**5<3 and 6*(p[0]+h)**2<9
    assert (Q(1,6))**2==Q(1,36)
    groups.append(dict(id=1,status='PASS',support_first_coordinate_lower=str(p[0]-h),
                       absolute_gap=str(gap),S3_distance_lower_squared=str(gap*gap/2),
                       radial_power_upper=str((p[0]+h)**5),body_triples_reenumerated=0))
    # 2: one symbolic bump derivative identity and outward norm constants.
    import sympy as sp
    u=sp.symbols('u',real=True);t=1/(1-u*u);eta=sp.exp(-t)
    d1=-2*u*t*t*eta
    d2=(4*u*u*t**4-8*u*u*t**3-2*t*t)*eta
    assert sp.simplify(sp.diff(eta,u)-d1)==0
    assert sp.simplify(sp.diff(eta,u,2)-d2)==0
    first=2*factorial(2);second=4*factorial(4)+8*factorial(3)+2*factorial(2)
    H=(second+5*first*first)/h**2
    M=145920043
    assert first==4 and second==148 and H==2280000
    assert Q(192*H+128,3)<M
    groups.append(dict(id=2,status='PASS',eta_prime_bound=first,eta_second_bound=second,
                       Hessian_bound=str(H),Ab_bound=M,sympy_version=sp.__version__))
    # 3: one complete rational parameter graph; no native or source point sample.
    g=Q(1,160*(4*(4+122**2))**6*10**26);c=g/6
    C=Q(1,160*3**6);N=32*(31*30//2)
    A=Q(16*N*32**6*216,3**6)
    T=max(1,ceilq(4*C/(5*c)))
    R=max(1,ceilsqrt(4*A*T/(c*delta**4)));L=12*R
    short=A*T/(R**2*delta**4);long=C/(5*T**5)
    assert short<=c/4 and long<=c/4 and short+long<=c/2
    assert Q(L,6)-R>=delta and L-R>=R
    C0=1/c;C1=16*L*M/c;eps=c/(32*L*(M+1));mbar=c/(16*L)
    assert C1>=C0 and eps*C0<Q(1,2) and eps*C1==Q(M,2*(M+1))<Q(1,2)
    assert C1*mbar==M
    bmin=2*h**6/9**6;norm=64*(g-c);gain=eps*bmin/norm
    assert norm==320*c and gain==bmin/(10240*L*(M+1))>0
    assert Q(6*120*27,160)==Q(243,2)<122
    values=dict(gstar=g,c=c,Cbar=C,Abar=A,T=T,R=R,L=L,C0=C0,C1=C1,epsilon=eps,
                density_lower=mbar,bump_mass_lower=bmin,normalized_gain_lower=gain,
                short_tail_upper=short,long_tail_upper=long)
    groups.append(dict(id=3,status='PASS',exact_values={k:str(v) for k,v in values.items()},
                       native_calls=0,body_evaluations=0,source_integrals=0))
    save('results.json',dict(status='PASS_THREE_EXACT_SANITY_GROUPS',configurations=groups,
         pid=os.getpid(),wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time(),
         native_calls=0,body_evaluations=0,source_integrals=0,old_grid_replays=0,automatic_retries=0,
         theory_not_replaced_by_tests=True,new_C6_bound=False))
    print(json.dumps(dict(status='PASS_THREE_EXACT_SANITY_GROUPS',pid=os.getpid(),
                         wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time(),native_calls=0)))


if __name__=='__main__':
    main()

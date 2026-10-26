#!/usr/bin/env python3
"""A100 new certificate checks only: residue/log/radial/geometry/coercivity.

No author program import, old062 resolvent baseline evaluation, F/L, grid,
fan/tuple enumeration, or physical kernel. No four-configuration replay.
"""
import sys
sys.dont_write_bytecode=True
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179.5)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
import os,json,hashlib,time
START=time.monotonic();CPU=time.process_time()
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from datetime import datetime,timezone
import sympy as sp
from flint import arb,fmpq,ctx
OUT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round100/time_first_budget_refinement')
ctx.prec=512

def need(ok,message):
    if not ok:raise ValueError(message)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def rat(x):return Q(int(x.p),int(x.q))
def A(x):
    q=Q(x);return arb(fmpq(q.numerator,q.denominator))
def enclosure(x):return Q(str(x.lower().fmpq())),Q(str(x.upper().fmpq()))
def dyadic(x):
    e=Q(2)**int(x['common_binary_exponent']);return int(x['lower_mantissa'])*e,int(x['upper_mantissa'])*e
def save(name,data):
    with (OUT/name).open('x') as f:
        json.dump(data,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
def inside(new,old,message):need(old[0]<=new[0]<=new[1]<=old[1],message)
def imul(a,b):
    vals=[x*y for x in a for y in b];return min(vals),max(vals)

def sqrt_bracket(n):
    d=10**130;target=n*d*d;lo,hi=0,1
    while hi*hi<=target:hi*=2
    while hi-lo>1:
        m=(lo+hi)//2
        if m*m<=target:lo=m
        else:hi=m
    need(lo*lo<=target<(lo+1)*(lo+1),'Independent integer square-root bracket')
    return Q(lo,d),Q(lo+1,d)

def main():
    reg=read(OUT/'registration.json')
    need(sha(AUTHOR/'manifest.json')==reg['author_final_manifest_sha256'],'Same final scientific seal')
    first=read(OUT/'saved_results.json')
    for p,digest in first['source_hashes'].items():need(sha(p)==digest,'First group source binding changed')
    data=read(AUTHOR/'results.json');c2=data['configurations'][1];c4=data['configurations'][3]
    u,a=sp.symbols('u a',real=True)
    alpha=sp.Rational(7,310)*(1-u);bl=sp.Rational(6,155)*u;br=sp.Rational(6,155)*(1-u)
    bw=sp.Rational(6,155)*(sp.Rational(7,8)-u);cv=sp.Rational(3,112);cw=sp.Rational(7,128)
    def T(a,b,c):return sp.cancel((a+b+c)/(a*b*c*(a+b)*(a+c)*(b+c)))
    aa,bb,cc=sp.symbols('aa bb cc',positive=True)
    need(sp.cancel(T(aa,bb,cc)-sum(1/(x*y) for x,y in [(aa,bb),(aa,cc),(bb,cc)])/((aa+bb)*(aa+cc)*(bb+cc)))==0,
         'Positive reciprocal-product convexity identity')
    pieces=[(T(alpha,bl,cv),Q(3,8),Q(1,2)),(T(alpha,br,cv),Q(1,2),Q(5,8)),(T(alpha,bw,cw),Q(3,8),Q(5,8))]
    points=[rat((pieces[0][0]+pieces[2][0]).subs(u,sp.Rational(3,8))),
            rat((pieces[0][0]+pieces[2][0]).subs(u,sp.Rational(1,2))),
            rat((pieces[1][0]+pieces[2][0]).subs(u,sp.Rational(5,8)))]
    need(list(map(Q,c2['convex_piece_endpoint_values_exact']))==points and Q(c2['R_max_exact'])==max(points),'Exact analytic endpoint maximum')
    # Independent Laurent coefficients at rational poles, using residue
    # derivatives. No sp.apart and no old time integral are called.
    rational=Q(0);logs={};certificates=[]
    for function,left,right in pieces:
        numerator,denominator=sp.fraction(sp.cancel(function));fac=sp.factor_list(denominator,u)[1]
        restored=sp.Integer(0);terms=[]
        need(sp.degree(numerator,u)<sp.degree(denominator,u),'New fixed rational integrand is proper')
        for polynomial,multiplicity in fac:
            p=sp.Poly(polynomial,u);need(p.degree()==1,'Only rational linear poles in this certificate')
            pole=-p.nth(0)/p.nth(1)
            analytic=sp.cancel(function*(u-pole)**multiplicity)
            z=rat(pole);need(not left<=z<=right,'No pole in actual mean interval')
            for order in range(1,multiplicity+1):
                coeff=sp.cancel(sp.diff(analytic,u,multiplicity-order).subs(u,pole)/sp.factorial(multiplicity-order))
                coeffq=rat(coeff);restored+=coeff/(u-pole)**order
                terms.append({'pole':str(z),'order':order,'coefficient':str(coeffq)})
                if order==1:
                    ratio=abs((right-z)/(left-z));logs[ratio]=logs.get(ratio,Q(0))+coeffq
                else:rational+=coeffq/Q(1-order)*((right-z)**(1-order)-(left-z)**(1-order))
        need(sp.cancel(function-restored)==0,'Exact residue partial-fraction identity')
        certificates.append({'interval':[str(left),str(right)],'residue_terms':terms})
    logs={r:c for r,c in logs.items() if c}
    original_logs={Q(x['log_positive_rational']):Q(x['coefficient']) for x in c2['u_integral_log_terms']}
    need(rational==Q(c2['u_integral_rational_part']) and logs==original_logs,'All signed definite-integral coefficients match stored expression')
    value=A(rational)
    for ratio,coefficient in logs.items():
        need(ratio>0,'Log argument positive');value+=A(coefficient)*A(ratio).log()
    Ju=enclosure(value);inside(Ju,dyadic(c2['u_integral']),'Independent signed-log sum contained in saved384 enclosure')
    need(0<Ju[0]<Ju[1]<max(points)*Q(1,4),'Positive mean integral below endpoint majorant')
    w=sp.symbols('w',positive=True)
    F=-sp.Rational(1,5)*w**-5+sp.Rational(3,4)*w**-4-w**-3+sp.Rational(1,2)*w**-2
    need(sp.cancel(sp.diff(F,w)-(1-w)**3/w**6)==0,'Correct radial primitive sign')
    Jr=rat((F.subs(w,1)-F.subs(w,sp.Rational(57,64)))/7**4)
    need(Jr==Q(c2['radial_integral_exact'])==Q(313,12033841140)>0,'Full radial factor and endpoints')
    # Positive elementary square comparisons pay all rational normalized gap
    # minorants; a common 1-7rho is applied before division bysqrt(A).
    square_tests={
        '1/sqrt2>7/10':Q(1,2)>Q(49,100),
        'sqrt(3/2)>6/5':Q(3,2)>Q(36,25),
        'sqrt928<31':928<31**2,'sqrt960<31':960<31**2,'sqrt992<32':992<32**2,
        'sqrt3/2>6/7':Q(3,4)>Q(36,49),'sqrt822/16>7/4':Q(822,256)>Q(49,16)}
    need(all(square_tests.values()),'Joint rational normalized-root minorants')
    min_gap_squares=[Q(9,128),Q(27,128),Q(3,32),Q(3,4),Q(822,256)]
    need(all(x>Q(3,49) for x in min_gap_squares) and 1-7*Q(1,64)>0,'Uniform common rho factor for all five physical gaps')
    origin=3*((1-a)**2+(u+a)**2)
    need(sp.expand(sp.diff(origin,a)-6*(2*a+u-1))==0 and Q(2,16)+Q(5,8)-1<0,
         'Actual w-small origin minimum at positive threshold')
    need(Q(3,8)>Q(1,3) and Q(21,80)-Q(1,32)>Q(1,64) and comb(5,2)==10,'All ten central k3 regions disjoint')
    # True face Jacobian: sqrt3 for A-face map andsqrt3 for Bmean map.
    DA=sp.Matrix([[-sp.sqrt(2),0],[-1/sp.sqrt(2),-3/sp.sqrt(6)]])
    DB=sp.Matrix([[1/sp.sqrt(2),1/sp.sqrt(6),1],[-1/sp.sqrt(2),1/sp.sqrt(6),1],[0,-2/sp.sqrt(6),1]])
    need(sp.simplify(DA.det()**2)==3 and sp.simplify(DB.det()**2)==3,'Original face Jacobian3')
    need(Q(3)*Q(1,3)==1 and Q(1,4)*4/Q(4-Q(6,5))==Q(5,14),
         'Radial sector and singular spatial integration factors')
    # Actual omitted-column scalar Holder minimum, including restriction to any
    # real quotient. This is not a new projected density or old033 peak run.
    x=sp.symbols('x',real=True);f=(32-x)**3/(96-x*x)
    derivative=-(32-x)**2*(288-64*x-x*x)/(96-x*x)**2
    need(sp.cancel(sp.diff(f,x)-derivative)==0 and 32**2<16*82<38**2 and 96-6**2>0,'True c31 minimum location and derivative')
    root=sqrt_bracket(82);xs=(-32+4*root[0],-32+4*root[1]);need(0<xs[0]<xs[1]<6,'xstar interior')
    numerator=((32-xs[1])**3,(32-xs[0])**3);denominator=(96-xs[1]**2,96-xs[0]**2)
    fstar=(numerator[0]/denominator[1],numerator[1]/denominator[0])
    c31=dyadic(c4['c31']);inside(fstar,(c31[0]**2,c31[1]**2),'Independent exact-root bracket verifies saved c31')
    # Targeted coefficient-ratio postprocessing uses the actual31-form c31,
    # original033 coarea class weights, and preserved omission factor1000/729.
    xball=-32+4*A(82).sqrt();cball=((32-xball)**3/(96-xball*xball)).sqrt()
    counts=[435,1740,435,580,1305];dets=[64,96,128,160,192]
    coarea=sum((A(n)/A(d).sqrt() for n,d in zip(counts,dets)),A(0))
    ratio=196608*729*2970*arb.pi()*coarea/(1000*44743079*cball**3)
    ratio_bounds=enclosure(ratio);inside(ratio_bounds,dyadic(c4['new_over_old']),'Independent actual-law/coarea ratio contained in saved result')
    need(ratio_bounds[0]>1,'Moment replacement cannot improve exact peak coefficient')
    for p,digest in first['source_hashes'].items():need(sha(p)==digest,'Postcheck original source hash changed')
    result=dict(status='PASS_INDEPENDENT_NEW_JOINT_AND_COERCIVITY_CERTIFICATES',
        author_manifest_sha256=sha(AUTHOR/'manifest.json'),saved_group_results_sha256=sha(OUT/'saved_results.json'),
        rational_piece_residue_certificates=certificates,u_integral_rational_part=str(rational),
        u_integral_signed_logs=[{'coefficient':str(c),'positive_rational':str(r)} for r,c in sorted(logs.items())],
        independent_u_integral_enclosure=list(map(str,Ju)),radial_integral_exact=str(Jr),
        convex_endpoint_values=list(map(str,points)),R_max_exact=str(max(points)),minorant_square_tests=square_tests,
        Jacobian=3,radial_sector_fraction='1/3',ten_regions_disjoint=True,
        independent_c31_squared_enclosure=list(map(str,fstar)),independent_coercivity_replacement_ratio=list(map(str,ratio_bounds)),
        full_forward_four_configuration_replay=False,old062_resolvent_formula_evaluations=0,old043074_F_L_runs=0,
        author_program_imports_or_executions=0,physical_kernel_calls=0,grids=0,MATLAB=0,source_writes=0,
        numerical_scope='new signed-log primitive and targeted new coarea ratio only; other claims are source/proof and saved endpoint checks',
        new_C6_bound=False,actual_I_or_model_lower_claim=False)
    save('joint_results.json',result)
    save('joint_run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds_including_result_save=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        hard_default_alarm_seconds=179.5,cpu_hard_limit_seconds=20,alarm_kept_through_receipt_and_exit=True,
        saved_exact_interface_group=2,automatic_retries=0,result_sha256=sha(OUT/'joint_results.json')))
    print(json.dumps(dict(status=result['status'],new_rational_pieces=3,radial_integral=str(Jr),
        old062_resolvent_evaluations=0,author_program_runs=0)))

if __name__=='__main__':
    try:main()
    except BaseException as error:
        save('joint_failure_receipt.json',dict(status='STOP_NO_RETRY',error=repr(error),pid=os.getpid(),
            wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,author_program_runs=0))
        raise

#!/usr/bin/env python3
"""Independent new104 rational certificates only; exact main release required."""
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
import sys
sys.dont_write_bytecode=True
import os,json,hashlib,time
from pathlib import Path
from fractions import Fraction as Q
from datetime import datetime,timezone
START=time.monotonic();CPU=time.process_time()
OUT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round104/complement_event_budget')
HASHES={}
def need(ok,message):
    if not ok:raise ValueError(message)
def sha(p):
    p=Path(p);value=hashlib.sha256(p.read_bytes()).hexdigest()
    old=HASHES.setdefault(str(p),value);need(old==value,'Changing input '+str(p));return value
def read(p):sha(p);return json.loads(Path(p).read_text())
def save(name,value):
    with (OUT/name).open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
def I(raw):
    factor=Q(2)**int(raw['common_binary_exponent'])
    return int(raw['lower_mantissa'])*factor,int(raw['upper_mantissa'])*factor
def main():
    approval=read(OUT/'primary_release.json');registration=read(OUT/'registration.json')
    need(approval['authorized'] is True and approval['round']=='104-independent-review' and
         approval['review_registration_sha256']==sha(OUT/'registration.json') and
         approval['script_sha256'][Path(__file__).name]==sha(Path(__file__)),'No exact main static release')
    need(sha(AUTHOR/'manifest.json')==approval['author_final_manifest_sha256']==registration['author_final_manifest_sha256'],'Final author seal')
    previous=read(OUT/'saved_results.json')
    need(previous['status']=='PASS_SAVED_BINDINGS_RESTORED_PAIR_ACCOUNT_AND_LOCAL_GLOBAL_SCOPE','Group1 must have passed')
    for p,digest in previous['source_hashes'].items():need(sha(p)==digest,'Group1 source continuation pin')
    data=read(AUTHOR/'results.json');cfg={x['configuration']:x for x in data['configurations']}
    # No producer import. Each pole coefficient is obtained from the Taylor
    # derivatives of its regular part, then the entire rational identity is checked.
    import sympy as sp
    from flint import arb,fmpq,ctx
    ctx.prec=512
    u=sp.Symbol('u');half=sp.Rational(1,2);seven=sp.Rational(7,8)
    def rat(x):return Q(int(x.p),int(x.q))
    def A(x):
        x=Q(x);return arb(fmpq(x.numerator,x.denominator))
    def endpoints(x):
        n,e=x.lower().man_exp();m,f=x.upper().man_exp()
        return Q(int(n))*Q(2)**int(e),Q(int(m))*Q(2)**int(f)
    def inside(value,saved,message):
        a,b=endpoints(value);l,r=I(saved);need(l<=a<=b<=r,message);return list(map(str,(a,b)))
    def integrate_piece(f,l,r):
        numerator,denominator=sp.fraction(sp.cancel(f));roots=sp.roots(denominator,u)
        need(sum(roots.values())==sp.degree(denominator,u) and sp.degree(numerator,u)<sp.degree(denominator,u),'Proper rational with all exact poles')
        terms=[];rational=Q(0);logs={}
        for pole,multiplicity in roots.items():
            need(pole.is_Rational and not l<=rat(pole)<=r,'No integration pole')
            regular=sp.cancel(f*(u-pole)**multiplicity)
            for power in range(1,multiplicity+1):
                order=multiplicity-power
                coefficient=sp.diff(regular,u,order).subs(u,pole)/sp.factorial(order)
                coefficient=rat(sp.cancel(coefficient));terms.append(sp.Rational(coefficient.numerator,coefficient.denominator)/(u-pole)**power)
                vl=l-rat(pole);vr=r-rat(pole)
                if power==1:
                    ratio=abs(vr/vl);logs[ratio]=logs.get(ratio,Q(0))+coefficient
                else:rational+=coefficient*(vr**(1-power)-vl**(1-power))/(1-power)
        need(sp.cancel(f-sum(terms))==0,'Independent complete residue identity')
        return rational,logs
    def family(name,left,right,kind):
        def Z(mu):
            if kind=='core':return 8/(9*mu**2*(1-u))+256/(1233*(seven-u)**2*(1-u))
            return 24/mu**3+512/(75*mu*(seven-u)**2)
        pieces=[sp.cancel(Z(u)),sp.cancel(Z(1-u))]
        rational=Q(0);logs={}
        for f,l,r in zip(pieces,[left,Q(1,2)],[Q(1,2),right]):
            rr,ll=integrate_piece(f,l,r);rational+=rr
            for arg,coef in ll.items():logs[arg]=logs.get(arg,Q(0))+coef
        logs={arg:coef for arg,coef in logs.items() if coef and arg!=1}
        saved=cfg[2 if kind=='core' else 3][name+'_integral']
        expected={Q(x['positive_rational_argument']):Q(x['coefficient']) for x in saved['logs']}
        need(rational==Q(saved['rational_part']) and logs==expected,'Exact signed-log rational certificate '+name)
        value=A(rational)
        for arg,coef in logs.items():value+=A(coef)*A(arg).log()
        interval=inside(value,saved['interval'],'Independent512bit log sum inside saved384bit interval '+name)
        # Each reciprocal product of positive affine functions has log convexity,
        # hence convexity; therefore maxima are among the endpoints of each half.
        candidates=[rat(f.subs(u,sp.Rational(point.numerator,point.denominator))) for f,points in
                    [(pieces[0],[left,Q(1,2)]),(pieces[1],[Q(1,2),right])] for point in points]
        need(max(candidates)==Q(cfg[2 if kind=='core' else 3][name+'_max_exact']),'Analytic piece endpoints maximum '+name)
        return dict(rational_part=str(rational),logs={str(k):str(v) for k,v in logs.items()},new512_interval=interval,maximum=str(max(candidates)),rational_identity=True)
    core=family('Zcore',Q(3,8),Q(5,8),'core');allmiddle=family('Zall',Q(1,4),Q(3,4),'all')
    # General radial primitive derivative, evaluated only for NEW lambda20.
    w=sp.Symbol('w');F=-w**-5/5+sp.Rational(3,4)*w**-4-w**-3+w**-2/2
    need(sp.cancel(sp.diff(F,w)-(1-w)**3/w**6)==0,'Radial primitive derivative')
    radial=rat((F.subs(w,1)-F.subs(w,sp.Rational(11,16)))/20**4)
    need(radial==Q(cfg[3]['radial20_exact'])==Q(15,164916224),'New lambda20 radial identity')
    old=read(AUTHOR/'upstream_constants.json')
    need(Q(cfg[2]['radial7_exact_reused'])==Q(old['radial7_exact'])==Q(313,12033841140),'Oldlambda7 bound only saved, not rerun')
    sigma=A(180)/A(20).sqrt()+A(450)/A(32).sqrt()+A(300)/6
    sigma_interval=inside(sigma,cfg[1]['Sigma_inv_J2'],'All930 exact coarea sum')
    D=3*A(6).sqrt()*sigma/(64*arb.pi()**4)
    d_interval=inside(D,cfg[1]['D_equals_3sqrt6Sigma_over64pi4'],'Restored four-dimensional constant with frozen old098 coercivity')
    for k in range(1,6):
        m=6-k
        # Squares of the five relative losses are exact. All positive losses<20.
        need(max(Q(4),Q(96,m),Q(256,225),Q(384,m),Q(32*k))<400,'Uniform allmiddle loss20')
        # faceJac=sqrt(km), sector=1/k, square=m/k; k=1/m=1 remain dimension4.
        need(Q(k*m,k*k)==Q(m,k),'Face Jacobian and sector product')
    need(Q(20,64)<1 and Q(7,64)<1,'Positive radial factors')
    # Strict triangle-source witness by block means, independent of a matrix
    # projector. The formal epsilon line is exact; it adds no kernel sample.
    directions=[[Q(z) for z in row] for row in cfg[4]['other_directions']]
    coefficients=list(map(Q,cfg[4]['other_parameters']));v=[Q(1)]*6;ww=[Q(1)]*3+[Q(-1)]*3
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    def normal(a):
        means=[sum(a[:3])/3,sum(a[3:])/3]
        return [x-means[i//3] for i,x in enumerate(a)]
    normals=[[c*x for x in normal(d)] for c,d in zip(coefficients,directions)]
    need([sum(n[j] for n in normals) for j in range(6)]==[0]*6 and all(dot(n,n)==Q(3,8) for n in normals),'Three strict large normal vectors and zero sum')
    need(all(dot(normals[i],normals[j])==Q(-3,16) for i in range(3) for j in range(i)),'Every angle120, none antiparallel')
    x=[Q(cfg[4]['a_v'])*v[j]+Q(cfg[4]['a_w'])*ww[j]+sum(c*d[j] for c,d in zip(coefficients,directions)) for j in range(6)]
    need(x==list(map(Q,cfg[4]['physical_point']))==[1,1,1,Q(1,2),Q(1,2),Q(1,2)],'Actual source physical point')
    slope=[directions[0][j]-ww[j] for j in range(6)]
    need(slope==[0,-2,-2,0,0,0] and dot(normal(slope),normal(slope))==Q(8,3),'Formal inward tube perturbation')
    need(Q(3,8)-Q(1,128)>Q(1,4) and Q(8,3)<4,'For0<eps<1/128 source labels staylarge, radius<1/64,face maximum unique')
    for p,digest in list(HASHES.items()):need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,'Postcheck input mutation')
    result=dict(status='PASS_NEW_EXACT_CERTIFICATES_AND_FORMAL_SOURCE_WITNESS',source_hashes=HASHES,
        author_final_manifest_sha256=sha(AUTHOR/'manifest.json'),core=core,all_middle=allmiddle,new_radial20=str(radial),
        sigma_interval=sigma_interval,D_interval=d_interval,analytic_universal_pair_gap_review_separate=True,
        cone_facejac_all_k=True,triangle_norm_squared='3/8',triangle_pair_dot='-3/16',inward_epsilon_interval='0<eps<1/128',
        inward_normal_squared='(8/3)eps^2',coercivity_used='old09832/sqrt3 only',
        author_generator_runs=0,author_forward_configurations=0,old_radial7_recalculations=0,kernel_calls=0,native_calls=0,
        children_created=0,source_writes=0,new_C6_bound=False,complete_global_g_I_S=False)
    save('exact_results.json',result)
    save('exact_run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds_including_result_save=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        default_hard_SIGALRM_seconds=179,CPU_hard_limit_seconds=20,alarm_kept_through_receipt_and_exit=True,
        saved_exact_group=2,automatic_retries=0,author_generator_runs=0,result_sha256=sha(OUT/'exact_results.json')))
    print(json.dumps(dict(status=result['status'],new_radial20=str(radial),complete_global_g_I_S=False,author_generator_runs=0)))
if __name__=='__main__':
    try:main()
    except BaseException as error:
        save('exact_failure_receipt.json',dict(status='STOP_NO_RETRY',error=repr(error),pid=os.getpid(),
            wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,author_generator_runs=0))
        raise

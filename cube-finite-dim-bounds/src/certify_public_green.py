"""Certified fixed-parameter cube Green model costs via Arb analytic integration.
No optimization. Tail estimates are the previously exact U02 estimates.
"""
from pathlib import Path
import hashlib,json,time
from fractions import Fraction as F
import flint
from flint import arb,acb,ctx
ctx.prec=96
P=arb.pi();ZP=acb(P);DELTA=arb(1)/1000;SPLIT=arb(3)/20;END=arb(32)
MIRROR=8;SPECTRAL=40

def early(t,analytic,a):
    rt=t.sqrt(analytic=analytic)
    sv=acb(1)
    for k in range(9):sv-=2*((-1)**k)*((2*k+1)/(2*rt)).erfc()
    pv=acb(0)
    for m in range(-MIRROR,MIRROR+1):
        pv+=(-(a-4*m)**2/(4*t)).exp()-(-(a-(4*m+2))**2/(4*t)).exp()
    pv/=(4*ZP*t).sqrt(analytic=analytic)
    return sv,pv

def late(t,analytic,a):
    sv=acb(0);pv=acb(0)
    for k in range(SPECTRAL):
        odd=2*k+1;e=(-ZP**2*odd**2*t/4).exp()
        sv+=4/ZP*((-1)**k)/odd*e
        pv+=(ZP*odd*a/2).cos()*e
    return sv,pv

def exact_interval(x):
    return [str(x.lower().fmpq()),str(x.upper().fmpq())]

def run():
    rows=[]
    prior=json.loads((Path(__file__).resolve().parents[1]/'upper/u02_results.json').read_text())
    cases=[(3,'1/2'),(4,'1/2')]+[(r['n'],r['a_decimal_input']) for r in prior['rows']]
    for n,astr in cases:
        exact_a=F(astr); space_a=acb(exact_a.numerator)/exact_a.denominator
        start=time.monotonic();parts=[]
        for fn,time_left,time_right in [(early,DELTA,SPLIT),(late,SPLIT,END)]:
            def integrand(t,analytic):
                sv,pv=fn(t,analytic,space_a)
                return sv**n+acb(0,1)*pv**n
            z=acb.integral(integrand,acb(time_left),acb(time_right),rel_tol=arb('1e-16'),abs_tol=arb('1e-18'),eval_limit=100000,depth_limit=40,deg_limit=80)
            assert z.is_finite(),(n,fn.__name__,z)
            parts.append(z)
        total=sum(parts,acb(DELTA))
        # Uniform safe envelope larger than all exact U02 omitted-tail bounds.
        omitted=arb(0,arb('1e-40'))
        mass=total.real+omitted;green=total.imag+omitted
        assert green>0
        cost=mass/(2**n*space_a.real**n*green)
        radial=27*P*arb(3).sqrt()/16 if n==3 else 2*P**2
        assert cost<radial
        bound=F(str(cost.upper().fmpq())); denominator=10**6
        simple_num=(bound.numerator*denominator)//bound.denominator+1
        simple=F(simple_num,denominator)
        rows.append({'n':n,'a':str(exact_a),'R':str(1/exact_a),'mass_interval':exact_interval(mass),'corner_green_interval':exact_interval(green),'model_cost_interval':exact_interval(cost),'model_cost_ball':cost.str(30),'exact_simple_model_upper':str(simple),'radial_model_cost_ball':radial.str(30),'strictly_below_radial_model':True,'strictly_below_general_target':bool(cost<2**n),'seconds':time.monotonic()-start,'certified_finite_integral_balls':[z.str(30) for z in parts]})
    return {'date':'2026-10-25','status':'certified model costs, not improvement of current global upper targets','precision_bits':ctx.prec,'flint_version':flint.__version__,'tail_proof_sha256':hashlib.sha256((Path(__file__).resolve().parents[1]/'upper/u02_tail_bounds.py').read_bytes()).hexdigest(),'cached_candidates_sha256':hashlib.sha256((Path(__file__).resolve().parents[1]/'upper/u02_results.json').read_bytes()).hexdigest(),'integrator':'python-flint acb.integral with analytic sqrt flag','parameters':{'delta':'1/1000','split':'3/20','L':32,'mirror_M':8,'spectral_N':40},'omitted_tail_absolute_envelope':'1e-40 each mass/Green; exact U02 tail proof reused','rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
    out=run();dst=Path(__file__).resolve().parents[1]/'results/round008_public_green.json';dst.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

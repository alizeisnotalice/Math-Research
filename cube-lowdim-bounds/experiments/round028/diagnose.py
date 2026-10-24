#!/usr/bin/env python3
"""Fixed one Fraction K point, finite symbolic checks, 18 numerical A values; no import experiment."""
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
from time import monotonic
import warnings

ROOT=Path(__file__).resolve().parent
SECONDS=180.
PER_VALUE_CALLBACKS,TOTAL_CALLBACKS=10000,180000
ALPHA=Q(47,40)
BETAS=(ALPHA,Q(6,5),Q(5,4))
SHAPES=((Q(1,2),Q(1,2)),(Q(1,4),Q(3,4)),(Q(1,16),Q(15,16)))
CANDIDATE_SHA='52a02acb669dcd33c05794abccab83658d02cf723a4a76cf3eadd2b38203fb97'
SERIES_TERMS,MP_DIGITS=48,60
QUAD_ABS,QUAD_REL,QUAD_LIMIT=2e-10,2e-10,200


class BudgetStop(Exception):pass
class NumericalFailure(Exception):pass


def dump(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temporary.replace(path)


def exact_B(value):
    x=abs(value)
    if x>=2:return Q(0)
    if x<=1:return Q(2,3)-x*x+x**3/2
    return (2-x)**3/6


def orbit(i,j):
    return {(x,y) for x in (i,-i) for y in (j,-j)}|{(y,x) for x in (i,-i) for y in (j,-j)}


def fraction_point(guard):
    path=ROOT.parent/'round08/candidate.json'
    if sha256(path.read_bytes()).hexdigest()!=CANDIDATE_SHA:raise NumericalFailure('Frozen candidate hash mismatch.')
    candidate=json.loads(path.read_text());u,v=Q(49,24),Q(43,96)
    points_path=ROOT.parent/'round06/point_manifest.json'
    if json.loads(points_path.read_text())['K'][22610]!=[str(u),str(v)]:
        raise NumericalFailure('Frozen point index/coordinate source mismatch.')
    if json.loads((ROOT.parent/'round08/result.json').read_text())['candidate_sha256']!=CANDIDATE_SHA:
        raise NumericalFailure('Frozen candidate/run provenance mismatch.')
    h,R=Q(candidate['h']),Q(candidate['R'])
    if (Q(candidate['alpha']),R,h,int(candidate['N']))!=(ALPHA,Q(33,20),Q(1,32),80):
        raise NumericalFailure('Frozen Fraction input parameters changed.')
    if not u+v>R:raise NumericalFailure('The single selected point must be outside the singular truncation.')
    value=Q(0);active=0
    for key,coefficient in candidate['a'].items():
        guard();i,j=map(int,key.split(','))
        for x,y in orbit(i,j):
            bx,by=exact_B(u/h-x),exact_B(v/h-y)
            if bx and by:value+=Q(coefficient)*bx*by;active+=1
    return {'u':str(u),'v':str(v),'r':str(u+v),'K_index':22610,'candidate_sha256':CANDIDATE_SHA,
        'point_manifest_sha256':sha256(points_path.read_bytes()).hexdigest(),
        'singular_part_exact':'0','spline_and_K_exact':str(value),'active_terms':active,
        'strict_sign':'negative' if value<0 else 'positive' if value>0 else 'zero',
        'C1_at_this_point':'violated' if value<0 else 'nonnegative',
        'scope':'One exact rational point only; no AK, repair or global proof.'}


def symbolic_checks(s,guard,announce):
    b,R,r,x,z,t=s.symbols('b R r x z t',positive=True)
    primitive=2*r**(2-b)/(2-b)-R**(-b)*r*r
    mass=b*R**(2-b)/(2-b)
    F=b*R**2/((2-b)*(R**b-1))
    Fx=b/(2-b)*x**(2/b)/(x-1)
    checks=[]
    tasks=[
        ('diamond_half_density',lambda:s.simplify(s.diff(primitive,r)-2*r*(r**(-b)-R**(-b)))),
        ('lower_endpoint_beta_below_2',lambda:s.limit(primitive.subs(b,2-z),r,0,dir='+')),
        ('half_mass_upper_endpoint',lambda:s.simplify(primitive.subs(r,R)-mass)),
        ('dominating_coefficient',lambda:s.simplify((1/(1-R**(-b)))*(1-R**(-b))-1)),
        ('dominated_mass',lambda:s.simplify(mass/(1-R**(-b))-F)),
        ('R_log_derivative',lambda:s.simplify(s.diff(F,R)/F-((2-b)*R**b-2)/(R*(R**b-1)))),
        ('stationary_x_log_derivative',lambda:s.simplify(s.diff(Fx,x)/Fx-((2-b)*x-2)/(b*x*(x-1)))),
        ('stationary_cost',lambda:s.simplify(Fx.subs({b:2-z,x:2/z},simultaneous=True)-(2/z)**(2/(2-z)))),
        ('monotone_integrand_derivative',lambda:s.simplify(s.diff(1/(1-z*t),z)-t/(1-z*t)**2)),
        ('beta1_cost_4',lambda:s.simplify(((2/(2-b))**(2/b)).subs(b,1)-4))]
    for name,task in tasks:
        guard();announce(name,'started');tic=monotonic();residual=task();guard()
        row={'name':name,'residual':str(residual),'passed':residual==0,'seconds':monotonic()-tic}
        checks.append(row);announce(name,'completed',row)
        if residual!=0:raise NumericalFailure('Fixed symbolic identity failed: '+name)
    return {'checks':checks,'scope':'Finite identities only; >=4 additionally uses analytic monotonicity for beta>=1, not a numerical grid.'}


def scipy_direct(alpha,beta,u,v,quad,callback,row):
    """Original symmetric differences, with an analytic even-series prefix."""
    r=u+v;directions=[];row['directions']=directions;total=0.;total_error=0.
    for a,b in ((u,v),(v,u)):
        cut=min(a/2,r/8)
        coefficient=1.;terms=[];next_term=None
        for k in range(1,SERIES_TERMS+2):
            coefficient*=((beta+2*k-2)*(beta+2*k-1))/((2*k-1)*(2*k))
            term=2*r**(-beta)*coefficient*r**(-2*k)*cut**(2*k-alpha)/(2*k-alpha)
            if k<=SERIES_TERMS:terms.append(term)
            else:next_term=term
        prefix=math.fsum(terms)
        # For beta<=5/4 and k>=1, coefficient ratio <=2. Integrated term
        # ratios are <=rho=2(cut/r)^2<=1/32. The exact analytic tail is
        # <= first omitted/(1-rho). Its FLOAT evaluation below is NOT an interval.
        rho=2*(cut/r)**2;tail_estimate=next_term/(1-rho)
        def integrand(t):
            if t<a:
                z=t/r
                difference=r**(-beta)*(math.expm1(-beta*math.log1p(z))+
                                        math.expm1(-beta*math.log1p(-z)))
            else:difference=(r+t)**(-beta)+(t-a+b)**(-beta)-2*r**(-beta)
            return difference*t**(-1-alpha)
        entry={'a':a,'b':b,'prefix_cut':cut,'series_terms':SERIES_TERMS,
            'prefix_value':prefix,'analytic_tail_estimate_float':tail_estimate,
            'tail_strategy':'Positive even binomial series, coefficient-ratio<=2; floating bound is not rigorous interval.'}
        entry['segments']=[];directions.append(entry);subtotal=prefix;error=tail_estimate
        for lower,upper in ((cut,a),(a,math.inf)):
            before=row['callbacks_started'];tic=monotonic()
            segment={'lower':str(lower),'upper':str(upper),'status':'started','callback_start_index':before}
            entry['segments'].append(segment);row['current_segment']=segment
            with warnings.catch_warnings(record=True) as seen:
                warnings.simplefilter('always')
                output=quad(callback(integrand),lower,upper,epsabs=QUAD_ABS,epsrel=QUAD_REL,
                            limit=QUAD_LIMIT,full_output=1)
            value,err,info=output[:3]
            segment.update(status='completed',value=float(value) if math.isfinite(value) else None,
                estimated_error=float(err) if math.isfinite(err) else None,
                quad_neval=int(info['neval']),callbacks=row['callbacks_started']-before,
                seconds=monotonic()-tic,warnings=[str(w.message) for w in seen],
                message=str(output[3]) if len(output)>3 else None)
            row['current_segment']=None
            if len(output)>3 or seen or not math.isfinite(value) or not math.isfinite(err):
                raise NumericalFailure('SciPy quad did not return clean finite convergence.')
            if segment['callbacks']!=segment['quad_neval']:raise NumericalFailure('Quad callback/neval ledger differs.')
            subtotal+=float(value);error+=float(err)
        entry.update(integral_without_factor_2=subtotal,estimated_error_without_factor_2=error)
        total+=2*subtotal;total_error+=2*error
    return total,total_error


def mpmath_second(alpha_q,beta_q,u_q,v_q,mp,callback,row):
    """Distributional g'' integral, including its negative atom; different expression."""
    aalpha=mp.mpf(alpha_q.numerator)/alpha_q.denominator
    beta=mp.mpf(beta_q.numerator)/beta_q.denominator
    u=mp.mpf(u_q.numerator)/u_q.denominator;v=mp.mpf(v_q.numerator)/v_q.denominator
    exponent=1/(2-aalpha);total=mp.mpf(0);error=mp.mpf(0);directions=[];row['directions']=directions
    for a,b in ((u,v),(v,u)):
        r=a+b
        def left(t):
            if t==1:return mp.mpf(0)  # analytic endpoint limit, not 0*infinity
            s=t/(1-t)
            return (s+b)**(-beta-2)*(s+a)**(1-aalpha)/(1-t)**2
        def middle(t):return a**(2-aalpha)/(2-aalpha)*(r-a*t**exponent)**(-beta-2)
        def right(t):
            if t==1:return mp.mpf(0)
            w=t/(1-t)
            return (r+w**exponent)**(-beta-2)/((2-aalpha)*(1-t)**2)
        segments=[];entry={'a':mp.nstr(a,40),'b':mp.nstr(b,40),'segments':segments};directions.append(entry)
        regular=mp.mpf(0);regular_error=mp.mpf(0)
        for name,function in (('s_negative_compactified',left),('s_0_to_a_desingularized',middle),
                              ('s_above_a_desingularized_compactified',right)):
            before=row['callbacks_started'];tic=monotonic()
            segment={'name':name,'status':'started','callback_start_index':before}
            segments.append(segment);row['current_segment']=segment
            value,err=mp.quad(callback(function),[0,1],method='tanh-sinh',error=True,maxdegree=7)
            segment.update(status='completed',value=mp.nstr(value,65),estimated_error=mp.nstr(err,12),
                callbacks=row['callbacks_started']-before,seconds=monotonic()-tic)
            row['current_segment']=None
            if not mp.isfinite(value) or not mp.isfinite(err):raise NumericalFailure('mpmath integral nonfinite.')
            regular+=value;regular_error+=abs(err)
        atom=-2*beta*b**(-beta-1)*a**(1-aalpha)
        factor=2/(aalpha*(aalpha-1))
        value=factor*(beta*(beta+1)*regular+atom)
        err=factor*beta*(beta+1)*regular_error
        entry.update(regular_integral=mp.nstr(regular,65),negative_atom=mp.nstr(atom,65),
            normalization_factor=mp.nstr(factor,65),coordinate_A=mp.nstr(value,65))
        total+=value;error+=err
    return total,error


def main():
    if (ROOT/'result.json').exists() or (ROOT/'calls.jsonl').exists():raise SystemExit('Existing diagnostic evidence: no repeat.')
    start=monotonic();trace=(ROOT/'calls.jsonl').open('x')
    counts={name:0 for name in ('Fraction_started','Fraction_completed','symbolic_started','symbolic_completed',
        'A_started','A_completed','A_failed','A_interrupted','callbacks_started','callbacks_completed')}
    result={'status':'started','counts':counts,'A_rows':[],'symbolic_rows':[],'comparisons':[],
        'limits':{'seconds':180,'Fraction_points':1,'symbolic_identities':10,'A_values':18,
                  'callbacks_per_A':10000,'callbacks_total':180000},
        'scope':'Exact Fraction point + symbolic identities + nonrigorous finite quadrature; no LP/full checker/candidate repair.'}
    def save():
        result['A_not_started']=18-counts['A_started'];result['elapsed_seconds']=monotonic()-start
        dump(ROOT/'result.json',result)
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('180s_internal_deadline')
    def event(value):trace.write(json.dumps(value,ensure_ascii=False)+'\n');trace.flush()
    save()
    try:
        result['phase']='Fraction_point';counts['Fraction_started']=1;save();tic=monotonic()
        fraction=fraction_point(guard);fraction['seconds']=monotonic()-tic
        result['Fraction_point']=fraction;counts['Fraction_completed']=1;save()
        print(json.dumps({'Fraction_point_completed':fraction},ensure_ascii=False),flush=True)
        # A negative exact point is a successful diagnostic; only computational
        # failures stop the predeclared beta investigation.
        guard();import sympy as sympy
        result['phase']='fixed_symbolic_identities';result['sympy_version']=sympy.__version__;save()
        def announce(name,status,data=None):
            if status=='started':
                counts['symbolic_started']+=1;result['symbolic_in_progress']=name
            else:
                counts['symbolic_completed']+=1;result['symbolic_rows'].append(data);result['symbolic_in_progress']=None
            save()
        result['symbolic']=symbolic_checks(sympy,guard,announce)
        guard();import scipy
        from scipy.integrate import quad
        import mpmath as mp
        mp.mp.dps=MP_DIGITS
        alpha_mp=mp.mpf(ALPHA.numerator)/ALPHA.denominator
        reference_C=-8*mp.gamma(1-alpha_mp)*mp.gamma(2*alpha_mp)*mp.cos(mp.pi*alpha_mp/2)**2/mp.gamma(1+alpha_mp)
        if not mp.isfinite(reference_C) or reference_C<=0:raise NumericalFailure('Gamma normalization reference invalid.')
        result.update(scipy_version=scipy.__version__,mpmath_version=mp.__version__,mp_digits=MP_DIGITS,
            Gamma_C_alpha=mp.nstr(reference_C,65),
            frozen_points=[{'beta':str(beta),'u':str(u),'v':str(v)} for beta in BETAS for u,v in SHAPES],
            numerical_policy={'scipy_epsabs':QUAD_ABS,'scipy_epsrel':QUAD_REL,'quad_limit':QUAD_LIMIT,
                'small_t_series_terms':SERIES_TERMS,'mp_maxdegree':7,
                'agreement_relative_tolerance':'2e-7','mp_error_and_control_relative_tolerance':'1e-35',
                'no_precision_escalation_or_retry':True,'errors_are_not_outward_intervals':True})
        completed={};case=0
        for beta in BETAS:
            for u,v in SHAPES:
                case+=1
                if u+v!=1 or not 0<u<1 or not ALPHA<=beta<2:raise NumericalFailure('Fixed beta/shape domain changed.')
                beta_mp=mp.mpf(beta.numerator)/beta.denominator;convex_lower=(beta_mp/alpha_mp)*reference_C
                for method in ('scipy_direct_symmetric_difference','mpmath_distributional_second_derivative'):
                    guard()
                    if counts['A_started']>=18:raise BudgetStop('18_A_value_cap')
                    row={'case':case,'alpha':str(ALPHA),'beta':str(beta),'u':str(u),'v':str(v),
                        'method':method,'status':'started','callbacks_started':0,'callbacks_completed':0,
                        'started_seconds':monotonic()-start}
                    result['A_rows'].append(row);counts['A_started']+=1;result['phase']='finite_quadrature';event(row.copy());save();tic=monotonic()
                    def callback(function):
                        def bounded(value):
                            guard();row['callbacks_started']+=1;counts['callbacks_started']+=1
                            if row.get('current_segment') is not None:
                                row['current_segment']['callbacks_started_so_far']=row['callbacks_started']-row['current_segment']['callback_start_index']
                            # The cap-th invocation is marked started/interrupted before
                            # evaluating its integrand, so the invocation cap is not exceeded.
                            if row['callbacks_started']>=PER_VALUE_CALLBACKS or counts['callbacks_started']>=TOTAL_CALLBACKS:
                                raise BudgetStop('fixed_integrand_callback_cap')
                            out=function(value)
                            if not mp.isfinite(out):raise NumericalFailure('Nonfinite integrand callback.')
                            row['callbacks_completed']+=1;counts['callbacks_completed']+=1
                            if row['callbacks_started']%128==0:save()
                            return out
                        return bounded
                    try:
                        if method.startswith('scipy'):
                            value,error=scipy_direct(float(ALPHA),float(beta),float(u),float(v),quad,callback,row)
                            value_mp=mp.mpf(repr(value));error_mp=mp.mpf(repr(error))
                            tolerance=mp.mpf('2e-7')*(1+abs(value_mp))
                        else:
                            value_mp,error_mp=mpmath_second(ALPHA,beta,u,v,mp,callback,row)
                            tolerance=mp.mpf('1e-35')*(1+abs(value_mp))
                        if not mp.isfinite(value_mp) or not mp.isfinite(error_mp) or error_mp<0:
                            raise NumericalFailure('Nonfinite numerical A/error.')
                        if error_mp>tolerance:raise NumericalFailure('Fixed estimated-error threshold not met; no retry.')
                        if beta==ALPHA and abs(value_mp-reference_C)>tolerance:
                            raise NumericalFailure('Beta=alpha Gamma normalization control failed.')
                        if value_mp<convex_lower-tolerance:
                            raise NumericalFailure('Numerical convexity control failed; not a strict mathematical counterexample.')
                        row.update(status='completed',value=mp.nstr(value_mp,65),estimated_error=mp.nstr(error_mp,16),
                            seconds=monotonic()-tic,convex_lower_reference=mp.nstr(convex_lower,65),
                            numerical_difference_from_convex_lower=mp.nstr(value_mp-convex_lower,65),
                            Gamma_control_difference=mp.nstr(value_mp-reference_C,65) if beta==ALPHA else None,
                            numerical_only=True)
                        counts['A_completed']+=1;completed[(case,method)]=value_mp
                    except BudgetStop as exc:
                        row.update(status='interrupted',reason=str(exc),seconds=monotonic()-tic)
                        counts['A_interrupted']+=1;event(row.copy());save();raise
                    except Exception as exc:
                        row.update(status='failed_unknown',reason=f'{type(exc).__name__}: {exc}',seconds=monotonic()-tic)
                        counts['A_failed']+=1;event(row.copy());save();raise
                    event(row.copy());save();guard()
                first=completed[(case,'scipy_direct_symmetric_difference')]
                second=completed[(case,'mpmath_distributional_second_derivative')]
                difference=abs(first-second);tolerance=mp.mpf('2e-7')*(1+abs(second))
                comparison={'case':case,'beta':str(beta),'u':str(u),'v':str(v),
                    'difference':mp.nstr(difference,65),'fixed_tolerance':mp.nstr(tolerance,30),
                    'numerically_consistent':difference<=tolerance,'strict_interval_agreement':False}
                result['comparisons'].append(comparison);save()
                if difference>tolerance:raise NumericalFailure('Two numerical expressions disagree; stop without extra points.')
        result['status']='fixed_Fraction_symbolic_and_18_A_values_closed'
    except BudgetStop as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_saved_checkpoint',reason=f'{type(exc).__name__}: {exc}')
    finally:trace.close()
    assert counts['Fraction_started']<=1 and counts['symbolic_started']<=10 and counts['A_started']<=18
    assert counts['A_started']==counts['A_completed']+counts['A_failed']+counts['A_interrupted']
    assert counts['callbacks_started']<=TOTAL_CALLBACKS
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        trace_sha256=sha256((ROOT/'calls.jsonl').read_bytes()).hexdigest(),python=sys.version)
    save();print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Same positive L1 input: bounded-entropy physical-scale variation pressure.

Certified rational moment/boundary lower bounds are separate from MC profiles.
"""
from pathlib import Path
from fractions import Fraction as F
import datetime,hashlib,json,math,time
import numpy as np

HERE=Path(__file__).resolve().parent
PREFIX='physical_scale_bounded_entropy_variation_20261007'
EPS=F(1,8);LAMBDA=F(1);DEGREE=10


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compact(v,precision=20):
    den=10**precision;num=v.numerator*den//v.denominator
    return dict(lower=str(F(num,den)),upper=str(F(num+1,den)),
                exact_value_binary_sha256=hashlib.sha256((str(v.numerator.bit_length())+':'+hex(v.numerator)+':'+hex(v.denominator)).encode()).hexdigest())
def atan_interval(x,N=24):
    total=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(N)),F(0))
    other=total+(-1)**N*x**(2*N+1)/F(2*N+1)
    return min(total,other),max(total,other)
def pi_interval():
    l5,u5=atan_interval(F(1,5));l239,u239=atan_interval(F(1,239))
    return 16*l5-4*u239,16*u5-4*l239


def certify_peaks(n,M,D,plo,phi):
    result=[];total=F(0)
    for k in range(M,2*M):
        h=F(2*k+1,2)
        slo=EPS**2/(phi*h)**2;shi=EPS**2/(plo*h)**2
        variance_lo=(1+slo/2)**n-1
        fourth_hi=(1+3*shi+3*shi**2/8)**n-4*(1+3*slo/2)**n+6*(1+shi/2)**n-3
        assert fourth_hi>=0 and variance_lo>fourth_hi
        G_lower=(variance_lo-fourth_hi)/18
        boundary_each=F(8*n,D)
        real_rise_lower=G_lower-2*boundary_each
        assert real_rise_lower>0
        total+=real_rise_lower
        result.append(dict(k=k,peak_L=str(h/M),preceding_valley_L=str(F(k,M)),
            amplitude_squared_interval=[compact(slo),compact(shi)],
            variance_lower=compact(variance_lo),fourth_central_moment_upper=compact(fourth_hi),
            periodic_entropy_excess_lower=compact(G_lower),
            actual_L1_input_rise_lower_per_W=compact(real_rise_lower)))
    universal=F(M,589824)
    assert total>universal
    return result,total


def one_dim_original(x,L,M,D):
    lo=max(x-L/2,-F(D,2));hi=min(x+L/2,F(D,2))
    if hi<=lo:return 0.
    phase_lo=float((M*lo)%1);phase_hi=float((M*hi)%1)
    integral=float(hi-lo)+float(EPS)/(2*math.pi*M)*(math.sin(2*math.pi*phase_hi)-math.sin(2*math.pi*phase_lo))
    return integral/float(L)


def input_validation(n,M,D):
    records=[];maxerr=0.
    for L in [F(1),F(2*M+1,2*M),F(3,2),F(2)]:
        for x in [F(0),F(3,32*M),F(-5,64*M)]:
            direct=one_dim_original(x,L,M,D)
            expected=1+float(EPS)*np.sinc(M*float(L))*math.cos(2*math.pi*float(M*x))
            error=abs(direct-expected);maxerr=max(maxerr,error)
            assert error<3e-14
            records.append(dict(L=str(L),x=str(x),direct_original_integral=direct,
                periodic_formula=expected,error=error,kind='interior_original_input'))
        for x in [F(D,2)-F(1,4),F(D,2)+F(1,4),-F(D,2)+F(1,4)]:
            direct=one_dim_original(x,L,M,D);assert 0<=direct<=1+float(EPS)
            records.append(dict(L=str(L),x=str(x),direct_original_integral=direct,
                kind='boundary_clipped_original_input'))
    assert M*D==int(M*D) and (1-EPS)>0
    return dict(source_density='f(x)=1_{[-D/2,D/2]^n}(x) product_i[1+epsilon cos(2pi M x_i)]',
        n=n,M=M,D=D,epsilon=str(EPS),lambda_=str(LAMBDA),
        full_source_mass_formula='D^n',full_mass_integer_bit_length=(D**n).bit_length(),
        full_mass_binary_sha256=hashlib.sha256(hex(D**n).encode()).hexdigest(),
        input_positivity_certificate='each factor>=7/8>0 on source box',
        exact_integer_period_certificate=M*D,exact_mass_certificate='Each one-dimensional source factor has integral D since MD is integer; product source has W=D^n.',
        full_space_vs_periodic_entropy_error_per_W_upper=str(F(8*n,D)),
        boundary_error_scope='Analytic probability/mass union bound, not a Monte Carlo boundary approximation.',
        direct_formula_validation=records,max_float_validation_error=maxerr,
        numerical_validation_scope='Fraction integration limits and stable phase reduction; trigonometric values floating, not interval certificates.')


def profile(n,M,samples,nodes,seed):
    rng=np.random.default_rng(seed)
    cosines=np.cos(2*np.pi*rng.random((samples,n)))
    sums=np.empty((samples,DEGREE));power=np.ones_like(cosines)
    for j in range(DEGREE):
        power*=cosines;sums[:,j]=power.sum(axis=1)
    L=np.linspace(1,2,nodes);a=float(EPS)*np.sinc(M*L)
    coeff=np.stack([(-1)**j*a**(j+1)/(j+1) for j in range(DEGREE)])
    logu=sums@coeff
    u=np.exp(logu)
    # E u=1 analytically. Subtract the exact linear term as a control variate.
    bregman=(u-1)/2-np.log1p((u-1)/2)
    mean=bregman.mean(axis=0);se=bregman.std(axis=0,ddof=1)/math.sqrt(samples)
    fine=float(np.maximum(np.diff(mean),0).sum())
    coarse=float(np.maximum(np.diff(mean[::2]),0).sum())
    amax=float(EPS)/(math.pi*M)
    log_taylor_error=n*amax**(DEGREE+1)/((DEGREE+1)*(1-amax))
    direct_errors=[]
    for L0 in [F(1),F(2*M+1,2*M),F(3,2),F(2)]:
        aa=float(EPS)*np.sinc(M*float(L0))
        exact_log=np.log1p(aa*cosines[:64]).sum(axis=1)
        poly=sum(((-1)**j*aa**(j+1)/(j+1)*sums[:64,j] for j in range(DEGREE)))
        direct_errors.append(float(np.max(abs(exact_log-poly))))
    assert max(direct_errors)<3e-14
    return dict(seed=seed,samples=samples,nodes=nodes,L=L.tolist(),
        periodic_entropy_excess_profile=mean.tolist(),iid_sample_standard_errors=se.tolist(),
        empirical_positive_variation=fine,nested_coarse_positive_variation=coarse,
        quadrature_difference=fine-coarse,
        arithmetic_log_Taylor_remainder_upper_diagnostic=log_taylor_error,
        direct_product_log_float_recheck_errors=direct_errors,
        sample_power_sums=sums,
        statistical_scope='iid phase sample standard errors, not certified confidence intervals. Common random phases across all L; nested mesh difference is not certified continuous total-variation error. Positive variation uses empirical mean profile; nonlinear estimation bias is not removed.',
        deterministic_profile_scope='Degree10 log Taylor remainder has an analytic bound; all profile/trigonometric arithmetic remains floating, not outward-rounded. Rigorous positive rises come only from separate rational moment certificates.')


def main():
    start=time.perf_counter();plo,phi=pi_interval();assert 3<plo<phi<4
    registration=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        fixed_lambda=str(LAMBDA),epsilon=str(EPS),rounds=[dict(n=n,M=m,samples=s,nodes=q,seeds=[202610071101+n,202610071901+n]) for n,m,s,q in [(16,4,8192,129),(64,8,16384,257),(256,16,32768,513)]],
        input_box_side='D=10^6 n^2',source='same full positive tensor cosine input in every L; no signed independent layers',
        certificate='Machin alternating rational pi interval + exact positive-source 2/3/4 moments + bounded-entropy Taylor lower + full L1 boundary mass error',
        numeric_scope='Exploratory profiles and nested mesh diagnostics only; no fit or general upper-bound claim.')
    reg=HERE/(PREFIX+'_registration.json')
    with reg.open('x') as f:json.dump(registration,f,ensure_ascii=False,indent=2)
    rounds=[]
    for rnd,item in enumerate(registration['rounds'],1):
        n,M=item['n'],item['M'];D=10**6*n*n
        cert,total=certify_peaks(n,M,D,plo,phi)
        original=input_validation(n,M,D)
        reps=[];arrays={}
        for j,seed in enumerate(item['seeds']):
            prof=profile(n,M,item['samples'],item['nodes'],seed)
            arrays['power_sums_rep'+str(j)]=prof.pop('sample_power_sums')
            reps.append(prof)
        ap=HERE/(PREFIX+f'_round{rnd}_power_sums.npz')
        with ap.open('xb') as f:np.savez_compressed(f,**arrays)
        core=dict(round=rnd,original_input=original,certified_peak_rises=cert,
            certified_sum_of_actual_positive_rises_per_W_lower=compact(total),
            uniform_general_lower_per_W=str(F(M,589824)),replicates=reps,
            power_sums_file=str(ap.resolve()),power_sums_sha256=sha(ap),
            input_sha256=hashlib.sha256(json.dumps(dict(n=n,M=M,D=D,epsilon=str(EPS),lambda_=str(LAMBDA)),sort_keys=True).encode()).hexdigest())
        rounds.append(core)
        print('ROUND',rnd,'n',n,'CERT_RISE_LOWER',float(total),'NUMERIC_VPLUS',[r['empirical_positive_variation'] for r in reps],flush=True)
    data=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),rounds=rounds,
        pi_interval_exact=[str(plo),str(phi)],rational_certificate_scope='Outward intervals stored by exact integer floor/ceil; underlying moment signs compared with Fraction. No floating value decides certificate validity.',
        main_scope='A necessary Omega(sqrt(n)) positive-variation lower bound on complete positive L1 inputs. No general O(sqrt(n)polylog) upper bound, FIRST/geom qualification or dimension fit is proved.',
        endpoints='u_2a is a 2^-n average of the 2^n translates of u_a at a*epsilon/2, so convex entropy endpoints decrease for arbitrary complete source; periodic fixture endpoints are equal, finite L1 source obeys the general Jensen inequality.',
        script_sha256=sha(Path(__file__)),registration_sha256=sha(reg),runtime_seconds=time.perf_counter()-start)
    out=HERE/(PREFIX+'_results.json')
    with out.open('x') as f:json.dump(data,f,ensure_ascii=False,indent=2)
    print('OUTPUT',out,'RUNTIME',data['runtime_seconds'],flush=True)


if __name__=='__main__':main()

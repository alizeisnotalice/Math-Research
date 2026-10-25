"""Certified alpha=3/2 body Green lower bound from tilted Gaussian mixing."""
from pathlib import Path
from fractions import Fraction
from itertools import product
import hashlib, json, sys, time
from flint import arb, ctx

ROOT=Path(__file__).resolve().parent

def main():
    began=time.monotonic()
    ctx.prec=256
    directions=[(1,)+v for v in product((-1,1),repeat=4)]
    frame=[[sum(v[i]*v[j] for v in directions) for j in range(5)] for i in range(5)]
    assert frame==[[16*int(i==j) for j in range(5)] for i in range(5)]
    beta=arb(3)/4
    u=arb(1)/10
    a=u**beta
    mean=beta*u**(beta-1)
    inverse=a.exp()*a.gamma_upper(1/beta)/beta
    entropy=16*(1-beta)*a
    prefactor=(4*arb.pi())**(-arb(5)/2)*(16*mean)**(-arb(5)/2)*(-entropy).exp()
    exponent=5*inverse/64
    # Integral_0^infty t^-10/3 exp(-exponent*t^-4/3) dt.
    time_integral=(arb(3)/4)*(arb(7)/4).gamma()*exponent**(-arb(7)/4)
    lower=prefactor*time_integral
    rational_target=arb(1)/65000
    assert lower>rational_target
    assert lower<arb(1)/10000
    out={
        'status':'certified_global_green_lower_bound_alpha_3_over_2',
        'model':'16 unnormalized body directions; independent standard symmetric alpha-stable coordinates',
        'alpha':'3/2','subordinator_beta':'3/4','dimension':5,
        'domain':'all x in [-1,1]^5, including all boundary strata and the origin (Green may be infinite)',
        'method':'positive conditional Gaussian density; exponential subordinator tilt; entropy Jensen; tight-frame inverse bound',
        'fixed_dimensionless_tilt_u':'1/10','lambda_at_time_t':'u*t^(-4/3)',
        'tight_frame_VV_transpose':frame,
        'constants_at_t_one':{
            'u_power_beta':str(a),'tilted_mean_S':str(mean),
            'tilted_mean_inverse_S':str(inverse),'relative_entropy_16_coordinates':str(entropy),
            'uniform_heat_lower_prefactor':str(prefactor),
            'uniform_heat_lower_exponent_a':str(exponent),
            'full_time_integral_factor':str(time_integral)},
        'heat_bound':'p_t(x)>=prefactor*t^(-10/3)*exp(-a*t^(-4/3)) for every x in cube',
        'green_bound_arb':str(lower),
        'certified_rational_g':'1/65000',
        'conditional_target_g':'1/10000',
        'conditional_target_achieved':False,
        'new_C5_upper_bound':False,
        'runtime':{'python':sys.executable,'precision_bits':ctx.prec},
        'seconds':time.monotonic()-began,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (ROOT/'tilted_gaussian_bound.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()

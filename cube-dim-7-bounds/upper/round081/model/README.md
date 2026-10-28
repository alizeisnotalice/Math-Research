# R081 one-variable model: analytic optimum

This is a model diagnostic only. It does not establish a C_n bound or transfer the original s=3/2 geometry, typicality events, local face bounds, weighted cap lemma, finite guards, or operator argument to variable s.

For s>1 put a=(s+1)/(2s), mu=(s-1)+(s+1)log a, and cmodel=log[s/(1+2mu)]. The coefficient audit is as follows. In the limiting full-target Fubini profile, the first moment per free coordinate is (2R/rho)(1-a+a log a). With R=sq, rho=tq and n/t=x, its n-coordinate coefficient relative to x is 2s(1-a+a log a)=(s-1)+(s+1)log a. This confirms both terms in the assumed ideal profile; it is not a proof that this profile controls all the actual finite-s geometry.

For chi>mu and sigma/2>chi(1-sigma), solving the latter gives sigma>2chi/(1+2chi), hence -log(1-sigma)>log(1+2chi). In the ideal limit chi downward to mu, sigma downward to 2mu/(1+2mu), and target mass upward to log s, their difference approaches log s-log(1+2mu). Strict inequalities mean this is a limiting model value rather than a finite feasible allocation.

Differentiating mu gives mu'=1-1/s+log a. Consequently

    cmodel'(s)=[1+2 log a]/[s(1+2mu)].

Here mu(1)=0 and mu'>0 for s>1, because log(2s/(s+1))<(s-1)/(s+1)<(s-1)/s, so the denominator is positive. The numerator decreases strictly from 1 to 1-2log 2<0. The unique global maximizer is therefore

    s*=1/(2 exp(-1/2)-1) = 4.6934844987231898...
    cmodel(s*) = 0.5553395620244103...

At this point mu=(s*-3)/2. Decimal output is 80-digit arithmetic, not an interval certificate for the transcendental optimum. Independent rational enclosures of the model and derivative numerator at five fixed rational points are in `result.json`; their displayed midpoint values are:

| s | cmodel(s) | derivative sign |
|---|---:|---|
| 4 | 0.553401017639 | positive |
| 9/2 | 0.555219319500 | positive |
| 5 | 0.555093648057 | negative |
| 11/2 | 0.553926173921 | negative |
| 6 | 0.552212306123 | negative |

Exact logarithm enclosures use 64 positive terms of the atanh series and a rational geometric tail bound; increasing-log monotonicity handles the interval for 1+2mu. The registered single batch used 20 log-enclosure calls and actual process PID 28078, CPU 1.42 s. Reproduce using `python3 model.py`; this overwrites the PID-bearing result. No MATLAB, scan beyond the registered points, commit or push was used.

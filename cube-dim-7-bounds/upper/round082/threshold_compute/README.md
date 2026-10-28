# R082: finite arithmetic threshold check for R081 s=4

Scope: the single base L0=2^24 was registered before execution. Exact arithmetic checks the explicit envelopes used in the R081 proposed proof, conditioned on its mathematical assertions. It does not audit the five-event construction, weighted-cap lemma, same-history operator bridge or whole proof, and does not alter the original proof or assert that its unquantified threshold was previously 2^24. The frozen proof SHA-256 is 1081f7052cd0fdea55186a71d062dd09137d4c9f6eaf31fe82c68871db64516d.

Set N=L^16/16! <= e^L=n, x=2L/5, h=L^-2, and all e=1/10000. The one actually executed Fraction run PID 28988 evaluates 54 gates, all PASS, with 64-term atanh log enclosures and positive S16 exponential prefixes. No nearby threshold or parameter search was run. Full rational values, names and booleans are in `result.json`; displays below are approximations:

| Guard | Certified upper/lower display | Required |
|---|---:|---:|
| N | 1.8832e102 lower | n>=N |
| source rho | 7.0155e94 lower | >=12 and >4rhi |
| baseline log-sum mean | 1.3775e-86 upper | <=e/2 |
| baseline Chebyshev failure | 2.2040e-77 upper | <=e/8 |
| each count failure | 2.1241e-94 upper | <=e/8 |
| first-profile failure | 3.3985e-97 upper | <=e/8 |
| source-face Chebyshev failure | 0.00035436951 upper | <=0.001 |
| weighted-cap payment | 5.9751e-50 upper | <=e |
| tau payment | 3.6991e-57 upper | <=e |
| root remainder | 1.0106e-90 upper | <1/200 |
| complete deficit budget | 1.366988355 lower | >13/10 |
| target mass log lower | 1.376231500 | (1-e) times it >137/100 |

The gates additionally cover the prescribed active geometry, epsilon qmax<=h/4, 2h<delta, log8/5<1/2, 80/L<=e, individual J, 80h<=1/100, the pin and five-event error payments, the target cap fraction and probability, count calibration, Chebyshev gap, radial width, B positivity, log(5/2)<23/25, and gamma-price-c=1/200. The two finite-history arithmetic guards check the size inequalities and sufficient displacement margin; they are not a verification of the history argument itself. For h=L^-2 the tau envelope uses rhi L^(5/2)/S16(L/50), with sqrt(L0)=4096 exactly.

For extending these *arithmetic envelopes* to larger L, N is proportional to L^16 and x proportional to L: all upper bounds expressed as positive sums of L^a/N with a<=2 decrease, while rho=N/(4x) increases. The individual logarithmic and width errors decrease; the count, first-profile and baseline failure bounds do likewise. The face Chebyshev bound 2x/(a x-1/2)^2 decreases when a x>1/2. The weighted-cap exponent (nu-1/N)eta x/(1+e) increases, hence its S16 reciprocal decreases. For the tau expression, yS16'(y)-(5/2)S16(y)>0 for y>=4: terms of order >=3 have positive coefficients and the cubic contribution dominates the negative degrees 0,1,2. The analogous yS16'(y)-(1/2)S16(y)>0 for y>=2 gives decrease of the root sqrt(L)/S16 term. The remaining positive-prefix and N reciprocals decrease. These are explicit tail comparisons, conditional on the proof's bounds and on n=e^L>=N.

Reproduce with `python3 verify.py`; it writes a fresh PID-bearing result. Python standard library only, no MATLAB, no commit or push, no modification of R081 source or other work.

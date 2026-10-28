# R084 single-threshold arithmetic for R083 s=9/2

Scope: fixed L0=2^50, N=L^16/16! <= n=e^L, x=(387/1000)L, h=L^-2, and all errors e=1/100000. The proof's SHA-256 is 6d2dff6fedf7cde05d4823830bca674235bd84be4608fc35ea129aa44d649c0e. Registration preceded one actual Python Fraction run (PID 29861) yielding 55 base gates PASS. These are arithmetic checks conditional on the cited proof's bounds; they do not independently audit the geometry, event construction, source-face laws, weighted-cap lemma, same-history bridge, or C_n theorem. No search or other threshold was run.

| Guard | Approximate certified envelope | Required |
|---|---:|---:|
| N | 3.1869e227 lower | n>=N |
| source rho | 1.6253e212 lower | >4rhi |
| baseline mean | 5.9111e-196 upper | <=e/2 |
| baseline failure | 9.4578e-185 upper | <=e/8 |
| each count failure | 1.2552e-217 upper | <=e/8 |
| first-profile failure | 2.0105e-217 upper | <=e/8 |
| source-face Chebyshev | 5.10008e-6 upper | <=e |
| weighted source-cap | 4.5914e-142 upper | <=e |
| tau | 2.1717e-130 upper | <=e |
| full root remainder | 3.4715e-212 upper | <1/1000 |
| (1-e)Gamma | 1.50280462097... | >1.501 |

Exact rational bounds and all gate names are in `result.json`. The calculation uses the R083 constants 7/rho; mean 49n/rho²+14/rho and variance 196n/rho²; pin 81/(16n); profile bad bound 64072265625/n; whole-band log shift 180xh=(3483/50)/L; profile transport 80h<=1/100000; source face mean .6+.79147x and variance <=2x with positive gap (3/100000)x-3/5; nu=22221/100000 and eta=1/19998; tau rate 379/2000000; and the exact source and target radius checks 4rhi>6>1+rhi. The numerator of the target pin bound is derived from rho>=n/(sx). At L0, sqrt L0=2^25 exactly. The scalar allowance is gamma-price-c=1/1000, and the complete root cross-term is retained.

For arithmetic tail extension, N is proportional to L^16 and x to L: the baseline, count and profile upper envelopes decrease, while the source-radius lower bound increases. The whole-matching shifts and radial errors decrease. The source-face 2x/(a x-p)^2 decreases when a x>p; the lower bound on the weighted exponent (nu-1/N)eta x/(1+e) increases. The tau prefix tail uses the *correct* four positive/negative low-degree terms:

    y S16'(y)-(5/2)S16(y)
      >= -5/2 - 3y/2 - y²/2 + y³/12 + y⁴/16.

The exact quadratic coefficient in yS16'-(5/2)S16 is -1/4; the displayed -1/2 is a deliberately weaker lower bound, not the exact coefficient. This lower polynomial equals 29/6 at y=4. Its first and second derivatives are positive for y>=4 (at 4 they are 29/2 and 13, respectively), so it stays positive. All omitted degree >=5 terms have positive coefficients. Therefore L^(5/2)/S16(cT L) decreases after the certified domain. Similarly, yS16'(y)-(1/2)S16(y)>0 for y>=2 ensures decrease of sqrt(L)/S16(sigma L/2); the other scalar prefix and N reciprocals decrease. These comparisons certify propagation of the *displayed arithmetic envelopes* if the corresponding mathematical estimates hold, not the proof of those estimates.

The code uses 64-term atanh log intervals and the positive S16 exponential prefix, all in exact Fraction arithmetic, with decimal values for display only. Run `python3 verify.py` to reproduce; the PID-bearing result is overwritten. No shared MATLAB, original-source modification, commit or push was used.

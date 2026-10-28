# R091 separately registered exact candidate

Scope: conditional fixed-tuple arithmetic, pending independent acceptance of frozen R090 proof SHA-256 `12605f79fc17c83516fd4423c69fd58b1340a2a780588490591795ee3ec59631`. It is not an independent proof of that geometry, source laws, weighted cap, original eligibility history or C_n theorem.

The candidate was created from the full `parameter_decimals` fields of worker 9 in the locked `summary.json`: for each of s, chi, sigma, delta, qmax, gamma, m, ret, e and m_face, Python's shortest round-trip `str(float)` was quantized to eight decimal places with Decimal ROUND_HALF_EVEN and converted to an exact Fraction. This explicitly defines a NEW rational tuple; it is not a reconstruction of a prior exact certificate. The original decimal values, complete exact fractions and process are preserved in `exact_registration.json`. Separately fixed are c=11/10, trim=1/10000000 and conservative price 6581/2500=2.6324.

The pre-registered plan was first L0=2^33 and, only on failure, one attempt at 2^34. At L0=2^33, actual run PID 33459 returned PASS for all 58 exact Fraction gates, so no 2^34 evaluation was performed. All gates, exact rational margins and machine-readable inputs are in `exact_result_L8589934592.json`.

| Quantity | Fraction-enclosure display | Required |
|---|---:|---:|
| d_mu | 0.1300827432887... | >0 |
| source-face failure | <0.000347214 | <=1-Ret=0.00060738 |
| target cap failure | <1.327e-19 | <=e/8=1.25e-6 |
| weighted cap/k | <4.706e-17 | <=e=1e-5 |
| tau/k | <1.434e-72 | <=e |
| complete H_tau remainder/n | <1.195e-80 | <0.00079131 |
| (1-e)Gamma | 3.73956140127... | >gamma=3.73319131 |

The scalar allowance gamma-price-c=79131/100000000 exactly, with log(1/(1-sigma)) below the separately fixed price. The code uses the frozen arbitrary-s constants D_s=4s²+6s+1, 2s/rho, s/rho, 5s/rho, 4s² and 16s²; source variance 2s²x and pin 1; rho>4(s+1) and target radius <s+1. It retains mu and d_mu as rational upper and lower logarithm enclosures, not decimal substitutes. It pays the complete non-square-root scalar cost 1/N+1/S24(delta L)+1/S24(sigma L)+L/S24((sigma-delta)L), with N=L^16/16!, 96-term exact atanh log bounds and S24 positive exponential prefixes.

## All-tail envelopes, conditional on R090 estimates

For L>=2^33, x=(1-sigma)L and N=L^16/16!; every tuple constant remains fixed. The following are consequences of the displayed arithmetic envelopes, not independent proofs of the R090 mathematical estimates.

* **Geometry and local bounds:** h=L^-2, x/N and x²/N decrease; N/(sx) grows. The base checks persist for 2h<trim, R>1, g>0, epsilon<min(1,h), epsilon qmax<=h/s (the ratio is proportional to L^-31), beta<2, kappa>1, t>=2rhi, rho>4(s+1), and individual J<=e. Thus active bands remain in the local domain and target radii remain <s+1.
* **Five full-target failure budgets:** the baseline mean <=4s^4x²/N+4s²x/N stays <=e/2, and Chebyshev 4(16s^4x²/N)/e² stays <=e/8. Each of the two count events costs <=4/(e²N)<=e/8. For the target cap, t>=2rhi gives fraction <=1-1/(4s); N>=10001 and (1-v)^(n-1)<=exp[-10000/(4s)]<=1/S24(10000/(4s)) pay <=e/8. The actual first-profile target pin is paid by s²/N<=d_mu/2, and its failure by 4s^4/(d_mu²N)<=e/8. No independence between these five events is assumed.
* **Whole matching and retention:** 2D_s xh and 10s h decrease, so the log-product and first-profile transports stay within e and d_mu. The count check 2h+1/N<=3e/4 and radial upper 1/[1-(h/rlo)²]<=1+e improve; the radial lower side is >=1. The R090 inverse-face mean and variance give failure <=2s²x/[(chi-m_face)x-1]^2, which decreases for positive gap. Its fixed Ret remains paid, provided that inverse-face law has been proved.
* **Weighted cap:** nu and eta are fixed positive, N grows and (nu-1/N)eta x/(1+e) increases. Thus M/S24 of this exponent decreases and stays <=e; the underlying weighted-cap lemma is an external proof dependency.
* **Tau, scalar cost and root positivity:** for y>=4, yS24'(y)-2S24(y)>=-2-y+y³/6>0, so L²/S24((delta-chi(1-sigma))L) decreases. The complete scalar remainder includes both exponential terms and the L cross term: their inverse prefixes decrease, and yS24'(y)-S24(y)>=-1+y²/2>0 for y>=2 gives monotonicity of L/S24((sigma-delta)L). Positive-prefix S24(L/2)>=x remains true by its quadratic term; hence the cited B>=nL/4>gamma n remains positive.
* **Original finite history:** for a_alg>=8n^4 and n>=N>=2, x>=2 implies a_alg>2t; beta<2 gives n^-3-(beta+2)/(a_alg+1)>1/n³-4/(8n^4)>0. The target item radius is <s+1, whereas source rho>4(s+1), uniformly on the tail. The actual row, envelope domination, energy identity and finite-cover theorem still require the cited proof, not just these inequalities.

To reproduce the original single threshold, `python3 exact_verify.py 8589934592`; it overwrites the PID-bearing receipt. No prior search receipt, paper, other round, commit or push was changed.

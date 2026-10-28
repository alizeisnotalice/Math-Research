# R087 fixed tuple and L0 candidate

This is a conditional arithmetic certificate for the frozen R085 parameter-family proof, SHA-256 `ded160a902cdd543a7be9ebbb9c1f4d076c33f866d8aa96e9481f92d69bb41fa`. It does not independently prove or accept its local profiles, probability laws, weighted cap, history or operator transitions, and does not assert a new C_n theorem. Both registrations preceded their respective computations. No grid, second threshold, MATLAB, commit or push was used.

Primary parameters s=47/10, chi=42431/50000, sigma=62927/100000, qmax=99999/100000, gamma=619/400, c=111/200 pass eight strict Fraction/96-term rational-log gates (PID 31183): chi-mu>8.032282e-6; sigma/2-chi(1-sigma)=0.0000261074; log(s qmax)-gamma>0.00005250866; gamma+log(1-sigma)-c>0.0002187558. These are interval lower bounds, shown with conservative decimal truncation.

The separately registered fixed choices are trim=e=1/10000000, Ret=999999/1000000, m=154753/100000, m_face=424309/500000, pin=3/5, price upper=99229/100000, and L0=2^80. Here d_mu=(chi-mu)/4 is not replaced by a floating decimal: exact log intervals bound it above and below. The base Fraction run PID 31252 returns 56 PASS, with exact bounds in `threshold_result.json` and numerical displays:

| Gate | Certified display | Requirement |
|---|---:|---:|
| d_mu | 2.0080705e-6 lower | >0 |
| transported profile upper/(1-e) | <0.848618 | m_face=0.848618 |
| source-face failure | <5.5781e-12 | <=1-Ret=1e-6 |
| weighted cap | <5.8657e-270 | <=e |
| tau | <1.6297e-238 | <=e |
| full scalar remainder | <2.3964e-352 | <0.00021 |
| (1-e)Gamma | 1.5475265954... | >gamma=1.5475 |

The mass and price are checked by strict rational atanh-log enclosures: (1-e)log(rhi/rlo)>m, and -log(1-sigma)<0.99229. The complete deficit uses P0, Ulo/Uhi, M and b exactly as R085; gamma-price-c=21/100000 exactly. The source face is paid by the actual pre-cap facewise Chebyshev bound 10x/[(chi-m_face)x-3/5]^2, not a normal approximation. Every target good-event probability, including 400/(d_mu² n), is separately paid <=e/8. Target cap uses (19/20)^399 and N>=400; matching uses 240/L<=e, 120/L²<=d_mu, the original count guard and 9/rho individual bound. New source-cap eta=(1-qmax)/(2qmax), nu=1/s-(1+1/s)e, and the finite-N correction (nu-1/N) are all retained. Tau uses rhi L^(5/2)/S16(cT L); root includes 1/N, A and its positive square-root cross term. Since sqrt L0=2^40, all are Fractions.

For the full-tail *arithmetic envelope* propagation, N=L^16/16! and x=(1-sigma)L imply x/N and x²/N decrease, while rho=N/(sx) increases. Trim, epsilon, individual J, count, target cap and profile-transport guards improve; the five failure estimates decrease. For face gap ax-p>0, 10x/(ax-p)^2 decreases. (nu-1/N)eta x/(1+e) increases, so M/S16 of it decreases. With y=cT L>=4, yS16'(y)-(5/2)S16(y) exceeds the weaker polynomial -5/2-3y/2-y²/2+y³/12+y⁴/16, positive and increasing on y>=4. The exact quadratic coefficient is -1/4; -1/2 is a conservative lower bound. The root sqrt(L)/S16 term decreases because yS16'(y)-(1/2)S16(y)>0 for y>=2, and its other two terms decrease. Scalar root positivity and finite-history size/displacement gates then persist. These arguments still assume the R085 estimates to which the envelopes apply.

With the archived primary receipt present, reproduce the threshold check using `python3 threshold.py` (which overwrites its PID-bearing receipt). Re-running `primary.py` overwrites the primary receipt and changes its PID/hash; preserve it or run it in a copy before invoking the threshold script. No original proof file was altered.

## Full-tail checklist

For every L>=L0 set N(L)=L^16/16!, x=(1-sigma)L, h=L^-2 and n=e^L>=N(L). All constants e, trim, d_mu, face, Ret, m, gamma and sigma stay fixed. The following checks are implications for the displayed envelopes, not independent proof of the R085 estimates:

* **Geometry and local domain.** The base gates give 2h<trim, rlo-2h>1, sqmin-h>0, (s-1)qmin-h>0, x/N<1, 1/(xN²)<1 and qmax/(xN²)<=h/s. Their left sides or relative-error ratios decrease thereafter: qmax/(xN²h) is proportional to L^-31. The bounds N/x>2rhi, N/(sx)>20 and >4rhi strengthen. Thus beta<2, kappa>1, t>2rhi, source separation, and the stated active relative-gap guard persist. Since a0<qmax<1<rlo, target pins stay outside the count threshold. The bounds 9/rho<=e and <=1/2, 240/L<=e and 120/L²<=d_mu improve.
* **Five full-target failures.** For the baseline event, absolute mean <=81s²x²/N+18sx/N remains <=e/2 and the failure bound 4(324s²x²/N)/e² remains <=e/8. The two distinct count events each have failure <=4/(e²N)<=e/8. The positive-face target-cap fraction is <=19/20 because t>=2rhi, and N>=400 gives failure <=(19/20)^399<=e/8; negative faces have no such failure under the cited geometry. The target-pin first-profile contribution satisfies 10/N<=d_mu/2 and the first-profile failure <=400/(d_mu²N)<=e/8. Thus all five events are separately paid, without assuming mutual independence.
* **Whole matching and radial calibration.** The bound 2h+1/N<=3e/4 decreases, so the cited endpoint-to-whole-interval count estimate persists. The radial upper envelope 1/[1-(h/rlo)²] decreases and stays <=1+e; the lower side is at least 1 and hence >=1-e. The independent-face mean bound uses (mu+2d_mu)/(1-e)<m_face; its variance is <=s x/(1-e)<10x. Once ax-p>0, a=chi-m_face, the actual upper envelope 10x/(ax-p)^2 decreases, so Ret stays paid.
* **Weighted source cap.** nu and eta are fixed positive; N grows and (nu-1/N)x increases. Consequently S16 at its lower exponent increases and M/S16 decreases. The estimated cap remains <=e, subject to the cited source-face lemma.
* **Global tau and scalar root.** S16(L/2)>=x persists because the positive quadratic prefix dominates linear x; consequently the cited B>=nL/4>gamma n and tau is positive. For y=cT L>=4, the polynomial lower bound on yS16'(y)-(5/2)S16(y) above proves rhi L^(5/2)/S16(cT L) decreases. For the root remainder, 1/N and 1/S16(sigma L) decrease, while sqrt(L)/S16(sigma L/2) decreases by the y>=2 derivative bound. All remain below the fixed slack 21/100000.
* **Finite history arithmetic.** For any a_alg>=8n^4, x>=2 and n>=N>=400 imply a_alg>2t=2n/x and ccap>1. The target radius is <1+rhi<6 while a source radius is >=rho>20>6. Since beta<2 and n>=2, zeta-(beta+2)/(a_alg+1)>1/n³-4/(8n⁴)>0. These arithmetic inequalities persist; the domination and energy identity are separate mathematical claims in the frozen proof.

The full list of 56 base gate booleans is the machine-readable checklist. All source probabilities, especially the retained source-face probability, refer to the R085 ambient inverse law; the actual increment measure is not assigned such a law.

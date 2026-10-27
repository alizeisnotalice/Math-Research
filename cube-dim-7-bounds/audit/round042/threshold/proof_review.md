# R042 explicit threshold: independent second audit

2026-10-28. **ACCEPT at s=(log n)^(1/3)>=2^36**, equivalently
log n>=2^108. The final locked proof quantifies all eventual conditions
of the accepted second-moment operator theorem and establishes

```
C_n<=F(A_n(t),B_n(t)-2^(-31)n/s)
   <n log n-n log log n-2^(-32)n/s,
t=n/(log n-(log n)^(1/3)).
```

This is a sufficient explicit starting point, not a claim of the smallest
one. It is not obtained by extrapolating a finite numerical batch.

## 1. Final source and frozen dependencies

The source read in full was
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/explicit_threshold/THRESHOLD_PROOF.md`.
Its final checked SHA256 is
`a73bae47dbef7b56b95dcc28bf5dc368828265e91c0de49d03a3dcf518e8c0c5`.
The initially submitted hash was
`76f0b794ba4fc3e7319bb87bf955a2627e18ffc5ed42bc62836bcbc975e7b19f`.
During review the author clarified only the deleted-coordinate upper-log
argument in Section 3.3; the parameters, estimates, script and threshold
were unchanged. The final text explicitly avoids the false inference
-log(1-v_j)<=v_j and uses log product_except<=-T_w+v_j instead.

The frozen full operator proof remains
`722764654cd15a8bf7fb18a799c298b54ced2f0188f5554ccae173240ddf90a2`;
the frozen explicit moment proof remains
`425a383952f341f42ff0dee9651924b00fff61e958c2957f355e8cc9fdc72d00`.
Their paths and actual hashes are retained in `exact_receipt.json`.
Neither frozen file was modified. The full operator proof was previously
accepted only sufficiently far out; this review checks the additional
explicit quantification, rather than accepting a threshold in advance.

Only `audit/round042/threshold/` is written. There was no production
replay, Arb work, numerical search, n=7 calculation, MATLAB or Git operation.

## 2. Basic geometry, matching and the actual kappa perturbation

Let S=2^36, s>=S, L=s^3, d=1/s and h=1/(100s^2). Already for s>=8,
x=s^3-s lies in [L/2,L], x>1, d<=1/8, q>=1/2 and r<2.
The exponential positive-term bound n=e^L>=L^2/2=s^6/2 yields
t>=s^3/2 and rho>=s^3/4>4. Hence target/source supports are disjoint,
rho>r, 2d<r and 2d<rho. The fixed-face moments and binomial formulas
therefore apply on their genuine domains.

The real finite perturbation satisfies

```
zeta*t=1/(x n^2)<=8/s^15<=h/2,
```

since the final comparison requires only s^13>=1600. Thus
h+zeta*t<=3h/2<d/2, and kappa>1 follows from zeta*t<1.
Also h/d=1/(100s)<1/2. The matching center remains in I at both
endpoints, so the shortest partner interval has length h. The source
logarithmic length is between h and 4h and the target one at most 2h.
These statements hold uniformly over all matching radii, not only
at interior radial points.

## 3. Target mass and row products are explicitly uniform

The fixed-face target mean is the exact m_w. Its Chebyshev bad rate
is at most 20000/(3s^3), which is below 1/4 for s>=S and decreases.
The inequality nd>=s^5/2>2 gives the K_w Markov bound at most 1/4.
Thus at least half every target face is good. The total target annulus
has mass at least nd/4. The expression

```
1/(1+3/(4s))-2/s^3
```

increases with s and is above 1/2 at s=4, so target corner lengths
are below 3/4. Because n>=6, the captured mass is at most nd/16.
Therefore nd/16<=lambda(A)<=nd, with no remaining eventual premise.

For row products, the center's error from s/4 is at most
(1+d)/[4(1-d)]<=9/28<1/3 for s>=8. The tolerance and shift error is
bounded, in the correct direction, by

```
1/100+18Ldh+4d/L=19/100+4/s^4<1/5.
```

The bound used q>=1/2, h+zeta*t<=3h/2, x<=L and t>=L/2;
it retains the nonvanishing constant compensation at this scale.
Consequently |T_w-s/4|<8/15 and T_w<=s/2. Every coordinate loss
is at most V=4/s^4<1/2. Its quadratic log correction is at most
[V/(1-V)]T_w<=4/s^3.

For the upper log of a product deleting coordinate j, use directly

```
log product_{i!=j}(1-v_i)<=-sum_{i!=j}v_i=-T_w+v_j.
```

The lower log is no smaller than the full product's lower bound.
Thus adding V+4/s^3 is sufficient, and
4/s^4+4/s^3<1/10 at s=8 and thereafter. The total log error is
less than 1/3+1/5+1/10<2, proving every required full and deleted
row product e^(+/-2)p uniformly. No assertion that the actual log
increment on deletion is at most v_j is needed.

## 4. Source corner deletion has an all-tail power bound

The explicit radius guards give corner exponent at least
3Ld/64>=s^2/32. Its complete source-band mass divided by M0 is
bounded by the safe larger envelope

```
100s exp(-s^2/32+s/4).
```

For s>=16 its exponent is at most -s^2/64. The single positive
Taylor term exp(u)>=u^2/2 then gives

```
Corner/M0<=819200/s^3<1/36,  s>=S.
```

The last power is decreasing. Since e^2<9, this removal is strictly
smaller than e^(-2)/4, leaving raw row lower bound e^(-2)M0/4.
Together with lambda(A)>=nd/16, this gives raw energy above
(1/1024)ndM0. The complete raw row upper bound remains 36M0.
This controls only the cap complement; no unavailable global
source-good-statistics mass argument has reappeared.

## 5. Moment reuse and energy constants do not hide an asymptotic step

The independently accepted moment proof is explicit for L>=64 with
C=48 and the same Y, kappa, h and matching interval. It includes the
pinned maximal source coordinate and the positive cubic moment term.
Its domain contains all s>=S. It can therefore support the frozen
main proof's weaker constant C2=144 without relying on that main
proof's formerly eventual moment estimates. This is a noncircular
use of an independent one-dimensional sublemma.

The same cutoff Y>HM0, H=2^20, removes energy at most
(144/H)ndM0. The exact comparison H>=4*144/cE and H>=36 gives
retained energy at least (3cE/4)ndM0 and both degrees at most HM0.
The cutoff removes energy, not source mass. Its measurable definition
requires neither a selector nor an assumption on the actual mu.

## 6. Tau, finite a and the original all-scale operator

The old scalar estimates with x in [L/2,L] give B0>8n for L>=64,
and delta=2^(-31)nd<n, so the root budget is positive. Also
A0<=ne^(-s), B0<=nL. Thus

```
tau/M0<=100s^(7/2)e^(-s/4)
       <=100*4^8*8!/s^(9/2)
       <=100*4^8*8!/s^4<2^38/s^4
       <=2^(-106)<cE/4.
```

This uses only Taylor's eighth positive term and decreasing powers,
not a numerical exponential sample. Subtracting tau*lambda(A)
uses the correct upper bound lambda(A)<=nd. The deficit is therefore
cE/(2H)nd=2^(-31)nd at this explicit threshold.

For every real a>=8n^4, a>2t, c>1 and
alpha*beta-kappa-2/(a+1)=zeta-(beta+2)/(a+1)>0.
Target actual radii are below 3 and source actual radii above 4,
so the same selected target's eligibility pays the relevant earlier
sources. The exact ambient alpha^(-n) scaling on both supports
gives the same finite-a deficit; energy scales quadratically while
degrees and target mass scale linearly. No TV convergence is used.

All source/target measures are still obtained from actual increments
of one history. The original witness cubes remain aQ_i with arbitrary
positive common radii and arbitrary L1 input. The existing compact
strict-level-set covering and inner-regularity bridge therefore applies.
For each fixed n, only after uniform finite-a bounds are established
is a taken to infinity in their scalar constants. No new input class,
selector limit, or exchange with n-to-infinity is required.

## 7. The new scalar cost is paid at the explicit starting point

The exact relation x=L(1-d^2) and d^2<=1/2 give
B0<=T+2nd^2+1. Using A0<=ne^(-s), B0<=nL and Taylor-eight bounds,
the total positive remainder after subtracting delta, normalized by
nd=n/s, is at most

```
R(s)=2/s+2/s^5+8!/s^7+(2^9*8!)/s^5.
```

The four terms pay, respectively, the log-x change, constant one,
A0, and the square-root cross term. The last exponential initially
gives an even smaller s^(-11/2) term; replacing it by s^(-5) is
an upper bound for s>=1. All terms decrease, and exact rational
comparison gives R(2^36)<2^(-32). Hence

```
F(A0,B0-delta)<T-delta+2^(-32)nd=T-2^(-32)nd.
```

The scalar calculation uses the current x=L-s throughout. No W_n
bound from the different x=L-3log L scale is imported.

## 8. Actual independent checks and final status

`python3 check_exact.py` was actually run and returned PASS, PID 74948.
It used only stdlib Fraction and integer powers, imported no producer
code and performed no search. Its 17 minimal exact gates independently
verify the key rational bases, Taylor coefficient, cap loss, energy
constants, tau and scalar remainder. The other simple domain statements
are proved uniformly in Sections 2-7 rather than inferred from finite
samples. The exact inputs, results and all source hashes are retained
in `exact_receipt.json`.

The source author's separate 21-check receipt was inspected, not
substituted for this independent proof. Every infinite tail is covered
by a stated monotone rational or inverse-power envelope. The full and
moment frozen hashes were checked unchanged. The final explicit theorem
is accepted, without a threshold-optimality claim or any C_n lower bound.

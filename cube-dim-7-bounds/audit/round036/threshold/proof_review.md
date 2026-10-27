# R036 explicit R035 threshold: final independent audit

2026-10-28. **ACCEPT for log L>=128, where L=log n.** The repaired
source proof quantifies every implicit R035 condition and proves

```
C_n<=F(A_n(t),B_n(t)-delta)<W_n-delta<T_n-delta/2,
delta=(n/65536)sqrt(log L/L),
t=n/(L-3log L), T_n=nL-n log L.
```

In particular C_n<T_n-(n/131072)sqrt(log L/L) for every integer n with
log log n>=128. This is a sufficient threshold, not a minimality claim.
The proof covers the entire infinite tail; it does not infer a tail theorem
from finitely many Arb boxes.

## 1. Scope, frozen source and a corrected harmless endpoint

The author source read in full is
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round036/threshold_theory/THRESHOLD_PROOF.md`.
The locked final SHA256 is
`064e3e6d697220817642250961180fcf2f747dc36446eb693344fb74f719bed1`.
The initial draft had the strict inequality e^(-m)<(3/8)^m also at m=0;
this review reported that harmless endpoint, and the author changed it to
strictness for m>0 and equality at zero. All actual exponents are positive.
No substantive repair of the threshold argument was needed beyond the
source's already explicit correction of its earlier tail interpretation.

The accepted R035 file was not changed. Its checked SHA256 remains
`f7d98eb21b89d83264053fda2c594387387facfdcbdbec7cf83a41e5dd55f0bd`.
This audit quantifies its previously accepted full-scope theorem; it does
not repeat baseline production calculations or certify an old Arb run.
Only `audit/round036/threshold/` is written. No n=7 work, parameter search,
Arb replay, shared MATLAB, or Git operation was performed.

## 2. Independent reconstruction of the envelope directions

Set z=log L>=128, d=sqrt(z)exp(-z/2), h=exp(-z/2)/(100sqrt z),
q0=1-3d/4 and x=L-3z. The inequality n=e^L>=L^2/2 and x>=3L/4 give

```
zeta*t=1/(x n^2)<=16/(3L^5)<8/L^5.
```

Thus the source's shift condition controls the real kappa perturbation,
not just a nominal bandwidth. Since h/d=1/(100z)<1/2, matching intervals
still have length at least h even at the endpoints of I.

For r=1+u,q=1-v and u,v in [d/4,3d/4],
|rq-1|<=d/2+9d^2/16<=3d/4 and rq>=q0. The resulting common Gaussian
mean error is bounded by the stated E_m. This uses the actual positive
radial intervals and has the correct denominator q0.

For rows, dividing the true tolerance/shift/mean-correction bound by
2rho=2qt yields respectively 1/100, 12/100, 48d/L^4 and 4d/L in its
numerator. This is precisely E_r=(13/100+48d/L^4+4d/L)/(2q0).
For columns, r>=1, Ld*h=1/100 and the shift estimate give
E_c=17/200+d+64d/L^4. Both are upper bounds, not asymptotic substitutions.

Every side loss is at most V=(d+h+8/L^5)/2. With
E=E_m+E_r+E_c, the actual total loss satisfies |T-z/4|<=E. The logarithmic
product lower error is at most E+[V/(1-V)](z/4+E). A deleted coordinate
can increase the upper log product by at most V. Therefore the source's

```
Perr=E+V+[V/(1-V)](z/4+E)
```

simultaneously controls full and every deleted-coordinate product on both
sides. In particular Perr<2 establishes the exact uniform products used
by the energy argument; no deleted-coordinate positive error is omitted.

## 3. Statistics, mass removal and the repaired exponential tail

For the target statistics, n/t^2=x^2/n<=2, so the bad-S probability is at
most (20000/3)d^3. Together with nd>2 the K Markov bound gives at least
half good points on each target face. The whole annulus has mass at least
nd/4. The explicit corner-face gap implies ell_t(r)<3/4. Its captured
mass is at most (nd/4)(3/4)^(n-1), which is at most nd/16 once n>=6.
Subtracting from the half-good mass gives lambda(A)>=nd/16. The hidden
n>=6 follows from n>=L^2/2 and the already implied L>=64.

For sources, q>=q0>1/2 and Ld>1 justify the fixed-face variance and
binomial estimates. The complete bad-statistics mass divided by k equals
at most the first two Bad terms

```
(2000000/3)z^(5/2)L^(-1/4)+(400/121)sqrt z L^(-1/4).
```

These bounds concern the whole source band, so they may be subtracted
from every capture row with kernel bounded by one. They are not
conditional-good-set probabilities.

The old explanation exp(-u)<=u^(-2) alone gives only an L^(-3/4)
corner envelope at this point and cannot justify the claimed L^(-5/4).
That explanation must stay withdrawn. The replacement is correct:

```
L=e^z>=z^4/24>=4096z,
Ld/32=sqrt(Lz)/32>=2z,
sqrt L=e^(z/2)>=z^2/8>=16z.
```

The two polynomial comparisons hold at z=128 and thereafter, since z^3
and z increase. Consequently

```
Corner/k<=50z L^(1/4) exp(-Ld/32)
        <=50z L^(-7/4)
        <=(25/8)L^(-5/4)<51200L^(-5/4).
```

This establishes the script's original third Bad term with a valid and
stronger argument over the *entire* stated domain.

## 4. Infinite-tail monotonicity, without finite-box extrapolation

Each domain test and envelope has the favorable direction for all z>=128:

- d,h,zd,zh and all required positive exponential decays decrease;
  q0 increases. The x/L, q0L/2 and target corner-face gap increase.
- The shift test is treated as a ratio:
  (h+8/L^5)/d=1/(100z)+8z^(-1/2)exp(-9z/2), which decreases.
  Merely subtracting two decreasing functions would not have sufficed.
- Ld and L^2d/2 increase; the target Chebyshev envelope decreases.
- Each numerator term of E_m and E_r decreases while q0 increases;
  E_c decreases. V and zV both decrease.
- Perr can be written
  E+V+(zV)/[4(1-V)]+VE/(1-V). Each positive factor has a nonincreasing
  direction, including (1-V)^(-1), so the whole Perr decreases.
- The first two Bad terms have polynomial powers 5/2 and 1/2 against
  exp(-z/4), and are decreasing above 10 and 2 respectively. The third
  term decreases. Tau follows the same half-power calculation.
- Rem splits into four positive terms proportional to z^(1/2)exp(-z/2),
  z^(-1/2)exp(-z/2), z^(3/2)exp(-3z/2) and z^(-1/2)exp(-3z/2).
  All decrease on this domain.

All denominator positivity needed here is established at the base and
preserved by these directions. Therefore the one exact base certificate
below proves every condition for every z in [128,infinity).

## 5. Independent coarse rational base certificate

This audit used a separate, coarser certificate rather than importing or
rerunning the author's script. Let p=3/8. The positive exponential series
gives e>8/3, and 1<=sqrt128<12. Exact integer arithmetic gives

```
p^8=6561/16777216<1/2500,
p^32<1/2500^4=:T0<10^(-9).
```

All larger powers p^m are also below T0. In particular at z=128,
d<12p^64<10^(-6), q0>99/100, and h<p^64/100. Replacing these and all
positive inverse-L powers in the upper envelopes by the coarser rational
bounds gives

```
E_m<1/1000, E_r<33/500, E_c<43/500,
E<153/1000, V<10^(-6), Perr<1/5.
```

For Bad, sqrt128<12 yields its first coefficient 131072000000, its
second coefficient below 40, and its third coefficient 51200. Thus

```
Bad<(131072000000+40+51200)T0<1/250<1/36,
Tau<86400T0<1/100,
Rem<131072(385+147457)T0<1/1000.
```

The exact independently computed rational upper values are

| Envelope | Coarse rational upper bound |
|---|---|
| Perr | 376002457/2457000000 |
| Bad | 3276801281/976562500000 |
| Tau | 27/12207031250 |
| Rem | 75695104/152587890625 |

All ten remaining domain inequalities follow using the same d<10^(-6),
q0>99/100 and inverse-L-power bound 10^(-9); their exact tests are retained
in `check_base.py` and `base_receipt.json`. These include the shift/d ratio,
source radius, target face gap, target Chebyshev bound and both lower
growth conditions. The hidden h/d<1/2 and n>=6 interfaces, and the two
polynomial comparisons for the repaired corner tail, are covered as well.

`python3 check_base.py` was actually run and returned PASS, auditor PID
66281. It used only stdlib Fraction and small integer powers. It did not
run an Arb box, a production batch, a numerical optimization, or any
n-dimensional integration. The source author's separate verification
receipt was inspected but not used in place of this independent base
certificate or the analytical tail proof.

## 6. Threshold closes the original energy and operator interfaces

Bad<1/36, together with e^2<9, gives source removal strictly smaller
than e^(-2)k/4. The exact matched-band lengths and product bounds therefore
give the accepted R035 row lower bound r_star=e^(-2)k/4 and ambient row/
column upper bounds 4e^2k and 2e^2k.

The positive root budget is legitimate: B>=n[L-log L-4/L]>8n, while
delta=nd/65536<n. Also A<=n/L^3 and B<=nL, so tau=sqrt(A(B-delta))<=n/L.
Tau<1 gives

```
tau/k<=100sqrt z L^(-1/4)<1/72<e^(-2)/8=r_star/(2k).
```

The ambient-degree energy lemma then supplies total deficit at least
e^(-4)nd/512>nd/65536. The strict comparison follows from e<3, hence
e^4<81<128. This uses actual increment domination, not concentration of
the unknown actual measure.

For every a>=8n^4 the explicit kappa guard and the source/target radius
gap prove the actual same-history eligibility bound. On their supports
the finite ambient measure is exactly alpha^(-n) times the limiting one;
row degrees, column degrees, energy and target mass scale in precisely
the powers required by the accepted energy lemma. Thus the same deficit
holds for every such finite-a original instance. The existing finite
cover and compact-level-set inner-regularity bridge remains valid for
general L1 inputs and all positive common cube radii.

Only afterwards is a taken to infinity in uniform scalar operator bounds.
No selected histories or selectors need converge. This yields the claimed
F(A,B-delta) bound. Rem is exactly 131072(W_n-T_n)/(nd); Rem<1 makes the
remaining positive scalar error smaller than delta/2. This closes the
stated strict negative correction at the explicit threshold.

## 7. Final status and limitations

The final locked R036 proof is accepted. The R035 frozen source remains
unchanged. No old erroneous exponential interpretation is retained as a
proof, and no assertion about a finite Arb sweep is needed for this result.
The explicit threshold quantifies an operator upper bound; it proves no
lower bound for C_n, no optimal threshold, no globally optimal covering
parameters, and no leading n log n coefficient below 1.

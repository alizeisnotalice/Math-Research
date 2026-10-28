# R070 second profile: independent LOCAL audit

**LOCAL ACCEPT.** The frozen source is
`/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/upper/round070/second_profile/PROOF.md`,
SHA256 `4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8`.
No correction is requested. This decision covers the stated H2 identities,
relative multiplicative domain, finite target means/concentration,
a.e./AC derivative and matching bridge, and the **per captured source
face** centered variance interval. It does not certify a Cantelli/BE
application, cap-conditioned independence, target count/column budgets,
actual-history/general-L1 transfer or any C_n theorem.

## Source and computation identities

The complete frozen proof was actually read. Its code and historical
receipt were identified at hashes
`0f28087a2600ac84d7dfd0f1cf8b35439097be4251e0e8ab8e7ff2e5c08857ee`
and `8d1e1d722b2a773bda7a42db4be259ea1c0d8b8dbf1b8448c5016273742be97c`.
Neither author code nor an old baseline program was executed/imported.

Two distinct small, single-process independent runs are preserved:

1. `check_independent.py`: 29 SymPy/Fraction checks PASS. This run
   occurred **before** the author's complete source was frozen/read and
   records independent formula checks only. It was not relabeled afterward
   as a source-bound run. Its log enclosure for the new H2 rate is a
   different recipe, even12/odd13 alternating log1p sums; no author's
   positive-atanh code or old baseline calculation was reused.
2. `check_source_binding.py`: after the frozen proof was read, 22 new
   scope/finite-guard checks PASS. This run binds all three final author
   hashes, reads the prior independent rational enclosure, checks the
   **actual n=2000** endpoint and source variance bounds, and checks
   the wide-h guards using weaker Taylor denominators 8 and 32.
   It does not repeat the first run's symbolic calculations or log sums.

Their original receipts, actual times/PIDs, interpreters, exact rational
values and code hashes are retained in `verification.json` and
`source_binding_verification.json`. Final identities are in
`source_lock.json`. Neither receipt is rewritten to invent prior source
knowledge. No large computation, shared MATLAB, Git, n=7 work or edits
outside this audit directory were performed.

## 1. Correct inverse power, closed form and Fubini law

The definition is precisely

`B=H2=(2rho)^-1 integral_capture phi^2/ell_R dz`.

It contains **one** inverse ell. Substituting
`Q=-log(1-u/(2R))` on a tail gives the primitive
`(2R/3)Q^3`. After the `1/(2rho)` normalization and two retained
tail lengths g and v, the exact result is

`B(w;q)=R[Q(R,g)^3+Q(R,v)^3]/(3rho)`.

If m>g, the extra deleted central segment has phi=0, so there is
no -m term. The relative domain and rho>=4r ensure at most one source
endpoint is deleted and the full other tail remains. Strict capture
endpoints have zero integral measure and do not change the formula.

At the baseline R=r, the target W capture probability of each fixed
source z is exactly ell_R(z). Fubini therefore cancels that **one**
inverse ell and gives `E_W B=E_fullsource phi^2`, not
`E_fullsource(phi^2/ell)`. The exact second-moment primitive is

`K2(a)=integral_a^1(log y)^2 dy=2-a(Q_a^2+2Q_a+2)`.

Consequently `E_W B=(2R/rho)K2(a)`, and the limiting coefficient
is `mu2=3K2(5/6)=1-5log(6/5)-(5/2)log^2(6/5)`.
Both raw integration identities and the derivative `K2'=-log^2 a`
were independently checked symbolically. All target laws here are
full-target laws before D restriction; they are not laws for actual mu.

## 2. Positive lower mean requires the true multiplicative chain

The proof correctly does more than impose g/R<=1/3. Its exact relations

`a-5/6=(2h-3epsilon q)/(6R)` and `R/q=3/2-h/q`

give `0<=a-5/6<=h/(3R)` and `R/q>=3/2-3h/2`.
Thus `K2(a)>=K2(5/6)-h/75`, leading to the stated
`mu2-h(mu2+1/25)` free-mean lower bound. The product lower factors
are positive throughout the displayed h range. Epsilon's positive
contribution to g is retained and bounded in the correct direction.
Under only the relative upper gap constraint, g/R could tend to zero;
no uniform positive H2 mean would follow. The proof explicitly avoids
that unsupported extension.

For an independent different enclosure, set
`cmin=(3/2)-(3/2)2^-20` and `amax=(1+1/cmin)/2`.
Monotonicity gives the lower free coefficient
`2cmin K2(amax)` and upper `3K2(5/6)`.
The new alternating rational log bounds enclose these by approximately

`0.005289304047020864 < n E_free(B)/x < 0.005289341165014467`.

The actual target profile has **n-1** free coordinates plus a fixed pin.
The pin is nonnegative and is below `1/(250rho_b)`. At n>=2000,
paying the omitted free coordinate below and the pin above using the
independent enclosures yields

`0.005286659394997354 < E S_b/x < 0.005292341165014466`.

These bounds imply the source's wider strict interval
`0.00528x<E S_b<0.005295x`. No pin is treated as a free Uniform
variable. Both finite corrections improve as n increases.

## 3. Strong derivative, genuine kink and continuous matching

The proof's derivative chain is the new multiplicative one:
`R'=3/2`, `g'=1/2+epsilon<1`. For u<=g<=R/3,
`|Q_R|<=1/(5R)`, `Q_u<=3/(5R)`, `Q<=1/5`.
The complete and residual Q derivatives are bounded respectively
by `9/(10R)` and `6/(5R)`. Therefore for their cubed sum S3,
`S3<=2/125` and `|S3'|<=63/(250R)`.

The prefactor derivative is essential and included:
`(R/(3rho))'=h/(3q rho)`. With q>=2/3 it contributes at most
`2h/(250rho)`. The combined bound is exactly
`|partial_q B|<=(21+2h)/(250rho)<1/(10rho)`.
At maximal h, its coefficient is `0.08400000762939454`, well below .1.

At the capture onset `|w|=a0`, the one-sided derivative jump is
`(1-epsilon)Qpin^2/(2rho ell_pin)>0`; this is generally a real kink.
At `|w|=R`, the residual Q^3 derivative vanishes and the first
derivative joins. Finite continuous smooth branches and the uniform
derivative bound give local Lipschitz continuity and absolute continuity
on the whole matching interval. Thus integration uses a.e. derivatives,
not a false everywhere-smoothness assertion.

Using `n/rho<=3x/2` and q-width `4h/3` gives the uniform matching
shift at most `xh/5`. For h<=2^-20, the baseline interval
`[.00527,.00530]x` is therefore contained after matching transport in
the claimed safe interval `[.005269,.005301]x`. This stronger statement
also implies the coarser `[.0052,.00531]x` interface mentioned in the
task. The first-profile derivative below `1/(2rho)` gives its separate
shift at most xh, within the original `10^-5 x` budget here.
No finite-difference or mesh argument is used.

## 4. Target concentration and the two different pin bounds

For any allowed target coordinate, R<3/2 and Qpin<1/5 give
`B(w;q)<1/(125rho)`. At the **baseline target's truly pinned**
coordinate w=+-r, only one Q^3 term remains, giving the half bound
`B(+-r;q_b)<1/(250rho_b)`. The latter bound cannot be applied to
an arbitrary coordinate index omitted by a source-face pin.

On each full target face the n-1 free coordinates are independent and
the pin has zero variance. Thus
`Var(S_b)<=n/(125^2 rho_b^2)`. The two-sided event endpoints
are at least `x/200000` from the full mean interval above.
Chebyshev gives the correct common failure coefficient
`5760000/n`. This is an upper bound even where it exceeds one;
it becomes useful in the stated large-n/finite-L scopes. It is uniform
over the actual local target radius range and signed face. No independence
between different good events or after D restriction is claimed.

## 5. Every inverse-source face: raw second moments are centered

For a fixed good target and captured signed source face j, the free
coordinate inverse law is `capture/(ell_R J_i)` relative to fullsource
Uniform[-rho,rho]. Therefore its phi mean is exactly `A_i/J_i`,
its raw second moment exactly `B_i/J_i`, and its variance is
`B_i/J_i-(A_i/J_i)^2`. Free coordinates are independent in this
law; the source pin contributes deterministic Qpin and zero variance.
The exact sum in the source is correct. It is neither a sum of raw
moments nor the raw moment of the full V without its cross terms.

Let `|J_i-1|<=eta<=1/1024`. The source's definite pointwise bounds
`A_i<3/(50rho)` and `B_i<1/(125rho)` hold for all target coordinates.
Removing source pin j must use the latter **full** bound; division
by J remains in both the second-moment and squared-mean terms.
For S2 in `[l2,u2]x`, the exact directed bounds are

`Var(V|j)/x >= [l2-C2/(nq)]/(1+eta)`
`                 - C1^2 (x/n)/[q^2(1-eta)^2]`,
`Var(V|j)/x <= u2/(1-eta)`.

Under n>=2000, q>=2/3, x/n<=1/1000 and eta<=1/1024, independent
exact arithmetic confirms the source's fractions:

`lower=9776676208/1862312890625 > 13/2500`,
`upper=912/171875 < 531/100000`.

Numerically these are about .0052497495 and .0053061818. These
remain **per-face** bounds. Total mixture variance includes variance
of face means and cannot inherit this upper interval without another
argument. Actual mu and cap/soft-conditioned measures are not covered.
No first-profile good event is needed for this variance lower bound:
the squared means are paid with their pointwise O(1/rho) bound.

## 6. New wide-h local finite guards at L>=4096

The proof's L>=4096 checks are correct in their explicitly limited
scope: x in [L/2,L], hc in [1/512,1/256], inward trim delta=2^-19,
source qmin=2/3, qmax<=1 and target upper radius <3/2. At the widest
endpoint h=2^-20=delta/2, the closed lower reference may reach R=1,
which is allowed by the local R>=1 domain; strict actual matching
remains above the endpoint. No assertion needs an old dimension threshold.

Taylor2 gives n>=L^2/2, hence n>=2000, x/n<=2/L<.001 and
rho>=L/3>6>=4r. The epsilon upper bound `8/L^5` is already
below `(2/3)/(512L)` because L^4>6144. This certifies the true
relative gap even for the narrow h choice. All associated ratios improve
as L increases.

The independent extension used weaker bounds `exp L>=L^3/8` and
`exp L>=L^4/32`, rather than the author's denominators 6 and 24.
They give `|J-1|<4/L^2<2^-20` and target failure at most
`5760000*32/L^4`, about `6.54836e-7<2^-20` at 4096.
Their negative powers prove the infinite tail. Thus these guards also
validate the stronger author's explicit figures. A complete mechanism
still must check its own trim, source/radius ranges, epsA allocation,
other kernel/degree/eligibility guards and operator passage.

## Reproduction and scope boundary

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python /Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/audit/round070/second_profile/check_independent.py
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python /Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/audit/round070/second_profile/check_source_binding.py
```

The original two actual receipts are the audit record; rerunning naturally
creates new time/PID identities. No Cantelli mean condition, BE third
moment/error bound, face-mixture independence, cap-conditioned law or
complete C_n theorem is added by this local ACCEPT. Those future uses
must preserve the fixed-face law and all displayed hypotheses.

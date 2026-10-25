# B044: the sparse4 frame is legal; all 21 old cones require refinement

At beta=1/80, adjoining the 120 unoriented support-four directions gives
a genuine conservative Cauchy convolution semigroup. The exact frame
Gram is (2561/80)I and each coordinate has Cauchy rate 33. However, every
one of the frozen 21 cones is split by at least one new linear form.
Changing only the old six ray rates q in each cone is invalid.

This is a frame/symbol preflight, not a new fan or a kernel certificate.
There are zero new G/H/K/p/native evaluations, zero time parameters,
zero grids or optimizers, and no new C6 bound. No old tube, point value,
Green floor, or all-time obstruction is transferred to this frame.

## Exact frame and normalization

Let V={(1,s2,...,s6):si=+/-1}. Let W consist of all vectors in {0,+/-1}^6
with exactly four nonzero coordinates and first nonzero coordinate +1.
There are binom(6,4)=15 supports and 2^3=8 unoriented sign patterns per
support, giving 120 distinct columns. The new matrix has columns

    A_beta=[V, beta W],       beta=1/80.

Every column carries an independent standard Cauchy(t) scalar, all at the
same time t. The canonical exact 152-column input is `frame.json` with
common denominator 80 and SHA256

    67a12cd67961baf5ee554c138c82d6b4b627701c66dc691cb6f8b2a823c9be1c.

This hash identifies a new frame; the old body32 kernel memo entries
do not apply. The sign convention removes a global +/- duplication,
not a factor of two in the symbol.

A fixed coordinate occurs in binom(5,3)*8=80 new columns. For two distinct
coordinates, sign reversal in one coordinate pairs the corresponding
products with their negatives; reorienting the whole vector to restore
its first positive entry does not change that product. Hence

    sum_(w in W) ww^T=80I,   sum_(v in V) vv^T=32I,
    A_beta A_beta^T=(32+80 beta^2)I=(2561/80)I.

The checker independently sums the exact matrices and coordinate
incidences. A coordinate marginal of a sum of independent Cauchy
variables has rate equal to the sum of the absolute coefficients, so

    gamma_i=32+80 beta=33,
    marginal density at a = 33t/[pi(a^2+(33t)^2)].

The rate depends linearly on beta; the frame Gram depends on beta^2.
The Gram is a geometric right-inverse identity, not a covariance matrix
of Cauchy variables and not the Fourier symbol of a diffusion. The full
unoriented orbit is preserved by signed coordinate permutations, so the
new kernel retains the required full cube symmetry.

## The true positive conservative semigroup and generator

Put Z_a(t) independent standard symmetric Cauchy scalars. The probability
law of X_t=sum_(a column of A_beta) a Z_a(t) has characteristic function

    exp(-t psi_beta(xi)),
    psi_beta(xi)=sum_(v in V)|v.xi|+beta sum_(w in W)|w.xi|.

This is a positive probability pushforward of a product of probability
densities. Equivalently it is the convolution, at the same time t, of
the original body32 law with the independent sparse4 law. Since the
body frame spans R6, psi_beta>=psi_body>=c|xi| for some c>0. Fourier
inversion therefore gives a smooth density p_t^beta for each t>0, with
integral one. Independent Cauchy addition proves
p_(s+t)^beta=p_s^beta*p_t^beta. Approximation to the point mass as t goes
to zero gives a strongly continuous positive contraction semigroup
P_t on L1(Lebesgue), and integration is preserved.

To make the generator sign unambiguous, let A be the nonnegative Fourier
multiplier psi_beta. The Markov generator is L=-A, and P_t=exp(tL).
For a smooth compactly supported f,

    Lf(x)=(1/pi) sum_a integral_0^infinity
       [f(x+s a)+f(x-s a)-2f(x)] ds/s^2.

The symmetric second difference makes the near-zero integral finite,
and the bounded numerator makes its far tail integrable. The Fourier
identity integral_0^infinity(1-cos(us))ds/s^2=pi|u|/2 gives exactly
symbol -psi_beta. The generator kills constants and obeys the positive
maximum principle (Lf<=0 at a global maximum); its semigroup preserves
positivity. The generator itself is not asserted to map every positive
function to a positive function. Scaling the column beta w produces
the factor beta in the symbol, exactly as in the marginal calculation.

The degree-one homogeneity of psi_beta implies

    p_t^beta(x)=t^-6 p_1^beta(x/t),
    p_t^beta(r x)=r^-6 p_(t/r)^beta(x),
    K_T^beta(r x)=r^-5 K_(T/r)^beta(x),
    K_T^beta=integral_0^T p_t^beta dt>=0,  integral K_T^beta=T.

The last mass identity follows from Tonelli; K_T may have infinite
point values. In distributions, LK_T=p_T-delta_0. No old frame's values
or singular-neighborhood constants appear in these identities.

For completeness, this frame has a valid direct time-average upper-bound
interface with one common spatial radius:

    r^-6 K_T^beta(x/r)/T
       = (1/(rT)) integral_0^(rT) p_s^beta(x) ds.

The positive conservative semigroup maximal average has weak(1,1)
constant one. A self-contained route uses, for P=P_h and f>=0,
R_0=0 and R_n=(f-lambda+P R_(n-1))_+. Positivity gives increasing R_n,
finite-measure D_n={R_n>0}, and

    integral_Dn(f-lambda)
      = integral R_n - integral_Dn P R_(n-1)
      >= integral R_n - integral R_(n-1) >= 0.

Also R_n>=sum_(j=0)^(n-1)P^j f-n lambda, so all discrete average bad sets
lie in the increasing union D_n. Strong L1 continuity passes finitely
many rational times through Riemann sums, then rational exhaustion and
the almost-everywhere local absolute continuity of time integrals give
all times and arbitrary nonnegative L1 data. This argument uses no
compact support, boundedness, special measure, or direction-specific
kernel estimate in the final interface.

Thus if a future certificate proved
m_T=ess inf_[-1,1]^6 K_T^beta>0, cube domination would give the conditional
cost T/(64m_T). A strict cost below 64 requires m_T>T/4096. This batch
does not evaluate m_T or any T. The old body32 all-time obstruction is
a statement for a different symbol and cannot be inherited here.

There is also a cheap symbolic necessary condition, not a new numerical
budget. Fubini on almost every slice x1=a in (-1,1), together with the
coordinate marginal, gives

    32 m_T <= integral_R5 K_T^beta(a,y)dy
             = log(1+(33T/a)^2)/(66pi).

Taking a sequence of such slices increasing to 1 avoids any boundary
regularity assumption and yields

    m_T <= log(1+1089T^2)/(2112pi),
    T/(64m_T) >= 33pi T/log(1+1089T^2).

In particular log(1+1089T^2)>33pi T/64 is necessary for this model's
strict-64 target. No time was substituted and no optimization of this
condition was performed. It alone proves no model exclusion or C6 bound.

## Exact failure of the old fan interface

The frozen old fan input SHA256 is

    7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812.

It contains 21 full-dimensional simplicial cones in the ordered Fourier
chamber, with 18 distinct integer rays. We reused its previously proved
coverage and linearity, and independently checked the selected witness
cone determinant. We did not regenerate the old fan or run its writer.

The preflight computes exactly 18*120=2160 distinct ray/form dot products.
All 21*120=2520 cone/form sign tests then use cached values on the six
rays. A form splits a cone exactly when those values contain both a
strictly positive and a strictly negative value. Zeros are retained as
boundary data, never perturbed or discarded. All 21 cones split. The
split-form counts in the frozen cone order are

    5,11,9,13,11,5,19,17,21,19,13,18,20,24,22,16,22,16,14,10,10,

for 315 cone/form incidences. These are not counts of refined cells or
distinct hyperplanes globally. The complete per-cone form indices are
saved in `results.json`; indices use lexicographic supports and then
lexicographic signs of the last three nonzero entries.

The first certificate is old cone1, whose six integer rays are

    r1=(2,1,1,1,1,0), r2=(3,3,2,2,2,0), r3=(2,2,2,1,1,0),
    r4=(1,1,1,1,0,0), r5=(1,1,1,1,1,0), r6=(1,1,1,1,1,1).

Their determinant is exactly 1. The new support-four form
ell(xi)=xi1-xi2-xi3+xi4 has ray values (1,0,-1,0,0,0). Let S=sum_i r_i.
The two exact integer points

    xi_plus =S+r1=(12,10,9,8,7,1),
    xi_minus=S+r3=(12,11,10,8,7,1)

have respectively coefficients (2,1,1,1,1,1) and (1,1,2,1,1,1) in the
cone basis. All coefficients are positive, the determinant is nonzero,
and both coordinate chains are strictly decreasing and positive. Thus
these are strict interior points of that same old cone and chamber.
Their ell values are +1 and -1, and the midpoint has ell=0.

On this cone the old symbol is linear. Each new absolute linear form is
convex, and the selected term alone contributes a strict midpoint defect
beta=1/80. Therefore

    (psi_beta(xi_plus)+psi_beta(xi_minus))/2
                   -psi_beta((xi_plus+xi_minus)/2) >= 1/80 > 0.

No other new term can cancel this defect because its coefficient is
positive. This proves that the new symbol is not linear on old cone1,
independently of merely counting sign changes. The same cached sign
criterion supplies a split for every other old cone.

Assigning the true new symbol values at the six old rays would give only
a linear interpolant on each old simplex. By subadditivity it is an
upper envelope of the new symbol there, generally strict. Replacing the
symbol by that interpolant in an oscillatory Fourier inverse is not an
identity or a justified spatial-density inequality. The old C013/C014
integration engines cannot receive this 21-cone q-only update as a true
new-frame input.

## Smallest proposed refinement acceptance object — not executed

The first local test should be only cone1 and ell above. In old nonnegative
cone coordinates the cut is lambda1-lambda3=0. The two half-cones have
generator lists

    (r1,r1+r3,r2,r4,r5,r6),
    (r3,r1+r3,r2,r4,r5,r6).

The union identity follows by subtracting the smaller of lambda1 and
lambda3 and assigning it to r1+r3; their interiors are disjoint and the
common cut is retained. This is an analytic proposed object, not a
generated or executed fan. It still need not be a linearity partition
for the other four splitting forms.

If separately authorized, the next minimal complete local object is
cone1 cut by its five cached crossing forms, indices 2,10,26,50,82.
At most 32 closed sign patterns need be considered. A valid local
deliverable must certify each feasible polyhedral cone, retain its
boundary/zero forms, give exact simplicial subdivisions, prove exact
union and disjoint interiors, and verify all 120 signs on every child
ray. The other 115 forms already have a fixed weak sign on the parent.
This local test would not certify the remaining 20 parents or create
a usable full Fourier fan.

For an eventual full fan, each old parent must be completely refined;
shared faces must be represented consistently, or overlaps must be
proved lower-dimensional. Every child must have an exact determinant,
primitive rational/integer rays, and the true rational rates
q_beta(r)=sum_V|v.r|+beta sum_W|w.r|. Coverage, sign constancy, scaling
normalization and full signed-permutation multiplicities must all be
checked. Only then could a kernel engine be assessed for zero/repeated
Fourier nodes and exact cancellation. No such refinement, kernel call,
precision test or new candidate point is authorized or performed here.

## Runtime, scope and reproduction

Registration and source hashes were written before arithmetic. The one
recording run completed as PID47900, exit0, tool chunk45ca04. The exact
checks and result recording took 0.007798 seconds, and the complete
printed process work took 0.008186 seconds. There were 2160 distinct
ray/form dot products, 2520 cached cone/form configurations, and ten
registered orbit/Gram/normalization relation groups. No additional beta,
spatial kernel point, time, scalar budget, or fan enumeration was tried.

`check_preflight.py` uses Python3.15.0 and standard-library Fraction only.
It imports no author code; runtime/source/input hashes are preserved in
`run_receipt.json`. The inherited old fan and bridge/audit records are
read-only and hash-locked in `input_manifest.json`. Every process exited;
there are no child handles. No old artifact was modified and no commit,
push, expansion of agents or MATLAB action occurred.

Read-only replay of the same saved checks, with no kernel or new fan:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round044/sparse4_frame_preflight/check_preflight.py
```

The recording mode refuses to overwrite an existing result. The batch
has stopped at the exact fan incompatibility. New frame viability and
any competitive weak-type constant remain unresolved; public/internal
C6 upper bounds remain 64. The parent retains final proof integration.
The read-only replay of these same configurations passed as PID48205,
exit0, in 0.008062 seconds (tool chunkb9c102); it added no configuration.

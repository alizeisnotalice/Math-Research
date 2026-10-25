# A043: a true short-time middle-tube bound from two parameters

The new positive interface closes for the single registered middle tube.
Let v=(1,1,1,1,1,1), w=(1,1,1,-1,-1,-1), U0=span(v,w), and

    D={x=(1,y), y in[0,1]^5: dist(x,U0)<1/64,
                                mean_(coordinates4,5,6)(x) in[1/4,3/4]}.

For the actual original body32 semigroup, with the true positive-face
coordinate measure, the certificate proves

    ||integral_0^(1/1000) p_s ds||_L^(6/5)(D)
       <615635661087/125000000000000000
       =4.925085288696e-6.

All complementary parameter events are included. No full-time Green
estimate, Green asymptotic, whole-face I, global g or C6 improvement is
claimed. Times above1/1000, the spatial endpoints and other plane types
remain unpaid. C6 remains64.

## Why the third distance can now be paid

A uniform positive third distance along the retained parameter is false.
The saved exact counterexample is u=a_v=1/4:

    x(u)-a_v v=(3/4)(1,1,1,0,0,0)
       =(3/8)[(1,1,1,1,1,-1)+(1,1,1,-1,-1,1)].

The residual lies in this non-v pair plane, although the original x has
squared distance1/6 to it. This input is an exact structural certificate,
not a spatial kernel evaluation.

The new step keeps the actual Cauchy density in the third denominator.
For every s>0, lambda>0 and z!=0,

    integral_R [s/(pi(a^2+s^2))]/(|z-a b|^2+lambda s^2) da
       <=(1+|b|/sqrt(lambda))/|z|^2.

Theorem section4 derives the integral exactly from the one-dimensional
Cauchy convolution formula before applying the inequality. It remains
valid when |z-a b| has zeros; no infimum over a replaces that distance.
Here lambda=928 and |b|<=sqrt(6). The denominator on the right is the
distance of the original x to each relevant pair plane, bounded below
on D by g=1/(4sqrt(2))-sqrt(3)/64>0. It is not the false distance of the
shifted residual uniformly over a.

## Positive partition, quotient factors and the full remainder

Use the three actual disjoint parameter events

    E:   |a_v|>=1/2 and |a_w|>=1/16;
    E_v: |a_v|<1/2;
    E_w: |a_v|>=1/2 and |a_w|<1/16.

The E branch removes just those two parameters, with the true coarea
factor1/J2=1/6. Its quotient is the original type-k3 Q4, and its density
is at most C2 s^-2, C2=64R2/pi^6, R2=2051/47038464.

Each remaining branch retains the actual density of its small specified
parameter, and performs three positive deletions in the other31
parameters. Their projected Gram is32I-p p^T, between26I and32I, rather
than a tight32I frame. The actual counts31,30,29 give992,960,928 in the
time denominators. The first two distances have proved positive gaps;
the third is integrated by the displayed Cauchy identity.

A quotient with one parameter omitted cannot use the original Q3 without
adjustment. Convexity, coordinate symmetry and psi(ones)=60 prove
psi(xi)>=10||xi||1, hence psi_without_b>=(9/10)psi. Fourier inversion
therefore gives the rigorous dimension3 peak loss(10/9)^3. This proves
the needed comparison; no new quotient peak was computed or guessed.

For each retained parameter, all26970 ordered triples avoiding it are
included, with the previously certified coefficient
C_rem=44743079/(2970pi^6). This source value is reused explicitly; it is
not reported as a new041 run or a new projection-kernel computation.
The resulting density bounds are

    B_v,s<=18894.419782083,
    B_w,s<=25585.702063042,
    B_s=B_v,s+B_w,s<=M, M approximately44480.121845124.

The second event may be enlarged to |a_w|<1/16 for its L-infinity bound.
For the L1 step we return to the actual disjoint complement B_s, which is
at most p_s. Thus it gets one whole-density slice budget, not two. All
other large parameters, cancellation events and the isolated distance
zeros are included. No Q2, fourth deletion or unequal time is used.

## The single new numerical budget

Interpolation uses the true positive-face slice bound integral p_s dy<=s/pi.
The two contributions are

    (2/3) C2^(1/6)pi^(-5/6)epsilon^(3/2),
    (6/11) M^(1/6)pi^(-5/6)epsilon^(11/6), epsilon=1/1000.

| Paid contribution | L6/5 short-potential upper budget, display |
|---|---:|
| Both specified parameters large | 9.699500983943049e-7 |
| Entire actual complementary event | 3.955135190301248e-6 |
| Sum | 4.925085288695553e-6 |

The complement accounts for about80.3059% of this upper budget; it has
not disappeared or been treated as negligible. The sum is approximately
0.060024477 of the conditional whole-face threshold4/48750. This comparison
does not verify the separate g or Hnorm goals or pay the rest of the face.
The reported bound is not a recomputation of C037's full-time M envelope
or its old uniform990-term background budget.

## Actual computation and verification

The six configurations were registered before execution:

1. Exact32-direction symbol sum and the missing-one peak-loss constants.
2. The explicit residual pair zero, its original positive distance,
   the two-parameter Gram determinant, and the symbolic denominator identity.
3. The single tube's positive geometric gaps and Cauchy average factor.
4. The Q4/J2 double-large branch and its one short-time budget.
5. The two fully summed, Cauchy-averaged complementary-event bounds.
6. Their combined short-time budget and conditional threshold comparison.

The actual new driver PID48192 completed all six in0.011203 seconds,
exit0 (chunk62d113). It did not run a kernel or another round's driver.
Full512-bit Arb lower/upper dyadic endpoints are preserved in results.json.
There were no arithmetic failures, retries or additional configurations.

The separate saved-data checker PID48361 passed in0.001459 seconds,
exit0 (chunk95f64b). It verifies the exact zero geometry, source hashes,
positive gap certificates, the missing-peak loss and the final rational
upper bound using exact sixth-power inequalities. It calls neither Arb
nor a spatial kernel and introduces no new parameter configuration.

There were zero G/H/K/p/native calls, grids, optimizers or MATLAB calls.
No old source was edited, and no stage/commit/push or delegation occurred.
Sources are locked in registration.json, and the final file hashes and
total wall time are in manifest.json. The original033 peaks and041
restricted-triple coefficient are read-only dependencies; only the043
calculations and receipts count as new work here.

## Scope of the next acceptance object

The new useful object is the positive Cauchy-averaged resolvent estimate,
together with the proved missing-parameter quotient comparison. The
remaining tasks would be the longer-time contribution on this same D,
spatial endpoints, other plane types and legitimate union/partition
accounting. None is evaluated in this batch. The present result must not
be summed over overlapping regions as though they were disjoint.

Review files: registration.json, theorem.md, evaluate.py, results.json,
run_receipt.json, check_saved.py, checked_results.json, this report and
manifest.json. The writer refuses a second batch. The parent retains
final proof review and integration responsibility.

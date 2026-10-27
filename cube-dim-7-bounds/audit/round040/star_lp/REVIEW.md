# Independent Sol audit of the R040 tree/containment interface

Verdict: **PASS for the stated finite LP obstruction and containment-payment
interface. No mathematical correction is required.** This audit does not
endorse a new maximal-operator bound; the frozen note expressly makes none.

Audited file:
`/Users/zhengzhihao/Desktop/T/cube-r040-star-lp/cube-dim-7-bounds/upper/round040/star_lp/PROOF.md`

Frozen SHA256:
`6aece373926c7dfc65e7373b7245d117ba3fe4ae8f17fcd09492d1acc55a2a95`.
The actual file hash was checked before the independent computation. The
R039 source was also read, including its representative-preserving
arrangement-cell enlargement, to verify the exact LP being obstructed.

## 1. Full common-star LP and the lower certificate

The laminar assumptions imply that any nonempty common-star family is an
ancestor chain: two incomparable nodes have disjoint stars. Conversely every
ancestor chain has the open lowest star as a common intersection. Hence no
additional non-chain common-star columns escape the argument. Strict radius
decrease identifies the highest node as the unique maximum-radius
representative in every non-singleton chain.

The leaf probability nu is absolutely continuous with respect to Lebesgue
measure. Thus integrating an almost-everywhere representative-load bound
against nu is legitimate. Its mass on a depth-k star is exactly 2^-k,
because that star contains exactly its 2^(d-k) descendant leaf stars and
meets no other leaf star.

The proposed p_i are nonnegative. For any column with highest depth k, the
sum is bounded by the full tail

    sum_{l=k}^{d-1} 2^(-l-1) + 2^(-d) = 2^(-k).

Missing levels or termination above leaf depth only reduce this sum. The
item-cover inequalities therefore give a lower bound sum_i p_i=1+d/2 on
every feasible fractional certificate, including all possible pivots and
all nonconsecutive chains. This is an explicit feasible dual argument,
not an assumption that an optimal certificate is symmetric.

## 2. Upper certificate and all item/load inequalities

For an ancestor v at depth k and descendant leaf w, the full v-to-w path
has representative v and coefficient q_v 2^(-(d-k)). Summing all columns
anchored at v gives precisely q_v. A node at depth l receives from its
ancestor at depth k the amount q_v 2^(k-l). Consequently its exact total
coverage is

    2^(-l) + (1/2) sum_{k=1}^l 2^(k-l) = 1.

At any point, the active representative stars form one ancestor chain,
so their summed coefficients are at most 1+d/2; throughout a leaf star the
bound is attained. The d=0 case separately gives value one. Both sides
therefore establish the exact minimum, not merely one bound direction.

The same arithmetic with q_v=1-1/m for nonroots proves the stated m-ary
extension, though that extension is not needed for the obstruction.

## 3. Fixed-dimensional cube realization and beta separation

The construction stays in centered common-radius axis-parallel cubes in
every fixed n. It uses only first-coordinate translations, which suffice
to establish disjointness of the full n-dimensional cubes.

Child-star containment follows from 1/2+1/8<1. Sibling inflated stars have
center distance R_parent and radius sum beta R_parent/4<=R_parent/2,
so they are strictly disjoint. Each inflated child lies strictly in its
inflated parent because

    1/2 + beta/8 < beta,    1<=beta<=2.

Iterating this containment reduces any incomparable-node pair to their
two distinct children below their lowest common ancestor. Their inflated
stars are thus disjoint, at arbitrary depths; checking only immediate
siblings would not have sufficed without this recursive containment.

For an ancestor small cube and its child star, the first-coordinate gap
is positive exactly under the sufficient condition

    1/a+1/8<1/2, equivalently a>8/3.

Every proper descendant lies in one child star. This proves separation
of all ancestor/descendant small-cube pairs. Incomparable small cubes are
separated by their branch stars. Thus all actual original small Q_i are
pairwise disjoint and all actual increments equal Q_i with positive volume.

## 4. Eligibility and one actual common witness

At a depth-k center, all k ancestors' inflated stars contain it. No
incomparable node can contain it, because its inflated star is disjoint
from the node's own inflated star. Radius order excludes earlier
descendants, and equal-depth ties are incomparable. Accordingly the
actual previous-only eligibility sum is exactly k a^-n.

Under d a^-n<=tau, every item remains eligible, including equality at
the deepest level. Induction therefore makes this a legal largest-eligible-
radius history, rather than just a geometric family satisfying a final
load constraint. The supplied original parameterization also respects
t>1, beta=1+1/t and a>2t when its additional stated restrictions hold.

The one function (1+epsilon)lambda 1_(S_root) genuinely witnesses every
star. No independent witness masses are used. Its integral over a node
star is exactly the claimed multiple of that star's volume.

## 5. The 2N containment LP and its direction

For a general actual history, the increments I_i are pairwise disjoint.
If Q_i is contained in the anchor star S_j then so is I_i. Hence for each
containment column

    sum_{i in J} v_i <= |S_j|,

which proves its actual cost is at most one. For a singleton,
v_j<=|Q_j| gives cost at most a^-n. The coordinate test
||y_i-y_j||_infinity+r_i<=a r_j is the correct open-cube containment
criterion, including its equality case. Boundary equality does not remove
any positive-volume set or invalidate containment of open cubes.

The proof of lambda U<=K||f||_1 uses the inequalities in the correct order:
item coverage increases the sum of actual volumes; the cost definition
turns that sum into c_h |S_j|; the original witness increases it to an
integral; and the arrangement-cell load bound bounds the combined integral
by K||f||_1. It uses nonnegative coefficients at every step. A feasible
primal certificate is sufficient for the claimed upper payment bound.

The N containment and N singleton columns really give at most 2N columns;
the representative-load constraints are finite arrangement-cell
constraints. Duplicate removal is understood for identical anchored
columns, not merely identical index sets with different anchors. This is
the ordinary meaning of the note's h=(J,j) column definition.

Taking all singleton coefficients equal to one recovers the exact old W.
The containment LP is therefore no worse than that certificate on each
fixed history. The proof does not assume arbitrary enlargement preserves
an actual-volume cost; the frozen note explicitly warns that it need not.

For appended common-star columns, the R039 bound
sum_{i in J}v_i<=((a+1)/a)^n |S_j| was separately checked in the R039 source.
Thus the proposed coarse cost is large enough; using the exact ratio is
also valid. No reverse LP/primal-dual inference is hidden here.

## 6. Exact tree cost and scope of the obstruction

At depth l there are 2^l disjoint small cubes, each with volume
(2/a)^n 8^(-nl). Therefore

    U=(2/a)^n A_(n,d),       A_(n,d)=sum_{l=0}^d (2*8^(-n))^l.

The root witness constraint forces common witness mass at least 2^n,
and g=1_(S_root) realizes that weak mass cost. Thus the exact selected-union
ratio is a^-n A_(n,d). This quantity is strictly below (4/3)a^-n for finite
d and every n>=1.

The root containment column includes the full family and attains this
ratio. Conversely the same uniform root function, either as a weak witness
or after strict multiplication and epsilon tending to zero, forces any
feasible containment certificate to have at least this value. The asserted
exact optimum is therefore established on both sides.

For fixed n,tau, varying a allows d=floor(tau a^n) and hence makes the
unmodified unit-index common-star LP grow as tau a^n. The legal-depth
restriction is explicit; the proof does not incorrectly let d grow at
fixed a while keeping its constant increment weights. The stronger
n log n obstruction is expressly conditional on tau a^n/(n log n)
tending to infinity, which is sufficient for the claimed conclusion.

All these lower statements concern the specified certificate LP, not C_n.
The selected union can be paid cheaply by its root witness, and the note
does not identify it with a whole maximal level set. The obstruction does
not rule out another selector, actual-volume columns, full mass incidence,
or a uniform containment LP theorem. Such a theorem remains unproved.

## 7. One new minimal exact computation

The author's eight-case batch was not rerun. A single independent
near-boundary case was run with standard-library Fraction arithmetic:

    n=2, d=3, a=267/100, beta=199/100,
    t=100/99, tau=3a^-2=10000/23763.

This has a only 1/300 above 8/3, beta near 2, and an exact equality
eligibility test at every leaf. The check enumerated 15 original cubes,
all 85 possible nonempty ancestor-chain columns for the dual inequality,
and all 32 path columns in the primal upper certificate. It checked all
small-cube pairs, all incomparable inflated-star pairs, recursive inflated
containment, and every actual previous-center load. Results:

    K_star=5/2,
    K_contain=7046875/48666624,
    minimum small-cube gap=7/5696>0,
    minimum incomparable inflated-star gap=201/25600>0,
    maximum actual eligibility=tau=10000/23763.

These are exact checks of one new case. The general conclusions above are
supported by the analytic audit, not inferred from this finite case.

Reproduction: `python3 audit/round040/star_lp/boundary_check.py` from this
dedicated clone. The new check and its output are confined to this audit
directory. No original proof modification, MATLAB call, production scan,
n=7 work, Git operation, or change outside this audit ownership occurred.

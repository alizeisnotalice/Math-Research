# F015 independent review: positive-band Taylor certificates

**PASS for the local theorem, the eight saved calibrations, and the stated
finite-cover acceptance contract. No global lower bound for K_T, g, or new
C6 upper bound is certified.** Fixed body32 T=4/125 and T=6/125 are already
excluded by separately certified point obstructions. The local inequalities
remain valid at these times; they cannot remove those obstructions.

This review owns only this directory. A015 and every earlier F directory
were read-only. The finite batch used eight existing calibrations and two
geometry toys, with no new kernel, gradient, Monte Carlo, MATLAB, agent,
commit, or push. No C016 review is included or initiated.

## 1. Positive pushforward gives the actual negative Hessian bound

Write f_t(z)=t/[pi(t^2+z^2)], F_t(z)=product_j f_t(z_j), and let p_t be
the positive pushforward of F_t(z) dz by the spanning frame matrix A.
For the frames in scope, AA^T=lambda I, with lambda=32 for body32 and
lambda=32+beta^2 for [body32,beta e_1,...,beta e_6], beta>0. The right
inverse R=A^T/lambda satisfies AR=I and ||Rh||^2=||h||^2/lambda.

The elementary identity

    (log f_t)''(z)+2/t^2
      =2 z^2(3t^2+z^2)/[t^2(t^2+z^2)^2] >= 0

implies, with w=Rh,

    D_w^2 F_t / F_t
      =(D_w log F_t)^2 + sum_j w_j^2 (log f_t)''(z_j)
      >= -2 ||h||^2/(lambda t^2).

Here is a direct justification of the pushforward step. For nonnegative
smooth compactly supported phi on R^6, D_w phi(Az)=D_h phi(Az).
Integration by parts in parameter space therefore gives

    integral p_t(x) D_h^2 phi(x) dx
      = integral F_t(z) D_w^2(phi(Az)) dz
      = integral D_w^2 F_t(z) phi(Az) dz
      >= -2||h||^2/(lambda t^2) integral p_t(x) phi(x) dx.

Although phi(Az) need not have compact support in parameter space, this is
legitimate by inserting parameter cutoffs and passing to the limit:
|f_t'/f_t|<=1/t and |f_t''/f_t|<=4/t^2 give integrable first and second
product derivatives. The cutoff error terms tend to zero. The Fourier
symbol has coercivity sum_j |v_j.xi| >= sqrt(lambda)||xi||, so p_t is smooth
for each t>0. The distribution inequality is thus a pointwise inequality:

    D_h^2 p_t(x) >= -2||h||^2 p_t(x)/(lambda t^2).

This proof uses positive measures and integration by parts, with no
conditioning on a zero-probability fiber, omitted coarea determinant, or
assumption that a lower curvature bound follows from an upper one.

For 0<a<T, the bounded time interval permits integration and differentiation:
B=integral_a^T p_t=H_a-H_T and J2_B=integral_a^T t^-2 p_t satisfy

    D_h^2 B >= -2||h||^2 J2_B/lambda.

All these band quantities are finite and smooth even on S2. This is not
a subtraction of infinite Green values.

## 2. The finite time mesh and the Taylor sign are correct

The previously audited sparse real Cauchy shift bound is

    p_t(c+h) <= Phi(||h||inf/t) p_t(c),
    Phi(u)=exp(12 asinh(u/12))
          =(u/12+sqrt(1+(u/12)^2))^12.

One factor has ratio at most exp(2 asinh(|b|/(2t))). A point on a cube
facet is a convex combination of at most six vertices, so a spatial shift
h can be produced by at most six signed body-column parameter shifts
with absolute sum ||h||inf. Concavity of asinh gives Phi. For the augmented
frame the extra parameter shifts can be zero; the same Phi is valid.

For C=conv{v_j}, c in C, define r=max_j||v_j-c||inf and
D=max_j||v_j-c||^2. The whole cell obeys these bounds by convexity.
On each positive time band [t_i,t_(i+1)], the weight
t^-2 Phi(r/t) is decreasing. Hence, using certified upper mass endpoints,

    J2_B(x) <= U2
      = sum_i t_i^-2 Phi(r/t_i) upper(H_ti(c)-H_t(i+1)(c)), x in C.

There is no infinite-time tail. Nonnegative true mass permits truncating
the lower mass endpoint at zero, but never permits lowering its upper
endpoint arbitrarily. A negative upper endpoint is inconsistent input.

For h=x-c, the segment c+s h remains inside C. Taylor's integral remainder
pays integral_0^1(1-s) ds=1/2, and yields

    K_T(x) >= B(x)
       >= lower(B(c)) + min_j lower(grad B(c).(v_j-c)) - D U2_upper/lambda.

The linear minimum and the squared-norm maximum may be bounded at vertices
independently; they need not occur at the same vertex. The gradient must be
grad H_a-grad H_T at the same exact center, frame, and endpoints. No
convexity of B or K_T is assumed. The coefficients 1/32 and 1/(32+beta^2)
are correct. Body32 values or the body32 singular barrier do not become
augmented-frame data merely because this formula has an augmented version.

## 3. Independent replay of the eight saved calibrations

`verify_saved.py` does not import or run A's producer. It reads the frozen
A007/B010 H and gradient records, verifies their hashes and center/time
identities, parses their 384-bit balls, and extracts exact rational endpoints.
Band subtraction, gradient subtraction, all vertex dot lower bounds,
barycentric cell witnesses, final Taylor arithmetic, and threshold tests
are checked using rational arithmetic. It then bounds Phi independently
at 768 bits by the algebraic square-root formula above, rather than A's
exp/asinh expression. Each saved rational U2 upper endpoint is at least
this independently recomputed upper endpoint. Every saved final lower
endpoint is reconstructed exactly from its saved, directionally verified
components.

The fixed center is

    (1000000007,998836847,892738545,807464540,467291939,96968636)
      /1000000007,

T=4/125, and the only cases are four a values and the two original radii.
The clipped boxes have 32 vertices, are contained in the ordered positive
face, contain c with an exact convex-combination witness, and have D=5r^2.
The clipped upper x2 endpoint is used, rather than silently replacing the
cell by a symmetric box for the linear calculation.

|a|r=1/250: certified formula|r=1/125: certified formula|
|---|---:|---:|
|1/1000|8.70143157169e-6 PASS|-1.13496088053e-5 insufficient|
|1/500|8.80590055727e-6 PASS|3.63629634628e-6 insufficient|
|1/250|8.88274132543e-6 PASS|6.32049223733e-6 insufficient|
|1/125|8.71277463254e-6 PASS|7.44759742192e-6 insufficient|

The displayed decimals are summaries only. Acceptance uses the exact
rationals saved in `results.json`. Four cases pass both 1/125000 and the
weaker bare threshold 1/128000. All four passing cases certify the same
r=1/250 geometric cell. An insufficient lower formula is not an actual
counterexample to K_T.

The matched older target was also 1/125000, with pure-Phi radius
991/500000. The passing tested radius ratio is exactly 2000/991. No optimal
radius is claimed, and no new time, radius, center, or frame was evaluated.

The actual independent replay exited 0: PID 35277, script time
0.08714129199506715 seconds, Python 3.15.0, python-flint 0.9.0,
SymPy 1.14.0, tool chunk df9e0c. It executed exactly eight saved calibrations
and two geometry toys, with zero kernel or gradient calls. Reproduce with:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_saved.py

## 4. Same-plane tubes and exact tree coverage

A valid TUBE leaf must have one fixed subspace U containing every vertex
within delta. The projection onto U-perp is linear and its norm is convex,
so the whole convex cell then has distance to U at most delta. Consequently
distance to the union S2 is at most delta. Different closest planes at
different vertices do not imply this statement.

Toy 1 independently computed distances to all 496 body-pair planes for
the six root vertices and their barycenter. Every vertex has distance
zero to S2, while the barycenter has squared distance exactly 1/9. No one
of the five canonical planes covers the root vertices at delta=1/100000.
This rejects the incorrect vertexwise-closest-plane rule without assuming
A012's general closest-plane formula.

The already audited A009 barrier gives, for body32 and T>=4/125,

    dist(x,S2)<=1/100000 => K_T(x)>=3/250000.

At T=6/125 this is strictly above T/4096=3/256000. This is a local tube
statement only, compatible with the B015 obstruction outside that tube.

Toy 2 verifies that the root's longest edge is (0,5), its exact midpoint
gives child volumes 1/240 from root volume 1/120, and—more importantly—the
two barycentric replacement identities hold exactly. If lambda_i>=lambda_j,
the first child uses lambda_i-lambda_j on v_i and 2lambda_j on the midpoint;
the opposite inequality uses the other child. The weights are nonnegative
and sum to the original total. These identities prove the exact parent
union; a volume sum by itself would not do so.

The schema requires a finite valid rooted tree, exact child identities,
no orphan, missing child, cycle, or duplicated ID, and a rigorous common
threshold on every terminal leaf. Every unresolved leaf remains unpaid,
regardless of volume. The schema is a correct mathematical contract;
neither A015 nor this audit implements or executes a complete global
kernel-proof reader and cover checker. The eight calibrations are not a
root partition.

## 5. Transfer scope and excluded model/time pairs

For these signed-permutation-invariant order-one frames,

    K_T(r y)=r^-5 K_(T/r)(y)>=K_T(y), 0<r<=1.

Thus a genuinely complete ordered-face lower certificate would extend to
the whole cube. The separately audited positive Markov maximal theorem
would then give C6<=T/(64m_star) if m_star>T/4096. The relevant dominating
kernel is the full short potential K_T, with mass T. The positive band is
used only to lower-bound it; this argument does not give the band a
mass-(T-a) maximal theorem or allow replacing T by T-a in the final cost.

The source contract correctly checks approved obstructions before a global
claim. B010 excludes body32 T=4/125, and B015/F015 excludes body32 T=6/125.
Those obstruction records were hash-checked, with no counterexample kernel
replay. Consequently neither time is a viable full-domain target here.
No new time or augmented-frame result is inferred. The exact remaining
gap is a complete, model-matched global lower certificate at a nonexcluded
time/frame; the current local method and one passing geometric cell do
not provide it.

All source identities, including A's documented post-run context-only
theorem/schema update, are bound in `results.json` and `manifest.json`.
This review used the L03 rational-and-interval-certificate skill only for
directed certificate workflow; it does not rely on a claimed general API
verification by that skill or on any unread literature theorem.

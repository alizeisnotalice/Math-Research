# F020 independent review: body32 plus the pair orbit

**PASS for the mathematical interface and the frozen finite structure.**
The mathematical interface is valid for the real positive 62-column
Cauchy frame A_beta=[V,beta P], beta>0, where V contains the 32 unoriented
body representatives and P={e_i+e_j,e_i-e_j:i<j}. No new kernel minimum,
quantitative singular barrier, or C6 improvement follows from the frame
construction. The finite structural replay and its receipt are recorded
at the end of this report.

## 1. Symbol, Gram, and the old fan

For real a,b, |a+b|+|a-b|=2max(|a|,|b|). Therefore, with z the decreasing
rearrangement of |xi|,

    psi_pair(xi)=2 sum_(i=1)^6 (6-i) z_i.

On the entire closed ordered chamber xi1>=...>=xi6>=0 this is
(10,8,6,4,2,0).xi, including ties and zeros. A sum hyperplane intersects
that chamber only where its two nonnegative coordinates vanish; a
difference hyperplane is a chamber boundary. No new pair hyperplane
cuts through the interior of an old canonical cone. This is a proof on
the entire cone, stronger than a test of the ray values alone.

Each unordered support pair contributes
(e_i+e_j)(e_i+e_j)^T+(e_i-e_j)(e_i-e_j)^T=2(e_i e_i^T+e_j e_j^T).
Every coordinate occurs in five supports, giving

    A_beta A_beta^T=(32+10 beta^2)I,
    psi_beta=psi_body+beta psi_pair.

The coordinate marginal Cauchy rate is 32+10 beta, not the Gram scalar.
At beta=1 both happen to equal 42. The proposed marginal normalization
42*(32/875)=32*(6/125)=192/125 is correct algebraically and is not an
evaluation or prediction for the new short potential.

On each old cone, preserve its integer rays and determinant and update

    q_new(r)=q_old(r)+beta(10r1+8r2+6r3+4r4+2r5).

At beta=1 the new rates are positive integers. The cone weight becomes
det(rays)/product q_new. The full signed-permutation group preserves
both direction orbits, so the 46080 old group images still give the
whole-space geometry. This reuses the frozen exact geometry and bridge
certificates; it does not re-enumerate or numerically validate an old fan.

A new fan hash denotes a different kernel. Its rates change node
denominators and divided-difference collision equations. Old value,
gradient, overflow, genericity, or quantitative S2-distance certificates
cannot simply be relabeled. The raw zero forms x.r=0 retain the same
rays, but this alone is not the whole evaluator input contract.

## 2. Actual S2 plane classes and dependent triples

Define S2 as the union over all 1891 distinct column pairs. Every such
pair is independent: there are no repeated unoriented directions among
these 62 columns. Positive beta rescales pair columns and does not change
these subspaces. The following count is of plane objects, never disjoint
tubes or disjoint volume.

Two body directions split coordinates into a nonempty agreement block
and a nonempty disagreement block. A body vector in their plane has
block multipliers of absolute value one, so that plane contains exactly
two body lines. Thus all 496 old body-pair planes are distinct.

For p=e_i+s e_j, half of the 32 body lines have v_j=s v_i. Then p is,
up to sign, the restriction of v to that two-coordinate support. Flipping
those two entries gives another body line, so the body/pair plane is an
old size-2/4-block body plane. The 480 matching pairs duplicate 240 old
planes, with exactly two body lines and one pair line on each plane.

The other 480 body/pair pairs satisfy v.p=0 and create distinct new
planes. A body vector in such a plane must equal +/-v outside the two
pair coordinates; the two remaining sign equations force its pair
coefficient to vanish. A pair vector in the plane must have zero v
coefficient, since otherwise all four coordinates outside the pair
support are nonzero. Hence these planes contain exactly one body line
and one pair line and cannot coincide with old body planes.

For two pair directions there are three classes:

|Support relation|Direction pairs|Distinct planes|Frame lines per plane|
|---|---:|---:|---:|
|Same support|15|15|2|
|Exactly one shared coordinate|240|80|3|
|Disjoint supports|180|180|2|

The same-support plane is a coordinate plane. For a shared coordinate,
the two pair vectors and their difference form the three lines in a
three-coordinate plane whose normal has three entries +/-1. There are
20 coordinate triples and four unoriented normal signs, hence 80 planes.
For disjoint supports, choose four coordinates, one of three pairings,
and four relative signs: 15*3*4=180. A combination using both disjoint
lines has four nonzero coordinates, so no third pair line is in the plane.
These classes have support unions of sizes 2, 3, 4 and cannot coincide;
none contains any body direction with six nonzero coordinates.

Consequently the exact totals are

    1891 direction pairs =496+480+480+15+240+180,
    1251 distinct planes =496+480+15+80+180,
    755 new planes =1251-496.

Exactly 240 old matching planes and 80 triangle planes contain three
frame lines. The other 931 planes contain two. Every dependent triple
of distinct frame directions spans one of these planes; there are
therefore exactly 320 such triples. In particular,

    (e1+e2)-(e1+e3)-(e2-e3)=0

is a new dependency. Any-three-independent reasoning is invalid here.

A useful sufficient strict-point criterion needs no distance computation:
if all six coordinates are nonzero with pairwise distinct absolute values,
the point is off S2. A body/body plane has at most two absolute levels;
a body/pair plane has the same absolute value on four coordinates outside
the pair support; a pair/pair plane is supported on at most four coordinates.
This criterion proves nonmembership only, not a quantitative distance.

## 3. True singular set and continuity, allowing all dependencies

The general finite positive-frame argument already proved in F012 applies
without requiring any triple to be independent. Here is a direct check
with the current 62-column physical matrix.

For any independent six-column basis B, positive change of variables gives

    p_t(x)=1/|det B| integral
       product_(i in B) f_t((B^-1(x-A_R u))_i)
       product_(j in R) f_t(u_j) du.

This is a positive pushforward formula with the single Jacobian
1/|det B|, not a conditional law given a zero-probability fiber. It defines
a continuous fixed-time density. No extra full-frame coarea Jacobian
may be multiplied into this already solved-basis formula.

For x in S2 choose its two independent columns and extend them to B.
Then B^-1 x has its last four entries zero. Restrict all 56 remaining
parameters to |u_j|<=t, an event with exact positive mass 2^-56. For
0<t<=1 the first two solved parameters remain bounded and their two
Cauchy factors are bounded below by c t^2. The other four solved
parameters are O(t) and contribute at least c' t^-4. Thus
p_t(x)>=c_x t^-2 and every K_T(x) and G(x) is infinite. The proof also
covers a zero pair coefficient, a line, the origin, and intersecting
planes. It chooses one pair and never multiplies by the plane count.

For dist(x,S2)>=delta>0, set M=sum_j||a_j|| and eta=delta/(2M).
If the columns with |z_j|>eta had rank at most two, their sum would lie
in a column-pair plane while the remaining sum has norm at most delta/2.
Hence a parameter vector mapping to this region contains an independent
large triple. Sum over the finitely many independent triples, allowing
overlap, and extend each to a six-column basis B_I. The three restricted
densities are at most t/(pi eta^2), the other three basis densities are
at most 1/(pi t), and the remaining probability integral is one. Each
event therefore has density at most 1/(|det B_I| pi^6 eta^6).

This gives p_t=O(1) uniformly away from S2. The measure inequality is
first proved on an open neighborhood, then upgraded from almost
everywhere to every point by continuity. The fixed-basis peak C t^-6
controls large time. Dominated convergence gives finite, continuous G
and K_T off the true S2, and joint continuity for positive cutoff times.

This also validates the substance of A's alternative quotient proof.
After removing any two columns at least 30 body columns remain, with
coercivity at least sqrt(32)-2sqrt(6)>0. In an orthogonal quotient the
deleted selected columns vanish, so the remaining quotient symbol is
the full symbol restricted to that quotient. The quotient density and
Jacobian can be justified with an independent projected basis and positive
Tonelli integration. The six-basis proof above avoids any marginal-value
ambiguity at quotient zero. Neither proof gives a usable new quantitative
singular tube barrier by itself.

At an off-S2 cube point a strict new-model K_T<T/4096 would extend to a
positive-volume cube neighborhood by continuity, even for a face point.
Finite point values above threshold still cannot certify a global minimum.

## 4. Positive Markov interface, domain, and all scales

Independent common-scale Cauchy variables give the characteristic function
exp(-t psi_beta), positive probability convolution kernels, P_t1=1,
L1 and L-infinity contractions, and symmetry. The estimates

    sqrt(32+10beta^2)||xi|| <= psi_beta(xi)
      <=(32+beta sqrt(220))||xi||

follow respectively from l1>=l2 and Cauchy-Schwarz for the body sum and
the chamber vector ell. Hence the L2 generator -psi_beta(D) has domain
H1 and form domain H^(1/2). Strict positivity follows from one invertible
six-column product density convolved with the other probability factors.
Order-one scaling supplies a probability approximate identity and strong
L1 continuity. None of these statements assumes finite ambient measure.

With K_T=integral_0^T p_t, Tonelli gives mass T. Testing the semigroup
identity gives L_beta K_T=p_T-delta_0 in distributions. This does not put
the singular K_T in the L2 generator domain. For mass-preserving dilation,

    (K_T)_r=r^-6 K_T(x/r)=r^-1 K_(rT),
    L_beta(K_T)_r=(p_(rT)-delta_0)/r.

The previously independently proved positive Markov maximal inequality
therefore applies. Its discrete recursion has finite-measure positivity
sets because D_n is contained in {f+P R_(n-1)>a}; positivity, P1=1 and
L1 contraction give a|D_n|<=integral_(D_n) f. Finite rational-time common
meshes, strong L1 convergence and Chebyshev pass to continuous time;
Tonelli supplies a common almost-everywhere time-continuity set before
exhausting all positive times. Thus the weak constant is one for general
L1 inputs after taking absolute values, not just for finitely many scales.

If a future matching certificate gives K_T>=m>T/4096 almost everywhere
on [-1,1]^6, kernel domination consequently gives

    C6<=T/(64m)<64.

Every dilation uses the same 62-column frame and time rT. This is a
common-radius cube result, not six independently chosen coordinate times.
Signed-permutation symmetry and K_T(r y)=r^-5 K_(T/r)(y)>=K_T(y) for
0<r<=1 also justify an ordered-face reduction, if that whole face is
certified. The required new-model global m is absent in this audit.

## 5. Taylor/Phi extension and hard non-inheritance boundary

R=A_beta^T/(32+10beta^2) is a right inverse with
||Rh||^2=||h||^2/(32+10beta^2). The product Cauchy logarithmic Hessian
bound gives, for a positive band,

    D_h^2 B>=-2||h||^2 J2_B/(32+10beta^2).

Taylor's half-weight gives loss D U2/(32+10beta^2). Sparse shifting of
at most six original body parameters, leaving all thirty pair parameters
fixed, preserves Phi(u)=exp(12 asinh(u/12)). Both facts are theoretical
extensions of the previously audited positive-pushforward argument.

The frozen A015 checker accepts only body and body-plus-equal-axes frames;
it does not now validate this pair frame. New source records, values,
gradients, quantitative singular barriers and a matching checker interface
remain necessary. Increasing the Gram scalar does not prove an ordering
of old and new kernels or improve any numerical global minimum.

## 6. Independent finite structural replay

The independent checker in this directory never imports or runs A's
writer. It uses exactly twelve fixed structural configurations: the one
beta1 frame with its marginal units; the symbolic chamber identity; the
existing 21-cone rate artifact; the five canonical plane classes; a full
1891-pair plane classification of that same frame; and the three already
registered symbol vectors. No new beta, spatial kernel point or cutoff
was tested. The time fraction is algebraic normalization only.

For plane classification the checker uses the fifteen 2x2 minors of two
integer columns, divided by their gcd and normalized by sign. These
primitive Pluecker coordinates identify their two-dimensional span.
Thus the independent deduplication examines every actual column pair,
instead of merely summing proposed class counts. Dependent triples are
derived from the resulting plane memberships; no triple enumeration
is needed. The fan check reads the frozen old coverage/bridge evidence,
compares all 126 existing ray occurrences, recomputes their true 32+30
symbol rates and new weights, and confirms unchanged rays/determinants.
No old fan baseline or high-dimensional kernel is replayed.

The actual independent run exited 0: Python PID **37385**, script time
**0.020548332977341488 seconds**, tool chunk **f800bd**, Python 3.15.0 and
SymPy 1.14.0. All twelve configurations passed. The independent 1891-pair
deduplication returned 1251 planes with 931 two-line planes and 320
three-line planes; the latter account for exactly 320 dependent triples.
All 126 rate entries passed, with new rates ranging from 42 to 230.
There were zero new kernel, Monte Carlo or MATLAB calls.

The author had announced the final report/results/input-manifest hashes
and stopped its writer before this run. Those hashes were explicitly
verified. The new integer fan hash is
`4202b6c685567e51fbd53275add62b6ed978176c834b9fb3f3aaac1f1e5f8a0e`;
the descriptive JSON fan hash is
`913f40e375e95dcd57ae94ea6b30ea8f68c0e53450d49484f3ff5f6d61ff9ed7`.
Reproduce this independent finite review with:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_structure.py

The exact result, source hashes, and receipt are in `results.json`
and `manifest.json`. Only this review directory is writable. No author
source, earlier audit, checker or artifact is modified; no agents,
MATLAB, stage, commit, or push are used. All new-model numerical minima
remain outside the proved scope.

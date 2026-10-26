# A080 independent review of the specified Stein–Stromberg allocation

The reviewed claim is precisely the formula

    K(beta,m)=(1+((beta+2)/beta)^6)
              (1+m*((beta+1)/beta)^6*beta^(6/m)),

for beta>1 and integers m>=1, together with its validity as an allocated
weak-type bound for the centered maximal operator of a convex symmetric
body in R^6, and the diagnostic K(beta,m)>1602/25>64. This is not a lower
bound on the best weak constant or on C6. It does not exclude other
covering selections, estimates, or allocations.

The primary manifest, theorem, code, five saved exact facts and execution
receipt are bound by registration.json. The local original PDF was read
only at PDF pages3–5, print261–263. No PDF was downloaded or read in full
again. The primary's finite five-fact program and MATLAB display were not
run here. The original source supplies the greedy and shell argument;
the all-beta finite-six diagnostic is internal work in 080, not a claim
copied from the source or a claim of publication/novelty.

## 1. Complete greedy selection

For a finite family of open norm balls, define I_j to be the jth selected
ball minus the union of all earlier selected balls. These increments are
measurable, pairwise disjoint, and their measures are nonnegative. They
are fixed once their balls are selected; later choices do not alter old
weights w_j=|I_j|/|B_j^*|. An empty increment contributes zero and causes
no difficulty.

At each step select a largest-radius remaining ball whose center has
sum_(j<k) w_j 1_(B_j^**)(center)<=1. Initially every center qualifies.
Scores can only increase when another nonnegative term is appended, so
an ineligible ball never becomes eligible later. It follows that selected
radii do not increase. Finite termination must mean no remaining eligible
ball exists, not that an arbitrarily ordered candidate was rejected.

Fix a ball U that is never selected, with radius r. If any selected radius
is less than r, consider the first such selection. U is still available
and would have been chosen before this smaller ball if it were eligible.
Therefore its center score is already >1, and every selected ball in that
score has radius >=r. If no selected radius is less than r, use the same
ineligibility argument at final termination. This proves the required
strict score >1 using only earlier selected radii >=r. It deals explicitly
with the possible equal-radius ties.

If the center of U belongs to B_j^**, and r<=r_j, then every x in U satisfies

    |x-y_j|_B < r+(beta+1)r_j <=(beta+2)r_j.

Hence U is contained in B_j^***. The sum of triple-expanded indicators with
these positive weights is >1 on every unselected ball, thus on their union.
Integrate over that union, rather than summing its possibly overlapping
volumes. Its measure is at most

    sum_j w_j |B_j^***|=((beta+2)/beta)^6 sum_j |I_j|.

The selected union has measure sum_j|I_j|. Adding the two union bounds proves
the first factor 1+((beta+2)/beta)^6. No disjointness of the original balls
or of the unselected balls was assumed.

## 2. The outer radii really use an earlier center condition

Fix x in at least one selected beta-expansion. Choose one containing
selected ball k of smallest radius, translate x to zero, and divide all
lengths by that radius. Then r_k=1, |y_k|_B<beta, and all contributing radii
are >=1. For a contributing j with r_j>=beta,

    |y_j-y_k|_B <= |y_j|_B+|y_k|_B
                 < beta*r_j+beta <=(beta+1)r_j.

The last inequality is exactly r_j>=beta, including equality. In addition,
beta>1 implies r_j>r_k, so the nonincreasing greedy order places every such
j before k. They are therefore a subset of the terms controlled by the
selection condition at y_k. Their weighted sum is <=1.

Both features are necessary: mere membership of y_k in B_j^** would not
justify using the center condition for a ball selected later, and merely
r_j>=1 would not justify beta*r_j+beta<=(beta+1)r_j. The 080 proof uses
r_j>=beta and beta>1, so both issues are resolved. No beta=1 extension is
accepted here; equal-radius contributions would require a new argument.

## 3. Half-open shells and disjoint increment volume

For a shell [a,b), 1<=a<b<=beta, only selected balls with x in B_j^* enter
the sum. Such a ball has |y_j|_B<beta*r_j and r_j<b. Its original ball,
and therefore I_j, lies inside the norm ball about x with radius
(beta+1)b. The increments, not the original balls, are disjoint. Hence
the sum of their volumes is at most |B|((beta+1)b)^6. Meanwhile
|B_j^*|=|B|(beta*r_j)^6 >= |B|(beta*a)^6. Dividing gives

    sum_shell w_j 1_(B_j^*)(x)
      <=((beta+1)/beta)^6 (b/a)^6.

Set a_i=beta^(i/m) for i=0,...,m. The m shells [a_i,a_(i+1)) partition
[1,beta); their common ratio is beta^(1/m). The outer part already includes
r_j=beta. All shell boundaries, tied radii and the smallest radius1 are
therefore assigned. Combining the shell sums with the outer bound1 gives
the second factor of K. This argument uses only norm triangle inequality,
volume dilation and disjointness of I_j, so it applies to cubes as well
as other convex symmetric bodies. No Euclidean-volume substitution occurs.

## 4. The common-scale witness and general L1 reduction

Define the centered maximal operator by averages over x+R B for all R>0.
For f>=0 in L1 and a point x with M_B f(x)>lambda, select a positive finite
R_x with average over x+R_x B strictly greater than lambda. Use the small
ball x+(R_x/beta)B in the finite selection. Its star expansion is exactly
the witness x+R_x B. Thus the witness time/scale is not changed inside
the mass comparison.

The small balls are open and contain their respective x. On each compact
subset K of the level set, a finite subfamily covers K. Apply the finite
selection to that subfamily. Each selected star retains its original
witness inequality, so

    lambda sum_j |I_j|
       <= sum_j w_j integral_(B_j^*) f
        = integral f(y) sum_j w_j 1_(B_j^*)(y) dy
       <= (1+m*((beta+1)/beta)^6*beta^(6/m)) ||f||_1.

Combining with the covering factor proves lambda|K|<=K(beta,m)||f||_1.
The equality uses a finite sum of nonnegative integrals, not a pointwise
assumption about f. Vanishing increments do not affect it.

For fixed R, the average is continuous in x: its change under translating
x is bounded by ||f(.-h)-f||_1/|R B|, which tends to zero. Therefore the
strict maximal level set is open, hence measurable; Lebesgue inner regularity
permits taking the supremum over compact subsets. This works even before
one knows the level set has finite measure. Replacing f by |f| gives the
usual signed-input statement. Every positive R_x was allowed, and the
compactness reduction imposes no finite global list of radii or restriction
to a fixed scale. Boundary changes are Lebesgue null for these bodies.

The original print261 proof uses a slightly enlarged auxiliary maximal
operator. The direct centered-witness reduction in 080 is sufficient for
the stated centered operator; it need not introduce any extra enlargement
constant or claim an equality of the two operators.

## 5. Uniform analytic exclusion of this allocated expression

For u>0,

    d(exp(u)/u)/du=exp(u)(u-1)/u^2.

The function decreases up to1 and increases after1, so its global minimum
is e. Set a=6log(beta)>0 and u=a/m. Then
m*exp(a/m)=a*exp(u)/u>=e*a. This is valid for every positive real m and
therefore for the required positive integer m. Discreteness cannot evade
the lower estimate. The first factor of K is decreasing in beta, with
derivative -12(beta+2)^5/beta^7; the additional ratio ((beta+1)/beta)^6 is >1.

The five saved exact facts are read-only dependencies of this audit. Static
review confirms their formulas: exp(1)'s first four terms equal8/3 with
positive remainder; log2=2atanh(1/3) has first term2/3 with positive remainder;
log7=2atanh(3/4) has the stated eight-term positive prefix. The saved prefix
5223903223377/2687038914560 exceeds97/50. The saved middle product is
3245450/50421, and the stated large threshold is1602/25. The prefix sums
and these five rational comparisons were not executed again.

The domain partition has no gap: (1,2], (2,7], [7,infinity). On the first
piece the first factor is >=65 and the second >1. On the second, the first
is >=1+(9/7)^6 and the second >35/3, so the saved middle product suffices.
On the last, the first is strictly >2 and the second strictly >801/25,
giving K>1602/25. At beta2 use the first piece, and at beta7 either of the
last two pieces applies. The strict inequalities follow from positive
tails and finite beta, not from a strict version of the u=1 minimizer.
There is no conclusion at beta1 or at an infinite-valued parameter.

## 6. Independent checks and their limited scope

One new one-dimensional fixture consists of five intervals with (center,
radius) (-4,4), (4,4), (59/10,3), (8,1), (10,1), and beta3/2. It is frozen
in the registration. It tests the dimension-independent selection logic:
the radius3 candidate becomes ineligible before radius1 choices, and a
later admissible center has score exactly1. Nontrivial overlapping original
balls produce genuine increments. At x9 the fixture also tests the outer
part and the one shell [1,3/2). A step function indicator[-10,12] tests the
star witness, weighted mass identity and weak-bound bookkeeping.

These finite exact checks do not prove a universal theorem by examples and
do not provide six-dimensional body32 numerical evidence. The universal
argument is the proof in sections1–5. Independent symbolic checks verify
the derivative, the outer threshold identity, the shell-volume factor,
and assembly of the six-dimensional K. No beta/m search, main five-fact
rerun, old baseline, kernel, MATLAB, or recursion is involved.

No substantive gap has been found in the specified claim. The final status
depends on the fixed checks; a failure would be saved without rerun.
Acceptance is limited to the allocated K-family diagnostic. The parent
retains final mathematical responsibility and no C6 ledger threshold changes.

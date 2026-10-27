# C160: normalize the joint count model with a genuine beta2 upper cap

Status: SOURCE PREPARED, UNRUN. Sol has not executed this packet, a
mathematical preflight, optimization, root or saved numerical checker.
Terra prepares it in an explicit independent root. Science requires the
independent saved audit of the complete C156 result, primary acceptance
of that audit, and separate exact-hash gates for each phase. Preparation
does not authorize either phase. No automatic phase transition exists.

## 1. Exact input and sole changed feature

Read the immutable C156 result whose SHA-256 is
`f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f`.
Let B be its complete rational value at
`budgets.finite_best.valid_true_beta2_upper`. This is the full positive
cone-weighted expression, not a displayed decimal, a single pair, or a
lower endpoint. Its validity is a dependency to be independently audited;
this packet cannot supply that audit itself.

Set D=2^96, N=ceil(D B), betaD=N/D. Integer arithmetic verifies

    (N-1) denominator(B) < D numerator(B) <= N denominator(B).

Consequently B <= betaD < B+1/D. The exact original B string, field
selector, result hash, N, D and inequalities are saved in the preparation
certificate. The model uses betaD; it never embeds B's huge denominator
in the feature coefficients. Upward quantization is conservative.

Use the complete accepted round128 `upstream.json` cap table unchanged
except beta_by_h[5][1]. The h at index5 is exactly 1/128. Feature rank2
is exactly Phi_2(n)=n(n-1)/2 for every integer n=0,...,32, so this is
E[d binom(N_(1/128),2)], precisely the C130/C156 true beta2. Its flattened
49-feature index is 1+6*5+(2-1)=32. Preparation and both phase preflights
verify this identity and compare all other entries against the pinned
round128 input. The existing beta0, all true beta1 caps, all other rank
caps, eight h values, Phi table and CDF intervals are unchanged. The
primary release must explicitly accept this complete base cap table.
No C5 upper bound, optimistic C130 lower endpoint or displayed decimal
is used. The original complete round128 simplex is the only warm start.

## 2. The cost is at most one

Keep the original body32 covariance law of rounds123--128:
Sigma=sum_i S_i v_i v_i^T, Y_i=1/(4 S_i) iid Gamma(1/2,1),
d=det(Sigma)^(-1/2), M=E d^(1/6). The nested counts are
N_j=#{i:Y_i<=h_j}. Let F0=1 and F_(j,k)=Phi_k(N_j). The accepted caps
satisfy E[d F_i]<=cap_i, including the new upward-rounded cap after C156
saved acceptance. For x_i>=0, sum_i x_i=1, x0>0, define the exact rational
coefficients a_i=x_i/cap_i and W=sum_i a_i F_i. Then

    E[d W] <= sum_i a_i cap_i = sum_i x_i = 1.

Holder with exponents6 and6/5, applied to
d^(1/6)=(dW)^(1/6)(W^(-1/5))^(5/6), yields

    M^6 <= E[dW] (E W^(-1/5))^5 <= H(a)^5.

This is a newly normalized genuine-upper-cap model. There is no C146
cost correction and no cost-one assumption based on an optimistic bL.
The x weights are dyadic; a=x/cap is generally rational, and is NOT
quantized. Rounding the evaluation of W outward on the 2^-192 grid
changes no coefficient and no cost accounting.

## 3. Complete joint law and positive tails

The complete round125/128 argument is reused. For each increasing CDF
vector p, the nine categories have probabilities p1,p2-p1,...,p8-p7,1-p8.
Their 32-input multinomial law gives the same nested counts jointly.
Nonnegative coefficients and nondecreasing Phi imply that W^(-1/5)
decreases coordinatewise. Common-uniform coupling of the exact complete
CDF laws therefore gives H(pU)<=H(ptrue)<=H(pL). Endpoint categories are
differences within each whole CDF vector, not independent intervals.

The prefix contains all 735471 tuples 0<=N1<=...<=N8<=16, at their
original subprobabilities. At K=N8>16 keep the positive intercept and
last-h six features to obtain Wtail(K)<=W(N). All16 binomial tail masses
at K17,...,32 are retained in the complete upper bound. The prefix at pU
is a lower bound on H only; prefix at pL plus the16 positive tail terms
is an upper bound on H. Coupling is applied before this tail relaxation;
no monotonicity across the prefix/tail seam is assumed. Both endpoint
laws have exact unit mass. No tail loss or probability renormalization
is permitted.

## 4. Phase A: one finite numerical selector

There is one SLSQP invocation, starting at the pinned full128 simplex.
The new midpoint-law objective uses all735471 prefix rows plus16 true
tail rows and49 features. There are at most149 optimizer objective plus
analytic-gradient calls and one final call at the new exact2^-44 simplex
candidate, total<=150. Iterations<=120, ftol=1e-12, numerical x0>=2^-24;
no wall/CPU cap, restart, second candidate search or parameter escalation.
One process, BLAS threads1, nice>=10. Float probability mass is checked
and never normalized. The positive model uses the genuine new cap, but
all phaseA objective values and comparisons remain numerical diagnostics.
They are NOT a certified M upper, global optimum or exclusion. A local
failure or exhausted budget excludes no parameter family.

## 5. Phase B: separately selected exact new-H certificate

Only a completed phaseA with a numerical improvement and an explicit
primary selection may be prepared for phaseB. A separate main plus
independent exact-hash release must then authorize its single actual run.
PhaseA never launches phaseB. The candidate is frozen before science.

The pinned128 Joint evaluator is reused as source, not as cached root
values for a different a. It forms exact rational a and outward W group
tables at S=2^192. For each new candidate state the saved integer root
certificate checks vL^5*wU <= S^6 <= vU^5*wL. Fifteen disjoint contiguous
prefix jobs cover exactly735471 states, with at most65536 cached roots
per process; the collector requests the16 new tail roots. There are
735487 new logical requests, no old functional evaluation. All integer
probabilities and positive sums have exact common denominators S^32 and
S^33 respectively. Both complete endpoint laws must have unit mass.
The table hash binds the SAME new rational a in every worker.

Save H_L,H_U and compare exact U=H_U^5 against the immutable old128
`M_power6_upper`, Theta(g0) and Theta(gstar), where

    Theta(1/q)=(8/9)^6 * 24576^5 / q^6,
    q=65000 or63500.

The q65000 value must equal the pinned128 exact threshold. The whole-I
interface remains I^5 <= (3/16)^6 U / 384^5. It has no extra cost factor.
H_L is not a lower bound on M. A strict new M-slot comparison remains
conditional on the unchanged determinant/whole-I interface; neither a
global g proof nor a new C6 bound follows from this packet alone. No new
C6 bound is registered even if the M-slot passes. An independent saved
certificate audit is a further requirement before accepted use.

## 6. Finite stop and source boundaries

Each phase has its own release, exclusive started marker, results and
receipt. First genuine error stops with files retained, no automatic
retry. Neither preparation nor either phase runs C156 again, compiles a
density, replays the old343456 functionals, evaluates heat kernels/erfc,
calls native/MC, reads live129, expands the76904685 full state law, or
changes the fixed h, prefix cutoff, arithmetic precision or caps. No
additional optimization or phaseB is authorized implicitly.

Sources read for the mathematical interface: round125 proof, round128
proof/input/source manifest/Joint/candidate/full certificate, C156 frozen
proof/result/inputs/release/receipts/source manifest, and round123
whole-I theorem. Source bytes are pinned by STATIC_SOURCE_PINS.json.
The mathematical proofs here are explicit; no unexamined literature or
generic interval-library guarantee is inferred from a skill name.

# C140: common entropy weight for two fixed positive bands

Scope: only the 11 C134 depth6 unresolved leaves inside its original
seven local boxes, at g*=1/63500. The 62 already paid leaves are only
referenced through the frozen C134 closed partitions. No original
source atlas is enlarged and no source point exponential is repeated.

## Fixed-source log concavity

For one authenticated source of positive band B_j, let c be its center,
0<L<=U its value interval, G its compact rectangular gradient enclosure,
and K=1/(32*a^2)>0 the already justified C134 frame32 constant. Define

 T(x)=inf_{f in [L,U],v in G}
       f exp(v.(x-c)/f-K*||x-c||^2).

The uncertainty set here is fixed independently of x. Each logarithm
log(f)+v.(x-c)/f-K||x-c||^2 is finite and concave in x. Compactness and
L>0 ensure that their infimum log T is finite; its hypograph is the
intersection of convex hypographs, so it is concave. C134's saved
formula is exactly this infimum: first minimize v.(x-c) at rectangular
endpoints to obtain b, then minimize f exp(b/f) at
f=min(U,max(L,b)). Its cached rational lower l(x)>0 therefore satisfies
0<l(x)<=T(x)<=B_j(x). These source inequalities are consumed from C134;
their old exponential or native certification is not rerun here.

## Self-contained finite Gibbs bound

Choose exactly one source in each of the two disjoint positive time
bands, giving T0,T1. For fixed w in [0,1], put

 H(w)=-w log(w)-(1-w)log(1-w), with 0log0=0,
 log K_w(x)=H(w)+w log T0(x)+(1-w)log T1(x).

The endpoints mean K_0=T1 and K_1=T0. For 0<w<1, weighted AM--GM on
T0/w and T1/(1-w) gives T0+T1>=K_w. Equivalently take the two-atom
reference Q=(1/2,1/2), P=(w,1-w) and V_j=log(Tj/g*). Everything is
finite, exp V_j>0, Z=(T0+T1)/(2g*) is positive and finite, P is
absolutely continuous with respect to Q, and KL(P||Q)<infinity.
Expanding KL(P||Q_V)>=0 gives

 log((T0+T1)/g*) >= H(w)+w log(T0/g*)+(1-w)log(T1/g*).

This is the finite-domain identity in the D02 skill; no general
Donsker--Varadhan theorem or unread physics paper is used. There is
no extra log2: it cancels between log Z and KL(P||Q).

For one fixed source pair and one fixed w, log K_w is concave. Every
point in the entire closed five-dimensional rectangle is a convex
combination of its 32 corners. Concavity implies that log(K_w/g*) on
the whole rectangle is at least the minimum of its corner values.
At corner i the saved positive rational lowers l0_i,l1_i yield the
reliable lower score

 s_i(w)=H(w)+w log(l0_i/g*)+(1-w)log(l1_i/g*).

If Arb512 gives a nonnegative exact rational lower endpoint for every
s_i, then B0+B1>=K_w>=g* on the entire leaf. No ordering/clipping of
the rectangle is needed. The other three positive bands may be
discarded. The same w must be used at every corner of this certificate.

This proof neither sums sources of the same band nor asserts that
their maximum has a concave logarithm. After each whole-leaf
certificate is independently valid, one may select the best certificate
among source pairs/weights. Source selection itself is never interpolated.

## Locked finite numerical packet (prepared, not run)

For each of the 11 leaves, use only its original parent's 16 source
keys per active band, hence exactly 256 pairs, 2816 pairs total. All
32 corner lowers must already occur in C134's pinned point cache;
any missing/nonpositive entry stops, rather than generating a new exp.

For each pair, a deterministic 70-call float objective search chooses
a candidate w: three endpoint/midpoint calls, two golden-section initial
calls, 64 further calls, and one final midpoint call. The objective is
H(w)+min_i[w log(l0_i/g*)+(1-w)log(l1_i/g*)], which is concave. The
search is solely a selector and makes no numerical-optimality claim.
Round its selected w to nearest n/2^24 (Python round ties-to-even),
with n clamped to [0,2^24]. Exactly that rational weight is certified.

At precision512, cache log(l/g*) once per unique source-point identity;
then cache the entropy once per exact weight. This entails at most
8544 Arb log evaluations for the exact selected unique source points
(the geometric upper count was11264) and 5632 for entropy, 14176 total.
Endpoint entropies require no log. There are at most
90112 new Arb corner score combinations, 2816 pair certificates,
197120 float objective calls, 382976 selector entropy log/log1p calls
and 8544 diagnostic float source logs, counted separately from Arb logs.
There are zero native, heat-kernel, source-point exp, final exp,
source-jet generation, refinement, or old paid-leaf score calls.

Both Arb and float operations start only after exact-hash main plus
independent-review release. Execution is a single process at nice>=10,
with one-thread environment variables and no wall-time cap. Saved-only
preflight performs only JSON/hash/Fraction identity checks and zero logs.
Any structural/hash/count failure stops; no automatic restart or budget
extension is permitted. Even 11 successes establish only these local
leaves and, using C134's remaining paid leaves, the original seven local
boxes. They do not establish global g*, an actual G upper, or a new C6
bound. A failure does not exclude other weights or source families.

Skill read: /Users/zhengzhihao/.codex/skills/math-d02-gibbs-entropy-duality/SKILL.md
and references/provenance.md. The required review/release is imposed by
the parent C140 assignment, not inferred from the skill.

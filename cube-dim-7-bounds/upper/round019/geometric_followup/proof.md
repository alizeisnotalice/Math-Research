# R019: three threshold-specific finite-parameter mechanisms fail

Date 2026-10-27. Native Sol, pure theory, zero numerical runs. Responsibility is this new directory only; the R014 files and other projects remain untouched. No commit or push. No new literature lookup was performed. C1=(11+sqrt(61))/12 is used from the already checked R014 source receipt.

## Result and operator

For Lebesgue L1 on R^n, M_n is the centered axis-parallel cube maximal operator, with one common radius r>0 in all coordinates and with all positive r allowed. We examine single-threshold mechanisms for

    lambda |{M_n f>lambda}| <= B ||f||_1,   B<2^n.

None of the following three mechanisms supplies such a bound. Each has an explicit exact counterexample. In particular no new C7 upper bound is claimed: the accepted bound is still 128. These counterexamples do not disprove the desired operator estimate C7<128, or the earlier target C7<=64C1.

This is not a replay of the excluded global L1 factorization or arbitrary-bad-cube-selector arguments: mechanism 1 is allowed to depend explicitly on lambda through the global mass radius; mechanisms 2 and 3 use threshold mass radii of source clusters and attempt a finite-cost coalescence cover. All their failures occur at one fixed threshold.

## Mechanism 1: finite global scale bands followed by vertical Melas

### Proposed rule and finite parameters

Write d=n-1 and x=(x',t). Fix a in (0,1] and K>0. At a threshold lambda>0 define

    R = (||f||_1/lambda)^(1/n)/2,
    h_(f,lambda)(x',t)
      = K sup_(aR<=r<=R) (2r)^(-d)
          integral_(||x'-y'||infinity<=r) f(y',t) dy'.         (1)

Here f>=0. The cap r<=R is legitimate for actual bad cube witnesses: if r>=R, their average is at most lambda. The proposal is to replace all smaller transverse scales by the lower scale aR, or a finite band of scales no smaller than aR, then use the vertical one-dimensional maximal operator M_v. If

    {M_n f>lambda} subset {M_v h_(f,lambda)>lambda} a.e.,     (2)

were true and ||h||_1<=B0||f||_1, Melas would yield B=C1 B0. Formula (1) is integrable: bounding every transverse average by the radius-R integral divided by (2aR)^d gives ||h||_1<=K a^(-d)||f||_1. In particular parameters with C1 K a^(-d)<128 would have been a concrete finite-parameter candidate in dimension 7.

### Counterexample for every a>0 and every finite K>0

Set lambda=1, choose p in (0,1/2) so small that

    p^(d/n) < a^d/(32K),

and put

    s=p^(1/n)/2,   epsilon=eta=s/16,
    b_epsilon=1_[-epsilon,epsilon]^n/(2epsilon)^n,
    f=2^n p b_epsilon + 2^n(1-p)b_epsilon(.-4e1).

The coordinate e1 is transverse, not vertical. Thus f is an explicit compactly supported L1 function, ||f||_1=2^n, and R=1. Consider the positive-volume receiver box

    V = {||x'||infinity<eta, |t-s|<eta}.

The genuine cube at each such point of common radius r0=s+epsilon+eta=9s/8 contains the small source bump. Therefore

    M_n f >= 2^n p/(2r0)^n = (16/9)^n > 1  on V.            (3)

The distant bump does not meet any transverse radius-r averaging box in (1), since r<=1. For every x' in the displayed receiver box,

    h_(f,1)(x',u)
      <= K [2^n p/(2a)^d] b_epsilon^1(u),

where b_epsilon^1 is the unit-mass interval bump. If |t-s|<eta, any centered interval meeting its support has radius at least t-epsilon>=7s/8. Its average is thus at most 1/s. Consequently

    M_v h_(f,1)(x',t)
       <= 4K a^(-d) p^(d/n) < 1/8  on V.                  (4)

Equations (3)-(4) disprove (2) with a positive-volume margin. The chosen witness cubes all have genuine common positive radii. The proof also excludes replacing the supremum in (1) by any finite subset of its bands. There is no issue of a weak-type composition estimate here: the required threshold inclusion itself is false.

**Diagnosis.** A positive lower transverse scale determined by total mass misses a sufficiently small isolated component when most mass is far away. Merely tuning the finite parameters a,K, even without a useful mass budget, cannot repair this rule. A successful threshold-specific construction must adapt to local mass across arbitrarily small scales.

## Mechanism 2: overlap coalescence of scalar threshold halos

### Proposed rule and potential cost

Fix alpha in [1,2). For a source cluster of mass m use its threshold mass radius

    rho(m) = (m/(2^n lambda))^(1/n)

and a halo cube of radius alpha rho(m) centered at that cluster. Start with the individual source atoms, merge clusters only when their current halos intersect, and stop when the halos are pairwise disjoint. Any mass-preserving output family of this form would have total halo volume

    sum (2 alpha rho(m_j))^n = alpha^n ||f||_1/lambda.

Thus a cover of the true level set by these terminal halos would directly give C_n<=alpha^n<2^n. The rule is threshold-specific and its source masses, scales and test for merging are explicit.

### Corner family defeats the rule before the first merge

Set lambda=1. Choose

    1/2 < rho < 1/alpha,

which is possible exactly because alpha<2. At each of the 2^n vertices y_sigma=sigma of [-1,1]^n put mass

    m=(2rho)^n.

Every initial halo has radius alpha rho<1. Distinct vertices differ by 2 in at least one coordinate, so these halos are pairwise disjoint with a positive gap. The rule therefore makes no merge, and its terminal halos miss the origin and a whole receiver neighborhood of it. Nevertheless the genuine radius-one cube centered at the origin contains all sources and has average

    (2^n m)/2^n = (2rho)^n > 1.                             (5)

There is a positive-volume failed receiver neighborhood: choose its radius eta smaller than both 1-alpha rho and 2rho-1. At a point in that neighborhood, a genuine common cube radius 1+eta still has average greater than 1, while all terminal halos miss that point after eta is made slightly smaller if needed.

This also has a direct L1 version. Replace each atom by a unit-mass bump of sufficiently small cube support, multiplied by m. Use a common witness radius 1+eta+epsilon. The strict margins alpha rho<1 and 2rho>1 preserve both the lack of any merger and the failed receiver set. No numerical atomic-to-L1 assertion is used.

**Diagnosis.** Source-cluster halos can be mutually disjoint while their masses cooperate inside a larger bad cube. Hence overlap-only coalescence with any scalar expansion alpha<2 cannot prove an improvement over 2^n. Allowing a new merger criterion beyond overlap is a new mechanism; it was not assumed here.

## Mechanism 3: mixed-aspect halos with the exact 64C1 budget

### Proposed rule

In n=7 let beta=C1 and set the fixed aspect vector

    b=(2,2,2,2,2,2,beta).

A mass-m cluster at threshold lambda is represented by the box

    H(c,m) = c + product_j [-b_j rho(m), b_j rho(m)].          (6)

This is a covering/auxiliary box, not a replacement of the target operator's common-radius cubes. Its volume is exactly

    |H(c,m)| = (64 beta) m/lambda.

Consequently, if intersecting source boxes could be merged into a single box of the form (6) with the summed mass, retaining all previously certified bad receivers, the final family would preserve the hoped-for budget 64C1. In particular this would give a literal finite-parameter implementation of the intended one-coordinate replacement, rather than an unsupported multiplication of weak norms.

### Two equal sources already invalidate the needed merger

Write q=2^(1/7). Take lambda=1, fix rho>0, and place two equal masses m=(2rho)^7 at vertical positions

    y_+=+(s/2)e7,   y_-=-(s/2)e7,   s=beta q rho.

Their initial boxes (6) intersect: their transverse centers coincide and

    s<2 beta rho,

since q<2. A coalescence rule for intersecting boxes must merge this pair. The combined mass 2m permits only radius rho(2m)=q rho in (6), so every possible new box, whatever its center, has vertical projection length

    2 beta q rho.                                           (7)

But the true threshold-one level set already contains each atom's open centered cube of radius rho: inside that cube the corresponding single atom alone has maximal average greater than 1. Since beta q<2, these two receiver cubes have overlapping vertical projections, and their union's vertical projection is the interval of length

    s+2rho = beta q rho+2rho > 2 beta q rho.                  (8)

Their horizontal receiver sections agree and have positive 6-dimensional volume. Hence no vertical interval of the permitted length (7), and therefore no new box (6) with the summed mass, can retain these genuine bad receivers up to null sets. The missed set has positive 7-dimensional volume. Choosing a different cluster center cannot solve a projection-length obstruction.

The inequality beta q<2 is exact, without decimal calculation: beta<8/5 follows from sqrt(61)<41/5, and q<5/4 follows from (5/4)^7>2. These two strict inequalities give beta q<2.

For a fully L1 counterexample, replace the two atoms by small equal cube bumps. The mass radii of the proposed source boxes are unchanged. Choose delta>0 with 2(1-delta)>beta q and then support half-width epsilon<delta rho. The genuine bad set contains the two receiver cubes of radius (1-delta)rho, using common witness radii (1-delta)rho+epsilon<rho. Their vertical projection remains longer than (7). This preserves the contradiction with a positive margin.

**Diagnosis.** The one-dimensional weak norm C1 is not a vertical geometric dilation factor for an n-dimensional source-mass merger. Already for one equal-mass pair, retention of the two genuine single-source halos would require beta>=2/q=2^(6/7). The proposed beta=C1 fails that necessary condition. This is a failure of this specific coalescence rule, not an impossibility theorem for every mixed-aspect method or every threshold slice construction.

## Minimal open gap and stopping condition

The three finite mechanisms have been fully tested analytically; all fail before any experiment is useful. No fourth mechanism or numerical batch is opened.

A viable threshold-specific slice proof must jointly supply the following missing operation: it must combine arbitrarily small local source scales with cooperative mass from separated sources, while preserving true receiver coverage and a uniform total payment below 128||f||_1/lambda. The local scale cutoff in mechanism 1 loses the small component; the overlap test in mechanism 2 misses cooperation; the additive-mass fixed-aspect replacement in mechanism 3 loses already genuine receivers.

The next finite acceptance should therefore be one explicit, threshold-dependent merger/transport rule with a proved receiver-retention invariant and a proved payment invariant. Before any integration with Melas, it must pass the three exact families above. An argument that simply assumes the level-set volume budget, uses an arbitrary supplied bad-radius selector, or constructs one global L1 majorant for all thresholds does not fill this gap. No such rule has been proved in this branch, and the C7 upper ledger remains 128.

## Execution record

Only reading the clone status, creating this owned directory, writing this report and receipt, and ordinary file-integrity checks were performed. No MATLAB, Python/SymPy/Arb calculation, numeric parameter scan, new literature search, code execution for mathematical results, commit or push. Every counterexample and inequality above is an analytic proof and is submitted for the parent's independent audit.

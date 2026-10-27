# R014: one-coordinate Melas replacement — two obstructions and an open local interface

Date: 2026-10-27. Independent Sol theory branch. No numerical calculation, MATLAB, experiment, commit, or push was run. The clone starts at `14f3507f406185667f84146328f1db9da3f9949a`; only this new directory is authored here. The parent retains proof and integration responsibility.

## 1. Scope and result

For n >= 2, write x=(x',t) in R^(n-1) x R and

    A_r^n f(x) = (2r)^(-n) integral_[x-r,x+r]^n |f(y)| dy,
    M_n f = sup_(r>0) A_r^n f.

The same r is used in every coordinate. The input class is general Lebesgue L1, and all positive scales remain available. Let C_n be the optimal weak (1,1) constant. The desired assertion is

    C_n <= 2^(n-1) C_1.                                      (T)

Melas's original Annals paper, Corollary 1, printed p.649, confirms

    C_1 = (11 + sqrt(61))/12.

Thus (T) for n=7 would give the genuine operator bound

    C_7 <= (16/3)(11 + sqrt(61)) < 128,

where the strict comparison follows exactly from sqrt(61)<13. No decimal evaluation is used. **This branch neither proves nor disproves (T).** It rigorously excludes two natural sufficient lemmas. The accepted C7 upper bound therefore stays 128; neither model integral costs nor the results below change it.

Aldaz arXiv:1811.10478v2, Lemma 2.6, describes the genuine descending-radius center-cover selection. Its final estimate integrates the selected overlap against the input. Theorem 2.5 and Corollary 3.5 identify the uniform-over-measures constant as 2^n for the sup norm. That sharpness over measures is not Lebesgue sharpness. The following obstruction is directly Lebesgue and does not rely on that distinction being blurred.

## 2. Route A: replacing the pointwise overlap load is false

### Finite family

Fix a=3/4. For each sign vector sigma in {-1,1}^n let

    c_sigma = a sigma,     Q_sigma = c_sigma + [-1,1]^n.

These 2^n cubes form an intersecting Besicovitch family: any two centers differ by 3/2 in at least one coordinate, so neither lies in the other's radius-one cube. Their common intersection is [-1/4,1/4]^n, of positive Lebesgue volume. No selected subfamily can cover all the centers unless it contains every cube. Even fractional weights w_sigma >= 0 with

    sum_sigma w_sigma 1_Qsigma(c_tau) >= 1  for every tau

must satisfy w_tau>=1. Consequently their overlap is at least 2^n on the common intersection. This disproves every universal replacement of the selected L-infinity overlap bound 2^n by 2^(n-1)C1, since C1<2.

This is compatible with bad averages for an actual L1 input: take f to be a nonnegative unit-mass bump supported in the common intersection, and any 0<lambda<2^(-n). Every displayed cube is a valid bad witness of the same radius, and

    integral f sum_sigma w_sigma 1_Qsigma >= 2^n ||f||_1.

Thus changing an unweighted multiplicity argument into a fractional or f-weighted estimate does not repair this particular lemma.

### Positive-volume receiver version, including arbitrary fractional selection

The obstruction is not an artifact of requiring cover only at finitely many centers. Fix 0<delta<1/8,

    E_sigma = c_sigma + [-delta,delta]^n,
    E = union_sigma E_sigma,
    K = [-1/8,1/8]^n,
    f = 1_K / |K|.

Supply, at every x in E, the prescribed witness Q_x=x+[-1,1]^n. Every Q_x contains K, because 3/4+delta+1/8<1. Hence A_1^n f(x)=2^(-n) for all x in E. For any two different sign clusters, Q_x with x in E_sigma misses E_tau: in an opposite coordinate the distance is at least 3/2-2delta>1.

Let nu be **any** nonnegative finite Borel measure on E (possibly chosen using f and lambda), and suppose its fractional cube load

    L_nu(z) = integral_E 1_Qx(z) dnu(x)

satisfies L_nu>=1 almost everywhere on E. On E_tau one has L_nu(z)<=nu(E_tau), so nu(E_tau)>=1 for every tau. Since all Q_x contain K,

    L_nu(y)=nu(E)>=2^n  for y in K,
    integral f L_nu >= 2^n ||f||_1.                           (A)

All objects are Lebesgue measurable; Tonelli justifies the load integrals. This is an exact analytic L1 counterexample, not a numerical atomic surrogate. It uses one common radius and therefore remains admissible when all scales are allowed.

**Precisely excluded:** a uniform theorem that starts from an arbitrary supplied bad-radius selector, covers an arbitrary bad subset by those same cubes (even fractionally), and bounds the paid source load by 2^(n-1)C1. A successful theorem could choose additional witness radii, use the geometry of the full level set, or retain a quantitative receiver-volume deficit. Assertion (A) does not exclude those possibilities and is not a counterexample to (T).

## 3. Route B: every global finite-mass vertical factorization is false

Let M_v be the one-dimensional centered interval maximal operator in the t variable, with x' fixed. The following seemingly powerful sufficient interface is impossible, even without asking that h depend linearly, positively, locally, or measurably as a map on f:

    For every f in L1_+(R^n), there is h in L1(R^n)
    such that M_n f <= B M_v h almost everywhere,              (B)

where B is any finite positive constant. In particular no norm bound such as ||h||_1<=2^(n-1)||f||_1 can make (B) true.

**Proof.** Take f=1_[-1,1]^n and put d=n-1. For x' with R=||x'||_infinity>=1 and |t|<=R, the genuine cube of common radius r=R+1 contains [-1,1]^n. Therefore

    M_n f(x',t) >= (R+1)^(-n).

Suppose h satisfying (B) exists. By Fubini, for almost every x', the section h_x'(t)=h(x',t) lies in L1(R), and (B) holds for almost every t. Set

    beta_R = [2 B (R+1)^n]^(-1).

Then {M_1 h_x'>beta_R} contains [-R,R] up to a null set. Applying the established one-dimensional weak inequality gives

    C_1 ||h_x'||_1 >= beta_R (2R)
                     = R/[B(R+1)^n]
                     >= [B 2^n R^d]^(-1).

But integral_(||x'||infinity>=1) ||x'||infinity^(-d) dx' is infinite. For example, every dyadic sup-norm shell contributes a fixed positive amount. Integrating the last estimate contradicts h in L1(R^n). QED.

The argument needs only some finite weak (1,1) bound for M_1; its role is structural rather than dependent on the sharp Melas value. It also rules out the stronger scale-by-scale inequality A_r^n f<=A_r^1 h. The usual H(x',t)=M_(n-1)^horizontal f(x',t) construction fails for the same reason: H has infinite total mass for this f. Consequently a pointwise composition with a weak-(1,1) operator cannot justify multiplying weak constants.

**Precisely excluded:** a single global finite-mass h for every f, valid simultaneously at all points and levels. This does not exclude h depending on lambda and needing to capture only one level set. It also does not exclude an exceptional-set/obstacle method. The proof never replaces common cube radii by independent coordinate radii in the target operator.

## 4. Remaining level-specific interface and finite acceptance condition

No third proof route or numerical batch is opened. To make the remaining gap explicit, the following lemma would suffice and escapes Route B's quantifiers:

    For every n>=2, f in L1_+(R^n), lambda>0,
    compact E contained in {M_n f>lambda}, and epsilon>0,
    there exists h_(f,lambda,E,epsilon) in L1_+(R^n) with
      ||h||_1 <= [2^(n-1)+epsilon] ||f||_1
    and
      E contained in {M_v h>lambda} up to Lebesgue null sets.   (L)

Fubini and Melas give lambda |E| <= C1 ||h||_1. Let epsilon decrease to zero, then use inner regularity on the measurable level set, to obtain (T), for general L1 and all positive scales. Nonnegative reduction is exact because M_n uses |f|. This is a sufficient interface only; it is **unproved** and is not being advertised as a new tensorization theorem.

The next finite acceptance task is analytic: produce one explicit construction of h in (L), prove both displayed requirements without assuming the desired bound on |E|, and check it on the unit-cube input and the separated 2^n receiver clusters above. Any construction that produces the same h for all levels contradicts Route B. Any construction based solely on fractional selection of arbitrary prescribed witness cubes contradicts Route A. No additional atomic search is needed to verify these two failures, and no new numerical experiment is requested by this branch.

## 5. Sources and public/internal boundary

* A. D. Melas, *The best constant for the centered Hardy-Littlewood maximal inequality*, Annals of Mathematics 157 (2003), 647-688, DOI 10.4007/annals.2003.157.647. Original publisher PDF: https://annals.math.princeton.edu/wp-content/uploads/annals-v157-n2-p08.pdf . Read abstract, definition (1.1)-(1.3), Theorem 1, Corollary 1, formula (1.8), Theorem 2, and the introductory explanation of the discrete covering proof. Corollary 1 is the exact general-L1 source; Theorem 2 also covers finite Borel input measures relative to Lebesgue receiver measure. The whole 42-page proof was not independently re-proved. PDF extraction was read through the web tool; binary retrieval timed out, so no PDF binary hash is claimed.
* J. M. Aldaz, *Kissing numbers and the centered maximal operator*, arXiv:1811.10478v2, 2020-03-10. Read definitions 1.1/2.2, Theorem 2.5 and full proof of Lemma 2.6, Corollary 3.5. Cached PDF SHA256: `d42a1d6518298227e2af47876761d47225e38b28e2f7726f91c0cc2498ab12e7`. The published source statement gives the uniform measure-class constant; the present Lebesgue obstructions are separate arguments.
* Lin author repository HEAD was checked as `72c022ba09032b522c878cca2f9b6c68217d77ea`. The original pinned README states c2<=3879/1000; its exact kernel expression is R^2/(R-1)-epsilon/6 with R=3797/2000 and epsilon=159411/200000. The same pinned `Challenge.lean` was successfully read in full through the web tool: it explicitly defines Lebesgue centered cubes with one common r and states `weakTypeConstant_two_le_upper`. That is a statement file using `sorry`, not this branch's proof verification. Its companion Solution/dependency proof was not audited or rebuilt. This is a public two-dimensional claim, not a discovered seven-dimensional bound; multiplying it by 32 has no justification from the two excluded routes. README SHA256: `387c56f90e42e7d6ffd31367f8daad770d973b2281b3bd8e2e013390d7f1f64a`.

No stronger public explicit finite-n7 bound was verified. The lookup was narrowly targeted and is not an exhaustive literature claim. The two obstruction proofs are internal research results of this branch and are not claimed as published or novel. The fixed accepted ledger remains strict C7>2.136578727666 and C7<=128; the lower branch was not touched.

# R048 manuscript: independent mathematical source review

2026-10-28. **ACCEPT. No blocking mathematical error, definition gap, or quantifier change was found in the locked source.** This accepts the fixed-d eventual theorem and its closed endpoint, while retaining the independent explicit R042 theorem. It does not certify any R049 threshold or public-priority claim.

Audited TeX: `/Users/zhengzhihao/Desktop/T/cube-cover-manuscript/cube-dim-7-bounds/paper/main.tex`, SHA-256 `ecd7fd7487083b7eb6318a0a94f969f9d291b9d4df443cc3374497a7f103b248` (1982 lines). The main coordinator retains final integration/proof responsibility. No author file was edited. This report and its small exact-check attachment are the only outputs; no Git, n=7 work, MATLAB, production computation, author verification script, 13-check transcription rerun, or compilation was performed.

## 1. Read scope and preservation of the accepted manuscript

The entire new mathematical content was read directly, including lines 566–902 (the six new subsections), the abstract/introduction/theorem statements, and lines 1644–1705 (complete-cost comparison). The shared actual-measure and energy definitions at lines 508–565, original finite-selection eligibility/cover interface, and cap/scalar bounds used in the new argument were also read directly. This is a mathematical audit of the actual new proof text, not acceptance merely because it cites an earlier reviewed mechanism.

The accepted 94c5 source was recovered from the writer's `upper/round048/manuscript_notes/private/main_94c5.tex` and its hash actually verified as `94c5e61d3bf86b16736e328ffe3e412d46df957dd2b8f8760939edc8831e8692`. A full line-by-line diff was inspected. It consists only of abstract/introduction synchronization, the fixed-d theorem/corollary, the six new proof subsections, the complete scalar-comparison subsection, and adding the new theorem to the minimum-of-applicable-bounds sentence.

In particular old lines 540–1282 equal new lines 901–1643 exactly: the entire R042 proof (row, column Y, fixed-face second moment, high-column truncation, finite-a bridge, scalar bound and explicit power/Taylor argument) and old energy/explicit-threshold proof remain intact. Old lines 1283–1383 equal new lines 1706–1806 exactly (R034 constant-family proof); old lines 1385–1559 equal new lines 1808–1982 exactly (old complete-root comparison, R037 actual-history restriction, fixed-level spike qualification, leading-coefficient restriction and bibliography). Thus no accepted R042/R034/R036/R037 argument was silently removed or migrated to the new eventual domain. The retained R042 theorem still explicitly requires L>=2^108. The historical fixed-level spike assertion still does not become a claim that the full weak norm tends to one.

## 2. Fixed-d domain, exact local integral, and pinned coordinates

Lines 105–124 correctly quantify: fix any real 0<d<=2^-24 first, then there exists N(d) working for every integer n>=N(d). The closed upper endpoint is a member of that domain. Neither a uniform-in-d start nor a varying-d statement is asserted. Abstract, main theorem, corollary and proof consistently give the endpoint saving 2^-41 n and an unspecified operator start.

The new local formulation in lines 608–706 is an exact rewriting of the R045/local-lemma integral, not an approximation. At a source radius rho=tq, set R=q+d-h, a=(1-epsilon)q and g=R-a. The two deleted-end lengths are m_+=(w-a)_+ and m_-=(-w-a)_+; at most one is positive. On a tail the change of variable s=g-(rho-|z|) makes the integral of ell_R^-1-1 equal Psi(R,s)=-2R log(1-s/(2R))-s. Deleting length m leaves Psi(R,(g-m)_+). Hence the stated exact formula for J retains both end corrections and the ordinary deleted length, including the case m>g. Since r-R<=2h, m<=g+2h; no deletion reaches the opposite tail.

The displayed Psi_R/Psi_s bounds give the stated a.e. derivative bound with constant 64. The proof expressly handles the real kink at |w|=a, uses absolute continuity across it, and identifies the C1 join at |w|=R. Consequently integrating the derivative gives a bound for **all** matched q, rather than only almost every q. Epsilon=zeta*t is fixed during this differentiation. The q0=r-d+h baseline may lie outside the source band, but the whole extended interval [q0-2h,q0] is defined and satisfies the local hypotheses eventually at each fixed d. The manuscript states that extension explicitly.

At R=r the captured fraction for a fixed source coordinate is exactly ell_r(z). Tonelli gives E J=1 with W uniform on [-r,r]. Central/tail estimates yield E(J-1)^2<=3g^3/rho0^2, and log inequalities/Jensen give the -6/12 mean/variance constants. The target maximum-face coordinate is instead fixed at +/-r and has its separate deterministic log bound [-g/rho0,0]. It is never incorrectly given a free-coordinate law.

On each target face the free log mean/variance sums tend to zero uniformly in r, since nd^3/rho0^2=O_d(L^2/n). Kw=1+Bin(n-1,2d/r) has a 4nd cutoff with failure tending to zero; no independence between the two good events is assumed. The full target annulus has mass at least nd/4 and its excluded original corner cap is exponentially small, so lambda(A)>=nd/16. Summing and integrating the a.e./AC derivative bound gives 1280Ldh=5/4096<1/8 and therefore product J in [1/2,2] throughout the matched interval.

## 3. Same weight in both degrees and inverse-tail payment

The pointwise column bound at lines 597–605 uses all 2n maximum faces, the factor 1/2 per face, deleted-product factor at most two and target log-width at most 2h/R. Normalized target overlap decreases with target radius because kappa*rho>|z_i|. This gives C_raw<=4nh G/R and b C_raw<=4k for b=min(1,theta/G), theta=exp(-12Ld^2), k=nh theta. The weight is Borel, positive, <=1, and remains in the post-selection kernel; it is not inserted into algorithm weights or an altered ambient source measure.

The exact inverse-row angular coefficient includes pinned source faces:

    U=(product J)/(2 ell_*) sum_i N_i/J_i,
    ell_*=1-g/(2R), N_i=1+1_{|w_i|<a}.

The strict inequality in N_i correctly counts captured faces at the endpoints |w_i|=a too. Bounds J in [1/2,2], product J in [1/2,2], ell_* in [3/4,1] give n/8<=U<=8n. The ordinary radial match then supplies the factor h; U itself is not a radial integral.

For the clipping lower bound the probability is normalized by capture/G on each captured source face. The free-coordinate density is capture/(ell_R J), with J>=1/2; the fixed coordinate contributes -log ell_*<=d. Because actual R>1 (stronger than the standalone local lemma's R>=1/2),

    E_inv[-log ell_R] <= 2 E_full[(-log ell_R)/ell_R]
                          <= 8g^2/(9R rho) <= d^2/rho.

The complete inverse expectation of -log G is at most d+2Ld^2<=3Ld^2, so Markov loses at most 1/4 at G<theta. This proves the **weighted** angular row lies between 3n theta/32 and 8n theta. The manuscript does not substitute a raw-source tail probability for this inverse probability.

After clipping, b<=1 permits deletion of the original corner-captured sources by their unweighted total mass. Relative to k this is bounded by a constant times (d/h) exp[-L(d/32-12d^2)], tending to zero for every fixed allowed d. Thus row>=3k/64, row<=32k, column<=4k hold under the same bounded kernel. There is no unpaid retention of old unweighted row energy after changing the column factor.

## 4. Actual history, finite a, and operator scope

The unchanged original measure is mu_a<=lambda_a. Both new supports are in D, are radially disjoint, and lambda_a=c_a lambda there exactly, c_a=alpha^-n>=1. A target actual point in increment j has cr_j<=||u_*-w||<alpha+r<3; a source in increment i has r_i>=rho>4. All source indices in the kernel are therefore earlier than **that same** j.

For every a>=8n^4, eventually a>2t and

    alpha beta-kappa-2/(a+1)=zeta-(beta+2)/(a+1)>0.

The source/target increment-to-centre distances are bounded by their respective r/(a+1); matching capture therefore implies y_j in the original bQ_i. Since the new b(z)<=1, its weighted source contribution is at most the original increment weight, so the original eligibility inequality pays row_mu<=tau. There is no theta^-1 payment and no nonconcentration assumption on mu_a. The two uses of b (original dilation and bounded source factor) are explicitly distinguished in the text.

The energy lemma's M>0 assumption and disjoint supports remain explicit. With M=32k, row lower 3k/64, tau<=3k/128 and target mass nd/16, the ambient deficit is at least 3nd/65536. Under finite-a scaling, energy scales by c_a^2 and degrees/target mass by c_a; the numerator is c_a^2 E-tau c_a a_*, so its divided lower bound can only increase. Keeping D_d=nd/65536 is conservative. This subtracts a deficit from the full actual uncaptured load, not just from a hypothetical continuum history.

The original finite witness integral, selected radii and common L1 input are retained. Adding captured/self terms gives the stated finite-a scalar bound. Taking only a scalar limit at fixed n,t,tau is legitimate even though the finite families may depend on a. The earlier strict-level compact cover and inner regularity then cover general L1 inputs, and all positive common radii are allowed. No common selector limit or family convergence is asserted.

## 5. Full scalar cost and closed endpoint

For x=(1-32d^2)L in [L/2,L], the accepted cap estimate gives B>8n eventually, so B-D_d>0. Bounds A<=n exp(-32d^2L), B<=nL imply tau/k<=(sqrt L/h) exp(-4d^2L)->0 at every fixed d. The required tau is consequently legitimate within the same N(d).

Expanding the complete root cost, not just its B term, gives

    F(A,B-D_d)<=T-D_d+64nd^2+o_d(n).

At all allowed d, 64d^2<=(d/65536)/4; equality holds at d=2^-24. The remaining o_d(n) can be made smaller than D_d/4. This proves the claimed saving D_d/2 and hence 2^-41 n at the closed endpoint. It does not need an unquantified limit uniform over d. The two uses of d in x and D are recomputed together at the endpoint.

Lines 1644–1705 compare full formula costs at their respective t values. Their lower envelope F(A,B-Delta n)>=T-n(4/L+Delta) has the correct direction and a positive root argument. Their upper envelope retains 2sigma n, 1, A, and the square-root cross term. Taylor positivity supplies the global decreasing tail E_sigma(L); at L0=2^240 the endpoint tail is exactly 2^-479+2^-393+2^-270, below its available margin 2^-41-2^-42. The R042 lower-envelope bracket uses its actual budget 2^-31 L^-1/3 (not its displayed half-saving), giving 2^-238+2^-111 at L0. The old d=2^-30 formula similarly uses its actual 2^-46 budget. Both strict full-cost orderings follow with all their positive terms accounted for.

Crucially L0 is explicitly a **formula-sorting threshold** only. The manuscript says it does not prove N(2^-24)<=exp(L0). Before its own new operator domain starts, R042 and other old applicable bounds retain their role. The minimum corollary still preserves each old theorem's domain and its complete functions; no new constant is promoted to global optimality, no lower bound for C_n is inferred from a formula lower envelope, and no R049 candidate threshold appears.

## 6. Source locks and actual local checks

All hashes below were recomputed from actual files, not copied from an author receipt:

| Dependency | SHA-256 |
|---|---|
| R048 TeX | `ecd7fd7487083b7eb6318a0a94f969f9d291b9d4df443cc3374497a7f103b248` |
| Accepted 94c5 manuscript snapshot | `94c5e61d3bf86b16736e328ffe3e412d46df957dd2b8f8760939edc8831e8692` |
| R045 complete geometry family | `add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9` |
| Independent one-coordinate lemma | `3fa91ce60a9ef05df270d1fa85a314de6a8573f8223fe1ea69f8433d6b6caa39` |
| Endpoint/complete-cost comparison | `ba7498d1ccdc2c9158b065ade728c60f7b83cc2706ce517580d866903451d2f9` |
| Independent comparison audit | `95e058c60d3856d41ab022d1cef414c7e6d0f197cb18446d78e759ddb96c8a86` |

`microcheck.py` was actually run locally with Python/Fraction; PID 82931, all eight new exact gates PASS. It checks the endpoint parameter/budget identities, scalar-quarter equality, bridge and cap margins, and three complete-cost dyadic gaps. It also checks the preserved old proof's exact 743-line equality and records the complete diff coordinates and actual source hashes in `microcheck.json`. The first preservation assertion used a mistaken line offset and was corrected to the exact diff mapping before the successful run; this was an audit-script indexing issue, not a manuscript or mathematical change. Finite arithmetic only corroborates the constants; the general-dimensional conclusions above use the analytic tail/AC/fixed-d arguments.

There are no requested author corrections for this locked version. PDF compilation/layout remains outside this review's scope.

# R054 manuscript: independent full mathematical review

2026-10-28. **ACCEPT for the final source and complete mathematical proof chain.** The manuscript establishes the simple-tuple original-operator bound for L=log n>=2^27 and the raw-tuple bound and strict full-cost ordering for L>=2^32. No theoretical, probability, threshold, finite-history, or operator-scope error was found. The one minor self-contained-definition clarification was incorporated by the writer and independently checked as a single local delta. Final locked TeX SHA-256: `17e3de43f0d841a285dbfcaf8ff12c57f9670270aab9ce95cf7fe5ce359b9554` (1340 lines).

The coordinator retains final integration/proof responsibility. This review is of mathematical TeX source. No PDF compilation/layout audit, author 51-check program, production computation, R052 search, old grid, n=7 work, shared MATLAB or Git operation was performed. Only this private `audit/round054/manuscript/` directory was written.

## 1. Source and read scope

The entire 1325-line source was read, including the definitions, preliminaries, original finite-cover proof, cap formula, actual-increment bridge, local J/Psi proof, generic budget theorem, both complete explicit threshold proofs, full-root comparisons, retained low-domain R034 theorem, real-history barrier and fixed-level spike qualification, and bibliography.

Initial locked source: `/Users/zhengzhihao/Desktop/T/cube-cover-manuscript/cube-dim-7-bounds/paper/main.tex`, SHA-256 `a05aecf23e8cc3c55e65afbe38845709e58f30fe61763d248ce8609159860cb0`.

The new proof was compared mathematically with the accepted R051 source `aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6` and its full audit `66355a1d58d3b6a832e844645b9e61ca9e3c5a2f27bc0958947054aa9144f960`, R053 threshold `26adb4e31f0026637ea17d8f7943bc2933dce97b8655e121f2604b6ddeea2775`, and R055 threshold `5d2c61adf0e341b2a11daabec61daedd6e009b65a0139bb468711175f0288e9f`. This is direct review of the manuscript's proofs, not acceptance merely because those sources were already audited. `sourcehash.json` records actual absolute paths and hashes.

## 2. Original finite cover and the actual missing measure

The object throughout is Lebesgue measure on R^n, centred axis-parallel cubes with one radius shared by their coordinates, arbitrary L1 inputs and all positive radii. The strict-level open-set/inner-regularity lemma is sound, including initially infinite level measure. The finite selection keeps actual increments and original weights fixed; increasing nonnegative eligibility sums imply selected radii are nonincreasing. Every unselected cube fails eligibility before any smaller selected radius, which gives the stated finite-cover factor.

For a point with positive load, the last selected positive-contribution item among the minimal-radius contributors is correctly chosen. Normalizing its radius to one preserves weights. Captured predecessors cost at most tau; the terminal item costs at most the original dilation to power -n. Every uncaptured contributing predecessor has radius in [1,t), and its actual increment is contained in the displayed normalized cap. Disjoint increments give the density domination with no independent-source assumption.

Coordinate monotonicity of the cap is proved on full faces: the positive face is always captured, the other face is counted by its strict indicator, deleted-coordinate products decrease with the corresponding centre coordinates, and the exceptional radial equalities have zero radial measure. The finite corner integral and uniform scalar limit therefore reproduce A_n(t),B_n(t), including n=1 and positivity of B. The witness integral uses the original selected a-dilates of the same general f; no support or scale restriction is introduced.

The actual uncaptured measure mu is dominated by lambda_a. Both new supports lie in D where lambda_a=alpha^-n lambda exactly. The energy identity uses one missing measure sigma=lambda_a-mu and disjoint supports, so its loss is at most M[sigma(A)+sigma(P)]<=M sigma(total), with M>0 explicit. The direction of the actual-energy payment and of the missing-mass lower bound is correct.

## 3. Generic fixed-tuple theorem and exact local formula

The full parameter domain is stated: a>3, s>2a, 0<d<=1/16, d<1/(32a), sigma=sd^2<1/2, 0<epsJ/epsA/epsD<1, epsS/epsT>0 and positive b_*. Thus mass, bracket, degree and deficit have the necessary signs. The letter a for the soft exponent is separated from the original dilation mathfrak a. The theorem fixes the complete tuple before n grows; its eventual N can depend on that tuple. It does not assert a common N for adaptive/varying tuples, and the scalar endpoint c=c_* is expressly excluded.

The Psi formula for J is exact. At most one end of the source interval is removed; a removed length m leaves a residual inverse-tail integral Psi(R,(g-m)_+). The ordinary removed length is separately subtracted. The two-tail expression also covers m>g, rather than extrapolating the wrong logarithmic branch. Psi_R/Psi_s and the lower J bound give the stated local derivative. Its true onset kink is treated by local Lipschitz continuity/absolute continuity and an a.e. derivative; at |w|=R the first derivatives match. Consequently the integrated estimate holds at every q, including kink points.

The baseline q0 has R=r and exact Fubini cancellation E J=1. Free log mean and variance bounds follow from the second moment and log inequalities; Jensen supplies the nonpositive log mean. The target maximum-face coordinate is pinned at +/-r and receives a separate deterministic estimate, not a free-coordinate law or a spurious variance. The whole extended interval [q0-2h,q0] is explicitly defined even when q0 leaves the source band. All its q,R,rho,g hypotheses are verified both eventually in the family and explicitly for the two tuples.

The target good set retains the requested fixed proportion 1-epsA of each face/radius uniformly; the good events need not be independent. Its mass uses the exact annular logarithm and the lower factor d/[2(1+3d/4)], yielding m_*. For an already good target the AC derivative integrates to 1280Ldh=5/4096, so the full log product is bounded by 517/4096 and the exact P0=3579/4096 enclosure follows. Individual J precision is proved separately by |J-1|<=2d/rho, not by a coordinate union probability or by replacing the product by (1+-epsJ)^n.

## 4. Same weighted rows and columns, then the actual history

The inverse shell coefficient correctly includes every pinned source face:

    U=(product J)/(2ell_*) sum_i Ni/Ji,
    ell_*=1-g/(2R), Ni=1+1_{|wi|<mathfrak-a0}.

Bounds for ell_*, individual J and product J yield precisely u_-n and u_+n. The clipping law is the normalized capture/G law on each captured source face. It has free-coordinate normalizers J>=1/2 and deterministic pin loss <=d. With actual R>1, E_inv[-log ell]<=8g^2/(9Rrho)<=d^2/rho, so E_inv[-log G]<=3Ld^2. Markov under this inverse law retains 1-3/a. A raw-source tail probability is not substituted for that law.

Matching log width lies in [h,2h/qmin]. The bounded source factor b<=1 stays in the post-selection kernel; it changes neither ambient/actual measures nor algorithmic weights. Thus the complete lower row, upper row and column all use the same k=nh exp(-aLd^2) scale. The pointwise column envelope is 4nhG/R, from face decomposition, deleted-coordinate product cost and target width; multiplying by b gives 4k. The exact maximum-degree constant is max(4,2u_+/qmin), rather than an independently averaged column budget.

Original corner-source removal is paid by unweighted full mass because b<=1. Its exponent d/32-ad^2 is positive by the stated strict d-domain. The separate tau/k exponent (s/2-a)d^2 is positive by s>2a. Once source loss<=epsS k and tau<=epsT k, target mass and degree yield Gamma n; keeping gamma=(1-epsD)Gamma preserves a strict deficit margin. Positive gamma<Gamma<d/16 makes the root valid once B>8n.

For every original dilation mathfrak a>=8n^4, the centre guard is exactly n^-3-(beta+2)/(mathfrak a+1)>0. Actual target items have radius below 3, whereas source items have radius at least rho>4, so every included source precedes the same target. The strict increment-to-centre bound gives that target's original eligibility condition. Multiplication by b<=1 preserves the payment tau with no inverse-theta factor and no regularity assumption on mu. Scaling ambient energy, degree and target mass by alpha^-2n, alpha^-n and alpha^-n can only improve the missing-mass lower bound. Only scalar inequalities are passed to the dilation limit; original finite histories/selector convergence are not assumed. The general L1 and all-positive-radius scope is retained.

## 5. Both explicit domains are fully paid

The manuscript reconstructs the same exact simple/raw coefficients as the standalone accepted threshold proofs: u_-=3579/8200, M_*=4; the listed m_*,b_*,Gamma fractions; gamma=(1023/1024)Gamma; and each exact c_*=gamma-sigma/(1-sigma)>2^-16. No displayed decimal, best-grid identity or top64 survival condition is a mathematical input. The finite search is explicitly only a way to suggest two fixed tuples.

At each own L0, n>=L^2/2 and n>=L^8/8! prove h/d<1/8, zeta*t<h/2, extended q>1/2, R>1, rho>4r, g<=d, T<2d, nd>=2, n>=33 and Ld>=1. These geometric margins improve on the whole tail. The complete free/pinned mean and variance bounds become decreasing sums of L^-6,L^-7. Their base values, the counting bound 8/(9nd), and corner rate (3/4)^32 each allocate at most epsA/3. Thus every face retains 63/64, not the old 5/8. The individual-J upper envelope decreases as L^-7.

The explicit cap coefficient has no hidden O constant: its original shell coefficient is <=(n/2)exp(-Ld/32), its log width <=d, and the row loss after division by k is bounded by 24*2^20*d^2/(eta^4 L^3). This decreases as L^-3. Simple eta=23/2^20 has eta L0=2944 and the displayed cap inequality is correct; the raw tuple's analogous exact comparison is also the accepted one.

Tau's global Taylor-64 envelope decreases as L^(-125/2). Simple nu=2^-21 has nu L0=64 and base prefactor 2^(101/2)<2^51. Raw nu=d^2/32 has nu L0=1046529/8192>64 and base prefactor 1023*2^48<2^58. The manuscript correctly says these coarse prefactor comparisons are only at the base; the actual decreasing power envelope covers the tail.

The complete scalar cost pays sigma/(1-sigma), A and the positive square-root cross term. Its three Taylor envelopes decrease as L^-2,L^-16,L^(-31/2). The displayed simple/raw base sums are respectively below their exact c_*-2^-16 margins; all exponents and prefactors match the standalone proofs. Thus L>=2^27 and L>=2^32 really close every operator guard, not merely a scalar comparison or an unquantified generic N. The generic theorem still has only a tuple-dependent eventual N for other tuples.

The rational base constants were already independently executed by this auditor in the recent R053/R055 standalone audits: 23 compressed Fraction certificates, PID 84467; and 20 certificates, PID 84659. Their actual receipts/hashes are locked in `sourcehash.json`. They were not rerun here and the writer's 51 Fraction/SymPy checks were not executed. The present review checks the manuscript's actual estimates and infinite tails directly.

## 6. Complete-root ordering, removed legacy claims and retained low-domain proof

For the main simple/raw comparison, the lower B envelope retains -log(1-sigma_s)>=sigma_s, so Phi_s>=T-n(gamma_s-sigma_s+4/L). The raw upper envelope retains all its positive terms. The accurate strict gap at Lr=2^32 compares c_*,r-(gamma_s-sigma_s) with the sum of raw error and 4/L; both positive losses decrease. Therefore Phi_r<Phi_s throughout the common valid domain. It is not a comparison of two c_* values or two displayed savings alone.

The old R042/R049 expressions Psi_1,Psi_2 are explicitly defined only as formulas. Their lower envelopes, positive root arguments and actual deficits are adequate to prove the displayed scalar sorting at L>=2^40. No remaining claim applies either old operator theorem without a proof or silently activates an omitted threshold. The older fixed-t family comparison likewise states a formula-cost lower envelope, not a C_n lower bound, and preserves each operator bound's own domain. Removed R042/R036 theorem/proof labels do not leave undefined references or unproved self-claims.

The retained R034 proof is complete: negative target faces lie in D, two Markov conditions retain the required mass, local shell overlap yields C0, the corner subtraction leaves positive D0, kappa shrinkage loss is paid with both products and face-availability terms, and all source centres precede one actual target. The variable P_kappa(w) deficit argument needs only one actual point, not a measurable selector. A larger interior q_* supplies the strict first root inequality. Its explicit Y=64 example proves q=1/20 by positive Taylor terms and exact fractions; W_n is defined only for this old parameter regime and used with its own t. Thus the lower-domain bound is self-contained.

## 7. Real-history barrier scope and public claims

The retained R037 proposition constructs real finite histories. The n-2 free-coordinate pairs give the full incoming-energy estimate, and Markov selects small full-source potential. Continuity follows from nonzero-normal boundary hyperplanes on maximum-coordinate cones. Compact block labels retain strict radial/corner gaps. Radius-label enlargement is paid by boundary layers relative to the original m(p), so the uniform mass error does not require a bound proportional to an uncontrolled block count. Fine grids work for every sufficiently large real original dilation; the complete incoming source bound, including future/self terms, keeps every item eligible. The final minimum-radius item has positive actual weight and remains uncaptured.

Its corollary fixes the nearlinear t-window exponent C and proves a leading pointwise W(0) obstruction for that algorithm. The text expressly excludes a C_n lower bound. The common-spike example fixes n,mass,level and takes support radius to zero; its ratio tends to one at that fixed level. It expressly does not claim the whole weak norm tends to one, which would conflict with dilation invariance. Joint input/coverage information remains outside the obstruction's scope.

The A/B/C known results remain distinct from numbered original statements. Nothing claims a public record, priority, continuous optimum, whole-grid completeness or a global best tuple. Sources omitted from the current manuscript are not being promoted to new public theorems by the formula comparisons.

## 8. Source/reference check and minor definition clarification

`check_references.py` was actually run, PID 85537: initial source has 77 labels and 63 references, no undefined or duplicate labels and no undefined citation keys. This is source/reference checking, not mathematical computation, compilation or an author-check rerun.

One P3 clarity suggestion was sent to writer and coordinator: at initial line 760 the bounded factor is said to multiply the “post-selection capture kernel,” but the applied K and P are not explicitly displayed. The existing matching/J formulas uniquely reconstruct the intended Borel kernel. For self-contained notation the suggested definition is

    P=D intersect {||z||_infty/t in 1-I},
    K(w,z)=1_A(w)1_P(z)1_{|r-q-d|<h}
           *1_{||z-w||_infty<kappa*rho} b(z),
    r=||w||_infty, rho=||z||_infty, q=rho/t.

Then 0<=K<=1 and A intersect P is empty by the already proved radial separation. This is definition clarification, not a changed estimate or a missing geometric payment. The writer was asked to make only that local clarification and return the exact new source/diff for incremental locking; no author source was edited by this auditor.

### Final incremental lock

The writer's corrected final source was actually inspected and independently diffed against the recovered a05 source. The only change replaces old line 760 by final lines 760–775: it defines P, r/rho/q, the full matching/capture kernel, 0<=K<=1, disjoint supports, and the R_lambda/C_lambda integrals. The factor b is extended by zero only **outside the complete source band**, preserving its value on cap-complement sources used in the full-row analysis. In particular it is not set to zero outside D before the cap loss is paid. The source's new definition is exactly the kernel audited above and has the required Borel/positivity/support properties. No estimate, parameter, fraction, inequality, theorem or threshold changed.

Actual final SHA-256 is `17e3de43f0d841a285dbfcaf8ff12c57f9670270aab9ce95cf7fe5ce359b9554`. The single-hunk patch `upper/round054/manuscript_notes/kernel_definition_delta.patch` has actual hash `ec87204f3d0602e312dca460598c3ab9c190e734adcf6f5fab2e0b035d8ab8ad`; the recovered original still hashes to `a05aecf23e8cc3c55e65afbe38845709e58f30fe61763d248ce8609159860cb0`. Final reference checking was rerun only because the definition edit changed the source: all 77 labels, 63 references and citation keys remain defined and unique. The resulting PID and final hashes are in `sourcehash.json`.

The full ACCEPT applies to this final lock. There are no outstanding requested mathematical or definition corrections. The writer/coordinator own compilation and PDF layout checks; none are asserted on the basis of this source audit.

# R059 wide-source-band independent full-proof audit

Conclusion: **ACCEPT the complete eventual original-operator theorem**, at final frozen source `fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f`: for all sufficiently large integer n,

    Cn<=F(A_n(t),B_n(t)-(1637/12500)n)<T-n/10.

No explicit finite n-domain is certified or imported. The initial frozen source `06a3ea9ac8e754ca0f3e63664aea99de20cd49a9c9080c49c8ceb3312e35f8ff` was fully audited, and its final definition/provenance-only patch was then read and independently reverse-checked. The author source is at

    /Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/PROOF.md

Only this private clone's `audit/round059/wide_source_band/` is edited. No author helper is imported, no old grid or passed baseline is rerun, and no MATLAB/large-CPU production or Git write operation is used. A new compressed Fraction/SymPy check of this fixed tuple is recorded separately.

## 1. Scope of the theorem under review

The proposed statement is eventual, for one fixed rational tuple:

    Cn<=F(A_n(t),B_n(t)-(1637/12500)n)<T-n/10,
    T=n log n-n log log n.

There is no explicit n threshold in this proof. None is manufactured by this audit or imported from R056, whose parameters and effective-domain proof differ. This is an original centered common-radius cube operator statement for Lebesgue measure, arbitrary L1 input and every positive scale, provided the full chain below is valid; it is not a conclusion from the R058 ideal model or its numeric ranking.

The domain values d=391/1905, a_band=283/2032 and b_band=63/64 have a new role: a_band/b_band are source-band endpoints, not the R051 soft exponent or the original algorithmic dilation. The proof correctly does not invoke R051's small-d/global parameter family by name.

## 2. New band geometry and the active-source distinction

The independent rational check verifies

    0<d<1/4<1/3,
    1/2<qmin<qmax<1,
    1<Rbase<rlo<rhi<5/4,
    2*inset<qmax-qmin.

Here vmin=a_band*d, vmax=b_band*d, qmin=1-vmax, qmax=1-vmin, Rbase=qmin+d, and target endpoints are inset by delta=2^-20. Since h->0, eventually h<delta/2. Every target's whole matching interval [r-d-h,r-d+h] is inside the source band, with exact ordinary width 2h; logarithmic width is in [2h/qmax,2h/qmin].

If a source has a matched target at all, then q>=qmin+delta-h and

    R=q+d-h>=Rbase+delta-2h>=Rbase>1.

This stronger lower bound is used only for active sources. On inactive sources R may be slightly below Rbase, but is still above one eventually; their columns are zero and the pin-moment bound is not applied to them. This is a real restriction, not an interchangeable bound on the full source band.

Over the whole active matching interval, R in [r-2h,r], q>1/2, rho>=4r eventually, zeta*t<h/2, and g=d-h+zeta*t*q lies in (0,d]. In particular g<R/2 and every guard of the accepted d<=1/3 local J lemma holds. The baseline is q_b=r-d+h, so R=r exactly; choosing r-d instead would lose its exact mean identity.

Target radii exceed one and source radii eventually exceed four; both annuli lie inside ||z||<t, with qmax<1. Thus the actual shell measures and disjoint target/source geometry are unchanged.

## 3. Five ambient bad events, including the new Kw event

At a fixed full target face/radius, one target coordinate is pinned and n-1 are independent uniforms on [-r,r]. The accepted baseline log-J moment estimates yield the displayed uniform mean and variance O(L^2/n), plus a pinned O(L/n) term. For the fixed tolerance e=2^-20, Chebyshev makes baseline failure o(1). A positive fixed tolerance may be small; it does not require a uniform threshold over a varying parameter family.

Unlike the d=0.3 R056 tuple, the present d<1/4 makes Kw<=4nd nonautomatic. The proof adds it explicitly. Since r>1>2d,

    Kw=1+Bin(n-1,2d/r),
    E Kw<=1+2nd, Var Kw<=n/4.

When nd>=1, the threshold-to-mean gap is at least nd, giving the correct `P(Kw>4nd)<=1/(4nd^2)=o(1)`. This is the extra event needed by the whole-q derivative bridge.

For each endpoint q_- and q_+, the pinned target coordinate contributes one captured source sign. The exact means are `n+(n-1)a0(q_+/-)/r`, with variance at most n/4. Chebyshev at epsN*n/4 gives failure at most 4/(epsN^2*n) for each endpoint. The n-1 term and its later 1/n drift have not been dropped.

The original target cap complement is bounded **per face** by (3/4)^(n-1) eventually; no averaged factor 1/2 is incorrectly used for a positive face. Thus baseline, Kw, count-minus, count-plus and target cap are FIVE bad events. Each has uniform o(1) probability, so their union is at most fixed epsA eventually. No independence between these events or after cap restriction is required. Consequently lambda(A)>=(1-epsA)n*log(rhi/rlo), with the new full target radial width, not an old nd or nd/2 shortcut.

The independent alternate logarithm enclosure verifies that this quantity is strictly larger than mL*n for mL=159477/1000000.

## 4. Whole matching J product and normalized source-face counts

The derivative bridge uses the accepted a.e./absolute-continuity local lemma, with the correct positive zeta*rho term in g. On a good target, Kw<=4nd and n/rho<=2L. Summation gives at most 640Ld, then length 2h gives 1280Ldh=1280h_c. Combining the baseline tolerance e with this bound gives the exact P0=1-e-1280h_c enclosure for the full product.

Separately the pointwise coordinate bound |J-1|<=5d/(8rho) tends to zero uniformly in all coordinate values and active matching q. This gives the fixed epsJ enclosure without a union bound over n. The proof does not replace the full product by `(1+-epsJ)^n`.

For the source pinned-face count, both Nsum(q) and ell_pin(q) are monotone increasing. Cross-endpoint calibration yields the displayed lower loss `2h+1/n+epsN/4` and upper gain `2hd+epsN/4`. The identity

    a_+/r-a_-/(r-2h)
      =2h(1-zeta*t)(d-h)/[r(r-2h)]

is exact. Its bound by 2hd uses r,r-2h>1 eventually. Both drifts tend to zero, so the endpoint events control the entire interval by [1-epsN,1+epsN]. There is no summation of independent worst positions.

The full signed source-face coefficient is

    U=(product J)/(2ell_pin)*sum_i N_i/J_i.

Combining the product, individual J and normalized count bounds gives the stated Ulo*n and Uhi*n pointwise for the same target. The pinned face and deleted-coordinate factors are included, rather than replaced by expected signs.

## 5. Active-source inverse-face moments and retention

For a fixed captured signed source face, the free coordinates have independent densities capture/(ell*J_i) under the full source law; the pinned contribution to Phi=-logG is deterministic vpin. The accepted exact integrals

    E_full(phi/ell)=R*vpin^2/rho,
    E_full(phi^2/ell)=2R*vpin^3/(3rho)

have the correct double-tail factors. Dropping capture in their nonnegative numerators and using J_i>=1-epsJ gives the face mean and variance estimates. Independence is used only inside this fixed, full captured source face.

For an active source, R>=Rbase and g<=d imply vpin<=-log(1-d/(2Rbase)). Also R/q=1+(d-h)/q<=Rbase/qmin. The new independent log enclosure and Fraction calculations verify vpin<Qbar<1/9, mean rate below 293/20000 and variance rate below 1/900. The active restriction is indispensable here and has been respected.

The threshold chi*x is at least x/40000 above every face's upper mean once x>=40000/9. Chebyshev gives at most 40000^2/(900*x), which tends to zero and is eventually below the fixed epsRet. Averaging these **face probabilities** preserves the bound even though the means can vary by face. There is no assertion that within-face variance alone is the mixture variance.

On G>=theta, soft weight equals theta/G; elsewhere its ratio to theta/G is in [0,1]. The full soft row therefore retains at least (1-epsRet) of the inverse row at scale theta, and is at most its full scale. This is established before source Dcap. Coordinates are not assumed independent after that deletion.

## 6. Wide source cap and sharp column budget

The wide full source-band mass is n*log(qmax/qmin), bounded by Csrc*n with Csrc=(qmax-qmin)/qmin. The independent exact check confirms Csrc>d for this tuple. Thus importing the old nd bound without a new logarithm calculation would be unjustified; the proof uses Csrc instead.

At source radius rho=tq, original cap complement free fraction is `1-(1-q)/(2rho)`, and its negative pinned face is excluded. Retaining n-1 in the exponent gives cap mass at most

    Csrc*n*exp[-(1-1/n)*eta*x],
    eta=vmin/(2qmax).

Since b_soft<=1 this full unweighted mass bounds each target's lost weighted row. The certified eta-chi>1/100000 makes the exponent coefficient eventually at least 1/200000. The prefactor Csrc/h is only linear in L. Thus the cap cost divided by k tends to zero uniformly in target, and any fixed epsS can be paid. No independence is reused after D restriction.

For target columns, kappa*rho>|z_i| and kappa*rho>r imply the coordinate capture interval contains zero. With fraction ell_i<1, exactly one target endpoint is captured and ell_i>1/2; with ell_i=1 the endpoint count is at most two. Therefore T_i/(2ell_i)<=1 pointwise. The full target shell coefficient is at most n*product ell_i(r)<=nG(z), with matched log-width at most 2h/R. After b_soft, active columns are at most 2k/R<=2k; inactive columns vanish. Since M=2Uhi/qmin>2, this and the full row bound share the same legitimate denominator Mk.

## 7. True eligibility, one actual history, and finite-dilation transfer

At the stated selection threshold, A_n<=n*exp(-sigma L) and B_n<=nL eventually give

    tau/k<=(sqrt L/h)*exp[-741L/1000000]->0.

The exponent is the complete scalar-buffer payment; it is not omitted in favor of a raw-kernel gain. Thus tau<=epsT*k eventually. The row excess is exactly b*k after epsS cap loss and epsT actual threshold are paid. Target mass then gives ambient missing mass at least Gamma_L*n; the independent rational check verifies chosen gamma is below both gamma_budget and Gamma_L.

For every original finite history, mu_a has its original actual increment density, not a new random input. Its domination by the accepted corner envelope is preserved. On A/P subset D, lambda_a=c_a*lambda exactly, with c_a=alpha^-n>=1, by the stated D-subset-H_a identity.

Actual uncaptured target items meeting A have radii below three. Source items have radii at least rho>4, hence precede the **same actual target** in the original nonincreasing-radius selection. The center-displacement guard is uniform for every algorithmic dilation mathfrak a>=8n^4 once n is large. Every included source piece therefore belongs to a genuinely captured earlier item. Its weighted included increment mass costs no more than its original item weight because b_soft<=1, giving row_mu<=tau on mu_a-a.e. target.

For nu_a=lambda_a-mu_a and disjoint target/source sets, expansion of lambda-mu gives

    E_lambda,lambda-E_mu,mu
      <=D_a*[nu_a(A)+nu_a(P)]<=D_a*Delta_a.

Energy, degrees and target mass scale by c_a^2,c_a,c_a. Thus the finite-a deficit lower bound is at least the infinity-envelope expression since c_a>=1. This is a uniform statement over all actual histories and all sufficiently large algorithmic dilations, not a claim requiring selector convergence or a distribution of mu.

## 8. Original operator and complete scalar cost

The original cover factor, captured budget tau and terminal self-weight remain present in the finite inequality. The improved deficit is inserted before passing algorithmic dilation to infinity in uniform scalar bounds. The chosen tau gives exactly F(A_n,B_n-gamma*n); it does not optimize only the overlap component.

For a general L1 input, strict maximal-level witnesses give open small cubes covering each compact subset, every selected dilate remains its actual witness, and the proof is invariant under the common terminal-radius normalization. Inner regularity supplies the full level set. No compact support, scale cutoff, radius lattice or independent witness replacement is introduced.

The scalar penalty is sigma/(1-sigma)=3/97. The independent exact check verifies

    gamma=1637/12500,
    gamma-3/97=121289/1212500,
    gamma-3/97-1/10=39/1212500>0.

The remaining complete positive root cost is bounded by n times `[1/n+exp(-sigma L)+2sqrt(L)exp(-sigma L/2)]`, which tends to zero. It can be absorbed into the strict fixed slack. The endpoint c_raw itself need not be claimed. All finitely many guard thresholds depend only on this one fixed tuple; each tends in the required uniform direction. This establishes an eventual statement only, not the explicit R056 threshold.

## 9. Independent new computation and resolved final identity

`check_new_margins.py` performs one compressed new fixed-tuple check, not a search. It uses 40-even/41-odd alternating log1p sums with an exact alternating enclosure, distinct from the author's 64-term positive atanh recipe. It checks the new target mass lower bound, pin-log/moment rates, chosen deficit below complete budget, cap/tau exponents, wide Csrc, full scalar slack, and three symbolic band/count identities. Eighteen Fraction checks and three SymPy checks actually passed. No author helper or production verifier was imported.

The arithmetic receipt records the initial source identity at the actual check and is retained unchanged. It is not relabeled as a fictitious later run. The author then supplied `REVISION.patch` and `DEPENDENCY_ADAPTATION.md`; the final source explicitly binds the accepted d<=1/3 lemma at SHA97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec, which this audit read in full and checked against the wide-band guards.

The historical author receipt's local-J hash 6d1742b7... was independently confirmed by that lemma's author to be a real intermediate identity. Only header/source-identity prose differs from final 97430; no mathematical lemma content changed. The adaptation honestly preserves the original certificate and does not pretend the old external path or a retained old byte snapshot exists.

The final proof also explicitly defines the FULL source band mathcal P_0, G/b on its whole pre-cap band, P=D intersect mathcal P_0, the final Borel kernel, both ambient degrees, disjoint A/P and inactive zero columns. These definitions agree with the law/order used in this audit. No parameter, estimate or computed scalar quantity changed. The new source-lock check reverses the exact patch in memory and reconstructs initial SHA06a3 exactly, while separately checking final SHAfe95, unchanged author code/receipt, the adaptation and the final local-J dependency. It performs no repeat of the arithmetic or baseline computations.

Final accepted source identity: `fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f`. The explicit kernel and dependency clarification closes the final documentation issue. All mathematical estimates and the original general-L1/all-positive-common-radius passage are accepted at that identity; the conclusion remains eventual only.

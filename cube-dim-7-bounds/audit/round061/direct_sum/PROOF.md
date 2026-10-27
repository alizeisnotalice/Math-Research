# R061 independent bounded audit: direct sums of layer energy deficits

**Decision: ACCEPT the general measure lemma under the stated disjoint point-support and common-history hypotheses.** An increment crossing several bands does not invalidate the payment. Overlapping projected supports do invalidate the same unmodified sum in general, even if the kernels have disjoint supports in the product space.

This is a measure-lemma audit only. It does not accept any layer's root/J estimates, full parameter gates, center-capture argument, or general-L1/all-scale operator bridge. Those are separate dependencies retained by the parent. No draft R061 geometry proof was available when this audit began; the proof below uses the parent's explicit abstract hypotheses directly.

## 1. Exact hypotheses and notation

Let lambda be a positive measure, let mu<=lambda be **one common actual measure**, and define the positive deficit measure nu=lambda-mu. It suffices to restrict to the finite-measure union B of the supports below; a global finite missing-mass quantity Delta_total must satisfy `Delta_total>=nu(B)`. More generally define nu as a measure, rather than subtracting two possibly infinite total masses.

Fix a finite number J of layers. For each j let A_j and P_j be measurable, with

    A_j disjoint from A_l for j!=l,
    P_j disjoint from P_l for j!=l,
    (union_j A_j) disjoint from (union_j P_j).                (1)

Thus all 2J point sets are pairwise disjoint. Disjointness up to lambda-null sets is sufficient because nu,mu<=lambda. Assume lambda(A_j),lambda(P_j) are finite.

Let `0<=K_j<=1` be a measurable kernel supported in `A_j x P_j`. Write

    R_lambda,j(w)=integral K_j(w,z) d lambda(z),
    C_lambda,j(z)=integral K_j(w,z) d lambda(w),
    E_eta,xi,j=double integral K_j(w,z) d eta(w) d xi(z).

Suppose constants `D_j=M_j k_j>0`, `r_j>=0` and `tau>=0` satisfy

    r_j<=R_lambda,j<=D_j       lambda-a.e. on A_j,
    C_lambda,j<=D_j            lambda-a.e. on P_j,
    R_mu,j<=tau                mu-a.e. on A_j.               (2)

The last condition is sufficient in its mu-a.e. form. If the geometry supplies it at every actual target point, that is stronger than required here. It cannot be replaced merely by a statement at item centers unless the separate geometric bridge has proved the required pointwise/mu-a.e. statement.

No hypothesis assigns a whole original item or increment to just one layer.

## 2. Single-layer identity and direct sum

**Lemma.** Under (1)--(2),

    Delta_total >= sum_j lambda(A_j)*(r_j-tau)/D_j.          (3)

Negative terms can be dropped, giving the stronger bound with `(r_j-tau)_+` in every summand.

**Proof.** Product measure expansion, with the same nu=lambda-mu, gives the exact identity

    E_lambda,lambda,j-E_mu,mu,j
      =E_nu,lambda,j+E_mu,nu,j.                              (4)

Since `C_mu,j<=C_lambda,j`, the upper degree bounds imply

    E_lambda,lambda,j-E_mu,mu,j
      <=D_j[nu(A_j)+nu(P_j)].                               (5)

The row lower bound and true actual-row budget give

    E_lambda,lambda,j>=r_j lambda(A_j),
    E_mu,mu,j<=tau mu(A_j)<=tau lambda(A_j).

Consequently

    D_j[nu(A_j)+nu(P_j)]
      >=lambda(A_j)(r_j-tau).                               (6)

These are finite integrals under the stated support/degree hypotheses. If the right side is negative, positivity of nu improves that individual lower bound to zero.

Finally (1) yields exactly

    sum_j[nu(A_j)+nu(P_j)]=nu(B)<=Delta_total.                (7)

Divide (6) by its own D_j and sum, proving (3). The D_j and k_j may differ between layers; there is no need to replace all of them by one worst degree before summing. All actual-row tests use the same tau, as specified. No exchange of optimizers, histories, witnesses, or kernels is made. QED.

The kernel bound by one is useful for the separate history-eligibility bridge, but the calculation (4)--(7) itself only needs nonnegativity and the displayed finite degree bounds.

## 3. Why an increment crossing bands is harmless

The accounting concerns **measure points**, not original item labels. If actual increments induce positive measures mu_i with `mu=sum_i mu_i`, then for any one item i,

    sum_j[mu_i(A_j)+mu_i(P_j)]=mu_i(B)<=mu_i(X).              (8)

An item crossing multiple target bands contributes different portions of its increment to different A_j. Each portion is counted once; the same is true of an item crossing source bands. Equality (7) is a statement about the single deficit measure and does not even require a deficit decomposition by item.

For instance an actual target item can have positive measure in both A_1 and A_2, while an actual source item has positive measure in both P_1 and P_2. Equations (4)--(8) still apply. Replacing each portion by the item's **entire** mass in each layer would be invalid, but the asserted lemma does not do that.

Two genuinely necessary integration conditions remain separate:

* Each layer's actual-row bound must hold for the restriction of the **same** actual mu. Bounds from independently optimized histories mu_j cannot be added as deficits of one history.
* Source precedes target and center capture must be established for every included actual pair/point in the layer kernel, including points of an increment lying in different bands. This audit assumes that row bound as given; it does not infer it from a one-item/one-layer assertion.

## 4. Finite-a ambient scaling

Let the finite-a ambient measure lambda_a satisfy `lambda_a=c lambda` on B, where c>0, and let mu_a<=lambda_a be the common actual measure. The kernels are unchanged, and the actual-row budget is still tau. Then

    lambda_a(A_j)=c lambda(A_j),
    r_j^a=c r_j, D_j^a=c D_j,
    E_lambda_a,lambda_a,j=c^2 E_lambda,lambda,j.

Apply the preceding lemma with lambda_a and mu_a. The exact resulting per-layer lower bound is

    (lambda_a-mu_a)(A_j)+(lambda_a-mu_a)(P_j)
      >=lambda(A_j)*(c r_j-tau)/D_j.                        (9)

In particular, **when c>=1**, (9) implies the base lower bound

    >=lambda(A_j)*(r_j-tau)/D_j.                            (10)

All the finite-a deficit supports are still disjoint, so (9), and hence (10), sum without loss. There is no residual factor c^(-1) multiplying the base deficit. The positive part version remains available.

In the R056 bridge used as the integration reference, `c=alpha^(-n)>=1`, with `alpha=mathfrak a/(mathfrak a+1)`. Thus the direction needed in (10) is the correct one. Equality `lambda_a=c lambda` is needed only on the target/source support sets, not on the whole cap domain.

If a future source calls c merely positive and omits c>=1, preservation of the **same** base constant is false. For one target and one source unit cell, let base lambda put mass one on each, c=1/2, `mu_a(target)=1/2`, `mu_a(source)=1/4`, K=1 and tau=1/4. Then r=D=1, the scaled deficit is 1/4, while the unscaled claimed lower bound would be 3/4. Formula (9) correctly gives 1/4. This is a guard on the scaling claim, not a counterexample to the R056/R061 setting c=alpha^(-n)>=1.

## 5. Explicit failure when the projected source supports overlap

Take three disjoint unit intervals A_1, A_2, P and let lambda be Lebesgue measure restricted to their union. Let mu have density one on both targets and density 1/4 on P. Define

    P_1=P_2=P,
    K_1=1_(A_1 x P), K_2=1_(A_2 x P),
    r_1=r_2=D_1=D_2=1, tau=1/4.

Each actual row equals tau at every target point. Both individual energy identities are exact:

    E_lambda,lambda,j=1,
    E_mu,mu,j=1/4,
    nu(A_j)+nu(P_j)=3/4.

But the **one common** global deficit is only `nu(P)=3/4`. Adding the asserted two layer lower bounds would yield `3/2`, which is false. The same source deficit was charged twice.

Notably `A_1 x P` and `A_2 x P` are disjoint product-space supports. Thus disjointness of the kernels as subsets of X x X is insufficient; the proof uses the projected point supports in (1). This example has strictly positive actual mass on every support and exact actual-row inequalities, so it does not rely on a vacuous bound on a mu-null target.

The example is a counterexample to an overlapping-support **measure-lemma extension**, not an original maximal-operator history or a C_n lower-bound construction. Replacing the intervals by three atoms of unit ambient mass gives the same exact arithmetic.

More generally an overlap multiplicity bound m would give only `sum_j nu(A_j)+nu(P_j)<=m nu(B)` and hence a 1/m loss in this proof. No multiplicity-one conclusion is asserted for overlapping bands.

## 6. Actual new exact checks and source identity

`verify.py` was actually run with Python 3.15.0 and standard-library Fraction arithmetic. It computes the kernels' complete finite row/column degrees and both energies, rather than assuming the desired deficit identity. The cases are:

1. Two fully disjoint layers, with one target-item measure crossing both target bands and one source-item measure crossing both source bands: lower sum = total deficit = 3/2.
2. The same disjoint construction at ambient scaling c=3/2: exact scaled lower sum = total deficit = 5/2, while the preserved base lower sum is only 3/2.
3. Disjoint target bands with one shared source band: claimed unmodified lower sum 3/2 exceeds the actual common deficit 3/4, confirming Section 5's obstruction.
4. One disjoint layer with c=1/2: the exact scaled lower bound is 1/4; preservation of the base bound 3/4 fails, confirming the necessity of the c>=1 guard for (10).

All exact checks passed. These are new measure-lemma examples, not reruns of any old LP, grid, geometry gate, root/J diagnostic, or MATLAB script. They do not certify an actual cube realization of the abstract examples.

The only file dependency read for the finite-a identification was `/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md`, specifically its Section 8 bridge/scaling statements. Its SHA-256 was `6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a`. The complete R056 proof and its tuple gates were not re-audited here. The parent's new R061 abstract hypotheses are reproduced in Section 1; no frozen R061 draft is asserted as a read source.

`verification.json` records the actual exact outputs, interpreter and source/artifact hashes. Reproduce with `python3 verify.py` from this owned audit directory. No fetch, rebase, commit, push, production rerun, shared MATLAB operation, or author/main source modification was performed. The parent retains the layer geometry and final operator proof responsibility.

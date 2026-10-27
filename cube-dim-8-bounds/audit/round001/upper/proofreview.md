# R001 independent audit of the dimension-eight cube-cap upper bound

Date: 2026-10-27. Verdict: **PASS for the complete Lebesgue centered common-radius cube operator chain and the stated finite scalar certificate.** The reviewed frozen theory file has SHA-256 `6d47b46bc873021b80c9ce92fbb401b70d060626b2dcf5a28132b73e4035b551`. Its exact source identity and the scalar package hashes are recorded in `manifest.json` and `result.json`.

The conclusion proved here is

\[
C_8\le F_8(3837029/524288)
 <\frac{56627260639}{2000000000}=28.3136303195,
\]

where

\[
F_8(t)=\big[(1+1/t)^4+\sqrt{B_8(t)}\big]^2,
\quad B_8(t)=1-[(1+1/t)/2]^8+8\log t
 -4\int_1^t[1+1/(2t)-1/(2\rho)]^7\frac{d\rho}{\rho}.
\]

Here \(M_8 f(x)=\sup_{r>0}(2r)^{-8}\int_{x+(-r,r)^8}|f|\), and \(C_8\) is the least constant for all \(L^1(\mathbb R^8)\) inputs and strict positive level sets. There is one common radius in all eight coordinates and **every positive radius** is included. The background is constant Lebesgue measure. No variable-density theorem, global minimizer, lower bound, or public novelty assertion is certified.

This audit read the R027 theory and earlier audit as references before reconstructing the argument. It is consequently a fresh mathematical recheck with independent computation, not a blind reconstruction. Prior PASS labels were not premises. Before the new author's files were sealed, only the reference mechanism was assessed; final acceptance here follows reading the sealed R001 proof and its review notes. The scalar author's own scope properly excluded a claim about \(C_8\); that claim becomes justified only after the operator chain below is combined with the independent scalar certificate. The parent remains responsible for project integration and the final record.

## 1. Finite selection really covers the finite supplied union

Fix \(t>1,a>2t,\tau>0\), and write

\[
q=1+1/t,\ b=aq,\ \alpha=a/(a+1),\ c=(b-1)/(a+1)>1.
\]

Given any finite family of open cubes, select a maximum-radius eligible unselected cube at each step. For selected \(Q_j=Q(y_j,r_j)\), freeze the actual unstarred increment

\[
D_j=Q_j\setminus\bigcup_{i<j}Q_i,\quad
v_j=|D_j|,\quad w_j=v_j/(2ar_j)^8.
\]

A center \(y\) is eligible when \(\sum_{j\text{ already selected}}w_j1_{bQ_j}(y)\le\tau\). The weights are nonnegative and frozen. Thus eligibility can only be lost, and selected radii are nonincreasing. Zero increments and ties cause no exception.

For any rejected \(Q_\beta\), the rejection sum using selected radii at least \(r_\beta\) exceeds \(\tau\): inspect the first selection of radius strictly smaller than \(r_\beta\). Immediately before that selection the still-unselected \(Q_\beta\) must already be ineligible, and all previous selected radii are at least \(r_\beta\). If there is no smaller selection, termination gives the same assertion. For any \(x\in Q_\beta\), each such summand obeys

\[
\|x-y_j\|_\infty<r_\beta+br_j\le(b+1)r_j.
\]

Consequently \(\sum_jw_j1_{(b+1)Q_j}>\tau\) on the rejected union. Integrating this charge and adding the selected union yields

\[
|\bigcup_{\rm supplied}Q|\le
\left[1+\frac{(q+1/a)^8}{\tau}\right]\sum_jv_j.\tag{A}
\]

This is an actual union-cover bound. It uses \(b+1\), not \(b\); possible overlap of selected and rejected unions merely makes the addition an inequality. Pairwise disjointness is asserted only for the \(D_j\).

## 2. Ties, captured terms, and the finite self term

At a fixed \(x\), let \(W(x)=\sum_jw_j1_{aQ_j}(x)\). If nonzero, choose the **last selected index** among its positive contributors having minimum radius. Call it \(k\). Translate \(x\) to zero and divide coordinates and radii by \(r_k\); each weight is unchanged because both its volume numerator and denominator scale by the same eighth power. Now \(r_k=1\). Nonincreasing radii imply every other positive contributor is earlier than \(k\); a later one would have radius at most 1 and therefore equal to 1, contradicting the last-index choice.

The eligibility test at \(y_k\) pays at most \(\tau\) for all captured previous contributors satisfying \(y_k\in bQ_j\). It does not contain the self term, which must separately satisfy

\[
w_k\le |Q_k|/|aQ_k|=a^{-8}.\tag{B}
\]

For an uncaptured contributor,

\[
\|y_k-y_j\|_\infty\ge br_j,
\quad\|y_k\|_\infty<a,\quad\|y_j\|_\infty<ar_j.
\]

If \(r_j\ge t\), these imply \(\|y_k-y_j\|_\infty<a+ar_j\le br_j\), a contradiction. Therefore every unpaid term has \(1\le r_j<t\), including earlier equal-radius terms. This closes both tie and large-radius payment issues.

## 3. The geometric cap domain really contains every unpaid increment

For \(y\in D_j\subset Q_j\), the reverse triangle inequality gives

\[
\|y_k-y\|_\infty>(b-1)r_j.
\]

Set \(z=y/(a+1),u=y_k/(a+1),m(z)=\max(1,\|z\|_\infty)\). The contributing-star condition gives

\[
\|u\|_\infty<\alpha,\quad\|z\|_\infty<r_j<t,
\quad r_j\ge m(z),\quad\|u-z\|_\infty>cr_j\ge cm(z).
\]

Thus the transformed disjoint increments lie in

\[
\Omega(u)=\{z:\|z\|_\infty<t,\ \|u-z\|_\infty\ge cm(z)\}.
\]

Replacing the strict inequality by a weak one enlarges the domain. Using the positive kernel \(r_j^{-8}\le m(z)^{-8}\), disjointness, and the Lebesgue dilation Jacobian gives

\[
\sum_{\rm unpaid}w_j\le\alpha^{-8}\Psi(u),\qquad
\Psi(u)=2^{-8}\int_{\Omega(u)}m(z)^{-8}\,dz.\tag{C}
\]

The coefficient is \((a+1)^8/(2^8a^8)\). There is no missing density, star disjointness, or multiplicity assumption. The \(b-1\) in the reverse triangle bound is necessary and differs correctly from the \(b+1\) in (A).

## 4. Corner monotonicity and the shell coefficient

Coordinate reflections reduce to \(h_i=|u_i|\in[0,\alpha]\). For \(\rho\ge1\), the captured interval length in coordinate \(i\) is exactly

\[
L_i(\rho)=\min(2\rho,(1+c)\rho-h_i)
 =2\rho-[h_i-(c-1)\rho]_+.
\]

It is positive and nonincreasing in \(h_i\). Inside the unit cube the captured fraction is \(2^{-8}\prod_iL_i(1)\); its complement is therefore coordinatewise nondecreasing.

On a shell of sup-norm radius \(\rho\), the face on the same-sign side as \(h_i\) is always captured in coordinate \(i\). The opposite face is captured precisely when \(h_i<(c-1)\rho\). Hence the captured face sum is

\[
S(\rho,h)=\sum_{i=1}^8
\left(1+1_{h_i<(c-1)\rho}\right)\prod_{\ell\ne i}L_\ell(\rho).
\]

In each summand, dependence on any particular \(h_j\) enters through one nonnegative nonincreasing factor: the face count if \(i=j\), the interval length otherwise. Each summand, and therefore \(S\), is nonincreasing. This proves the desired comparison shell by shell without differentiating a moving captured volume.

For almost every point outside the unit cube there is a unique coordinate attaining its largest absolute value. Parameterize that coordinate as \(\pm\rho\) and keep the other seven coordinates free in \((-\rho,\rho)\). This is a volume parameterization with determinant 1. Ties between largest coordinates lie in finitely many hyperplanes and have zero ambient Lebesgue measure. The full shell face sum is \(16(2\rho)^7\). Thus

\[
\Psi(u)=1-2^{-8}\prod_iL_i(1)
 +2^{-8}\int_1^t\rho^{-8}
 [16(2\rho)^7-S(\rho,h)]\,d\rho.
\]

Both inner and outer uncovered terms increase coordinatewise; consequently \(\Psi(u)\le\Psi(\alpha\mathbf1)\). The comparison corner need not be attained by an actual center. Strict cap inequalities agree with the complement of the weak cap domain. An opposite-face equality affects at most one \(\rho\) per coordinate because \(c>1\), so it changes no integral.

Since

\[
c-1=(a/t-2)/(a+1)>0,\quad
\alpha/(c-1)=at/(a-2t)>t,
\]

exactly one face per coordinate is captured at the corner throughout \([1,t]\). Also \(1+c=\alpha(2+1/t)\). Put \(g_t(\rho)=1+1/(2t)-1/(2\rho)\). Then

\[
\frac{L(\rho)}{2\rho}=\alpha g_t(\rho),\quad
\frac{L(1)}2=\alpha q/2,
\]

and direct substitution gives

\[
H_8(a,t):=\Psi(\alpha\mathbf1)
 =1-(\alpha q/2)^8+8\log t
 -4\alpha^7\int_1^t g_t(\rho)^7\frac{d\rho}{\rho}.
\]

The captured shell normalization is \(8/2^8\cdot(2\alpha g)^7=4\alpha^7g^7\); this is the source of the factor 4 (the general \(n/2\)). It is not 8. As \(0<g_t\le1\), the inner uncovered term is positive and the outer integrand is at least \(4/\rho\). In particular \(H_8>0\). Combining this with (B) and (C) proves

\[
W(x)\le\tau+R_8(a,t),\qquad
R_8(a,t)=a^{-8}+\alpha^{-8}H_8(a,t).\tag{D}
\]

The finite self term remains present.

## 5. Promotion to all L1 inputs and all positive common scales

For each fixed \(R>0\), translation continuity in \(L^1\) gives continuity of the centered average, with difference bounded by
\((2R)^{-8}\||f|(\cdot+h)-|f|\|_1\). Its supremum over **all** \(R>0\) is lower semicontinuous. Therefore \(E_\lambda=\{M_8f>\lambda\}\) is open.

For any compact \(K\subset E_\lambda\), choose at each \(x\in K\) an actual bad radius \(R_x>0\) and supply the small cube \(Q(x,R_x/a)\). These open cubes cover \(K\); compactness supplies a finite subcover. Its selected stars are exactly some of the original bad cubes. Thus

\[
\lambda\sum_jv_j\le\sum_jw_j\int_{aQ_j}|f|
 =\int|f|W\le(\tau+R_8(a,t))\|f\|_1.
\]

Apply (A) to bound \(|K|\), then take the supremum over compact \(K\subset E_\lambda\) using Lebesgue inner regularity. This does not presume the level set has finite measure. The result is the genuine finite-a full-operator bound

\[
C_8\le[1+(q+1/a)^8/\tau][\tau+R_8(a,t)].\tag{E}
\]

No atomic, product, boundedness, compact-support, or smoothness assumption was imposed on \(f\). No scale was truncated or discretized. Cubes with closed rather than open boundaries have the same Lebesgue averages.

## 6. Scalar limit and admissible threshold

Keep \(t>1,\tau>0\) fixed while finite \(a>2t\) tends to infinity. The **same** full operator is bounded in (E) for each such \(a\), so only its scalar right side is being limited. In the explicit formula, the integral is independent of \(a\), \(\alpha\to1\), and \(a^{-8}\to0\). It follows that

\[
C_8\le(1+q^8/\tau)(\tau+B_8(t)).
\]

Subcovers and selectors may depend on \(a\); none must converge. Moreover

\[
B_8(t)=1-(q/2)^8
 +8\int_1^t[1-g_t(\rho)^7/2]\frac{d\rho}{\rho}>0.
\]

Choose the valid positive threshold \(\tau=\sqrt{q^8B_8(t)}\). Expanding the preceding product gives \(q^8+B+\tau+q^8B/\tau=(q^4+\sqrt B)^2\). This establishes the exact \(F_8(t)\) operator bound used by the scalar certificate. No exchange of operator limits, suprema over inputs, or minimizers is needed.

## 7. Where constant Lebesgue background is essential

The geometric background is distinct from the input \(|f|\), which remains arbitrary. The following steps use its constant Lebesgue density:

- \(|Q(y,r)|=(2r)^8\), the ratios \(a^{-8}\), the \((b+1)^8\) integration in (A), and the matching normalization of every bad average.
- The frozen increment masses are Lebesgue volumes. Common translations and dilations preserve normalized weights, and the transformation in (C) has the constant determinant \((a+1)^8\).
- Coordinate reflections preserve background mass; captured boxes have product volume, and the shell parameterization integrates ordinary face lengths with determinant 1. These facts establish corner monotonicity and the factor 4.
- Translation continuity, zero cube-boundary measure, and inner regularity promote finite witnesses to the actual full level set.

For a variable ambient density, simple cube-volume ratios, reflected product volumes, and the constant transformed density are not available. The proof here does not infer any such extension.

## 8. Independent scalar computation and actual evidence

The change of variable \(x=1/\rho\) gives
\(I_8(t)=\int_{1/t}^1(p-x/2)^7dx/x\), \(p=1+1/(2t)\). Seven exact polynomial convolutions independently give the integrated coefficients

\[
-7/2,\quad21/8,\quad-35/24,\quad35/64,\quad-21/160,
\quad7/384,\quad-1/896.
\]

Specifically
\[
I_8=p^7\log t+\sum_{j=1}^7d_jp^{7-j}(1-t^{-j}),
\quad B_8=A(t)+(8-4p^7)\log t.
\]

The sign of the log coefficient is respected by sorting its two rational endpoint products. `symbolic_check.py`, actually executed with SymPy 1.14.0 in the existing mathematical venv, independently differentiates the primitive to the exact integrand and verifies five finite-a parameter identities. All six checks returned zero residual and exit code 0. The hand integration above and the Fraction convolution do not rely on SymPy.

`independent_check.py` imports no author code and uses only Python's standard library. All mathematical enclosures and comparisons use `Fraction`. It bounds \(\log t\) directly with \(z=(t-1)/(t+1)\) and 400 terms of the positive atanh series. The remainder is bounded by

\[
0\le\log t-2\sum_{k=0}^{399}\frac{z^{2k+1}}{2k+1}
 \le\frac{2z^{801}}{801(1-z^2)}.
\]

It uses rational square bisections for displayed enclosures and separately certifies the target by a positive squared residual. For \(s=q^4\), an outward rational upper endpoint \(B_+\), and \(U=56627260639/2000000000\), it proves

\[
d=U-s^2-B_+>0,\qquad d^2-4s^2B_+>0.
\]

Hence \(U>(s+\sqrt{B_+})^2\ge F_8(t)\), with no floating square root or interval-library assumption. The actual resulting enclosures are

\[
13.3365553664751128154150225918012012532221
 < B_8(t) <
13.3365553664751128154150225918012012532222,
\]
\[
28.3136303194220576219613586917743020205613
 < F_8(t) <
28.3136303194220576219613586917743020205614.
\]

The squared residual exceeds \(0.000000001900403275813764719335\), and \(U-F_8(t)\) exceeds \(0.000000000077942378038641308225\). The three grid candidates are independently enclosed, and the winner's upper endpoint is below both other candidates' lower endpoints. The source Arb intervals are contained in the independent rational enclosures; this checks agreement without claiming to inspect Arb's internals.

## 9. Registration, trace, receipts, and limitations of provenance

The scalar source code implements the written registration: seeds 6,7,8,9,10 in order; bounded SciPy minimization on [2,32] with `xatol=2**-24` and 80-callback cap; best observed point rounded exactly to the \(2^{-20}\) grid with ties to even; 512-bit Arb certificates at that point and its two neighbors; smallest certified upper endpoint; then the next strictly larger \(10^{-10}\)-grid rational.

The actual saved trace has 19 rows: five seeds and 14 optimizer callbacks. All optimizer inputs are in [2,32]; the saved best point gives grid index 7674058, so the exact winner is \(3837029/524288\). Every saved callback value agrees with an independently bounded true scalar to within \(10^{-11}\), sufficient for this ordinary float search's consistency check. This tolerance does not enter the strict proof. The chosen grid point, row identities, strict three-row winner, and final rational rounding are verified exactly.

The author run receipt reports success, 512 bits, 19 objective calls, Python 3.15.0, SciPy 1.18.1 and flint 0.9.0, with UTC start/finish and timing. Its registered-source hashes and all seven hashes in the author manifest match the actual files. The check receipt also reports PASS. These records are mutually consistent. Saved files and self-recorded timestamps do not independently establish immutable preregistration chronology; no stronger historical claim is made. I did not rerun the optimizer or rely on its success as a theorem of global optimization. The exact chosen-point proof is unaffected by that provenance limit.

The first local independent run passed its mathematical assertions but failed while serializing a Fraction with more than Python's default 4300 digits. The checker was changed to save outward 50-place rational endpoints instead of enormous raw fractions. The final complete run exited 0 in approximately 4.4 seconds, after checking all saved callback inputs, with its actual receipt in `actualreceipt.json`. This was a serialization correction, not a change to the mathematical model or chosen candidate. No MATLAB, subagent, new search, commit, push, or source-author modification was performed by this audit role.

## Reproduction and repository state

From this clone root:

```text
python3 cube-dim-8-bounds/audit/round001/upper/independent_check.py /path/to/scalar
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-8-bounds/audit/round001/upper/symbolic_check.py
```

The first command is sufficient for the independent exact scalar certificate and saved-package checks. The second requires the existing SymPy environment and repeats only the symbolic identities. The analytic operator proof above is a strict mathematical argument, not numerical sampling evidence.

The dedicated clone base is `79f8acb20c418e60c5196a3c830103253e34c900`. It was provided as a successful no-hardlinks/no-checkout clone. I initialized a cone sparse view and checked out only root instructions/README and the two R027 reference directories; those operations completed with exit 0. An attempted checkout of `cube-dim-8-bounds` reported an absent pathspec, as that project is new and absent from base HEAD. I did not test or repair unrelated repository blobs, so no full-repository integrity claim is made. All authored deliverables are confined to this assigned audit directory. Nothing was committed, staged, or pushed by this role.

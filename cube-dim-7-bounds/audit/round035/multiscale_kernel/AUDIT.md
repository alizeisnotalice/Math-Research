# R035 strip-degree audit (preliminary; exact model only)

Set I=[d/4,3d/4], w=d/2, d²=log(L)/L, h=1/(100Ld), and K(u,v)=1{|u+v-d|<h}. Assume L>1 and h<w/2, equivalently log L>1/25. Let λ be normalized Lebesgue measure on I. The reflection u↦d−u sends I to itself. Writing x=u−d/4, a row intersection has length

  ℓ(x) = min(w,d−u+h−d/4) − max(0,d−u−h−d/4)
       = 2h − (h−x)_+ − (h−(w−x))_+,  0≤x≤w.

Thus h≤ℓ≤2h, both endpoints sharp, and ℓ=2h outside the two boundary intervals of width h. Therefore normalized worst row and column degree equals 2h/w=4h/d=1/(25 log L); minimum is h/w=2h/d. The boundary-row λ-mass is 2h/w=4h/d. The exact product probability is (2hw−h²)/w² = 4h/d−4h²/d². The open inequality has no effect for Lebesgue measure, but matters if arbitrary measures have atoms on the boundary.

These formulas use normalized Lebesgue marginals and do not prove small degree for the actual source measures. Indeed both marginals can be uniform on intervals of width h/2 centered at u=v=d/2: every pair satisfies K=1; all mass lies on rows of degree 1. This example has absolutely continuous marginals. More generally if ν has density f relative to λ and ∫f²dλ≤M, Cauchy–Schwarz gives sup_u ∫K(u,v)dν(v)≤min{1, sqrt(M·4h/d)}. If f≤M∞ a.e., the sharper bound is min{1,M∞·4h/d}. These require verified uniform bounds on the true measure; finite energy without a quantitative bound is insufficient. The λ-small boundary set may contain all mass of an arbitrary measure. With density f and energy M, ν(boundary)≤sqrt(M·4h/d); with bounded density M∞, ν(boundary)≤M∞·4h/d. Apply separately to both marginals; union bound gives the sum, capped at 1. Neither type of estimate establishes the proposed operator gain without checking the actual row/column marginals, any normalizing constants and the link between degree, selection and deficit.

SymPy independently integrated the piecewise strip length in check.py and returned 2hw−h². This is a deterministic algebra check, not evidence for an operator estimate. Precise proof status remains open pending definitions and geometric derivation from the geometry team.

A sharper, optimal bound for a probability density is also available: for any set E with λ(E)=p and ∫f²dλ≤M, ν(E)≤min{1, p+sqrt((M−1)p(1−p))}. This follows from ∫(f−1)=0 and Cauchy–Schwarz applied to (f−1)(1_E−p). The bound is attained by a two-level nonnegative density whenever M≤1/p; when M≥1/p, ν(E)=1 is attainable. Thus as L→∞, degree tends uniformly to zero if the actual L² budget satisfies M=o(log L); mere O(log L) permits order-one concentration. Apply with p≤1/(25 log L) and, for the increasing small-p branch, use the corresponding supremum over 0≤p≤1/(25 log L). The two matched intervals above have M=w/(h/2)=d/h=100 log L, showing explicitly that a finite but growing L² bound cannot by itself prevent degree 1.

## Correction and concrete dominated-measure interface

The preceding arbitrary-marginal discussion is a warning about replacing the actual λ degree by uniform radial length; it does **not** require a separate non-concentration hypothesis on μ. Geometry specifies κ=β−n^−3=1+1/t−n^−3>1, not κ<1. My initial objection that a captured pre-shell might have zero faces in a coordinate assumed κ<1 and is withdrawn. For κ>1 each pre-shell coordinate has at least one captured face; both are captured exactly when |w_i|<(κ−1)ρ (up to null boundaries). The zero-face objection is inapplicable.

Let λ_A,λ_P be the restrictions to disjoint source/target sets, μ_A≤λ_A and μ_P≤λ_P, and let d_A=λ_A−μ_A, d_P=λ_P−μ_P. For nonnegative measurable K,

 E(λ_A,λ_P)−E(μ_A,μ_P)
 = ∬K d d_A dλ_P + ∬K dμ_A d d_P
 ≤ q_A d_A(A)+q_P d_P(P),

where q_A=ess sup_{w∈A} ∫K(w,z)dλ_P(z) and q_P=ess sup_{z∈P} ∫K(w,z)dλ_A(w). Since A and P are disjoint, the right side is at most max(q_A,q_P)(λ−μ)(A∪P). This identity and inequality require no density bound on μ. For a **lower** energy bound, replacing λ_A/P by their full shells goes in the wrong direction; one must subtract the actual excluded cap and bad-set contributions using appropriate row/column degree bounds. Upper degrees can use full shells.

Geometry's actual ambient measure is dλ=2^−n max(1,||z||∞)^−n 1_D∞ dz; on each signed face of the ρ shell the measure is dρ/(2ρ) times normalized uniform measure of the other coordinates. Hence the full-shell unconditioned pre degree is the integral of

 (1/2) Σ_i (1+1{|w_i|<(κ−1)ρ}) Π_{j≠i} min(1,(1+κ)/2−|w_j|/(2ρ)) dρ/ρ,

provided the interval intersection is positive, as holds when ρ>4 and ||w||∞<3. The radial-strip length alone cannot bound this degree at scale k; the product and multiplicity must be controlled from the actual w-good conditions. For post degree at a fixed z the corresponding product is Π_{j≠i} min(1,(r+κρ−|z_j|)/(2r)), and the face count is 1+1{|z_i|<κρ−r}, under κρ+r>ρ. In the matched strip the thresholds obey r−(κ−1)ρ= u+v+ρ/n³ and ρ−(κρ−r)=u+v+ρ/n³, so both lie within h+ρ/n³ of d. This is the bridge between the good-set truncated sums and the actual coordinate-product losses. The claimed uniform exponential constants and exceptional-set budgets remain to be audited against the full proof.

## Independent asymptotic check of degree mechanism

Use L=log n, ℓ=log L, t=n/(L−3ℓ), d=sqrt(ℓ/L), h=1/(100Ld), κ=1+1/t−n^−3, r=1+u and ρ=t(1−v). On the matched strip write η=u+v+ρ/n³ and ε=h+ρ/n³. Then |η−d|≤ε=o(d); the actual coordinate loss thresholds are r−η and ρ−η. For sufficiently large n all shifted-loss nonzero coordinates belong to the corresponding K_w or K_z counting sets. Since x↦(x−a)_+ is 1-Lipschitz in a,

 |Σ_j(|w_j|−(r−η))_+−S_w|≤12nd ε,
 |Σ_j(|z_j|−(ρ−η))_+−S_z|≤16Ld ε.

For the pre-shell coordinate product q_j=1−(|w_j|−(r−η))_+/(2ρ), the sums of losses satisfy Σ(1−q_j)=ℓ/4+O(1), uniformly on the specified post-good set and matched strip. More precisely the good-sum error contributes at most (1+o(1))/200, the threshold shift at most (6+o(1))Ldh=.06+o(1), and the deterministic mean is ℓ/4+o(1). For the post-shell product p_j=1−(|z_j|−(ρ−η))_+/(2r), the corresponding bounds are 1/200+o(1) and (8+o(1))Ldh=.08+o(1), with deterministic mean ℓ/4+o(1). In both cases max individual loss times the loss sum is o(1): post-shell uses dℓ→0; pre-shell has the smaller dℓ/t. Consequently Σ log q_j, Σ log p_j equal −ℓ/4+O(1) uniformly, and the O(1) absolute constants are eventually strictly below 1. Removing one factor changes either product by exp(o(1)). This substantiates bounds e^−2 L^−1/4 ≤ each full coordinate product and each all-but-one product ≤ e² L^−1/4 for sufficiently large n.

The full pre-shell face multiplicity lies in [n/2,n], and the full post-shell multiplicity is at most n. Each radial strip has length at most 2h; its intersection has length at least h at every point of I. Since dρ/ρ=dv/(1−v), dr/r=du/(1+u), these estimates yield the claimed full-shell row lower ≥ e^−2 nh L^−1/4/4, row upper ≤4e² nh L^−1/4 and full-shell column upper ≤2e² nh L^−1/4 for sufficiently large n. Restrictions to D∞ and the good set preserve upper estimates only. The row lower for the restricted pre region requires subtraction of all removed pre mass, uniformly bounded by total excluded pre mass. Thus proving excluded pre mass=o(nhL^−1/4) is indispensable.

For the claimed pre-bad rate, on each signed ρ face the remaining n−1 coordinates are independent uniforms, so E[(|Z|−ρ+d)_+]=d²/(2ρ), E[(|Z|−ρ+d)_+²]=d³/(3ρ). The fixed face contributes d to S_z. Therefore Var(S_z)≤(n−1)d³/(3ρ)=O(Ld³), and Chebyshev gives bad-S probability O(Ld³); the K_z tail also requires a separate binomial tail argument. Integrating over the full pre radial width O(d) gives O(nd·Ld³)=O(nℓ²/L). This is o(k), k=nhL^−1/4, since the ratio is O(ℓ^(5/2)L^−1/4). This check does not by itself control D∞ complement or establish selection/load hypothesis Eμμ≤τμ_A; those remain separate proof obligations.

The two remaining pre exclusions can be bounded from the same face model. For K_z, beyond the fixed face, each coordinate enters { |z_j|>ρ−2d } with probability 2d/ρ, independently; the mean is (2+o(1))Ld. A binomial Chernoff bound gives P(K_z>16Ld)≤exp(−cLd) for some absolute c>0 and large n. For D∞, 1−βρ=−ρ+v; conditional on any signed face, each remaining coordinate independently lies below this threshold with probability v/(2ρ). Since v≥d/4 and ρ≤t, the complement probability is at most exp(−(n−1)d/(8t))=exp(−(1+o(1))Ld/8). Multiplying either probability by the full pre-shell mass O(nd) gives o(k), because Ld=sqrt(Lℓ) and d/k times n grows only polynomially in L,ℓ. Boundaries and face overlaps have zero λ measure.

For the full post shell the corresponding coordinate-excess moments are d²/(2r), d³/(3r), plus fixed-face d. Chebyshev at tolerance t/100 gives bad-S probability O(nd³/t²)=o(1), and the K_w binomial mean is at most (2+o(1))nd, so Chernoff gives an exponentially small tail at 12nd. Also 1−βr=−u−r/t, and u≤3d/4=o(1), r/t=o(1), so each non-face coordinate has a uniformly positive probability (e.g. at least 1/3 eventually) of witnessing D∞; its complement is exponentially small in n. Hence post good-set λ-mass is (1−o(1)) n∫_I du/(1+u)=(1/2+o(1))nd. All these are asymptotic, with unspecified sufficiently-large thresholds.

## Review of producer draft

Reviewed `cube-r030-geometry/cube-dim-7-bounds/upper/round035/energy_theory/PROOF.md`, SHA-256 `f7d98eb21b89d83264053fda2c594387387facfdcbdbec7cf83a41e5dd55f0bd`. For the assigned scope, §§2–7 pass: face measures, matching, both face multiplicities, Gaussian product centering and error constants, bounds for λ degrees, both bad-set mass estimates, target good mass, and dominated-measure energy inequality. In particular the conditional free-coordinate law and the fixed-face d contribution match; the pre cap-complement comes only from positive faces and is bounded as claimed. My stronger Chernoff estimates above are optional; the draft's Chebyshev/Markov estimates suffice. Degree lower is established after the pre exclusions with total mass o(k), and the same ambient upper bounds are valid for the energy identity. Exact constants e^-4/512 >1/65536 follow from e<3.

This is a **scoped acceptance** of the ambient degree and energy estimates for sufficiently large n, not independent acceptance of the inherited R033 finite-history and scalar-budget interfaces in §8, nor an explicit finite threshold or global optimality claim. The original κ<1 concern was withdrawn and does not apply. The broad arbitrary-measure density discussion is illustrative only; no μ-density assumption is used in this proof.

# A01：继承规格引用的冻结来源定位

继承规格所指的旧审计文件未在输入中找到。这里逐卡列出冻结批次已记录的来源、假设和定位；本索引没有新增独立数学审核。未读外引、私有接口、被隔离原式和两轮缺失保持原状态。

## P-33450c5f81cf2615

[冻结阅读卡](papers/AD-P-33450c5f81cf2615.json)，输入SHA-256 `ddf6b9bad5dab454b88859ed8caf4bd0f8fb498b837f8eb89dab87e502f4da0e`。

Lemma 3.2 (relaxed CCT equality)

The exact CCT problem (P0) and the relaxed perspective problem (eP) have the same infimum; every exact-transform optimum is an optimum of the relaxed problem. Their optimizer sets need not be equal.

假设：["g>0 on K", "f≥0 on K"]

定位：{"start_line": 343, "end_line": 390, "pdf_pages": [7, 8]}

证明骨架：["K_=⊆K_≤ gives inf(P0)≥inf(eP).", "Every relaxed point is a radial scaling τ≥1 of an exact-transform point.", "The perspective objective scales by τ and is nonnegative, so the relaxed point cannot improve the exact optimum.", "The paper supplies a counterexample to equality of optimizer sets when f vanishes."]

缺口：["Proof uses g>0 to define the scaling; the conclusion is false without f≥0."]

Theorem 4.1 (attainment and equivalence)

The fractional optimum is attained; the convex CCT reformulation also attains its minimum, with t g(s/t)=1 at an optimizer and equal optimal values.

假设：["K nonempty", "f≥0 and g>0 on K", "f,−g,h_i are convex polynomials", "sup_K g<∞"]

定位：{"start_line": 525, "end_line": 573, "pdf_pages": [10, 11]}

证明骨架：["Set r*=inf(FP), and Φ=f−r*g; positivity of g makes Φ≥0 on K.", "A ratio-minimizing sequence and the bound g≤M show inf_K Φ=0.", "The cited minimizer theorem for convex polynomial optimization gives x* attaining Φ=0.", "Map x* to (s*,t*)=(x*/g(x*),1/g(x*)); Lemmas 3.1–3.2 transfer optimality to PCCT."]

缺口：["The upper bound on g is essential; Example 4.1 shows nonattainment without it."]

Theorem 4.2 (strong duality)

The CCT primal PCCT and its Lagrange dual LD have equal values, and the dual maximum is attained.

假设：["Slater point for the original constraints h_i(x)≤0", "polynomial-section assumptions and positivity conditions for FP"]

定位：{"start_line": 680, "end_line": 780, "pdf_pages": [13, 14, 15]}

证明骨架：["Construct strict feasibility in PCCT from a Slater point by choosing t=(1+ε)/g(x̂), s=tx̂.", "Separate the convex primal image set from points below the primal optimum.", "Show the multiplier on the objective coordinate is positive; normalize it to obtain feasible dual multipliers.", "Take t→0+ in the separating inequality to recover a dual value at least the primal value; weak duality gives equality."]

缺口：["Requires the finite-dimensional convex separation setup used in the paper; does not justify strong duality outside stated convex/Slater setting."]

Theorem 4.3 (asymptotic SOS hierarchy)

The SOS relaxation values sup(D_r) converge as r→∞ to max(LD)=inf(FP). Each relaxation value is a lower bound on the fractional minimum.

假设：["A1: K nonempty", "A2: f≥0 and g>0 on K", "A3: K bounded", "convex polynomial data", "Slater condition"]

定位：{"start_line": 815, "end_line": 881, "pdf_pages": [15, 16]}

证明骨架：["For feasible (λ,γ,η), the SOS polynomial is globally nonnegative; restricting to K and using h_i≤0 gives f/g≥γ−Cη.", "Compactness and g>0 give g_min>0; choose L above the radius of K and C>2/g_min.", "Apply the cited high-degree perturbation result to a nonnegative dual polynomial on the scaled cube; the SOS perturbation ε_r tends to zero.", "Use strong duality from Theorem 4.2 to identify the limiting value."]

缺口：["No explicit convergence rate or finite level is given for arbitrary convex polynomial data; numerical SDP output alone is not a proof."]

Theorem 4.4 (moment dual and rank-one extraction)

max(D_r)=inf(Q_r); under the stated rank-one and nonzero-mass conditions, x̄=(L_ȳ(x_i)/ȳ_0)_i is an optimal solution of FP.

假设：["A1–A3", "Slater condition", "optimal Q_r moment solution with ȳ_0≠0", "rank M_r(ȳ)=1"]

定位：{"start_line": 906, "end_line": 1055, "pdf_pages": [17, 18, 19]}

证明骨架：["Build a strictly feasible moment sequence by pushing forward positive measure from a neighborhood of a strictly feasible CCT point; this gives SDP strong duality.", "Rank one plus the Hankel moment structure gives ȳ_α=ȳ_0 x̄^α.", "The Q_r inequalities imply original feasibility and 1/ȳ_0≥g(x̄).", "Compare primal and dual objective values to force equality and optimality."]

缺口：["Rank-one condition is sufficient, not asserted necessary. The example later reports numerical rank one only approximately and separately verifies its candidate analytically."]

Theorem 4.5 (SOS-convex exact first level)

inf(FP)=max(D′)=inf(Q′); the first SOS relaxation is exact, and x̄=(L_ȳ(x_i)/ȳ_0)_i is optimal.

假设：["f,−g,h_i SOS-convex of degree ≤2d", "A1–A2", "Slater condition", "optimal Q′ solution ȳ with ȳ_0≠0"]

定位：{"start_line": 1222, "end_line": 1338, "pdf_pages": [22, 23, 24]}

证明骨架：["SOS-convexity plus a minimizer of the nonnegative dual polynomial shows that polynomial is SOS, making the Lagrange dual optimizer feasible for D′.", "Strict feasibility of Q′ gives SDP strong duality.", "Normalize ȳ by ȳ_0, apply the SOS-convex moment Jensen inequality to constraints, −g and f, then sandwich objective values."]

缺口：["The first-level exactness depends on SOS-convexity, stronger than ordinary convexity."]

常数：{'explicit': 'Choose L>max_{x∈K}||x|| and C>2/g_min, where g_min=min_K g>0. Set Θ_r(x)=Σ_j(x_j/L)^{2r}; the penalized dual objective is γ−Cη.', 'dependence': 'The valid penalty constant depends on the positive minimum of denominator g over compact K; scaling L depends on the radius of K. The paper proves perturbation error ε_r→0 but gives no general numerical rate. SDP size grows with the monomial/moment dimension s(n,2r).'}

端点与限制：{"denominator": "g must remain strictly positive on K; points with g=0 are outside the transformation and need separate modeling.", "relaxation": "f≥0 on K is used to equate exact and inequality CCT values.", "attainment": "sup_K g<∞ is required for the stated attainment theorem; bounded K is later imposed for the SOS hierarchy.", "convexity": "Ordinary convex polynomial data support asymptotic convergence; finite first-level exactness requires SOS-convexity.", "algorithm": "No claim that every SDP solver output is an exact certificate; numerical rank evidence in the example is approximate."}

## P-fdba03cbf10513db

[冻结阅读卡](papers/AD-P-fdba03cbf10513db.json)，输入SHA-256 `a801d2407c97c0ba4cad72a6976ae356d7e02e5be81b3e1408477e0da686cdff`。

Theorem 1 (exact deterministic reformulation)

The original random fractional chance constraint is equivalent to deterministic scenario inequalities involving Φ⁻¹(z_j) times the standard deviation of the Gaussian affine numerator-minus-benchmark expression, with Σ p_j z_j≥1−ε and 0≤z_j≤1. The objective and x∈X are unchanged.

假设：["a₁,b₁ jointly Gaussian", "a₂,b₂,γ jointly discrete with J outcomes and independent of (a₁,b₁)", "scenario denominator coefficients a₂^j≥0, b₂^j>0", "c(x)≥0 on X", "ε∈(0,1)"]

定位：{"start_line": 167, "end_line": 229, "pdf_pages": [4, 5]}

证明骨架：["Condition on each discrete scenario j and use independence to express the joint chance probability as a p_j-weighted sum.", "Introduce z_j to lower-bound each scenario’s conditional success probability.", "Because the scenario denominator is positive (a₂^j≥0,c≥0,b₂^j>0), multiply through without reversing the inequality.", "The remaining random scalar is Gaussian; convert its chance constraint to mean plus Φ⁻¹(z_j) times standard deviation ≤0."]

缺口：["The OCR source mangles the formulas; I checked the exact structure against rendered PDF pages 4–5. I did not independently validate every covariance expansion entry."]

Theorem 2 (claimed convex reformulation)

The authors claim a convex optimization reformulation using auxiliary Y variables and exponential-cone/norm-style inequalities for the covariance terms. Lemma 2 bounds each z_j below by Φ(1), so log Φ⁻¹(z_j) is defined and convex on the working interval.

假设：["Assumption 1 above", "μ₁−r_j a₂^j componentwise nonnegative", "ε≤min_j p_j(1−Φ(1))", "c₀(x) is stated to be convex in the PDF", "c_i(x)≥0 and log-convex for i=1,…,n"]

定位：{"start_line": 261, "end_line": 343, "pdf_pages": [6, 7, 8]}

证明骨架：["The probability budget plus the epsilon restriction implies z_j≥Φ(1)>0.", "The cited quantile lemma gives convexity of log Φ⁻¹ on this interval.", "Nonnegative log-convex c_i and the covariance-factor constraints provide convex lifted inequalities.", "The proof substitutes the lifted Y representation into the deterministic reformulation."]

缺口：["Material issue: as printed, the maximized objective is μ₀ᵀc₀(x), while Assumption 2 says c₀(x) is convex. Standard convex maximization requires a concave objective; the proof addresses constraint convexity but does not resolve this curvature mismatch. It may be a typo or omitted sign assumption, but the stated theorem is not established as written."]

Theorem 3 (piecewise chord quantile approximation; printed bound direction corrected)

作者印刷结论称(11)给max问题(8)的upper bound并随最大节点间距趋零而收敛。已核验局部修复：在所有z_j被限制到完整节点覆盖区间[Φ(1),ξ_{K+1}]⊂[Φ(1),1)、表达式有定义且F=max_k(u_k z+t_k)为logΦ^{-1}上包络时，(11)的可行集包含于(8)，故最大化最优值为lower bound。原文仍允许z_j>ξ_{K+1}或z_j=1，没有全域包络/端点认证，不能为印刷的原全域模型确认任一界方向或值收敛。

假设：["Assumptions 1–2", "piecewise linear interpolation of log Φ⁻¹ on nodes ξ₁=Φ(1)<…<ξ_{K+1}<1"]

定位：{"start_line": 345, "end_line": 396, "pdf_pages": [8, 9]}

证明骨架：["作者原文：凸性使相邻节点间的弦高于logΦ^{-1}，式(12)继而把不等式写到[Φ(1),1]；最后节点仅称close enough to1，未证尾区间（PDF8,lines362–370）。", "作者原文：式(13)实际写出(11)约束可行集⊆(8)约束可行集，却在其后称upper bound（PDF9,lines380–389）；这是原文方向错误，旧AD卡proof_steps此前也复述了错误。", "本次核验/局部修复：在被所有节点覆盖的[Φ(1),ξ_{K+1}]内，凸性给F≥logΦ^{-1}；指数函数单调，替换后各下限约束更强，其余目标与约束相同。最大化缩小可行集只能降低上确界，所以是lower bound，不是upper bound。", "本次核验：Section3.2切线包络满足F≤logΦ^{-1}（PDF9,lines398–407），在其有效且表达式有定义的域内使约束放松，最大化给upper bound；原文将其称lower同样方向反转。", "作者原文只用节点间距趋零及笼统set-distance断言值收敛（PDF9,lines390–395）；本次不把该断言认作已证明收敛，需要尾域控制、值稳定性/可行恢复等额外条件。"]

缺口：["原文max的弦upper/tangentlower标签与其可行集包含关系冲突。已局部修正方向，不能据原标签调用。", "最后节点ξ_{K+1}<1仅称close enough to1，但模型保留z_j≤1。凸弦外推在最后节点之外不一般高于函数，而logΦ^{-1}(z)随z↑1发散；需要显式z_j≤ξ_{K+1}的域限制或认证尾处理。", "max目标的曲率、Theorem1从min到max方向及零值log/端点有定义问题另见该卡其他项；本修复只是集合包含推出界方向，不补强凸性或等价性。", "节点mesh→0本身不保证全域一致逼近或最优值收敛，原文未提供必要的尾控制/可行集稳定性证明。"]

常数：{'explicit': 'Theorem 2 uses ε≤min_j p_j(1−Φ(1)); this is what ensures z_j≥Φ(1). Piecewise approximation error depends on max_k|ξ_{k+1}−ξ_k|. The paper gives a six-row experiment K=3…6 with observed upper/lower gaps 0.2066, 0.1263, 0.0782 and 0.0497.', 'dependence': 'No general rate is proved. The convex-lift dimension grows with the number n of random numerator coefficients and J discrete scenarios; the reported CPU results are for n=5 and J=2.'}

端点与限制：{"denominator": "Each scenario denominator must remain positive; theorem assumptions ensure it by a₂^j≥0,c(x)≥0,b₂^j>0. A denominator that can vanish cannot be multiplied through this way.", "quantile": "Finite interpolation nodes cannot include z=1 because Φ⁻¹(1)=+∞; the upper approximation must explicitly handle the tail above the final node.", "curvature": "The objective curvature condition is inconsistent as printed for a maximization model.", "probability": "Gaussian/discrete independence structure and finite scenario support are essential to the exact conditioning argument."}



## A01 Step 7 cross-model source locators (current, 2026-10-07)

The legacy frozen batch above contains only the primary CCT paper P-33450 and the chance-constraint paper P-fdba03. When Step 7 invokes a model below, use its separate complete reading card here; do not treat the two legacy cards above as exhaustive. Original PDFs and extraction caches are in the delivery root `evidence/papers/<paper_id>/` and are SHA-indexed in `audit_current/`.

- **Binary/projective multi-ratio convex hull:** [P-482f9333a90a4dc7 source card](papers/AD-P-482f9333a90a4dc7.json), PDF SHA `482f9333a90a4dc761a7ed1037645f9eb0c0e804b1e0fbd8389003dcdfae5967`, 38 pages incl. Appendix A.1–A.15. Projective hull Theorems 1–2, PDF pp.6–9; common binary/RLT monomial lift, pp.12–14; simultaneous multi-ratio hull Theorem 3 p.22 and Appendix A.10.
- **Robust fractional programming:** [P-42888155391c88e9 source card](papers/AD-P-42888155391c88e9.json), PDF SHA `42888155391c88e9a23f15a65b6a6eecc67f76b60bb993336f27338cc1ea56e0`, 21 pages incl. Appendices A/B. Independent numerator/denominator uncertainty special case, pp.6–7; x-independent denominator, pp.7–8; general root method pp.8–10; Sion exchange assumptions p.18.
- **Proximal fractional algorithms:** [P-436c21dd3224191c source card](papers/AD-P-436c21dd3224191c.json), PDF SHA `436c21dd3224191c0fe6f365271279b15eff055481eb3fb24c4aa52659348731`, 15 pages. Concave denominator Theorem 7, pp.3–7; convex denominator Theorem 10, pp.8–10; whole-sequence KL convergence Theorem 15, pp.11–13.
- **Binary multi-ratio submodularity:** [P-38deac7237303114 source card](papers/AD-P-38deac7237303114.json), PDF SHA `38deac7237303114fbd5b439961b710b4d8110f47dd10a2aa528390dc008807f`, 18 pages. Single-ratio theorem with positive intercept `b₀>0` and assumptions A1–A3, pp.3–6; `b₀=0` obstruction Proposition 3 p.8.
- **Fuzzy multiobjective fractional modeling:** [P-83d2d3b1c87ebaa1 source card](papers/AD-P-83d2d3b1c87ebaa1.json), PDF SHA `83d2d3b1c87ebaa1f4fa9b1fc7b933eb277dd2e05d45cd5d2f0eda3a067ea861`, 12 pages. The fuzzy/Taylor construction and local approximation scope appear in §2.2, PDF p.4; it is not evidence for a general CCT theorem.

P-fdba03 remains a local source-error case: Theorem 3, PDF pp.8–9, prints a maximization upper-bound label inconsistent with the inclusion in (13). On the covered node interval, a convex chord lies above `q(z)=log Φ⁻¹(z)`, tightening the maximization constraints and hence giving a lower bound; for `z` beyond the last node or `z=1`, no full-domain chord envelope is established. See the P-fdba card and current A01 correction log.

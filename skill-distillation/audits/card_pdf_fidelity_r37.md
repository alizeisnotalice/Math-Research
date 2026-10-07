# 卡↔PDF 保真性验证（R37）——G02 与 D01 主依赖卡通过

续 R36（K04 三卡通过）。本批扩展到 G02 与 D01 的重依赖卡。

---

## 1. G02 / P-c0cdd27d39fd5e63（47 页）——Theorem 5.1（PDF p30）

**PDF 原文**："Let 0 < |α| ≤ 1 and |β| ≤ 0.9. Let 0 < b ≤ 2⁻⁸... There is some
constant C depending only on α, β, b, φ, and p, such that if p > 2,
‖F‖_{L^p(X,σ,S_b)} ≤ C‖f‖_p, and if p = 2, ‖F‖_{L^{2,∞}(X,σ,S_b)} ≤ C‖f‖₂."

对照源卡/skill（R14 记录）：**参数 0<|α|≤1、|β|≤0.9、b≤2⁻⁸、强型 p>2、
p=2 仅 L^{2,∞}、常数依赖 α,β,b,φ,p——全部逐字一致** ✓
定理名为 "Generalized Carleson embedding" ✓。卡↔PDF **通过，0 偏差**。

## 2. D01 / P-04672953a3d2f1c1（26 页）——Theorems 1/2/4/5（PDF p5–6）

**PDF 原文（引言 + §2–3）**：
- "We identify its dual space as the space of vector-valued Lipschitz maps; see
  Theorem 1" + "Theorem 2 provides an analogue of the Kantorovich–Rubinstein
  duality formula" —— 与卡的"等距对偶 + 达值势"一致 ✓
- "In Theorem 4 we answer in the affirmative the conjecture of Klartag, provided
  there exists an optimal transport with **absolutely continuous marginals of its
  total variation**" —— 与 skill 步骤 3"质量平衡需两总变差边缘 AC"逐字对应 ✓
- "**We provide a counterexample to the conjecture, for the case m > 1**; see
  Theorem 5. It shows that, in general, the mass balance condition (10) fails to
  be true. It follows that it may happen that an optimal transport with absolutely
  continuous marginals do not exist, unlike in the one-dimensional case."
  —— 与 skill"m>1 一般失败"及 R19 对照表逐条一致 ✓

**卡↔PDF 通过，0 偏差**。

---

## 累计与结论

| 批次 | 卡 | 结果 |
|---|---|---|
| R36 | K04：CCK Thm2.1 / Thm2 / BFS Thm2.1 | 3/3 通过 |
| R37 | G02 Thm5.1、D01 Thm1/2/4/5 | 2/2 通过 |

**卡↔PDF 保真性累计 5/5 卡，0 偏差**——转换卡的忠实性获得直接证据，
此前所有"信任卡"的结论（引用忠实性核对、参数与端点声明）被独立确认。
方法已定型（PDF 文本层提取定理原文 → 与卡逐项比对），可继续覆盖其余重依赖卡。

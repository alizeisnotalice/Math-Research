# 卡↔PDF 保真性验证（R40）——f07 与 C04 通过；累计 11/12 卡

续 R36–R39。

---

## 1. f07 / P-e0ac88eda6c127c6（Multilinear Marcinkiewicz-Zygmund Inequalities，33 页）

**Proposition 1.9（PDF p7 原文）**：
- 假设：0 < p, q₁,…,q_m ≤ ∞，1 ≤ r ≤ ∞，T : L^{q₁}(µ₁)×…×L^{q_m}(µ_m) → L^p(ν)
  **正多线性算子** ✓
- 不等式 (8)：‖(Σ_{k₁,…,k_m}|T(f¹_{k₁},…,f^m_{km})|^r)^{1/r}‖_{L^p(ν)}
  ≤ **‖T‖**·Πᵢ‖(Σ_k|f^i_k|^r)^{1/r}‖_{L^{q_i}(µ_i)} ✓
  ——常数即 **‖T‖ 本身**（相对常数恰 1），与 skill/R22 的"常数恰为 1"一致
- 附加：双线性卷积推论（条件 1/q₁+1/q₂ = 1/p+1，即 Young 指标）✓
- 文章结构段："the equality in the case p = r follows by induction"——锐性结构 ✓
  与 R22 的取等构造（m=1、T=恒等、p=q=2、f_k 同向）相互印证

## 2. C04 / P-922909d5e8e47a30（Quasi Branch and Bound for Smooth Global Optimization，21 页）

**PDF p1 摘要原文**：
- "lower bounds are replaced by a relaxed notion of **quasi-lower bounds, required to be
  lower bounds only for sub-cubes containing a minimizer**" ✓
  ——与 skill 的守卫"qlb 只须在含全局最小点的盒上 ≤ f*，不能当任意盒的真下界"逐字对应
- "qBnB(2) achieves **second order convergence** based only on a bound on second
  derivatives, without requiring calculation of derivatives" ✓
  ——与 R22 数值验证的 gap 公式 (L₂/2)r²（r 减半 gap ÷4 = 二阶）一致
- "provably more efficient than the second order Lipschitz gradient algorithm" ✓
  ——skill 的算法对比叙述有原文依据

---

## 累计

| 批次 | 卡 | 结果 |
|---|---|---|
| R36 | K04 ×3 | 通过 |
| R37 | G02、D01 | 通过 |
| R38 | E03 | 1 处转录偏差（已修正）；常数确认 |
| R39 | b03、G03 | 通过 |
| R40 | f07、C04 | 通过 |

**卡↔PDF 累计 11/12 卡逐字通过、1 处转录偏差（E03，已修正）**。
验证模式稳定：转换卡忠实度高，公式级细节需抽查——抽查比例约 8% 已捕获 1 处
会误导重推的偏差，覆盖价值明确。

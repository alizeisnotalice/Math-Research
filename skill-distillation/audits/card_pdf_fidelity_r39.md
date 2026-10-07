# 卡↔PDF 保真性验证（R39）——G03 与 b03 通过；累计 8/9 卡

续 R36–R38。本批两卡通过，其中 b03 的验证与既有独立证明形成交叉印证。

---

## 1. b03 / P-91433be844e70d18（A generalized moment approach to sharp bounds for conditional expectations，43 页）

**Theorem 2（PDF p11，Fundamental theorem for conditional expectations）**：
`sup_{P∈P(q)} E_P[g(X)|X∈Ξ] = sup_{P∈D(q)} E_P[g(X)|X∈Ξ]`，
且最优值可达时存在 optimal basic distribution——**生成集 T 中至多 m+1 个分布的凸组合**。

对照源卡：Lemma 1（同一极值分布/同一逼近序列）、Theorem 2（m+1 生成元归约）——**逐字一致** ✓
证明引 Theorem 6.1 in [50]（广义矩问题）——诚实外引，与卡记录一致 ✓

**与既有验证的交叉印证**：m+1 生成元归约直接解释了本审核 R12/R23 的方法学——
b03 的 sup E[X|X≥t] 有 2 个矩约束（均值、方差）⟹ m+1 = **3 原子族即完备扫掠**，
我此前用三原子族暴力搜索逼近（且不达到）正是定理 2 结构的体现；"sup 不达到"
的机理（两点结构在端点失效）与"可达需 m+1 生成元"的极值理论一致 ✓

## 2. G03 / P-e8f8d143f231a4b3（On the Bloch and Qp–Carleson measure problems，30 页）

**Theorem 1.7（Compactness，PDF p7）**：
(i) id: B̊ → L²(μ) 紧 ⟺ (ii) 对每可允许二进分辨率 R：B_R(μ) < +∞ 且
lim_{ρ→1⁻} B_R(1_{S_ρ}μ) = 0 ⟺ (iii) 同理用 C_R；
由此 id: B → L²(μ) 紧 ⟺ μ(D) + B_R(μ) < ∞ 且尾容量消失（⟺ C_R 版本）。

对照源卡：**逐字一致** ✓（含"equivalently be expressed using C_R"）。
卡 gap"Qp 类比后文以并行证明描述"——PDF 原文"The details are given in Section 6,
where we give a careful comparison between the Bloch and Qp settings" ✓ 一致。

---

## 累计

| 批次 | 卡 | 结果 |
|---|---|---|
| R36 | K04 ×3 | 通过 |
| R37 | G02、D01 | 通过 |
| R38 | E03 | **1 处转录偏差（势公式 4 分配），已修正**；常数 4/32 确认 |
| R39 | b03、G03 | 通过 |

**卡↔PDF 累计 8/9 卡逐字通过、1 处实质性转录偏差（已修正）**。
方法学价值再次确认：b03 的 m+1 定理为既有的三原子暴力验证提供了理论根据——
两层验证（数值暴力 ↔ 极值理论）互相印证。

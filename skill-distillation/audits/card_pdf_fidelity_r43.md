# 卡↔PDF 保真性验证（R43）——A03 与 G05 通过；累计 15/16 卡

---

## 1. A03 / P-f683dde51b3674de（18 页）——Theorem 4.2（PDF p14）

**PDF 原文**：
- c 连续于 [0,1]²，二进方块 I^n_{i,j}，cₙ = 逐格 min、cⁿ = 逐格 max (12) ✓
- 收敛 (13)：lim∫cₙγ^max = lim∫cⁿγ^max = sup_{γ∈C}∫cγ ✓
- "the sequence of maximizers converges, **at least along some subsequence**, to a
  maximizer" ✓ ——子序列收敛与卡的标题声明一致
- **定义域 [0,1]²** ✓——skill 的范围限定（"仅是 [0,1]² 均匀边缘 copula 的源定理；
  任意固定边缘不是原假设"）有原文直接依据
- 证明依赖：(13) 已在 [25] 证明；其余按 Villani Thm 5.20 的思路 ✓——诚实外引

## 2. G05 / P-b7592e844240da1c（The Bellman Functions of the Carleson Embedding Theorem and the Doob's…，18 页）

**Theorem 1.1（PDF p1 原文）**：
- Carleson 条件：Σ_{J∈D, J⊆I} α_I ≤ C|I| ✓
- 嵌入：Σ_{I∈D} α_I|⟨f⟩_I|^p ≤ C_p·C‖f‖_p^p（p>1）✓
- **"the constant C_p = (p′)^p is sharp (cannot be replaced by a smaller one)"** ✓
  ——与卡的"常数 (p′)^p sharp"逐字一致
- Bellman 函数 B(F, f, M; C) 三变量定义 ✓——与卡的 Bellman 框架记录一致

**与既有验证交叉印证**：R24 的小树数值（p=2 时 4C）+ R31 的 Doob L² 证明
（‖max‖₂ ≤ 2‖f‖₂ ⟹ 2²=4）——原文的 (p′)^p 在 p=2 恰为 4，三个独立路径
（原文 Bellman / R24 数值 / R31 Doob）收敛到同一常数。

---

## 卡↔PDF 累计

```
R36–R43 八批：16 张卡核验，15/16 逐字通过、1 处转录偏差（E03，已修正）
G05 常数三路印证：原文 sharp 声明 + 数值 + Doob 独立证明
```

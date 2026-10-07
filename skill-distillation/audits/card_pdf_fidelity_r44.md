# 卡↔PDF 保真性验证（R44）——G01 与 C03 通过；累计 18/19 卡

---

## 1. G01 / P-e9e9211495f0b20e（Carleson Measure Estimates and ε-Approximation…，31 页）

**Theorems 1.1/1.2（PDF p4 原文）** vs 源卡——**逐项一致**：
- Theorem 1.1：A 向——存在 UR 边界超域 Ω̃（Ω̃⊂Ω、∂Ω⊂∂Ω̃）⟹ (a),(b) 成立；
  B 向——corkscrew + (1.4),(1.7) + (a) 或 (b) ⟹ 存在这样的 Ω̃ ✓
- Theorem 1.2：(a) 或 (b) ⟹ 调和测度 packing `Σ_j dist(p_j,∂Ω)^d ≤ C(ε)R^d`
  （ω(p_j,E_j) ≥ 1−ε、E_j 互斥）；反向同 ✓
- 证明深度：Part A = Whitney cube 简单应用 [HMM2]；Part B = corona 分解变体
  [GMT]、"deeper"——与卡的深度评估一致 ✓
- （G01 的另一主卡 P-c0cdd27d 已在 R37 随 G02 验证通过——同一 PDF。）

## 2. C03 / P-916bf62c1502285e（Best ℓ1-Approximation of Nonnegative Polynomials…，4 页）

**Theorem 2.1（PDF p3 原文）**：
- 最佳 ℓ1 逼近 `g(x) = f(x) + λ*₀ + Σᵢ λ*_i x_i^{2d}`（λ* ∈ R^{n+1}₊）✓
- `ρ_d = Σᵢ λ*_i` ✓
- λ* 是 SDP (2.5) `min_{λ≥0}{Σλᵢ : f+λ₀+Σλᵢxᵢ^{2d} ∈ Σ[x]_d}` 的最优解 ✓
- SDP 原对偶链 (2.6)→(2.7)→(2.8) 结构完整，最优性由紧可行集 + Gauss 权测度
  构造保证 ✓

**与 skill 的衔接**：C03 步骤 3 的证书方向（p∓γ 的 PSD）与步骤 5 的平坦延拓
条件——与该源文的 SDP 框架一致；卡↔PDF 通过，0 偏差。

---

## 卡↔PDF 累计

```
R36–R44 九批：19 张卡核验，18/19 逐字通过、1 处转录偏差（E03，已修正）
剩余重依赖卡约 4-5 张
```

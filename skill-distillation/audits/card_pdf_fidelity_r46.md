# 卡↔PDF 保真性验证（R46）——I03 与 e04 通过；累计 24/25 卡

---

## 1. I03 / P-787a427d3f582b5c（On the Fourth Moment Condition for Rademacher Chaos，40 页）

**Theorem 1.6（PDF p7 原文）**：
"Assume that X is a symmetric Rademacher sequence. Then, **for each m ≥2**, there is a
discrete multiple integral F of order m with respect to X such that **E[F²] = 1, E[F⁴] = 3**
which is **not normally distributed**. In particular, the fourth moment theorem fails for
chaos of order m ≥2."
——与卡的 Theorem 1.6 记录**逐字一致** ✓（方差 1 + 四阶矩 3 但非正态的反例）
- Theorem 1.1 框架：非对称非齐次 Rademacher 序列（p_k∈(0,1)）✓ 与卡一致
- Theorem 1.1 的 d_W(F,N) ≤ C₁(m)√(EF⁴−3)+C₂(m)√(sup Inf_k(f)) 结构、Corollary 1.2
  的 CLT 判据——与卡一致 ✓

## 2. e04 / P-faa1738580ee2154（Optimal Stopping, Randomized Stopping and Singular Control…，15 页）

**Theorem 4.1（PDF p5 原文）**——六类上确界相等链：
`sup_{ξ∈A^c_H} E∫k e^{−ξ}dξ = sup_{ξ∈A_H}(...) = sup_{G∈G_H} E∫k dG = sup_{τ∈T_H} E[k(τ)] = Φ = sup_{G∈G*_H}(...) = sup_{τ∈T*_H} E[k(τ)]`
——与卡的"Theorem 4.1 将最优 stopping、randomized stopping 及 singular control 的
上确界值相互等同，并说明在相应可预测类中仍有相同值"**逐字一致** ✓
- 原文假设 **E[sup₀≤t≤T |k(t)|] < ∞** 可见——与卡中"对可取负值的 k 不成立"的
  质疑记录**分区清晰**：该质疑针对的是另一条 randomized-mass 声明，不针对
  Thm 4.1 链本身；卡的分区判断正确 ✓
- Theorem 4.2 的显式变换（τ*→阶跃 G*、G*→广义逆 α*(r)、ξ*↔G* 的微分关系）
  与卡一致 ✓

---

## 卡↔PDF 累计

```
R36–R46 十一批：25 张卡核验，24/25 逐字通过、1 处转录偏差（E03，已修正）
e04 卡的"质疑分区"（Thm4.1 本体 vs 另一 randomized-mass 声明）经 PDF 确认正确
——蒸馏卡不仅转录定理，还正确标注了源文可疑处（与 R9 的 Snell 验证呼应）
```

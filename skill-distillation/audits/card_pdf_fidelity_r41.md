# 卡↔PDF 保真性验证（R41）+ 缺口清单清理

---

## 1. 缺口清单清理（账面维护）

`coverage_and_gaps.json` 中 3 条过时条目更新为已解决状态：
- b03 log 遗漏 →【已解决 R12】
- E01 Theorem 5.5 反例 →【已解决 R9】
- "0 个定理对着证明验证" →【大部分闭合 R31–R34】（源论文专有定理仍开放）

## 2. F05 / P-3da1b59567d35c36（Contributions of Issai Schur to Analysis，99 页）

**The weighted Schur test（PDF p17 原文）**：
- ζ_r(A) := sup_j (1/r_j)·Σ_k |a_jk|·r_k，κ_r(A) := sup_k (1/r_k)·Σ_j |a_jk|·r_j (3.14)
- **C_A ≤ √(ζ_r(A)·κ_r(A))** (3.15) —— 与源卡的 "C_A ≤ √(ζ(A)κ(A))" **逐字一致** ✓
- 插值形式 (3.16)：‖A‖_{l²→l²} ≤ √(‖A‖_{l^1,r→l^1,r}·‖A‖_{l^∞,r⁻¹→l^∞,r⁻¹})
  ——与 R20 的独立证明结构（Cauchy–Schwarz 拆分 + 行列界）完全对应：
  R20 证明的 R 即 ζ、C 即 κ，两条路径互相印证
- 附加确认：Toeplitz/Hankel 矩阵的 C_A ≤ Σ|w_l| 推论、Hilbert 矩阵反例（未加权
  Schur 失效的例子）——为"加权版本必要性"提供了原文实例

**卡↔PDF 通过，0 偏差。**

---

## 卡↔PDF 累计

```
R36: K04×3 ｜ R37: G02+D01 ｜ R38: E03（1 偏差已修正）
R39: b03+G03 ｜ R40: f07+C04 ｜ R41: F05
累计 12/13 卡逐字通过、1 处实质性转录偏差（E03，已修正）
```

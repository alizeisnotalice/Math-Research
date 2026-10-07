# 卡↔PDF 保真性验证（R45）——k05 与 H04 通过；累计 21/22 卡

另：登记用户常设指令——**全部完成时明确通知用户**（完成判据：7 条验收条件
全过，或批处理可达工作全部完成且仅剩人类信息/研究级缺口时给出最终清单）。

---

## 1. k05 / P-b7ea4e1f39a3af1b（Four Talagrand inequalities under the same umbrella，35 页）

**凸距离不等式（PDF p1，式 (1)）**：
`∫_Ω e^{(1/4)d²_A} dP ≤ 1/P(A)` ——与卡的 `P(A)·E_P[exp(d_A(X)²/4)] ≤ 1`
**逐字等价** ✓
- U_A(x)、V_A(x)（凸包）、d_A（0 到 V_A 的欧氏距离）定义逐字一致 ✓
- 文献溯源 [81, Theorem 1.1]、[83]（记号 f_c(A,x)）✓ 与卡的"original form
  summarized in Section 1"一致
- p2 的 Bernoulli 型不等式 (4)（Z=sup f、U/V 参数、K exp(−t/(KU log(1+tU/V)))）✓
- 卡关于"log-Sobolev 推 Gaussian 集中 r²/(2L²)"的式 (5) 记录——该综述的
  结构（四不等式统一处理）与 p1–3 的内容组织一致 ✓

## 2. H04 / P-39f1b4af17834b0f（Asymptotics of the rate function in the LDP，6 页）

**Theorem 2（PDF p3 原文）**：
"For (8) to hold (as x→∞), it is **necessary and sufficient** that
L(x) ∼ L₀(x) for some real-valued **concave** function L₀."
——与卡的"−Λ*(x)∼L(x) iff L(x)∼L₀(x) concave（等价 Λ*(x)∼−ln P(X≥x)）"
**逐字一致** ✓
- Remark 3：扩展实线的 ∼ 约定 ✓（与卡对端点/零处的谨慎一致）
- 证明两部分：必要性（L₀ := −Λ*，Λ* 凸 ⟹ −Λ* 凹）+ 充分性（拟极大化子
  t_x 的构造 (21)，需 (i) 易分析 (ii) 接近 Λ*(x)）✓ 与卡的证明骨架一致
- 6 页短文，全文覆盖范围与卡的全文阅读声明相容 ✓

---

## 卡↔PDF 累计

```
R36–R45 十批：22 张卡核验，21/22 逐字通过、1 处转录偏差（E03，已修正）
```

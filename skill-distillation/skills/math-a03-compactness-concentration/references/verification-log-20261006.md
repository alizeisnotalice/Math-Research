

---

## 引用忠实性核对（R22）

- 步骤 4：Theorem 4.2（dyadic cost approximation and subsequential optimizer convergence，
  PDF pp14–15）在源卡 AD-P-f683dde51b3674de 出现 5 次 ✓；skill 明确「仅是 [0,1]² 均匀边缘
  copula 的源定理；任意固定边缘不是该定理原假设」——范围限定与源卡一致 ✓
- 步骤 5：Theorem 1（≤k 原子分割）+ Corollary 4（≤k 原子代表保矩保积分）与 skill 的
  「原子代表仍可行；不能由原子数界推出近极值列收敛」范围声明一致 ✓
- 步骤 3：紧域耦合 tightness 由边缘尾概率并集界 —— 标准论证，方向正确 ✓

## 档位提升依据（R22）
待证据 → **部分可用**（引用与范围声明逐条一致 + 标准论证方向正确）。
仍待证据：copula 定理证明本身；分位变换迁移的 cost/边缘条件核验。


---

## 更正条目（R43）：卡↔PDF——A03/G05 通过

- **A03 Theorem 4.2**（PDF p14）：二进代价逼近 cₙ/cⁿ、收敛 (13)、子序列极大化
  收敛、**定义域 [0,1]²**（skill 范围限定有原文依据）、证明引 [25]+Villani 5.20 ✓
- **G05 Theorem 1.1**（PDF p1）：Carleson 条件、嵌入、**"(p′)^p is sharp"** 逐字一致 ✓
  **三路印证**：原文 Bellman sharp + R24 小树数值 + R31 Doob L² 独立证明
  （p=2 时三者均给 4C）
- 卡↔PDF 累计 16 张核验，15/16 逐字、1 偏差（E03）已修正

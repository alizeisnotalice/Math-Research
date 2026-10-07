

---
## 核心嵌入验证（R24）

- **步骤 3/5**：Theorem 1.1 的 (p′)^p sharp 常数，p=2 时恰为 **4C**——小树 50 例
  Σα_I|<f>_I|²≤4C‖f‖² 数值通过 ✓；skill 的弱构件（互不交孩子+统一 mass 流+全子树
  packing≤Cσ(J) ⟹ 蛋糕+Doob L² 嵌入 ≤4C‖f‖²）与源卡 Bellman 框架一致 ✓
- 步骤 1/8 守卫：TC-A4 缺定义只返回缺口；「不能当一般 Carleson 树常数 e」✓

## 档位提升依据（R24）
待证据 → **部分可用**（p=2 嵌入数值验证 + Theorem 1.1 引用一致 + 弱构件链正确）。
仍待证据：Bellman 全凹性/sharpness；TC-A4 私有接口定义。


---

## 更正条目（R31）：承重经典引理证明复现完成——依据类型升级

四条引理的完整证明已重写并数值确认（总表：`audit_r12/classical_proofs_verified_r31.md`）：

- **Hoeffding 引理**：凸性 + G(s) 二阶控制（G″≤1/4 by AM–GM）。数值 2000 组最差比值 0.999984（界紧）。
- **Bernstein 尾**：Bennett mgf + 子指数逐项系数（2·3^{k−2}≤k!）+ 最优 λ=t/(σ²+Mt/3) 代入化简。
- **Cantelli**：参数化 a>0 的 (σ²+a²)/(a+t)² 最小化（a*=σ²/t）。紧性此前已由两点分布证明。
- **Doob L^p 极大不等式**：Doob 极大引理（互不相交 A_k + 条件期望恒等）+ 层蛋糕 + Hölder。
  **p=2 给出 ‖max‖₂²≤4‖f‖₂²——本 skill 的 (p′)^p=4C 的 Doob 侧来源闭环。**
  （4C 的「sharp」仍以源卡 Bellman 框架为据，Doob 不等式本体不含锐性。）

本 skill 引用上述引理处的依据类型由 standard_theorem_citation 升级为 **source_proof_verified**
（不等式本体已重证；其各自上游——如 Bennett 界的一般推导、Doob 常数锐性——仍为引用，如实保留）。


---

## 更正条目（R43）：卡↔PDF——A03/G05 通过

- **A03 Theorem 4.2**（PDF p14）：二进代价逼近 cₙ/cⁿ、收敛 (13)、子序列极大化
  收敛、**定义域 [0,1]²**（skill 范围限定有原文依据）、证明引 [25]+Villani 5.20 ✓
- **G05 Theorem 1.1**（PDF p1）：Carleson 条件、嵌入、**"(p′)^p is sharp"** 逐字一致 ✓
  **三路印证**：原文 Bellman sharp + R24 小树数值 + R31 Doob L² 独立证明
  （p=2 时三者均给 4C）
- 卡↔PDF 累计 16 张核验，15/16 逐字、1 偏差（E03）已修正

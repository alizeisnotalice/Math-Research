---
name: math-b03-generalized-moments
description: "用于广义矩问题（Bertsimas–Popescu）的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **交接准备状态（结构性）**：本副本只整理技能结构与引用可移植性；未进行全文深读或数学验证。原评级、可用状态、假设、证明限制和未决条件均沿用原件。结构检查不构成科学验证。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# B03 · 广义矩问题（Bertsimas–Popescu）

当前证据状态：**部分可用：基础推导已检查，专题证据仍待完整验收**。全文转换、论文阅读和证明核验是三个独立状态。

## 输入与产出

支持集或结构类、测度质量、已知矩 E[h_j(X)]=q_j、待界定事件 Ξ 与目标 g。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. P-91433be844e70d18 (4)先把原未知分布P写成指定支撑的概率测度，明确h₀≡1,q₀=1及每个广义矩约束；变换后的Q则是有限非负测度，不能继续要求全域质量1。
2. P-91433be844e70d18 (5)将无条件事件概率/目标写成测度线性泛函；条件期望保分母P(Ξ)，逐可行P检P(Ξ)>0，单点正值不等于全族统一正下界。
3. P-91433be844e70d18 (6)条件比值设α=1/P(Ξ)>0,dQ=αdP：Q(Ω)=α且Q(Ξ)=1，Q不是一般概率；逆映射P=Q/α并检矩约束。Theorem1强对偶另需bounded g、uniform inf_P P(Ξ)≥ε>0及q在moment cone内点，外引conic双偶基础仍单列未核。
4. P-91433be844e70d18 (7)逐点majorant积分给弱对偶界；紧界须有可行极值分布或极值序列匹配。Prop1的均值方差例仅σ>0,t<μ给sup μ+σ²/(μ−t)，不达到；t≥μ为∞，不得称有限最大值。
5. P-91433be844e70d18 Theorem2仅对满足A1–A3和mixture-generated convex class的极端生成表示归约至≤m+1生成元（包含mass moment）；其Popescu外引cardinality证明仍未独立核，原子数、强对偶和达到分别验收。
6. P-91433be844e70d18 Theorem3的contextual DRO要保C1–C4及附加conic representability/Slater，非紧支持或事件概率趋零不可套紧性/正margin。若另用SOS/SDP，须另给半代数支撑、矩阶数与证书方向，不能由该篇直接认证通用SOS层级。

## 证据与失败处理

有限矩通常不唯一确定分布；条件比值在 P(Ξ)→0 时可能失控，且端点上确界可不达到或为无穷。弱对偶可用可行majorant，但强对偶需相应闭性/Slater条件；有限原子结论依赖具体 moment cone/生成类定理。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

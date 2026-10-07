---
name: math-b01-lagrange-farkas
description: "用于Lagrange 对偶与 Farkas 证书的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **交接准备状态（结构性）**：本副本只整理技能结构与引用可移植性；未进行全文深读或数学验证。原评级、可用状态、假设、证明限制和未决条件均沿用原件。结构检查不构成科学验证。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# B01 · Lagrange 对偶与 Farkas 证书

当前证据状态：**待证据验收：可执行方法工作流，不是已证明的课题结论**。全文转换、论文阅读和证明核验是三个独立状态。（此行声明已被下方【修订】标注替代：该 skill 经外部审核提升为「部分可用」，详见对应修订行与 references/verification-log-20261006.md）
  【修订 2026-10-06（外部审核 R26）】既有验证记录支撑档位提升：B1 案例通过（LP gap 恒等式一行解析证明，A 级）。证据状态提升为：部分可用：基础推导已检查，专题证据仍待完整验收。源定理证明本身未重证（守卫保持）。

## 输入与产出

优化方向、等式/不等式约束、凸性、乘子符号及可行性信息。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 最小化g≤0取λ≥0，L=f+Σλg+Σνh；先证明对每个原可行x及对偶可行乘子q(λ,ν)=inf_x L≤f(x)。
2. 有限实矩阵Farkas形式明确为：b=Ax,x≥0，或存在y使Aᵀy≥0且bᵀy<0，二者恰一成立；零矩阵另直接处理。用精确运算核验方向及严格负margin。
3. P-a9029e19a7586865用有限生成锥闭性及最近点构造y；证书conditioning依赖dist(b,AR₊ⁿ)，不提供统一正margin。
4. 混合LP来源P-0f4471d96cc82d4f采用A_Ix≥b_I、Aᵀλ=c、λ_I≥0：精确gap=cᵀx−λᵀb=λᵀ(Ax−b)。证书充分性不需rank/Slater；该文顶点与必要乘子证明另需可行、rank(A)=n及目标下有界。退化点只要求存在一个最优working set，不要求每个活动子集乘子皆非负。
5. 有限非负线性系统P-8b6157c4cbb03498：Ax=b,x≥0与Aᵀu≤0,bᵀu=ρ>0恰一可行；最小残差z=b−Ax*给bᵀz=||z||²，非零时ρz/||z||²是归一证书。nullspace参数化x=x̄−Kᵀy的显式双射另需rank(A)=m≤n；维数m≤n本身不足。近零残差及小分母不能作浮点可行性认证。
6. 锥/非线性/无限指标情形先定义对偶锥、拓扑和乘子空间；弱对偶、零间隙及对偶达到分别验收。
7. 锥LP来源P-08ccea72ef78920a用h_c(b)与h_c**(b)区分原/对偶值：有限值下零gap检h_c在b的下半连续，dual达到另需∂h_c(b)非空。保拓扑、连续A及正对偶锥；PDF的bar h_c proper不能误抄成h_c proper。Gale全无限指标的liminf尾与截距inf不能由有限采样认证。
8. DC-composite来源P-18e365da3a6d391f的sup_λ inf_μ必须保留一个共同λ。其Theorem3.6印刷∀μ∃λ不足以推出∃λ∀μ，缺量词交换证明时不调用等价。
9. 输出原与对偶可行候选及gap；浮点残差不冒充精确不可行证书。
10. 非闭锥来源P-5ae169908ed70fca限Hilbert与bounded linear A、bounded closed convex generator K∋0：cl A(coneK)的零support方向测试只证明∀ε>0近似可解；exact A(coneK)另需一个共同C使∀y ⟨b,y⟩≤Cσ_K(A*y)。ε>0 dual唯一不推ε=0 dual达到；general support-face inclusion需同时回验残差/normal，H/E恢复原非凸锥是额外条件。该文Remark1.7 upper semicont错误及Prop4.1平方根缩放错误按卡局部修复后才可调用，均尚待两轮审核。

## 证据与失败处理

有限维Farkas不证明一般锥/非凸强对偶。DC-composite原文存在已定位的量词缺口：逐μ选λ不是共同乘子，须修正或另证minimax/鞍点条件。 该有限LP扰动证明的共同小ε与固定working-set子序列已做局部补足，原文Step1计数须加终端检查；不得将rank假设移植为Farkas矩阵限制或无限维强对偶。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

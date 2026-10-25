---
name: math-b04-qualification-facial-minimax
description: "用于约束资格、面约化与极小极大交换的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> 外部文献、claims/cases 与 PDF 属于完整证据包，不随单个 Skill 安装。请使用[便携证据访问指南](references/evidence-guide.md)，按本 Skill 的[来源索引](references/handoff-evidence-index.csv)通过统一 resolver 定位。registry/显式包根路由均不依赖当前工作目录；SHA 检查只确认当前字节，不等于全文阅读或数学验收。


# B04 · 约束资格、面约化与极小极大交换

当前证据状态：**部分可用**。可在有限维闭凸锥LP且拆分锥/有限值/PPS前提均满足时核查资格与面约化；一般minimax、测度或无限维问题不直接继承强对偶。P-5ae 的局部缩放采用独立检查的比值比例，并已登记源平方根印刷错误及精确标量复核；源定理整体证明未独立重证。专题验收仍在进行。

## 输入与产出

有限维锥线性原/对偶问题或拟交换模型、闭凸锥分解K=K₁×K₂、线性子空间L、目标值有限性、候选相对内点slack，以及需要恢复强对偶/最小面的目标。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 先判断目标问题确属有限维闭凸锥线性规划；若是一般minimax/测度问题，只将本文当作邻近面约化工具，不直接套强对偶结论。
2. 将锥明确拆为K₁×K₂并指出K₂多面；逐侧检查Partial Polyhedral Slater：对偶slack须有s₁∈ri K₁、s₂∈K₂，原始侧则检查对应的ri(K₁*)点。
3. 只有一侧最优值有限且另一侧满足对应PPS时，才按已证 Proposition 2 推零间隙和对侧达到；仅有可行点或弱Slater不能替代这些前提。
4. 若无相对内点，按partial-polyhedral separation寻找x∈K*∩L⊥：区分<c,x><0的不可行证书与<c,x>=0、非正交块外分量的约化方向。每步验证包含原可行slack的face。
5. 需要复杂度界时先计算各块ℓpoly(Kᵢ)，使用FRA-Poly的≤1+Σℓpoly(Kᵢ)步界；需最小奇异度时检查相对内部方向可实现性。
6. 对sup/inf交换另列问题空间、凸凹性、闭性、紧性和适用的minimax定理；本论文的PPS/面约化结果自身不提供交换。
7. P-08ccea72ef78920a的四类资格条件要按Y0=span(A(P)−Q)及原拓扑逐一选择；Fréchet方案需P,Q,Y0闭且b∈icr像锥，有限维方案需dimY0<∞与icr。仍有−∞分支；仅有限值分支给零gap及dual达到，weak*紧另需Y0=Y。其Example4印刷正gap位置已修为第三坐标0、第二坐标正；非紧序列逃逸与值lsc另核。
8. Hilbert generator来源P-5ae169908ed70fca的exact dual达到以共享normal及⟨v,x*⟩>0认证，逆向另需Ran(A*)闭；有限维非零共享normal或closed P都不充分。Prop4.1的源平方根比例为印刷错误，PDF p.22推导和精确标量例支持局部修正比例λ*/σ_K(A*y)；须满足σ_K(A*y)>0及源条件，不能把此局部修补提升成原定理整体已证。Prop4.2另需∀x∈P Span{x}+P闭。详细记录见references/verification-log-20261007.md。

## 证据与失败处理

不能把可行性当作PPS或完整Slater；强对偶另需一侧有限最优值。定理只覆盖有限维闭凸锥线性规划及一块多面分解，步数界不等于运行时界。Minimax、非凸、无穷维和一般测度问题必须重新证明适用条件。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
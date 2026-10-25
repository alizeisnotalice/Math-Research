---
name: math-b01-lagrange-farkas
description: "用于Lagrange 对偶与 Farkas 证书的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> 外部文献、claims/cases 与 PDF 属于完整证据包，不随单个 Skill 安装。请使用[便携证据访问指南](references/evidence-guide.md)，按本 Skill 的[来源索引](references/handoff-evidence-index.csv)通过统一 resolver 定位。registry/显式包根路由均不依赖当前工作目录；SHA 检查只确认当前字节，不等于全文阅读或数学验收。


# B01 · Lagrange 对偶与 Farkas 证书

当前证据状态：**部分可用**。有限维 Farkas/弱对偶证书及列出的已核条件可用；一般锥、无限维或 DC-composite 强对偶需另证资格/拓扑/紧性。两项原文错误须按本入口更正后再调用：P-18e Theorem 3.6 的 `∀μ∃λ` 不推出共同乘子；P-5ae Remark 1.7 的上半连续说法被反例否定，Prop.4.1 的平方根比例修为 `λ*/σ_K(A*y)` 且仅在源条件下采用。新增相邻来源 E1258 的 arXiv v2 第2.2节 Lemma 2 局部反例显示其放松问题的对偶有限性条件错误；不得调用其 top-k 特征值刻画，后续近似结论需分别审查。其余一般定理证明未独立重证，专题验收未完。

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
10. 非闭锥来源P-5ae169908ed70fca限Hilbert与bounded linear A、bounded closed convex generator K∋0：cl A(coneK)的零support方向测试只证明∀ε>0近似可解；exact A(coneK)另需一个共同C使∀y ⟨b,y⟩≤Cσ_K(A*y)。ε>0 dual唯一不推ε=0 dual达到；general support-face inclusion需同时回验残差/normal，H/E恢复原非凸锥是额外条件。该文Remark1.7的上半连续声明由ℓ²反例否定；本例只否定该声明，不据此声称一般下半连续。Prop4.1由PDF p.22核得的校正比例为 λ*/σ_K(A*y)，且需要σ_K(A*y)>0和源中共享法向/稳定性假设；平方根比例是原文错误。精确标量例与适用范围见references/verification-log-20261007.md。
11. 对只含齐次项 `v_j^T A_jv_j`、各 `v_j` 相互独立的二次 Lagrangian，先把全部等式乘子代入后再计算 `g=inf L`，不能只靠一阶驻点筛选。若 `A_j=μ_jI−M`，在域 `v_j≠0` 上有限下界的充要条件是每个 `A_j⪰0`：若存在负方向，沿其任意非零倍数趋于无穷可令该项趋于 `−∞`；若 `A_j⪰0`，该项下确界为0，但当 `A_j` 正定时下确界不达到。不要把有限下确界、达到的极小点和特征值驻点混为一谈。E1258 的局部核验见 `references/verification-log-20261007.md`。

## 证据与失败处理

有限维Farkas不证明一般锥/非凸强对偶。DC-composite原文存在已定位的量词缺口：逐μ选λ不是共同乘子，须修正或另证minimax/鞍点条件。 该有限LP扰动证明的共同小ε与固定working-set子序列已做局部补足，原文Step1计数须加终端检查；不得将rank假设移植为Farkas矩阵限制或无限维强对偶。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
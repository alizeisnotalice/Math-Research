---
name: math-d04-entropy-chain-rule
description: "用于熵链式分解的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **本轮审核状态（2026-10-07）**：Vigneaux 与 Peres–Quas 全文已读；Nishiyama–Sason 与量子链式论文由独审全文阅读。现行经典链式结论以下述标准 Borel 条件核恒等式为准；邻近论文只支持其各自的有限混合或有限状态/量子结论。其余候选来源筛查仍在进行。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# D04 · 熵链式分解

当前证据状态：**标准 Borel 概率的经典 KL 链式公式由独立推导支持；来源迁移边界已收窄**。不把几何、有限混合或量子来源当作一般经典链式证明。

## 输入与产出

联合分布P_XY,Q_XY、边缘与正规条件分布、KL方向及有限性。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 先选标准Borel等具有正规条件分布的空间，固定联合概率P、Q与KL(P||Q)方向；有限值推导先检P≪Q及必要可积性。
2. 在标准 Borel 空间取正规条件核 `P_{Y|X=x}`、`Q_{Y|X=x}`。KL 链式恒等式为
   `D(P_XY‖Q_XY)=D(P_X‖Q_X)+∫D(P_{Y|X=x}‖Q_{Y|X=x})P_X(dx)`。
   两项均非负，右侧按扩展实数加法取值，故不出现 `∞−∞`。有限值情形可令 `r=dP_X/dQ_X` 并在 `P_X`-几乎处处有 `s_x=dP_{Y|x}/dQ_{Y|x}`，从 `dP_XY/dQ_XY=r(x)s_x(y)` 展开对数后用 Tonelli/Fubini 得式；若 `P_X` 非绝对连续于 `Q_X`，左侧与第一项均为 `∞`；若条件核在正 `P_X` 质量集合上不绝对连续，则联合绝对连续失败，条件 KL 积分为 `∞`。该说明使用标准 Borel 正规条件核的存在性。
3. 不把上述等式改写成熵差，也不对不可积的 `log r` 或 `log s_x` 做形式相减。若只能在更一般的可测空间上给出条件对象、或无法构造所需条件核，应退回已验证的子类，或另引适用版本并记录条件。
4. 多层重复用真实条件参考Q_{X_j|X_{<j}}及P前史权；Q相关时不能换成Q的独立边缘。无限层需说明一致的过程律、sigma代数和相对熵极限定理，不从形式有限求和推极限。
5. P-8692a72e4fbe27ed全文仅邻近证据：有限混合 `M=Σα_iP_i` 的熵凹性亏损等于 `Σα_iD(P_i‖M)`，即标签–样本互信息；不是一般条件链式法则的文献证明。
6. 自然对数nats版本为D(P||Q)=∫₀¹χ²(P||(1−s)P+sQ)ds/s及D(μ_C||μ)=ln(1/μ(C))。来源任意共同log底数公式的log_e指log(e)，nats时为1；调用保留μ(C)>0及有限/扩展值约定，含Q_min的收缩上界不称分布一致。
7. 新全文 P-a89b2a1ccfdb7c4f 式(11)在给定 `(T,ξ)`-disintegration `ν_t` 及 `ρ=rν` 下，先设 `M(t)=∫r dν_t`、`P_T=Mξ`、`ρ_t=(r/M)ν_t`；有限值推导需 `log r` 与 `log M(Tx)` 可积，再展开 `H_ν(ρ)=H_ξ(P_T)+E_{P_T}H_{ν_t}(ρ_t)`。取 `ν=Q`、`ξ=Q_T`、`ν_t=Q_t` 才转为 KL 链式。该式是已有 disintegration 后的熵分解，不是一般条件核存在性来源；零 `P_T` 纤维无贡献。
8. 其式(14)限有限个不相交 rectifiable carriers，`μ=Σμ_i`、`ρ=Σq_iρ_i` 且各熵项有限（或至少单边可积并保证右侧总和有定义），此时 `H_μ(ρ)=H(q)+Σq_iH_{μ_i}(ρ_i)`。标签唯一本身不排除某层微分熵 `+∞`、另一层 `−∞` 导致未定义的 `∞−∞`；来源不据此支持一般扩展值式。重叠混合不能默认样本决定唯一标签。Area/coarea 应用必须另检 Lipschitz、正切 Jacobian/rank 几乎处处为正及相应可积性；来源 p.5 的正 coarea 因子与 p.6 任意 Borel pushforward AC 均有常值映射反例。
9. 该来源p4曲线像长度式漏multiplicity；p8指定strong typical集合的dimension窗口漏stratum维数常数，m=(0,10)混合反驳系数1。修复窗口保nηΣm_i或已证更锐常数、0<ξ<1/2及外引AEP依赖，不凭full_read认证原陈述。
10. 邻近量子源P-c5a8ba8d74c8e1a5只有finite-dimensional channel链式不等式；Thm3.5需E TPCP、F CP及Dmax(E||F)<∞，加项为non-stabilized barDreg。来源base2转nats乘ln2。单letter替代、sup over任意finite R的达值和Prop3.1图像global上界均不默认；smooth proof调用限定ε<1且mε+sqrt(mε)+ε′<1。
11. 中心立方体捕获来源/层/赢家标签必须由同一联合律产生；返回逐条件熵账和零质量处理，不借邻近混合公式认证私有分布接口。

## 证据与失败处理

链式公式的权是P_X，参考是Q_{Y|X}；来源χ²论文只证明混合标签熵亏损及标量散度关系。无限层、非标准空间或零质量缺适用定理时留缺口，不将阅读全文计为一般链式全证明认证。 新disintegration来源的几何应用存在正Jacobian与Borel退化缺口，已用原PDF和显式反例确认；stratified AEP指定窗口漏维数常数。只使用已列有限值修复，不默认原公式字面全真。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

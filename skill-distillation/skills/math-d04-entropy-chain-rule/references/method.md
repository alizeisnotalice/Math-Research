# 熵链式分解：方法与验收边界

1. 先选标准Borel等具有正规条件分布的空间，固定联合概率P、Q与KL(P||Q)方向；有限值推导先检P≪Q及必要可积性。

2. 在第1步条件下写正规条件核 `P_{Y|x}`、`Q_{Y|x}`。有限值推导若 `P≪Q` 且 `D(P‖Q)<∞`，令 `r=dP_X/dQ_X` 并对 `P_X`-几乎处处的条件核取 `s_x=dP_{Y|x}/dQ_{Y|x}`；则联合 RN 密度为 `r(x)s_x(y)`，按 `P_X` 权重展开对数并积分，得
   `D(P_XY‖Q_XY)=D(P_X‖Q_X)+∫D(P_{Y|x}‖Q_{Y|x})P_X(dx)`。
   零 `P_X` 质量纤维不贡献，参考核必须是 `Q_{Y|x}`。

3. 扩展值版本按非负 KL 项作和，不对熵或 KL 作差。标准 Borel 条件下正规条件核存在，且 `Q_{Y|X}` 仅需按 `Q_X`-几乎处处定义；若 `P_X≪Q_X`，该版本也按 `P_X`-几乎处处确定。若 `P_X` 不绝对连续于 `Q_X`，则 `D(P_XY‖Q_XY)=D(P_X‖Q_X)=+∞`。若 `P_X≪Q_X` 但 `P_{Y|x}` 不绝对连续于 `Q_{Y|x}` 的集合有正 `P_X` 质量，则条件 KL 积分与联合 KL 均为 `+∞`。在剩下情形，由 disintegration 可知 `P_XY≪Q_XY`：对任意 `Q_XY`-零集，`Q_{Y|x}`-纤维零的 x 是 `Q_X`-零，亦为 `P_X`-零；其余纤维由 `P_{Y|x}≪Q_{Y|x}` 消去。取联合 RN 密度 `f=dP_XY/dQ_XY` 的可测版本，令 `r=dP_X/dQ_X`。再按 `Q_X` 与 `Q_{Y|x}` 分解 `Q_XY`，得 `r(x)=∫f(x,y)Q_{Y|x}(dy)` 对 `Q_X`-几乎处处成立；将 `P_X(dx)P_{Y|x}(dy)=f(x,y)Q_X(dx)Q_{Y|x}(dy)` 与 `P_X=rQ_X` 比较，并用标准 Borel 空间上正规条件分布的 `P_X`-几乎处处唯一性，在 `r>0` 的 `P_X`-几乎处处纤维上得到 `s_x(y)=f(x,y)/r(x)=dP_{Y|x}/dQ_{Y|x}`。因此联合密度为 `f=rs`，对数密度为 `log r+log s_x`。概率密度的负部积分由 `∫_{0<f<1}f|log f|dQ≤1/e` 控制；对条件核逐个用同一界再按 `P_X` 积分，故两个对数项各自的负部可积，正部允许为 `+∞`，扩展积分逐项有定义，且 `E_P(log r+log s_X)=E_{P_X}log r+E_{P_X}D(P_{Y|x}‖Q_{Y|x})`。这给出扩展链式式而不产生 `∞−∞`。任何不可积的熵差、非标准 Borel 空间或缺少可测核 RN 分解时，退回已验证有限域或另提供适用定理与条件。
4. 多层重复用真实条件参考Q_{X_j|X_{<j}}及P前史权；Q相关时不能换成Q的独立边缘。无限层需说明一致的过程律、sigma代数和相对熵极限定理，不从形式有限求和推极限。

5. P-8692a72e4fbe27ed全文仅邻近证据：有限混合M=Σα_iP_i的熵凹性亏损等于Σα_iD(P_i||M)，即标签–样本互信息；不是一般条件链式法则的文献证明。

6. 自然对数nats版本为D(P||Q)=∫₀¹χ²(P||(1−s)P+sQ)ds/s及D(μ_C||μ)=ln(1/μ(C))。来源任意共同log底数公式的log_e指log(e)，nats时为1；调用保留μ(C)>0及有限/扩展值约定，含Q_min的收缩上界不称分布一致。

7. 全文来源 P-a89b2a1ccfdb7c4f 的 Eq.(11) (PDF p.5) 在给定 `(T,ξ)`-disintegration `ν_t` 及 `ρ=rν` 下，先设 `M(t)=∫r dν_t`、`P_T=Mξ`、`ρ_t=(r/M)ν_t`；有限值推导检 `log r` 与 `log M(Tx)` 可积，再展开 `H_ν(ρ)=H_ξ(P_T)+E_{P_T}H_{ν_t}(ρ_t)`。取 `ν=Q`、`ξ=Q_T`、`ν_t=Q_t` 才转为 KL 链式。零 `P_T` 纤维无贡献。

8. 同一来源 Eq.(14) (PDF p.7) 限有限个互不相交 rectifiable carriers：`μ=Σμ_i`、`ρ=Σq_iρ_i` 且所有熵项有限（或至少单边可积并保证右端有定义）时，`H_μ(ρ)=H(q)+Σq_iH_{μ_i}(ρ_i)`。重叠混合不能默认样本决定唯一标签；仅有标签有限也不能排除 `+∞−∞`。Area/coarea 应用必须另检 Lipschitz、正 tangent Jacobian/rank 和相应可积性；来源 p.5 的正 coarea 因子与 p.6 任意 Borel pushforward AC 均有常值映射反例。

9. 该来源p4曲线像长度式漏multiplicity；p8指定strong typical集合的dimension窗口漏stratum维数常数，m=(0,10)混合反驳系数1。修复窗口保nηΣm_i或已证更锐常数、0<ξ<1/2及外引AEP依赖，不凭full_read认证原陈述。

10. 邻近量子源P-c5a8ba8d74c8e1a5只有finite-dimensional channel链式不等式；Thm3.5需E TPCP、F CP及Dmax(E||F)<∞，加项为non-stabilized barDreg。来源base2转nats乘ln2。单letter替代、sup over任意finite R的达值和Prop3.1图像global上界均不默认；smooth proof调用限定ε<1且mε+sqrt(mε)+ε′<1。

11. 中心立方体捕获来源/层/赢家标签必须由同一联合律产生；返回逐条件熵账和零质量处理，不借邻近混合公式认证私有分布接口。

## 不可省略的限制

链式公式的权是P_X，参考是Q_{Y|X}；来源χ²论文只证明混合标签熵亏损及标量散度关系。无限层、非标准空间或零质量缺适用定理时留缺口，不将阅读全文计为一般链式全证明认证。 新disintegration来源的几何应用存在正Jacobian与Borel退化缺口，已用原PDF和显式反例确认；stratified AEP指定窗口漏维数常数。只使用已列有限值修复，不默认原公式字面全真。

# 有限离散KL链式：已检查基础推导

设P,Q是有限X×Y上的概率，P≪Q。对P_X(x)>0的x，
p(x,y)/q(x,y)=[p_X(x)/q_X(x)]·[p(y|x)/q(y|x)]。
乘p(x,y)后对x,y求和，得到
KL(P_XY||Q_XY)=KL(P_X||Q_X)+Σ_x p_X(x)KL(P_{Y|x}||Q_{Y|x})。
P_X(x)=0的项按0计。条件期望权重是P_X，参考是Q_{Y|x}。
连续空间推广需要条件分布、Radon–Nikodym分解与积分条件，不能由有限证明自动完成。

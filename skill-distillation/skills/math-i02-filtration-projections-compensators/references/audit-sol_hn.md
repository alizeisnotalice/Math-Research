# I02 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-d69f9d48286d2a25 Thm3.1 lines297–517是single-jump filtration下compensator定义及判据；不是所有filtration的自动通式。源Thm4.1 numerator符号卡中有疑点保留。

## 证明骨架与已检查推导

给定filtration与predictable test H→验证E∫H dN=E∫H dA→局部鞅；projection、UI及stopping按各自conditions使用。

本次实际检查：γ~Exp1，Nt=1γ≤t，At=t∧γ，T1两期望同1−e^{−1}=.63212。

## 正反例与失败处理

正例：令 γ∼Exp(1)，N_t=1_{γ≤t} 为首次跳计数，A_t=t∧γ=∫₀ᵗ1_{s≤γ}ds 是其可预测补偿子；在 T=1 时，E N₁=P(γ≤1)=1−e⁻¹≈0.63212，E A₁=∫₀¹P(γ≥s)ds=∫₀¹e⁻ˢds=1−e⁻¹≈0.63212，可逐式核对。

条件缺失例：以放大后的滤过引入未来信息后仍沿用旧补偿 λdt，或对非 UI 鞅任意使用期望版停时定理。

## 中心立方体迁移推断

中心立方体迁移接口：先由原始样本及尺度历史定义 F_r，再明确哪些 cube-indexed 量是可选或可预测；若随机选中心/半径，须给出相对该滤过的补偿子，不能默认存在自然 Poisson 补偿。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

扩大filtration公开未来后compensator会改变；optional/current-jump integrand不自动predictable。

补偿子依赖滤过；可选投影与可预测投影不同，局部鞅也未必是真鞅。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。

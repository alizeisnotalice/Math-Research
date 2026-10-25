# 非局部 balayage 与退出测度：方法与验收边界

1. 对右连续过程与开集 D，定义 τ_D=inf{t≥0:X_t∉D}；有限退出时记录 X_{τ_D}∈D^c，允许它落在 D^c\∂D。

2. 给定寿命 ζ，定义存活到几何退出的测度 ω_D(x,A)=P_x(τ_D<ζ,X_{τ_D}∈A)（A⊂D^c 可测）；将无限不退出与先行/同时杀死质量分列。

3. 对调和/势函数应用 Dynkin 公式或非局部 Dirichlet 问题时核验测试域、可积性和边界数据在整个 D^c 上的定义。

4. 比较 balayage 前确认所用核、边界条件与概率退出测度一致；检查总质量仅在几乎必然退出且不被杀死时为1。

## 局部证明记录

右连续路径的退出时刻与 jump overshoot 示例、退出测度质量和判据见 [J03 steps 1–2 local derivation](evidence-guide.md#j03-steps-12-local-derivation)。Kuntz et al. 的下界、收敛和密度结论严格限于其可数状态CTMC假设。

## 局部证明记录

右连续路径的退出时刻与 jump overshoot 示例、退出测度质量和判据见 [J03 steps 1–2 local derivation](evidence-guide.md#j03-steps-12-local-derivation)。Kuntz et al. 的下界、收敛和密度结论严格限于其可数状态CTMC假设。

## 不可省略的限制

非局部跳跃会越过几何边界；仅给 ∂D 数据可能不足以定义非局部 Dirichlet 问题。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。

局部证据补充见 [SOL_HN记录](audit-sol_hn.md)；其审读者与范围以记录为准。

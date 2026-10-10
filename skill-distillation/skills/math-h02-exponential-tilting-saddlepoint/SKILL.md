---
name: math-h02-exponential-tilting-saddlepoint
description: "用于指数倾斜与鞍点尾近似：先核验模型、鞍点和误差范围，再选择精确倾斜恒等式、Gaussian Mills 界或有来源的渐近式。"
---

# H02 · 指数倾斜与鞍点

当前证据状态：**部分可用**。一般倾斜恒等式与标准 Gaussian Mills 界可直接推导。Niu–Ray Choudhury–Katsevich 的已深读结果是一个条件独立数组的 Lugannani–Rice 相对误差定理；在该定理之外，简化的 `u⁻¹` 前因子仍须单独核对模型和渐近区间。具体结论、条件和来源定位见[来源与阅读状态](references/provenance.md)及[方法说明](references/method.md)。旧入口已存档于 `references/history/SKILL.pre-evidence-revision-20261007.md`。

## 输入与产出

输入随机变量或和 `S`、对数矩母函数 `K(s)=log E exp(sS)`、尾方向与阈值 `x`，以及独立性/条件独立性、格性、矩母函数定义域和可用来源。产出必须区分恒等式、概率界、有限样本近似与相对误差渐近，并列出条件、常数、方向和缺口。

## 执行步骤

1. 确认 `K` 在含 0 与候选鞍点的开区间内有限且二阶可微。解 `K′(ŝ)=x`，核实解存在、唯一性、符号与尾方向；检查 `K″(ŝ)>0`。无根或退化时停止鞍点近似。
2. 先给出精确换测恒等式
   \[
   \mathbb P(S\ge x)=e^{K(\hat s)-\hat s x}\,\mathbb E_{\hat s}\!\left[e^{-\hat s(S-x)}\mathbf1_{\{S\ge x\}}\right]
   \]
   （上尾且 `ŝ>0`）。此式本身不是 Gaussian 局部近似，也没有推出 `u⁻¹` 前因子。
3. 若输入是 Niu 等人的条件独立数组，按[方法说明 §2](references/method.md#2-已读的条件-lugannanirice-定理原文原条件)逐项核对 Theorem 1：`W_{in}` 条件于 `F_n` 独立、条件均值为零；CSE 或 CCS 二选一；`n^{-1}Σ_i E(W_{in}²|F_n)=Ω_P(1)`；`F_n`-可测平均阈值 `w_n=o_P(1)`。记 `K_n(s)=n^{-1}Σ_i log E(e^{sW_{in}}|F_n)`、总和对数 MGF `Λ_n(s)=nK_n(s)` 和总和阈值 `x_n=nw_n`。先在来源保证的高概率唯一根事件上求 `K_n′(ŝ_n)=w_n`，再定义
   \[
   \lambda_n=\hat s_n\sqrt{\Lambda_n''(\hat s_n)},\qquad
   r_n=\operatorname{sgn}(\hat s_n)\sqrt{2\{\hat s_nx_n-\Lambda_n(\hat s_n)\}},
   \]
   并先报告完整条件 LR 式 `Q(r_n)+φ(r_n)(1/λ_n−1/r_n)` 的相对误差 `1+o_P(1)`；不能把它静默改成简单前因子。CSE/CCS 的精确条件见方法说明和来源卡，不用“存在指数矩”替换其条件。
4. 仅在上述 Theorem 1 条件下，再要求上尾 `w_n>0`（以概率趋于 1）及 `r_n→_P+∞`，才能使用该来源 Appendix H 的 `λ_n/r_n→_P1` 与 Mills 界，把完整 LR 式另行简化为 `exp{Λ_n(ŝ_n)−ŝ_nx_n}/[ŝ_n√(2πΛ_n″(ŝ_n))]` 的相对渐近式。其总误差率未量化；不覆盖 `r_n=O_P(1)`、固定非零阈值大偏差或任意依赖和。离散数组在该特定 shrinking-cutoff 定理中不需额外非格点假设，但这不产生有限样本格点无误差保证。
5. 标准正态上尾 `Q(u)=1−Φ(u)` 可单独用 Mills 界：对 `u>0`，`φ(u)u/(u²+1)<Q(u)<φ(u)/u`。所以 `A(u)=φ(u)/u` 是上界，满足 `0<(A−Q)/A<1/(u²+1)`，且 `0<A/Q−1<1/u²`；只有 `u→∞` 才有 `Q(u)∼A(u)`。这是 Gaussian 的精确界和渐近误差方向，不是一般倾斜分布的有限样本结论。下尾先对变量变号，并重新核验参数。

## 失败处理

鞍点不存在/不唯一、`K″=0`、矩母函数越界、方向符号不匹配、阈值不满足来源条件或 LR 与简化式的适用区间未核实时，报告精确恒等式或已证明的弱结论并指出缺项。数值案例只检查个别输入，不能代替渐近证明。用[适用案例](examples/positive.md)和[条件缺失案例](examples/negative.md)检查边界。

H–I 操作步骤与来源定理的逐项证据、证明状态及未覆盖接口见[H–I 操作证据与共享审计文件定位说明](references/portable-audit-access.md)。

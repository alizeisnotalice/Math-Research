# 独立复核记录（append-only，2026-10-06）

本文件由外部审核批次写入，**不修改**本 skill 的任何既有内容。
记录的是对该 skill 具体断言的独立复核：复核强度、精确定位、以及诚实的分级。
分级含义：

- **A 完整证明**：含取等/紧性/分支/退化/达到性的完整论证；
- **B 附加条件下成立**：在明确假设下推导完整，假设本身未独立证明；
- **C 仅数值/实例检查**：随机采样无违反、暴力未超界、误差很小——**不构成**一般命题的证明；
- **R 文献确认**：查证权威来源，非本审核的证明。


## 复核分级：C

- EN₁=EA₁=1−e⁻¹=0.63212：闭式精确。补偿子恒等式的一般理论未验（标准结果）；案例真正内容是 F₀ 揭示 γ 后 γ 变为可预测停时、补偿子不再是 t∧γ——该结构性陈述未独立证明。

## 复核者与范围声明

- 复核批次：R8–R11（2026-10-06）；本文件为 append-only 记录，原文件未改动。
- 分级修正：早前汇总中『全部成立』的表述过强，已按 A/B/C/R 重新分级；
  数值无违反 ≠ 一般命题已证。

---

## 更正条目（2026-10-06，append-only）

上述"复核分级：C"为**单档标记，已过时**。本 skill 的断言已拆分为带唯一 ID 的子断言
（登记簿：`audit_r12/assertion_register.json`），各子断言独立定级：

- **I02-CASE-01**（instance_check_only，numeric_experiment/standard_theorem_citation）：EN_1=EA_1=1-e^{-1}=0.63212 for A_t=min(t,gamma)

**统计口径更正**：此前汇总"累计20条全部成立"是错误计数——原始验证对象为 24 个 skill 级
对象（29 个验证动作），拆分后 **34 条子断言**，分布：independent_proof 10、
conditional_proof 7、instance_check_only 13、pending 4。随机采样无违反、暴力未超界、
LP 求解成功、误差极小均**不构成**一般命题的证明；"全部成立"的表述撤回。
本日志的分级字段以本更正条目及登记簿为准（验证状态 × 依据类型双维）。

---

## 档位提升（R26）

依据既有验证记录：EN₁=EA₁=1−e⁻¹ 闭式精确（instance_check_only）
待证据 → **部分可用**。源定理证明本身未重证（保持守卫）。

---

## 更正条目（R47）：P-2f81500e5dfa70c6 卡条件表述偏差（卡↔PDF 比对发现）

**发现**：本 skill 引用的源卡 HK-P-2f81500e5dfa70c6（Oertel, On Jump Measures of Optional
Processes with Regulated Trajectories, arXiv:1508.01973）经 PDF 逐字比对，定理陈述主体
（Prop 1 / Thm 4 / Thm 5 / Cor 1）逐字一致，但两处条件表述夹带了非原文内容：

1. Prop 1 assumptions 原写 "0∉A, hence is bounded away from zero"——PDF p4 仅假设
   0∉A，"hence" 推断不成立（A={1/n} 反例）；原文证明中 lim∈Ā⇒ΔX_{t*}≠0 一步实际
   需要 0∉closure(A)，系原文证明 gap，卡以此隐式修补但未标注。
2. 定义卡 "jump measure" 中 "A bounded away from zero" 的有限性条件——原文 Lemma 2
   仅写 A∈B*（0∉A）；该条件是必要的数学修补（0∉A 不足以保证 N^A_X(t)<∞），但
   原文 p5 "A⊆R\(−ε,ε) for all ε>0" 量词系原文疏漏，卡的表述应标注为卡的澄清。

**处置**（append-only）：三份卡副本的原文均已保留至 statement_original_r47 /
assumptions_original_r47 字段，新表述含【修订 2026-10-06 R47】内联标注。本 skill 的
SKILL.md 未转录这些假设表述（工作流性质），无需修订正文。使用该论文的迭代枚举
构造时，务必显式假设 0∉closure(A)。

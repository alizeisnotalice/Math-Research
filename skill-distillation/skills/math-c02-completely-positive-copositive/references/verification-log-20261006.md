# 独立复核记录（append-only，2026-10-06）

本文件由外部审核批次写入，**不修改**本 skill 的任何既有内容。
记录的是对该 skill 具体断言的独立复核：复核强度、精确定位、以及诚实的分级。
分级含义：

- **A 完整证明**：含取等/紧性/分支/退化/达到性的完整论证；
- **B 附加条件下成立**：在明确假设下推导完整，假设本身未独立证明；
- **C 仅数值/实例检查**：随机采样无违反、暴力未超界、误差很小——**不构成**一般命题的证明；
- **R 文献确认**：查证权威来源，非本审核的证明。


## 复核分级：B

- CP⊆COP*（⟨A,Σwwᵀ⟩=ΣwᵀAw≥0 对 A∈COP, w≥0）：一行代数，数值复核通过。
- **降级**：『CP=COP*』等式的反向 COP*⊆CP 是 Diananda 型深结果，属标准定理引用，本审核未验证。『∩C_d=COP』skill 自己声明依赖外引矩增长定理——与复核范围一致。

## 复核者与范围声明

- 复核批次：R8–R11（2026-10-06）；本文件为 append-only 记录，原文件未改动。
- 分级修正：早前汇总中『全部成立』的表述过强，已按 A/B/C/R 重新分级；
  数值无违反 ≠ 一般命题已证。

---

## 更正条目（2026-10-06，append-only）

上述"复核分级：B"为**单档标记，已过时**。本 skill 的断言已拆分为带唯一 ID 的子断言
（登记簿：`audit_r12/assertion_register.json`），各子断言独立定级：

- **C02-IDEN-01a**（independent_proof，own_derivation/numeric_experiment）：CP subset COP and CP subset COP* (one-directional inclusions)
- **C02-IDEN-01b**（pending，standard_theorem_citation）：CP = COP* (equality, reverse inclusion COP* subset CP)
- **C02-IDEN-01c**（independent_proof，own_derivation）：PSD basis size binom(n+d,d)

**统计口径更正**：此前汇总"累计20条全部成立"是错误计数——原始验证对象为 24 个 skill 级
对象（29 个验证动作），拆分后 **34 条子断言**，分布：independent_proof 10、
conditional_proof 7、instance_check_only 13、pending 4。随机采样无违反、暴力未超界、
LP 求解成功、误差极小均**不构成**一般命题的证明；"全部成立"的表述撤回。
本日志的分级字段以本更正条目及登记簿为准（验证状态 × 依据类型双维）。

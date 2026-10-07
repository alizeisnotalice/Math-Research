# 独立复核记录（append-only，2026-10-06）

本文件由外部审核批次写入，不修改本 skill 任何既有内容。
分级：A 完整证明 / B 附加条件下成立 / C 仅数值实例检查 / R 文献确认 / G 证据不足。

## 复核分级：C（案例级，R13 批量执行通过）

- **案例执行**：tower property E[E[X|G]|F]=E[X|F] for F⊂G
- 结果：通过
- **引用忠实性**：本 skill 无编号定理引用；引用为内容级（式/常数/假设），与源卡陈述方向一致，未发现矛盾

## 档位提升依据（R16）

待证据 → **部分可用**。依据：案例实际执行通过 + 定理引用忠实性核对（内容级一致性）。
仍待证据：源定理证明本身未重证；专题迁移主张保持原有守卫声明。
验证强度如实记录：案例为**实例级**检查，不构成一般命题证明。

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

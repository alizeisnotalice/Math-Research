---
name: math-g01-carleson-measure-box-and-tent-test
description: "用于Carleson 测度的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **交接准备状态（结构性）**：本副本只整理技能结构与引用可移植性；未进行全文深读或数学验证。原评级、可用状态、假设、证明限制和未决条件均沿用原件。结构检查不构成科学验证。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# G01 · Carleson 测度

当前证据状态：**待证据验收：可执行方法工作流，不是已证明的课题结论**。全文转换、论文阅读和证明核验是三个独立状态。（此行声明已被下方【修订】标注替代：该 skill 经外部审核提升为「部分可用」，详见对应修订行与 references/verification-log-20261006.md）
  【修订 2026-10-06（外部审核 R16）】案例已实际执行通过（实例级），定理引用忠实性核对无矛盾；证据状态提升为：部分可用：基础推导已检查，专题证据仍待完整验收。源定理证明本身未重证，验证强度为实例级。

## 输入与产出

底空间、上半空间/帐篷几何、正 Borel 测度 μ、盒/帐篷族和目标嵌入。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 路由适用性

本入口保留 partial：有限树路由只在显式树、底测度 σ、全部节点系数 α 与子分割给定后检逐节点 packing；当前叶例只产出 packing 常数 1，未提供 f、嵌入输出及匹配来源定理的全套假设，不能作为嵌入正例。UR/PDE 路由另需匹配的域 Ω、调和函数 u、实际局部能量/边界几何输入及指定方向的 UR 或 corkscrew+CDC 条件；缺这些输入时停止该路由，不从叶 packing 改称 PDE 定理适用。

## 执行步骤

1. 先定义实际T(Q)/S(Q)、底测度、尺度及边界，验证正Borel测度局部有限；写sup_Q μ(T(Q))/σ(Q)，不能只验根盒。
2. 有限过滤树的系数packing检每个J的Σ_{I⊆J}α_I≤Cσ(J)；有互不交子分割与统一mass流时可用已检查4C嵌入构件，树外多参数盒条件须另证。
3. 解析/调和空间的嵌入与反向kernel测试按指定空间定理，先明确核范数归一化；普通盒packing不自动证明所有reproducing-kernel thesis。
4. 邻近P-e9e9211495f0b20e限Ω⊂R^{d+1}有界调和函数局部能量r^{-d}∫_{B(x,r)∩Ω}|∇u|²dist(y,∂Ω)dy≤C||u||∞²；不是任意度量Carleson测试。
5. 其Theorem1.1正向需存在Ω̃⊂Ω、∂Ω⊂∂Ω̃且∂Ω̃为UR；逆向另需corkscrew+CDC。Theorem1.2需两两不交边界E_j及ω(p_j,E_j)≥1−ε，得Σdist(p_j,∂Ω)^d≤C(ε)R^d。
6. 作者Whitney近/远边界拆分→调和事件packing→corona/UR构造的常数依赖d、几何参数及ε；外引容量/UR理论未独立复证，不能称维数无关。
7. 中心立方体迁移须实际树化或帐篷化、共同输入及packing桥梁；区分有界、消失、加权与算子特定Carleson，不从邻近边界PDE定理得全尺度最大界。

## 证据与失败处理

根packing、抽象单树packing、解析kernel测试与边界调和函数Carleson估计有不同前提。邻近来源含UR/corkscrew/CDC与ε依赖，不能把缺Ahlfors假设读成无几何假设。

复用 P-e8f8d143f231a4b3 前，必须先读 [有效来源合同与函数空间纠正](references/source-contract-P-e8f8.md)，再按 active 来源卡核条件；不得用历史 little Bloch 误名替代原文 B˚={f∈B:f(0)=0}。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [有限树 packing 示例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

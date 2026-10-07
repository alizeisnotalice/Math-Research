---
name: math-g05-tree-carleson-square-function-tc-a4-audit
description: "用于TC-A4 Carleson 型平方函数界的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **交接准备状态（结构性）**：本副本只整理技能结构与引用可移植性；未进行全文深读或数学验证。原评级、可用状态、假设、证明限制和未决条件均沿用原件。结构检查不构成科学验证。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# G05 · TC-A4 Carleson 型平方函数界

当前证据状态：**待证据验收：可执行方法工作流，不是已证明的课题结论**。全文转换、论文阅读和证明核验是三个独立状态。（此行声明已被下方【修订】标注替代：该 skill 经外部审核提升为「部分可用」，详见对应修订行与 references/verification-log-20261006.md）
  【修订 2026-10-06（外部审核 R24）】p=2 嵌入常数 4C 数值验证（小树 50 例）；Theorem 1.1 引用一致、弱构件链正确；证据状态提升为：部分可用：基础推导已检查，专题证据仍待完整验收。Bellman sharpness 与 TC-A4 定义待补。

## 输入与产出

树 T、节点测度/权、系数 α_I、平均或投影算子、平方函数，以及 TC-A4 的正式接口定义。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. TC-A4的系数归一化、树序、测度与平方函数缺正式定义时，只返回缺口及标准模型，不声称私有结论。
2. 标准二进模型α_I≥0且每个J均Σ_{I⊆J}α_I≤C|J|；所有子树都需验证，根packing不足。
3. P-b7592e844240da1c的作者Theorem1.1：p>1有Σα_I|<f>_I|^p≤(p′)^p C||f||_p^p；p=2给4C，p=1强型不包含。
4. 作者sharp常数与Bellman remodeling源自其具体二进/infinitely-refining过滤模型；本次未独立证明全凹性或锐性，不将任意滤过都称同Bellman。
5. 已检查弱构件：一般测度有限过滤树有互不交孩子、统一mass流及所有子树packing≤Cσ(J)时，极大节点层蛋糕+Doob L²给嵌入≤4C||f||²。此上界不需Lebesgue或小边界；sharpness和私有树适用性另证。
6. 平方函数若定义使积分正好为上述非负嵌入和才用该界；多个算子的交叉项另证正交/几乎正交。
7. 有限深度、零测度节点和无限极限分开；中心立方体选族需先证明真实公共输入、树化和packing，不能直接替代任意中心上确界。
8. 新增复用P-6a781aff4000e334 stochastic常数e限连续平方可积、orthogonal/equal-bracket驱动及generalized CR系统，alpha条件tail mass≤1且其closed martingale在驱动stable subspace可表示；不能当一般Carleson树常数e，更未提供TC-A4。

## 证据与失败处理

目录未定义 TC-A4。标准二进 Carleson 嵌入不能冒充对私有系数规范或平方函数的精确结论。

复用 P-e8f8d143f231a4b3 前，必须先读 [有效来源合同与函数空间纠正](references/source-contract-P-e8f8.md)，再按 active 来源卡核条件；不得用历史 little Bloch 误名替代原文 B˚={f∈B:f(0)=0}。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

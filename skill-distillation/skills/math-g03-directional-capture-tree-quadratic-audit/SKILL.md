---
name: math-g03-directional-capture-tree-quadratic-audit
description: "用于方向性捕获树二次型与来源算子的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前状态（2026-10-07）**：仅可执行明确输入的有限 Hermitian/PSD 矩阵算术。私有捕获树、方向、可行叶、来源到 Gram 矩阵的映射及共享输入谓词未定义；缺任一项时返回 `evidencepending`，不声称项目捕获估计。12 项 G 组案例中本 Skill 两例的结果和源码散列见 `audit_current/g/cases.json`。
> 外部文献证据随完整证据包提供；包根须同时含 `workspace_revised/`、`audit_current/` 与 `evidence/`。`EVIDENCE_ROOT` 和 `MATH64_SOURCE_PACKAGE_ROOT` 均指包根。单独安装的 Skill 不含共享 PDF/claims/cases；按[本地访问说明](references/portable-audit-access.md)从 registry root `math64-20261007` 解析包根以定位绝对 resolver，或显式提供包根。用 `--skill-dir` 与来源 ID 核验来源 SHA、用 `--artifact-path` 定位 claims/cases 并核验当前 SHA。registry、root 或 resolver 缺失时报告配置缺口，不猜路径/内容。resolver 的 PASS 只验证路径和当前字节，不代表全文阅读、命题证明或数学验收。


# G03 · 方向性捕获树二次型与来源算子


## 输入与产出

有限根树、方向标签、节点/来源权重、捕获集合或算子、显式矩阵/核 A 和二次目标。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 先索取树、方向、可行捕获叶、来源到向量或Gram矩阵的映射及共用输入谓词的正式定义；任一缺失时停止捕获适用性判断和课题界的推导，标为 evidencepending。此时仅可执行明确给出的有限 Hermitian 矩阵算术。
2. 有限实例中逐项建立来源向量 x 和矩阵 A，列明对角项及不同来源的交叉项。
3. 计算 Q(x)=x*Ax；声称半正定时用精确主子式或认证特征值界核验 Hermitian 矩阵。
4. 已检查弱反例：N阶全1的PSD Gram矩阵每个对角1，最大特征值N、ones二次型N²；对角界不能控制总交叉项。这里只是一般矩阵反例，未声称它可由真实TC-A14几何实现。
5. 另查树约束、子节点覆盖、来源质量守恒和方向兼容是否由同一输入满足。
6. 在小树上与穷举对比；区分外放松上界与真实可实现的捕获构型。
7. 复用 [P-e8f8 冻结来源合同](references/source-contract-P-e8f8.md) 时，只限单位圆盘、有限正 Borel 测度及原点归零 Bloch 子空间 B˚={f∈B:f(0)=0} 到 L²(μ) 的来源结果；须核每层有限嵌套面积分割、所有 L¹ 条件期望极限、固定排序、实际 reduced Bergman 投影和有限 packet-Gram 表示。外引 Bergman 表示、SDP duality、complex Grothendieck 仍 open；来源结果不提供私有捕获映射。

## 证据与失败处理

来源树算子是未定义的项目私有接口。一般二次型或局部树估计不能定义捕获几何，也不能证明目标全局界。

复用 P-e8f8d143f231a4b3 前，必须先读 [有效来源合同与函数空间纠正](references/source-contract-P-e8f8.md)，再按 active 来源卡核条件；不得用历史 little Bloch 误名替代原文 B˚={f∈B:f(0)=0}。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [算术示例（捕获待证据）](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

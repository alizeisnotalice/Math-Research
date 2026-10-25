---
name: math-l04-mathematical-programming-workflow
description: "用于数学规划软件的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# L04 · 数学规划软件

> **当前核验范围（2026-10-07）：** 本 Skill 仍属于 L–N 组进行中的审核；逐条数学与能力状态记录在完整交付包的 `audit_current/ln/claims.json`，案例定义记录在 `audit_current/ln/cases.json`；单独安装目录不包含这两份包级账本。请按 [`references/portable-audit-access.md`](references/portable-audit-access.md) 通过 bundle-root resolver 访问。以下限定只陈述已绑定账本的证据边界，不代表整项验收完成。

> 便携文献定位见 [本 Skill 证据索引](references/handoff-evidence-index.csv)。若需读取包外 PDF，以数据包根目录为 `EVIDENCE_ROOT`，按索引中的相对路径定位。

> **本项边界：** 账本 L04-SOLVER-01 记录 SciPy/HiGHS 对单个两变量 LP 的真实 API 调用，并以 Fraction 检验原/对偶可行和零间隙。此结果不外推到一般 LP 误差界、MIP、非线性优化或其他求解器接口。

## 输入与产出

决策变量及其取值域、目标函数、完整约束、优化方向、求解器与容差、原/对偶或不可行证书。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 将数学问题逐项翻译成变量、目标和约束；记录每个放松及其界的方向。
2. 检查单位、定义域、凸性或整数性，并在手算可验证的小例子上检查可行性。
3. 运行适用的 LP、凸规划、混合整数或非线性求解器，保存状态、容差和见证。
4. 独立计算返回见证的目标值，并按问题类型核验对偶界、整数性或不可行证书。
5. 无法精确验证时，将结果称为数值候选，不要把求解器状态改写成定理。

## 证据与失败处理

求解器只处理输入的模型。约束放松、方向写反、容差或局部最优都可能改变数学结论。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

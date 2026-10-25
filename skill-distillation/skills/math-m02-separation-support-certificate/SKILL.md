---
name: math-m02-separation-support-certificate
description: "用于Hahn–Banach 分离与支撑函数证书的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# M02 · Hahn–Banach 分离与支撑函数证书

> **当前核验范围（2026-10-07）：** 本 Skill 仍属于 L–N 组进行中的审核；逐条数学与能力状态记录在完整交付包的 `audit_current/ln/claims.json`，案例定义记录在 `audit_current/ln/cases.json`；单独安装目录不包含这两份包级账本。请按 [`references/portable-audit-access.md`](references/portable-audit-access.md) 通过 bundle-root resolver 访问。以下限定只陈述已绑定账本的证据边界，不代表整项验收完成。

> 便携文献定位见 [本 Skill 证据索引](references/handoff-evidence-index.csv)。若需读取包外 PDF，以数据包根目录为 `EVIDENCE_ROOT`，按索引中的相对路径定位。

> **本项边界：** 账本 M02-SEP-01 覆盖单位圆盘中的一个精确分离例；M02-SEP-FD-01 仅支持 R^n 中非空闭凸集与外点的严格分离版本。局部凸空间、非闭集合或其他分离定理版本必须重新核对拓扑和闭性假设。

## 输入与产出

局部凸空间、凸集及拓扑、待分离点、闭性/紧性假设、候选连续线性泛函。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 检查集合凸且待分离点在集合之外，并选定满足假设的分离定理版本。
2. 构造连续线性泛函，计算其在整个集合上的上确界，明确要求严格分离还是弱分离。
3. 有限维紧凸集可用支撑函数 h_C(u)=sup_{x in C} u·x 比较候选点。
4. 将严格分离间隙作为显式证书，并对整个集合验证泛函，而非只抽查样本。
5. 若目标在闭包内或集合非凸，说明限制；只有在问题允许时才作凸化。

## 证据与失败处理

支撑超平面可能与集合接触而不能严格分离外点；无限维空间中泛函连续性和闭性假设不可省略。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

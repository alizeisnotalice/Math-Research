---
name: math-m01-convex-hull-extreme-points
description: "用于Carathéodory、Krein–Milman 极点的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# M01 · Carathéodory、Krein–Milman 极点

> **当前核验范围（2026-10-07）：** 本 Skill 仍属于 L–N 组进行中的审核；逐条数学与能力状态记录在完整交付包的 `audit_current/ln/claims.json`，案例定义记录在 `audit_current/ln/cases.json`；单独安装目录不包含这两份包级账本。请按 [`references/portable-audit-access.md`](references/portable-audit-access.md) 通过 bundle-root resolver 访问。以下限定只陈述已绑定账本的证据边界，不代表整项验收完成。

> 便携文献定位见 [本 Skill 证据索引](references/handoff-evidence-index.csv)。若需读取包外 PDF，以数据包根目录为 `EVIDENCE_ROOT`，按索引中的相对路径定位。

> **本项边界：** 账本 M01-BND-01a/b 已在有限维 R^d 证明 Carathéodory 的 d+1 上界及标准单纯形重心的锐性。Bachir 的结果只在论文明确的 Banach 函数、紧可度量 Φ-凸、连续性及分离条件下支持 exposed-point 闭凸包表示；一般 Krein–Milman 来源仍待核。有限维证明不推出一般 KM。

## 输入与产出

凸集及其环境拓扑、紧性/闭性、有限维或局部凸条件、极点定义、所需表示形式。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 核对凸性和具体拓扑；区分有限维 Carathéodory 条件与局部凸紧集的 Krein--Milman 条件。
2. 在 R^d 中对 conv(S) 内的点，用仿射相关性将有限凸组合化简到至多 d+1 个点。
   【修订 2026-10-06（外部审核）】三个不同陈述须区分：(a) 『至多 d+1』＝Carathéodory 定理（标准结果引用；LP 佐证：等式约束秩 d+1 ⟹ 基本解支撑 ≤d+1，未独立证明）；(b) 『锐性』＝存在需要 d+1 的实例（单纯形顶点本身即需 d+1 表示的极端情形；数值观察 150 例 max=d+1 与之一致）；(c) 『所有实例都须 d+1』＝不成立也不被声称（单纯形内点 1 个足够）。登记簿 ID：M01-BND-01a（pending）/ M01-BND-01b（instance_check_only）。
   【修订 2026-10-07（本轮独立证明）】上述旧状态由本轮断言卡覆盖：仿射相关消元证明 M01-BND-01a；标准 d-单纯形重心的唯一重心坐标证明 M01-BND-01b。150 次随机样本仅作演示，不再作为锐性的证据。当前假设、证明和边界见 `references/audit-current-lmn-20261007.md` 及 `references/method.md` 中的锐性证明。
3. 对 Hausdorff 局部凸空间中的**非空紧凸集**，才应用 Krein--Milman 得到它是极点的闭凸包；这不是有限凸组合结论。
4. 检查目标究竟需要有限凸组合、闭包还是积分表示。
5. 若文献只是提到 Krein 空间或 Krein 算子，须确认存在真正的凸极点论证。

## 证据与失败处理

无限维 Krein--Milman 一般只给闭凸包，不保证有限凸组合。Krein 空间是带不定内积的空间，不能凭同名认定为 Krein--Milman 结果。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

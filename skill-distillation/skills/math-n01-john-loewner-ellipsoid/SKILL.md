---
name: math-n01-john-loewner-ellipsoid
description: "用于John、Löwner 椭球与容量的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# N01 · John、Löwner 椭球与容量

> **当前核验范围（2026-10-07）：** 本 Skill 仍属于 L–N 组进行中的审核；逐条数学与能力状态记录在完整交付包的 `audit_current/ln/claims.json`，案例定义记录在 `audit_current/ln/cases.json`；单独安装目录不包含这两份包级账本。请按 [`references/portable-audit-access.md`](references/portable-audit-access.md) 通过 bundle-root resolver 访问。以下限定只陈述已绑定账本的证据边界，不代表整项验收完成。

> 便携文献定位见 [本 Skill 证据索引](references/handoff-evidence-index.csv)。若需读取包外 PDF，以数据包根目录为 `EVIDENCE_ROOT`，按索引中的相对路径定位。

> **本项边界：** 账本 N01-AFFINE-01 覆盖满维凸体在可逆仿射变换下的最大内接椭球对应关系；N01-JOHN-02 是非对称体不能误用 √n 因子的三角形反例。John 接触点矩条件及一般/中心对称因子定理的完整来源、双向条件和变形证明仍待核验；不得把反例当成 John 定理证明。

## 输入与产出

R^n 中满维凸体 K、对称性、待求内接/外接椭球、体积目标、仿射归一化。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 区分最大体积内接的 John 椭球与最小体积外接的 Löwner 椭球；只有核对原点和对偶体条件后才用极性联系两者。
2. 将 John 椭球平移到中心，再用可逆仿射映射化为 B^n。
3. 非对称情形用 K⊂nB^n；原点中心对称情形用 K⊂sqrt(n)B^n，并保留精确因子与唯一性条件。
4. 接触点证明需核验加权一阶/二阶矩恒等式；变形证明则分别核验包含关系和行列式增长。
5. 用单纯形或立方体验证尖锐性，并区分椭球问题与 John 域、John--Nirenberg 不等式及 Loewner 微分方程。

## 证据与失败处理

对称性和内接/外接约定会改变常数。John 域或 Loewner 演化的同名文献并未证明椭球结论。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

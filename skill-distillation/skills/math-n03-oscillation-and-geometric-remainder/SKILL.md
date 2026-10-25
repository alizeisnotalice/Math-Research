---
name: math-n03-oscillation-and-geometric-remainder
description: "用于振荡抵消与几何余项（van der Corput、半密度）的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# N03 · 振荡抵消与几何余项（van der Corput、半密度）

> **当前核验范围（2026-10-07）：** 本 Skill 仍属于 L–N 组进行中的审核；逐条数学与能力状态记录在完整交付包的 `audit_current/ln/claims.json`，案例定义记录在 `audit_current/ln/cases.json`；单独安装目录不包含这两份包级账本。请按 [`references/portable-audit-access.md`](references/portable-audit-access.md) 通过 bundle-root resolver 访问。以下限定只陈述已绑定账本的证据边界，不代表整项验收完成。

> 便携文献定位见 [本 Skill 证据索引](references/handoff-evidence-index.csv)。若需读取包外 PDF，以数据包根目录为 `EVIDENCE_ROOT`，按索引中的相对路径定位。

> **本项边界：** 账本 N03-JAC-01 与 N03-VDC-02 只支持列明的行列式恒等式和简单振荡例。Wu–Xi 的部分逐页记录只覆盖所列补充证明/公式，不解决其范围矛盾；Guo 的适用域、printed 主项、Thm 2.2 归一化及 Thm 2.1 引用链仍 pending。指数对优化/Lelechenko 程序主张也尚未全文和代码核验。

## 输入与产出

振荡相位、振幅正则性、导数/曲率下界、区间/区域几何、尺度参数、目标余项。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 先明确I(λ)=∫e^(iλφ)ψ、支撑、λ>0、导数阶数及全支撑下界；经典一维vdC k≥2可用非零k阶导数，k=1还需φ′单调等具体源条件。
2. 在驻点或退化Hessian区域分片，将非退化区估计、局部体积/切向消振与端点项分别列出，再按参数平衡；不能在零导数点套全局下界。
3. 每次换元写真实Jacobian与支撑。Guo Th2.6限d≥2，DγᵀADγ=2A给因子2^(d/2)；完整printed主项与由其定义的低阶R一并隔离。Th2.2的√(2π)及c_k归一化也禁止按printed精确主项调用；正确Th2.1及其明示p域独立保留。
4. 算术指数和使用Wu–Xi的squarefree/squarefull、模数、指数三元组与parent范围条件；addendum补证明不等于勘误，B1区间标签矛盾仍隔离，欧氏积分与算术和不混用。
5. 指数对搜索分别记录P与convP、分式分母正性、严格LC、可行incumbent与全部active-node下界；Lelechenko源Step1无可行性过滤的程序隔离。
6. 半密度/Kaufman/GP-BRANCH私有名需原始定义、变量映射及连接证明。缺任一项时仅记gap，不将消振机制改称该私有接口公开定理。

## 证据与失败处理

Guo Th2.6完整printed渐近及低阶R、Th2.2精确主项和Lelechenko约束搜索源程序已隔离；余项机制/修补不因反例而自动认证，局部修复仅按两审具体scope使用。相位导数、Jacobian、算术parent条件或严格可行性缺失均不允许调用更强结论。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

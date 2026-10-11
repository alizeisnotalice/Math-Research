---
name: math-n02-log-concave-marginal-inequalities
description: "用于Prékopa–Leindler、BBL 与对数凹边际的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# N02 · Prékopa–Leindler、BBL 与对数凹边际

> **当前核验范围（2026-10-07）：** 本 Skill 仍属于 L–N 组进行中的审核；逐条数学与能力状态记录在完整交付包的 `audit_current/ln/claims.json`，案例定义记录在 `audit_current/ln/cases.json`；单独安装目录不包含这两份包级账本。请按 [`references/portable-audit-access.md`](references/portable-audit-access.md) 通过 bundle-root resolver 访问。以下限定只陈述已绑定账本的证据边界，不代表整项验收完成。

> 便携文献定位见 [本 Skill 证据索引](references/handoff-evidence-index.csv)。若需读取包外 PDF，以数据包根目录为 `EVIDENCE_ROOT`，按索引中的相对路径定位。

> **本项边界：** 账本 N02-SLICE-PL-01 仅在逐点联合对数凹、固定切片可积并满足 PL 点态前提时推出边际不等式；Tonelli/Fubini 的 a.e. 结论不自动给每个切片或处处代表元。连续 PL/BBL 一般定理、Satomi 参数变换、逆 Minkowski 和零集代表元仍有明确来源缺口。

## 输入与产出

非负可测函数、维数与积分变量、插值参数、对数凹或 p-凹指数、可积性和支撑条件。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 连续Prékopa–Leindler：先固定0<λ<1和非负Lebesgue可测可积f,g,h，0<∫f,∫g<∞；核验h((1−λ)x+λy)≥f(x)^(1−λ)g(y)^λ在正支撑乘积上a.e.成立，再得∫h≥(∫f)^(1−λ)(∫g)^λ。
2. 连续BBL采用α≥−1/d约定，输出积分幂为α/(1+dα)；α=−1/d取−∞，α=∞取1/d。Satomi比值均值的p=−α，Q_d(p)=p/(1−dp)，两套记号不可混用。
3. 消去变量前固定两个外变量y0,y1，把联合对数凹条件应用到内部变量的两个切片，再用PL。需切片可积且边际有限；只由a.e.切片结论时，说明代表元与零集，不能暗称每点有限。
4. 多输出比值只在每个输入函数正支撑上定义essential infimum；slice递推保留Fubini零集与Q_e∘Q_d=Q_(d+e)适用范围。逆Minkowski只用均值幂≤1，低幂简化限p≤1/(d+1)。
5. 离散PL/BBL另核膨胀格点/输出函数约定；等号、稳定性和奇异Brascamp–Lieb各自保留源条件与尚未全读的外引。

## 证据与失败处理

正负指数约定与continuous/discrete模型不能混用。Satomi的多输出证明依赖外引二函数minimum不等式，当前保留外引缺口；只有两轮审核且绑定当前hash的scope可标checked。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

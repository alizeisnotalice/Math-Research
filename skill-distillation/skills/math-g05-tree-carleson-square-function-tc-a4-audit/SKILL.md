---
name: math-g05-tree-carleson-square-function-tc-a4-audit
description: "用于TC-A4 Carleson 型平方函数界的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前状态（2026-10-07）**：标准二进 Carleson 嵌入构件部分可用；TC-A4 仍未定义，专属平方函数界无可检对象。条件化有限过滤树的 4C 推导另列且不冒充 TC-A4。
> 外部文献证据随完整证据包提供；包根须同时含 `workspace_revised/`、`audit_current/` 与 `evidence/`。`EVIDENCE_ROOT` 和 `MATH64_SOURCE_PACKAGE_ROOT` 均指包根。单独安装的 Skill 不含共享 PDF/claims/cases；按[本地访问说明](references/portable-audit-access.md)从 registry root `math64-20261007` 解析包根以定位绝对 resolver，或显式提供包根。用 `--skill-dir` 与来源 ID 核验来源 SHA、用 `--artifact-path` 定位 claims/cases 并核验当前 SHA。registry、root 或 resolver 缺失时报告配置缺口，不猜路径/内容。resolver 的 PASS 只验证路径和当前字节，不代表全文阅读、命题证明或数学验收。


# G05 · TC-A4 Carleson 型平方函数界


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
8. 复用 P-6a781aff4000e334 Theorem 2.1 的随机常数 e 时，滤过满足 usual conditions（右连续且完备），整数 n≥1；驱动 X⁰,…,Xⁿ 为连续平方可积实鞅，彼此正交且二次变差相等。对每个 k=0,…,n，a^k 渐进可测且 E∫₀∞|a_s^k|²d⟨X^k⟩_s<∞；初值 u⁰₀,…,uⁿ₀∈L²，u^k 按 (2.2)–(2.3) 的广义 Cauchy–Riemann 积分式定义，因此有 L² 终值。α≥0 渐进可测且对所有 t≥0，M_t=E[∫_t^∞α_s ds|F_t]≤1 a.s.；总 α 质量闭鞅 N_t=E[∫₀∞α_s ds|F_t]须有表示 N_t=N₀+Σ_{k=0}ⁿ∫₀ᵗm_s^k dX_s^k，其中 N₀∈L¹、m^k 渐进可测且对每个有限 t 有 E∫₀ᵗ|m_s^k|²d⟨X^k⟩_s<∞。此定理仅给上述随机系统的嵌入≤e·EΣ_k|u^k_∞|²，不是一般确定性 Carleson 树常数 e，也未提供 TC-A4。

## 证据与失败处理

目录未定义 TC-A4。标准二进 Carleson 嵌入不能冒充对私有系数规范或平方函数的精确结论。

复用 P-e8f8d143f231a4b3 前，必须先读 [有效来源合同与函数空间纠正](references/source-contract-P-e8f8.md)，再按 active 来源卡核条件；不得用历史 little Bloch 误名替代原文 B˚={f∈B:f(0)=0}。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

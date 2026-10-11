---
name: math-g01-carleson-measure-box-and-tent-test
description: "用于Carleson 测度的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前状态（2026-10-07）**：部分可用。有限树逐节点 packing 例与来源条件化流程已核；Garnett Theorem 1.2 的全部先决条件已在步骤 5 列出。G01 的两个有限树案例由本地脚本实际运行，见 `audit_current/g/cases.json`；来源几何构造未独立重证，中心立方体迁移未建立。Holmes–Psaromiligkos–Volberg 的 Theorem 1.11 印出 bi-tree box/Carleson 分离，但显示的 `μ_A` 分量在顶盒 Q 的计数与 p.15 约 N 的说法不合；`μ_A` 只是总测度 `μ` 的一部分，其他分量及其对两侧的影响未说明。因此只记 source-proof gap，不把它当作已验证反例，也不判定定理真假，详见 `audit_current/independent-g/P-e38596-top-box-counting-proof-gap-20261007.json`。
> 外部文献证据随完整证据包提供；包根须同时含 `workspace_revised/`、`audit_current/` 与 `evidence/`。`EVIDENCE_ROOT` 和 `MATH64_SOURCE_PACKAGE_ROOT` 均指包根。单独安装的 Skill 不含共享 PDF/claims/cases；按[本地访问说明](references/portable-audit-access.md)从 registry root `math64-20261007` 解析包根以定位绝对 resolver，或显式提供包根。用 `--skill-dir` 与来源 ID 核验来源 SHA、用 `--artifact-path` 定位 claims/cases 并核验当前 SHA。registry、root 或 resolver 缺失时报告配置缺口，不猜路径/内容。resolver 的 PASS 只验证路径和当前字节，不代表全文阅读、命题证明或数学验收。


# G01 · Carleson 测度


## 输入与产出

底空间、上半空间/帐篷几何、正 Borel 测度 μ、盒/帐篷族和目标嵌入。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 路由适用性

本入口保留 partial：有限树路由只在显式树、底测度 σ、全部节点系数 α 与子分割给定后检逐节点 packing；当前叶例只产出 packing 常数 1，未提供 f、嵌入输出及匹配来源定理的全套假设，不能作为嵌入正例。UR/PDE 路由另需匹配的域 Ω、调和函数 u、实际局部能量/边界几何输入及指定方向的 UR 或 corkscrew+CDC 条件；缺这些输入时停止该路由，不从叶 packing 改称 PDE 定理适用。

## 执行步骤

1. 先定义实际T(Q)/S(Q)、底测度、尺度及边界，验证正Borel测度局部有限；写sup_Q μ(T(Q))/σ(Q)，不能只验根盒。
2. 有限过滤树的系数packing检每个J的Σ_{I⊆J}α_I≤Cσ(J)；有互不交子分割与统一mass流时可用已检查4C嵌入构件，树外多参数盒条件须另证。
3. 解析/调和空间的嵌入与反向kernel测试按指定空间定理，先明确核范数归一化；普通盒packing不自动证明所有reproducing-kernel thesis。
4. 邻近 P-e9e9211495f0b20e 限 domain Ω⊂R^{d+1}, d≥1。性质 (a) 是：对每个有界调和 u、x∈∂Ω、0<r<diam(Ω)，r^{-d}∫_{B(x,r)∩Ω}|∇u(y)|² dist(y,∂Ω)dy≤C||u||∞²。不是任意度量 Carleson 测试。
5. 同源性质 (b) 是：对每个有界调和 u 与 0<ε<1，存在 g∈W^{1,1}_{loc}(Ω)，||u−g||∞<ε，且存在 C=C(ε,Ω)，对所有 x∈∂Ω、r>0，r^{-d}∫_{B(x,r)∩Ω}|∇g(y)|dy≤C。Theorem 1.1 A 向只需存在 Ω̃⊂Ω、∂Ω⊂∂Ω̃ 且 ∂Ω̃ UR，即推出 (a),(b)；B 向还需 Ω 满足 corkscrew (1.4)+CDC (1.7)，并有 (a) 或 (b) 任一项，才推出这样的 UR 超域。Theorem 1.2 两向均要求 Ω 满足 (1.4),(1.7)，ε₀依赖两几何常数。A 向：若 (a) 或 (b) 成立，则对每个 0<ε<ε₀ 存在 C(ε)，对任意 x∈∂Ω、R>0、p_j∈Ω∩B(x,R) 及两两不交 E_j⊂∂Ω，只要 ω(p_j,E_j,Ω)≥1−ε，就有 Σ_j dist(p_j,∂Ω)^d≤C(ε)R^d。B 向：若对某个固定 0<ε<ε₀，存在有限常数 C 使上述 (1.10),(1.11)⇒(1.12) 对所有可行情形成立，则 (a),(b) 成立。不得简写成纯几何 packing 或一般盒测试。
6. 作者Whitney近/远边界拆分→调和事件packing→corona/UR构造的常数依赖d、几何参数及ε；外引容量/UR理论未独立复证，不能称维数无关。
7. 中心立方体迁移须实际树化或帐篷化、共同输入及 packing 桥梁；盒、帐篷与所有并集条件须分别定义并直接证明所需桥梁，不自动互换。Holmes–Psaromiligkos–Volberg Theorem 1.11 虽印出 product-bi-tree 分离，但显示的 `μ_A` 分量在顶盒 Q 有 `≳N log N` 项数，与 p.15 的约 N 说法不合；由于 `μ_A` 只是总测度 `μ` 的一部分，其他分量及其对两侧的影响未说明，这只构成来源证明缺口，不是已验证反例，也不判定定理真假（见 `audit_current/independent-g/P-e38596-top-box-counting-proof-gap-20261007.json`）。区分有界、消失、加权与算子特定 Carleson，不从邻近边界 PDE 定理得全尺度最大界。

## 证据与失败处理

根packing、抽象单树packing、解析kernel测试与边界调和函数Carleson估计有不同前提。邻近来源含UR/corkscrew/CDC与ε依赖，不能把缺Ahlfors假设读成无几何假设。

复用 P-e8f8d143f231a4b3 前，必须先读 [有效来源合同与函数空间纠正](references/source-contract-P-e8f8.md)，再按 active 来源卡核条件；不得用历史 little Bloch 误名替代原文 B˚={f∈B:f(0)=0}。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [有限树 packing 示例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

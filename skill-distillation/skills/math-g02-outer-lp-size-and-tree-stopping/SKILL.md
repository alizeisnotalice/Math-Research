---
name: math-g02-outer-lp-size-and-tree-stopping
description: "用于外测度 Lp、大小函数与树嵌入（Do–Thiele）的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前状态（2026-10-07）**：部分可用。外分布定义、有限 retained-set 反例和局部 Cauchy 目标已核；旧 R13 的 plain `supavg≤L2` 通过已撤回，历史脚本实际检查的是更弱的 `supavg≤L2·√(N/10)`。Do–Thiele 及 Fraccaroli 只在各自写明的 size、几何和参数范围内引用。一般中心立方体/树接口仍未证明。
> 外部文献证据随完整证据包提供；包根须同时含 `workspace_revised/`、`audit_current/` 与 `evidence/`。`EVIDENCE_ROOT` 和 `MATH64_SOURCE_PACKAGE_ROOT` 均指包根。单独安装的 Skill 不含共享 PDF/claims/cases；按[本地访问说明](references/portable-audit-access.md)从 registry root `math64-20261007` 解析包根以定位绝对 resolver，或显式提供包根。用 `--skill-dir` 与来源 ID 核验来源 SHA、用 `--artifact-path` 定位 claims/cases 并核验当前 SHA。registry、root 或 resolver 缺失时报告配置缺口，不猜路径/内容。resolver 的 PASS 只验证路径和当前字节，不代表全文阅读、命题证明或数学验收。


# G02 · 外测度 Lp、大小函数与树嵌入（Do–Thiele）


## 输入与产出

外测度空间、生成集、size 泛函 S、函数 F、指数 p，以及可选树结构。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 由生成集E、费用σ定义覆盖外测度μ，并完整给size的单调、齐次和准次可加公理。
2. 外分布μ(S(F)>λ)=inf{μ(B):outsup_{X\B}S(F)≤λ}，不是普通点态超水平集；outer Lp层析正是该分布的定义，无需假称外测度可加。
3. 树/tent停时同时证阈值捕获、异常集费用和packing；覆盖费用相加仅给上界，不等于最优覆盖外测度。
4. P-c0cdd27d39fd5e63 Theorem4.1在其R×(0,∞) tents与指定核/size下强型1<p≤∞；square size另需均值零等条件。
5. 其Theorem5.1的相位空间strong p>2、p=2仅weak，参数0<|α|≤1、|β|≤0.9、b≤2^{-8}及Fourier支持须匹配；常数允许依赖α,β,b,φ,p。
6. 外对偶、插值和配对依据具体size定理；不套普通Lp对偶，也不从参数依赖界推参数一致性。
7. 新增全文P-eb102a6795f6ce4e的Theorem1.3限σ-finite (X,μ,ω)、μ-null⇒ω-null及标准ellr size：p>a用sup局部La,q刻画，p≤a须sup-inf且B⊂A满足μ(A\B)≤μ(A)/2；参数常数C(a,p,q,r)不宣称端点一致。
8. 不能把被删部分条件换成μ(B)≥μ(A)/2：已检查有限μ(nonempty)=1、ω计数、f=1、r=2、q∞时，真weaknorm和正确K′均sqrt(m)，错误K″=1；比值无界。
9. 同源Theorem2.1对一般size须p有限且p>a并保分割常数K；p∞推广另需Remark4.3的size–mass不等式。新源原定理为sourceclaim，有限反例是checkedderivation，中心立方体接口仅推断。
10. 新全文P-aaa03e027f21d631的finite outerellr统一对偶限1<p≤∞、1≤r<∞或p=r=1/∞；p1,r>1有树反例非一致。指定R^d半空间box费用、ω=t^−1的对偶范围才完整1≤p,r≤∞，tent等价另保几何与维数常数。
11. 该v1 Thm6.1 weak-(1,q) 印刷式漏尺度因子 t^{d−d/q}；§6.3 控制的是尺度加权对象，缩放检验否定未加权字面式（q=1 不由该缩放反例否定）。Appendix Prop.8.6 p.31 印有“μ-null⇒ν-null 或逐集合局部积分–size 控制”两种替代假设。固定标准模型 X=N、μ=ω=counting、S=ℓ∞_ω、ν({n})=n 满足零集分支；对每个 n 取 f_n=δ_n，则 d_{f_n}(λ)=1 对 0≤λ<1、其后为0，outer L¹(S)=1 而 ∫|f_n|dν=n。故在同一个 σ-finite 模型和同一个 ν 下，n→∞ 给出比值无界，零集分支单独不足以推出 (8.1)，印刷的 OR 命题因该分支而有反例。这个序列论证是解析反例；有限案例只复算若干 n，不是一般证明。第二个局部积分–size 分支需按原文完整保留，本审计不否定该分支本身，也不猜写替代定理。
12. Do–Lewers 2020 Theorems 1–2（P-e448d2077c4eb584）仅支持其连续上半三维 wave-packet P(f)、定制 tent σ_w/S_w（含 sup 与 lacunary 加权 L² square size）、Fourier 紧支 Schwartz φ、δ=2^{-8b}、2<q<∞、w∈A_{q/2}及常数依赖范围；证明用三套 dyadic grids、well-separated partial families 与 good-λ。它不是抽象外测度或任意树定理，中心立方体接口另证。Fraccaroli 2021《Duality for double iterated outer Lp spaces》（P-6920b727cff4af3e，PDF 44页、含附录 A/B）研究三层迭代 outer Lp 的有限模型对偶、collapse 与特定上半三维 dyadic 几何；它仅作双迭代结构旁证，不替代本组单层 size/stopping 定理。P-aaa03e027f21d631 是 2020 年《Duality for outer L^p_μ(ℓ^r) spaces and relation to tent spaces》，正文第10步引用保留该 ID。

## 证据与失败处理

非加性不妨碍按定义做外分布层析，错误在将覆盖上界当等号或替换分布定义。相位空间p=2只弱型，几何/频率参数进入常数。

复用 P-e8f8d143f231a4b3 前，必须先读 [有效来源合同与函数空间纠正](references/source-contract-P-e8f8.md)，再按 active 来源卡核条件；不得用历史 little Bloch 误名替代原文 B˚={f∈B:f(0)=0}。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

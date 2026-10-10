# 2026-10-07 当前验证记录

本日志对应 `math-g01-carleson-measure-box-and-tent-test` 的当前工作副本。当前 SKILL SHA-256：`9d69373d3dd8c13ad36e9458bbde5a9381e656b6c01bf0cb8742c03358afe8b4`。

断言级记录：`audit_current/g/claims.json` 中 6 条（含条件、来源定位、证据强度、scope、未决缺口）。案例：`audit_current/g/cases.json` 中本项 2 条，均通过有限/结构门检查；见下表。

| case_id | kind | result | target claim | 本地验证范围 |
|---|---|---|---|
| `G01-positive-leaf-mass-packing` | positive | pass | G01.C01 | All-node packing constant equals 1. |
| `G01-negative-root-only-test` | negative | pass | G01.C01 | Root test is 1/2 while LL ratio is 2; all-node gate returns 2. |

执行方式：在当前工作副本上实际调用 `audit_current/g/cases/run_g_cases.py`；结果由精确有限计算与断言产生。它不是自动调用 LLM 执行 Skill，也不替代来源定理证明。

来源页码/全文读卡、PDF SHA、文本 SHA、逐页覆盖/附录状态与页内行范围见 `audit_current/g/papers-read.json`；附件身份与首页相关性筛查见 `candidate-screening.json`。

独立审核复核：Garnett 2020 PDF SHA `e9e9211495f0b20e966daccecfef4833a62cd3c233b88b7293b24b1eff385732`，31/31 页；Theorem 1.2 p.4 已逐字核对且视觉复看。A/B 两向均要 corkscrew (1.4)+CDC (1.7)，ε₀依赖这些常数；A 对每个 ε<ε₀ 得 C(ε)，B 只要求某一个固定 ε 的 packing 对所有合格配置成立，随后推出 (a),(b)，不得写成同一 ε 的等价。x、R、p_j 位置、E_j 互不交及调和测度条件已纳入 Skill step 5。来源几何证明仍属 partial check。

来源证明缺口复核：P-e38596bd68bb977e 的 exact-PDF/SHA 独审卡为 `audit_current/independent-g/P-e38596-top-box-counting-proof-gap-20261007.json`（SHA-256 `75650342f7a7d28bfd6d436f23547814cea3c9bb30560bf610a4d4cb72d70513`）。原文 (1.11) 对所有 `R0∈D` 测试，包含 Q；其显示的 `μ_A` 分量在 Q 给出至少 `cN log N` 个下/超曲线内钩形矩形，与 p.15 的约 N 说法不一致。但 `μ_A` 只是总测度 `μ` 的一部分，其他分量及其对两侧的影响没有说明，因此只记录 source-proof gap，不能认定定理为假或将其用作已验证反例。两条 G01 有限树案例与此来源构造无关，本轮未重跑；它们不验证 G01.C05。

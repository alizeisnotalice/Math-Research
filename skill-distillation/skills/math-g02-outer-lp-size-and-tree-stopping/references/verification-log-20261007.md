# 2026-10-07 当前验证记录

本日志对应 `math-g02-outer-lp-size-and-tree-stopping` 的当前工作副本。当前 SKILL SHA-256：`e2d1ea4d45b82fca783cb85b1effacf28b7a3800f9620257729ab62375403092`。

断言级记录：`audit_current/g/claims.json` 中 10 条（含条件、来源定位、证据强度、scope、未决缺口）。案例：`audit_current/g/cases.json` 中本项 5 条，均通过有限/结构门检查；其中2条属于主验收集，3条为补充案例；见下表。

| case_id | kind | acceptance_case | result | target claim | 本地验证范围 |
|---|---|---|---|
| `G02-positive-outer-distribution-p>a` | positive | yes | pass | G02.C01 | d_F(3)=1, d_F(4)=0 and outer L^{2,∞} norm is 4. |
| `G02-negative-retained-set-condition` | negative | yes | pass | G02.C03 | Correct condition yields 4 and wrong replacement yields 1. |
| `G02-positive-local-cauchy-target` | positive | no | pass | G02.C05 | All 6339 windows satisfy (average f)^2≤average(f^2). |
| `G02-negative-global-rms-bound` | negative | no | pass | G02.C06 | Window average 1000 exceeds global RMS since 1000²>667667/2. |
| `G02-negative-prop86-null-branch` | negative | no | pass | G02.C09 | In one fixed σ-finite counting model, outer L¹(δₙ)=1 while ν({n})=n, so no finite uniform C works. |

执行方式：在当前工作副本上实际调用 `audit_current/g/cases/run_g_cases.py`；结果由精确有限计算与断言产生。它不是自动调用 LLM 执行 Skill，也不替代来源定理证明。

来源页码/全文读卡、PDF SHA、文本 SHA、逐页覆盖/附录状态与页内行范围见 `audit_current/g/papers-read.json`；附件身份与首页相关性筛查见 `candidate-screening.json`。

R13 旧 plain `supavg≤L2` 通过已撤回。归档脚本真实断言含 `√(N/10)`；当前独立核验的是每个有限窗口的局部 Cauchy–Schwarz 与反例。

Prop. 8.6 p.31 的零集分支由固定标准 ℓ∞-size 计数模型反例核验。G02.C10 单列全 n 解析论证：在同一 ν({k})=k 下 outer L¹(δ_n)=1、ν 积分=n，故比值随 n→∞ 发散；有限套件只复算 n=2,10,1000，不是一般证明。该反例只否定第一分支单独充分性，未否定第二个局部积分–size 分支，也不构造替代定理。

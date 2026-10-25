# 2026-10-07 当前验证记录

本日志对应 `math-g03-directional-capture-tree-quadratic-audit` 的当前工作副本。当前 SKILL SHA-256：`1623de947e5e60280029a48decd8b1fc3491057601bf5d203d26429c222971ea`。

断言级记录：`audit_current/g/claims.json` 中 3 条（含条件、来源定位、证据强度、scope、未决缺口）。案例：`audit_current/g/cases.json` 中本项 2 条，均通过有限/结构门检查；见下表。

| case_id | kind | result | target claim | 本地验证范围 |
|---|---|---|---|
| `G03-positive-explicit-PSD-quadratic` | positive | pass | G03.C01 | x*Ax=16 and diagonal-only sum=4. |
| `G03-negative-private-interface-missing` | negative | pass | G03.C03 | Gate returns evidencepending and makes no project estimate. |

执行方式：在当前工作副本上实际调用 `audit_current/g/cases/run_g_cases.py`；结果由精确有限计算与断言产生。它不是自动调用 LLM 执行 Skill，也不替代来源定理证明。

来源页码/全文读卡、PDF SHA、文本 SHA、逐页覆盖/附录状态与页内行范围见 `audit_current/g/papers-read.json`；附件身份与首页相关性筛查见 `candidate-screening.json`。

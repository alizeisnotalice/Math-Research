# 2026-10-07 当前验证记录

本日志对应 `math-g04-rubio-de-francia-square-function-scope-check` 的当前工作副本。当前 SKILL SHA-256：`76b6ab79f69cbe6789f8b370c2d8753941bb10ecedf5b9036066bbf7e085360e`。

断言级记录：`audit_current/g/claims.json` 中 3 条（含条件、来源定位、证据强度、scope、未决缺口）。案例：`audit_current/g/cases.json` 中本项 2 条，均通过有限/结构门检查；见下表。

| case_id | kind | result | target claim | 本地验证范围 |
|---|---|---|---|
| `G04-positive-disjoint-projection-Parseval` | positive | pass | G04.C02 | Total energy equals sum of block energies, 30. |
| `G04-negative-overlapping-interval-counted-twice` | negative | pass | G04.C03 | Energy doubles 5→10; square norm factor is √2. |

执行方式：在当前工作副本上实际调用 `audit_current/g/cases/run_g_cases.py`；结果由精确有限计算与断言产生。它不是自动调用 LLM 执行 Skill，也不替代来源定理证明。

来源页码/全文读卡、PDF SHA、文本 SHA、逐页覆盖/附录状态与页内行范围见 `audit_current/g/papers-read.json`；附件身份与首页相关性筛查见 `candidate-screening.json`。

# 2026-10-07 当前验证记录

本日志对应 `math-g05-tree-carleson-square-function-tc-a4-audit` 的当前工作副本。当前 SKILL SHA-256：`3ac59ea662027598521020bd96921eb20c9b39205a293cff2622c7b98b2d7a9a`。

断言级记录：`audit_current/g/claims.json` 中 4 条（含条件、来源定位、证据强度、scope、未决缺口）。案例：`audit_current/g/cases.json` 中本项 2 条，均通过有限/结构门检查；见下表。

| case_id | kind | result | target claim | 本地验证范围 |
|---|---|---|---|
| `G05-positive-depth1-Carleson-embedding` | positive | pass | G05.C02 | Embedding sum 1≤4C‖f‖₂²=4. |
| `G05-negative-root-only-packing-gate` | negative | pass | G05.C02 | Root ratio 1 but left-child ratio 2; all-subtree constant is 2. |

执行方式：在当前工作副本上实际调用 `audit_current/g/cases/run_g_cases.py`；结果由精确有限计算与断言产生。它不是自动调用 LLM 执行 Skill，也不替代来源定理证明。

来源页码/全文读卡、PDF SHA、文本 SHA、逐页覆盖/附录状态与页内行范围见 `audit_current/g/papers-read.json`；附件身份与首页相关性筛查见 `candidate-screening.json`。

TC-A4 仍未定义；本日志不将标准二进模型结论转称为 TC-A4。

新增来源条件核对（G05.C04）：Domelevo et al. Theorem 2.1 固定 n≥1，假设滤过 usual conditions；对每个 k，渐进可测 a^k 的全时域积分平方可积且 u_0^k∈L²；α 的条件尾质量≤1；总质量鞅 N_0∈L¹，并有每个有限时间上平方可积的渐进可测 stable-subspace 系数 m^k。驱动的连续性、平方可积性、正交与等二次变差条件及定义 (2.2)–(2.4) 见 PDF p.3。该相邻随机结论不定义 TC-A4。

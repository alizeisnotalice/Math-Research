# 独立复核记录（append-only，2026-10-06）

本文件由外部审核批次写入，不修改本 skill 任何既有内容。
分级：A 完整证明 / B 附加条件下成立 / C 仅数值实例检查 / R 文献确认 / G 证据不足。

## 复核分级：C（案例级，见下）

- **K02-CASE-01**（instance_check_only，own_derivation+numeric_experiment）：P(G>2)=Ψ(2)=0.0227501、联合 p²、并 2p−p²=0.0449827 全部精确复核（erfc 闭式）。案例正确指出 G2=G1 时交叉项不可忽略。

## 更正条目 I（统计口径）

本 log 为 R13 末补建（此前批量写回遗漏了本 skill 的验证记录）。
验证状态与依据类型以 `audit_r12/assertion_register.json` 双维登记簿为准。

---

## 档位提升（R26）

依据既有验证记录：Ψ(2)=0.0227501 与 2p−p²=0.0449827 erfc 闭式精确（instance_check_only）
待证据 → **部分可用**。源定理证明本身未重证（保持守卫）。

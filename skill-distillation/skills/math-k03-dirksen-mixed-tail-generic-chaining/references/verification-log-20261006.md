# 独立复核记录（append-only，2026-10-06）

本文件由外部审核批次写入，不修改本 skill 任何既有内容。
分级：A 完整证明 / B 附加条件下成立 / C 仅数值实例检查 / R 文献确认 / G 证据不足。

## 复核分级：C（案例级）+ 引用忠实性已核（R14）

- **K03-CASE-01**（instance_check_only，numeric_experiment）：混合尾界
  P(|X|≥t) ≤ 2exp(−min(t²/K₁,t/K₂)) 的实例检查通过（X=G₁+G₂²−1, t=4）。
- **K03-REF-01**（source_comparison，R14）：步骤 3 引用 Dirksen Theorem 3.5
  (E sup_t‖X_t−X_t0‖^p)^{1/p} ≤ C(γ₂+γ₁)+2sup_t(E‖X_t−X_t0‖^p)^{1/p} 及尾形式
  exp(−u)——与源卡 HK-S-9ecc39a3d40ae7fe Theorem 3.5 **逐字一致**；
  步骤 5 的 Hu–Simchi-Levi Theorem 3.1（多体制逐点同时界）与
  HK-P-b19fb9f102d33f2e 一致（含有限 m 条联合可测半度量、可分性、锚点条件）。

## 档位提升依据（R14）

待证据 → **部分可用**（引用忠实性已核 + 案例实例通过）。
仍待证据：Dirksen/Hu–Simchi-Levi 定理证明本身；generic chaining 的 γ₁/γ₂ 泛函计算实例。

# 独立复核记录（append-only，2026-10-06）

## 档位提升（R26）

依据既有验证记录：B1 案例通过（LP gap 恒等式一行解析证明，A 级）
待证据 → **部分可用**。源定理证明本身未重证（保持守卫）。
验证详情见 audit_r12/assertion_register.json 对应子断言。


---

## 更正条目（R32）：三条承重定理证明复现——依据类型升级

总表：`audit_r12/classical_proofs_verified_r32.md`。

- **Carathéodory**（凸版 d+1 / 锥版 d / 有限生成锥闭性）：极小支撑 + 仿射相关扰动论证完整重证；
  锥版闭性 = 有限个 simplicial cone 之并——**本 skill 步骤 2 所依赖的定理本体已重证**，
  M01-BND-01a 的依据由"标准结果引用（LP 佐证）"升级为 source_proof_verified
  （LP 基本解佐证保留为独立第二路径）。
- **Farkas 引理**：不可行 ⟹ 分离证书的方向完整（锥闭性由 Carathéodory 锥版给出，
  分离定理的使用条件齐备）——本 skill 的证书逻辑两端闭环（对偶侧 gap 恒等式 R26 已证）。
- **不交谱投影 Parseval 分解**：P_I 正交投影代数（P_IP_J=0）+ 分辨率恒等 ⟹ p=2 常数
  恰 1 且**不可改进**（单频取等）——本 skill 的 p=2 常数机理完整。

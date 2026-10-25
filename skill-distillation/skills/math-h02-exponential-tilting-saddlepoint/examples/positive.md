# H02 适用示例：Gaussian 大 `r` 的 Mills 渐近

实际运行见 `audit_current/hi/run_h02_i01_cases.py`。取 `n=10^6` 个标准正态变量之和，平均阈值 `w=n^(−1/4)`，于是 `u=r=√n w=31.6228`。精确尾为 `8.97916×10⁻²²⁰`，Gaussian 首项 `φ(u)/u=8.98813×10⁻²²⁰`，首项/精确尾约 `1.000998`。该有限实例与 Mills 界及 `u→∞` 的相对等价方向相符；单个数值例不证明渐近结论。

同一脚本另算 Rademacher 格点三角阵，记录精确尾与简化项的有限样本偏差，避免把 Gaussian 结果误读为任意分布的精确近似。

Package-relative `audit_current/...` locators in this file are resolved with the package-root resolver described in [portable audit access](../references/portable-audit-access.md); they are not Skill-local file links.

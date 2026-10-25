# I01 适用示例：条件中心化 Rademacher 增量

脚本 `audit_current/hi/run_h02_i01_cases.py` 实际枚举/核算 4 个独立 Rademacher 增量，取 `c=1`、每步 `v_k=1`、总上界 `v=4`、阈值 `t=2`。每个整数 `p≥2` 的条件绝对矩为 1，满足 `1≤p!/2`；条件均值为零。脚本得 `P(ΣX_k≥2)=5/16=0.3125`，Bernstein 界为 `exp(−1/3)=0.71653`。另在 `θ=1/3` 核对乘积 MGF `cosh(1/3)^4=1.24387≤exp(1/3)=1.39561`；每步包络为 `exp(1/12)`，四步相乘才是 `exp(1/3)`。这是实例检查，普遍结论由 `references/method.md` 的推导支持。

Package-relative `audit_current/...` locators in this file are resolved with the package-root resolver described in [portable audit access](../references/portable-audit-access.md); they are not Skill-local file links.

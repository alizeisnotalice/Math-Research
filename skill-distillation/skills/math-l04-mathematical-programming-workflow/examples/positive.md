# 适用案例

最大化 x+y，约束 x+y<=1、x,y>=0。原问题点 (1,0) 的目标值为 1；对约束 x+y<=1 取对偶乘子 1，也给出上界 1。

本例现已额外通过真实 `scipy.optimize.linprog(method="highs")` 调用：求解器返回 optimal 状态、原问题点 `(1,0)` 和不等式边际 `−1`。将边际换成最大化对偶乘子 `λ=1` 后，用 `Fraction` 精确核对原问题可行、`λ≥0`、`λ≥1` 及原/对偶值相等。版本、输入模型、容差、原始 API 输出和逐项证书在 `audit_current/ln/l04_solver_probe.json`；可用 `audit_current/ln/run_l04_solver_probe.py` 与 `requirements-l04.txt` 重跑。

这一条只验收两变量 LP 的 API 调用、状态读取和返回证书检查，不代表混合整数、非线性或任意 LP 工作流已验收。求解器使用浮点计算；精确 Fraction 检查只核对本次返回的十进制见证，不给出一般误差界。

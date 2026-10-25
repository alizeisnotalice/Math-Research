# 当前反例：依赖输入必须拒绝 Hoeffding 调用

目标 claim: `H05-HOEFFDING-INDEPENDENT-SUM-01`。本例对应当前 case `H05-Hoeffding-dependent-counterexample`，用于测试独立性门控。

取 `X_1=…=X_20=Z`，其中 `Z~Bernoulli(1/2)`。每个变量单独仍在 `[0,1]`，但 20 个变量完全依赖，因此 claim 的相互独立假设不成立。令 `S=Σ_i(X_i−1/2)`、`t=10`，则

`P(S≥10)=1/2`，

若错误套用独立 Hoeffding 公式，会得到

`exp(−2·10²/20)=exp(−10)=0.000045399929762484854`，

明显小于真实尾概率。结论是必须拒绝对该输入调用独立 Hoeffding 定理；这不是该定理的反例，也不表示独立情形的结论为假。

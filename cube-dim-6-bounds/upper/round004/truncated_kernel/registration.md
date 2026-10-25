# 第004轮：截断 Hermite 核的有限探路

计算冻结 21-fan 与完整 signed-permutation 群的
\(H_\epsilon=(2\pi)^{-6}\sum\alpha\operatorname{DD}_5[\arctan(\cdot/\epsilon)]\)。
输入为三个有理代表（round-003 最小点附近、两层 \(S_2\) 点、已认证的三层
round-002 Green anchor）和 \(\epsilon\in\{1/100,1/200,1/500\}\)，共 9 次。

每个 cone 的节点均按有理交叉乘法排序及判等；相等节点采用完整 confluent
Hermite Newton 表，使用 \(F_\epsilon^{(k)}(b)/k!\)，不删项。C++ 的
`long double` 在本机实际是 53-bit、8-byte binary64，只作探路；至少一个
generic 点由 Arb 高精度独立校准。即使该校准通过，本批也不认证 \(g\)、\(I\)、
截断误差或 \(C_6\)。

# 第003轮：直接谱核有限诊断

目标是以冻结的 21 个排序 simplex 和完整 signed-permutation 群，直接诊断
generic face 点的 \(G\) 值量级。不是条件 Monte Carlo，不生成或认证连续
最小值与积分。

初始分母 257 的四个固定分层点均在完整群作用后出现精确零或重复节点；这些
点已登记为 `STRUCTURAL_FAILURE`，其删项和不作数值结论。修补子批使用分母
1000000007、固定 seed 6003004 的 32 个 Latin-hypercube face 点。C++ 在每个
cone 先以 `__int128` 比较 \(b_i=n_i/(Dq_i)\)，只有无零、互异节点才计算
DD5(sgn)；本子批全部通过该门控。

一个 generic 点另由 python-flint `fmpq` 对全部 967680 项直接作准确有理和，
再以 384-bit Arb 计算 \(G\)。其与 C++ long-double 的相对差为
\(3.534950879305865\cdot10^{-12}<10^{-6}\)，故允许执行注册的 32 点诊断。
包含 zero/repeated 节点的 v1 点另以四分量 Laurent C++ 诊断重跑；其负幂在
浮点输出中为零，但未取得逐点精确 cancellation 证书，因而不作为本批成本样本。

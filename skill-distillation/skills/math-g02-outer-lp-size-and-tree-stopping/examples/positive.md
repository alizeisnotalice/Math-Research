# 外分布定义的有限正例

令 X 有 16 个点，ω 为计数测度，μ(∅)=0 且 μ(A)=1 对所有非空 A 成立，F≡1、r=2。此时 size(F)=4，故 d_F(3)=1、d_F(4)=0，按外分布定义得到 outer L^{2,∞} 范数 4。这里检查的是定义本身的有限模型，不声称满足 Fraccaroli 或 Do–Thiele 定理的全部假设。

实际结果见 `audit_current/g/case_results_20261007.json`（G02-positive-outer-distribution-p>a）。另用 1,…,1000 的向量穷举所有指定长度连续窗口，逐个核验局部 Cauchy–Schwarz：`average(W)^2 ≤ average(|f|^2 over W)`；这才是样本平均所对应的正确独立目标。

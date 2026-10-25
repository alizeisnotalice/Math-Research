# 条件缺失案例

高维DNN矩阵就称CP：要求CP分解或适用的低维等价依据。

本例已由本轮独立执行；实际输入、条件、结果及当前 Skill 哈希见 [本地证据解析说明](../references/evidence-guide.md)；bundle-relative artifact：`audit_current/cd/cases.json`（由说明中的resolver命令解析）。该实例检查不等于一般定理证明。


## HIER-02：固定零阶不精确

取 `A=diag(1,−1)`、`K=R₊²` 与同一乘积 Exp(1) 测度。零阶局部化矩阵是标量
`M₀(f_Ay)=∫(x₁²−x₂²)dμ=2!−2!=0`，所以 `A∈C₀`；但取 `x=e₂` 得 `xᵀAx=−1`，因此 `A∉COP`。固定有限阶可行不能作为共正证书；本例只检验这个端点，不否定所有阶的交等式。

证据计算记录：`audit_current/cd/c02-hier02-instance-checks-20261010.json`。

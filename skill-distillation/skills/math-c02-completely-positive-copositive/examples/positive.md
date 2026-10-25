# 适用案例

v=(1,2)≥0，vvᵀ=[[1,2],[2,4]]是CP；逐项检验定义。

本例已由本轮独立执行；实际输入、条件、结果及当前 Skill 哈希见 [本地证据解析说明](../references/evidence-guide.md)；bundle-relative artifact：`audit_current/cd/cases.json`（由说明中的resolver命令解析）。该实例检查不等于一般定理证明。


## HIER-02：指数正交锥的条件迁移

取 `A=I₂`、`K=R₊²`、`dμ=e^(−x₁−x₂)dx₁dx₂`，则 `f_A=x₁²+x₂²≥0`。`K` 闭、`f_A` 连续、`supp μ=K`；任意多项式加权 `x^αf_A` 可积；对两个坐标的每个固定整数移位 `s≥0`，矩为 `(2(k+s))!`，Carleman 级数发散。对任意多项式 `g`，`∫g²f_A dμ≥0`，故每阶局部化矩阵均半正定且 `I₂∈COP`。这是实例条件与接口检查，不是一般层级定理的独立证明。

证据计算记录：`audit_current/cd/c02-hier02-instance-checks-20261010.json`。

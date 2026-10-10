# 条件缺失案例

二维取 K₀=K₁=I₂。分层范数平方和为 2+2=4，但 ||K₀+K₁||²_HS=||2I₂||²_HS=8；忽略交叉项会低估一倍。

已执行检查：`audit_current/f/cases.json` 中的 `F03-NEG`；逐项输出见 `audit_current/f/example_execution_20261007.txt`。这些有限计算只核验该例，不证明一般定理。

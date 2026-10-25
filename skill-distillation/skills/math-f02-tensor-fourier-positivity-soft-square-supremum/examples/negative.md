# 条件缺失案例

在 `(Z/2)^2` 取 `K(x,y)=1-2(-1)^(x+y)`，其 Walsh Fourier 系数在频率 `(1,1)` 为 `-8`，所以 K 不正定。任意交叉项都须单独检验。

已执行检查：`audit_current/f/cases.json` 中的 `F02-NEG`；逐项输出见 `audit_current/f/example_execution_20261007.txt`。这些有限计算只核验该例，不证明一般定理。

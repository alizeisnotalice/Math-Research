# 条件缺失案例

若在 K=0 处直接取 log K，就会遇到未定义值，可能错误排除真实的边界候选。

已执行检查：`audit_current/f/cases.json` 中的 `F01-NEG`；逐项输出见 `audit_current/f/example_execution_20261007.txt`。这些有限计算只核验该例，不证明一般定理。

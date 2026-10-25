# 条件缺失案例

本例测试最小化剪枝所需界的方向与缺项。目标为最小化；当前记录只有 `node_upper=3/10` 和 `incumbent_lower=1/4`，它们是最大化方向的界，不能充当最小化所需的认证节点下界 `node_lower` 与可行 incumbent 上界 `incumbent_upper`。因此应报告缺少这两个最小化界，不把最大化的 `U<L` 比较套用于该目标，也不据此剪枝。

历史附加例仍保留最大化算术：`U=3/10`、`L=1/4` 时 `U<L` 为假，所以不剪枝。它只说明最大化方向的有限比较，不测试或反驳最小化剪枝命题。

本例仅检查有限输入与目标方向是否匹配，不证明一般剪枝定理。实际输入、当前条件、结果及 Skill 哈希见 [本地证据解析说明](../references/evidence-guide.md)；bundle-relative artifact：`audit_current/cd/cases.json`（由说明中的 resolver 命令解析）。

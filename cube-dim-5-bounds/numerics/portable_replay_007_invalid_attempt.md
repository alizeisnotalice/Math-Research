# 可移植重放首次无效启动

首次启动针对主 clone 已清理过的 `body_fan_arb_anchor.py` 仍假定旧的
`body_fan_check` import，导致临时复制阶段的断言失败；未生成节点、未编译、
未运行数值计算。replay 脚本随后改为同时接受已经使用 `body_fan_geometry` 的
干净 helper，后续的实际重放以新 receipt 为准。

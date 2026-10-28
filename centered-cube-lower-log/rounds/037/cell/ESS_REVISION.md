# ESS 后处理修订记录

本修订是同一冻结协议的重放纠错，不是独立确认。原 `probe.py`、`registration.json`、`RESULTS.json` 均保留且 SHA256 与 `ESS_revision.json` 完全一致。

首次主实验耗时 26.1449329853 秒。它在 n=63、R=8 的极端尾前缀产生 ESS 的 0/0 警告：原实现先计算线性权重，再将其平方，部分平方在浮点中下溢为零。不能据此忽略这些前缀。首次未保存逐样本 logM/logD，因此无法只对原文件后处理。经主线同意，使用同版本解释器和 NumPy、原九个 seed、原 batch 顺序重放一次；重放耗时 25.8362388611 秒。主实验数据生成实际执行了两遍，总耗时 51.9811718464 秒，两遍不是两个独立数据集。小例/校准和首次宏观诊断另见原记录；宏观诊断未在纠错重放中重新生成。

修订只替换 `weak_scan` 中累计尾权重、累计权重平方、ESS 与系数的运算方法。用 `logaddexp.accumulate(logW)` 和 `logaddexp.accumulate(2 logW)` 计算它们；未过滤极端尾候选与 ESS≥64 的扫描均使用 log 域。最小变更见 `ESS_FIX.patch`。`replay_stable_ess.py` 导入原脚本并替换这个函数，还附加保存日志与计算 hash；来源、query grid、提议分布、RNG 调用和选择规则均未改。

重放逐 case 核对以下原有字段完全相等：metadata、权重质量均值/SE/max、uniform 分支数、各模型/标签次数、最大实际系数贡献及其帽、PB 质量残差、正 subnormal 个数和下溢量级预算。八个未受影响 case 的旧/新扫描 best 系数相对差最大为 1.8431e−14。此一致性支持同协议重放；首次没有保存 query hash，因此不能声称有两次完整逐坐标 hash 的比较。新重放的 query 流 hash、logM/logD 流 hash 和 NPZ 文件 hash 均保存在修订结果中。

纠错后的 n63/R8 joint 最佳值仍为 1.0014778442（ESS≈1750），极端尾没有成为新的最佳候选。只读日志复核发现原公式在该 case 的 joint 有 13 个非有限 ESS 前缀、prior 有 57 个，各 fixed 表为 1–12 个；其余八个 case 为零。新扫描所有模型、所有 case 从保存日志重算后与 `RESULTS_stable_ESS.json` 完全一致。详见 `saved_log_audit.json`。该复核不生成样本。

运行环境为 `/opt/homebrew/opt/python@3.15/bin/python3.15`，Python 3.15.0、NumPy 2.5.4。纠错脚本核验版本一致并禁止第二次纠错重放。原源码、注册、原结果、修补脚本和最小 patch 的 SHA256 已写入 `ESS_revision.json`。此后交付以 `RESULTS_stable_ESS.json` 为数值主表，原 `RESULTS.json` 仅保留为不可静默替换的历史记录。

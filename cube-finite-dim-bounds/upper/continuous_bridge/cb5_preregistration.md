# cb5 第5轮增量预登记

2026-10-25；原生 upper_continuous_bridge 分支，不重启第5轮。唯一写入当前目录；主 clone 第5轮文件只读。模型固定为 Lebesgue 背景、中心轴平行方体、共同半径、一般 L1 与全尺度。

## 输入和问题

输入绝对路径及内容哈希写入 cb5_sources.json。已读主 clone rounds/005-expansion.md、results/round005_expansion_snapshot.json、rounds/005-registration.md、STATE.md、results/bounds_ledger.json；这些状态只作上下文，当前已审 L3 以 auditclone audit/007b_results.json 完整准确端点为准；005c 未审。上界只研究有限差分核 LP 与连续分布 Delta 的桥接，不比较下界。

问题：对紧支撑张量 Q1（多线性帐篷）插值，准确计算分布 Delta 的面源，检验逐点有限差分 Delta 非负是否足够；另外比较二维标准对角三角剖分的 P1 插值，如能直接给精确公式则纳入，不能则停止在 Q1。至多这两种机制。

## 有限批和停止条件

1. 理论计算 Q1 面源，预期逐轴二阶差分，而不是各轴和。
2. 用 Python Fraction 对一个显式有理离散可行三维格点函数作完整有限差分检查，并提取离原点有正距离的负面源。首选离散 Green 函数的最小网格 [-2,2]^3，若其不显式可解只再试 [-3,3]^3；不扩大。
3. 若原点邻面失败可直接证明一般紧支撑非负 Q1 无法只有原点负 Delta；同时给有限元边界/外延公式。
4. 提出可验证接口：单元体源、内面跳、外边界跳、原点通量、障碍和质量，说明误差如何修补/为何尚不自动修补。
5. 交一个命题/反例、接口、最小准确 checker 和实际运行摘要，交 independent_audit。没有连续候选即不声称改善 U3/U4。

Python 实际运行：/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python。MATLAB 只经主代理串行，若需要提交前缀 cb5 脚本，不自行调用。不会重跑主 LP 或扩网格；准确 dual 只用标准库核一次用于区分范围。不自行提交推送/扩招。

## 证据状态

预期可严格证明 Q1 的局部公式和明确反例；数值计算与严格证明分开。任何自由 LP 数值解不升级为连续可行核；准确有限对偶 >8 不升级为连续障碍。主模型矩阵和 Fraction 算术共享须披露。新接口仅在所有列出条件被核后可称连续证书。

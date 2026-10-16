# 第241轮：共同软化参数的极值窗口几何

日期2026-10-17。一般n^(1/2+o(1))一致弱型上界尚未闭合。

新估计保留全部坐标共用的soft参数，以显式一维凹优化Ψ=max_s∏[(1−s)(1−g_i)+s(1−ε)]定义窗口亏损Gsharp=−logΨ，ε=(1/2)exp(−3)。它在旧距离E=0处仍可能正，并保留稀疏坐标的对数增益。真实adaptive接收列付原B菜单，联合FIRST后全部节点、端点及外系数核清。

单一来源门max(E,Gsharp)≥2log+ρD包含旧距离门，使用一次来源预算，费用为84Λe⁵[1+log(e ceilπ√n)]³||f||₁。原文量词核验允许从同一原分裂统一取K=1；不删除低pe等补项。新门只在几何参数层面可严格扩大，未声称实际新增FULL正交通非空或覆盖比例。两个亏损均低于阈值的近来源联合预算仍未支付。

## 验证与精简归档

[research.zip](research.zip) 保存13个必要成员，包含详细LaTeX、独立推导、K来源审计、独审、两组数值程序/结果、说明、报告和重放驱动。不编译PDF。

- ZIP 46,612 bytes；SHA256 `b184d16687b8c9896035a10e2aa05c88affbab5b3d869cc52a80698dc6670817`。
- 本机实际40个核/连续L1几何配置、25个共同soft极值配置、192项Fraction阈值检查通过。
- 包外完整驱动复跑通过；两份完整JSON逐字节一致，12项manifest在前后核验通过。普通高精计算不是区间证书，有限模型不认证一般阶数或FULL历史。
- 标量30/50位复核；产品40/70位仅指基于32位scalar输入的后续log算术，精度边界明确保留。

解压后执行 `python3 reproduce.py --output-dir <包外不存在的新目录>`，本机Python/mpmath路径见numeric_notes.md。完整checkpoint：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round241_20261017_final_checkpoint`；回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round241_20261017_final_replay/receipt.json`。

本轮仅新增README和ZIP两个tracked文件，另更新项目总览；本轮work经checkpoint/备份ZIP/归档逐成员一致核验后移除。下一步集中处理同时低E、低Gsharp下的真实接收与高来源质量联合几何。

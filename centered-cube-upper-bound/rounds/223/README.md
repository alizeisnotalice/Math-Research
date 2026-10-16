# 第223轮：有符号后验交通与负信用核验

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合，本轮没有新的根阶付款。恢复原 good−bad score 方向，证明同一原来源测度下完整端点、删除跳、外部来源正交通、suffix 恢复及正门损失恒等式。原 N_del 是连续反向交通，不能替代删除原子；N_Z 已消费，不能重复扣除。R222 正项放大不能自动移交给负项。

删除跳在该坏来源后验势中有利。真正未解决的是带真实接收测度的边界与 suffix 恢复费用的联合几何预算，逐行势≤1尚不足以给根阶。未将恒等式或有限模型称为一般闭合。

## 实际验证

- 完整 LaTeX、来源账、联合账及独立审核保留，不编译 PDF。
- 本机 Python Fraction：258组配置、90对阈值连续性检查，维数至10000；有理多项式求导及59/420连续证书通过。独审另重构全部原始记录。
- 本机 MATLAB R2026a Symbolic Math Toolbox：四个符号恒等式化简为0，删除跳12/95、连续预算59/285。MCP入口失败后实际采用本机batch；失败路径在报告中保留。
- 最终checkpoint迁移到包外重跑：Python约0.456秒，MATLAB约17.507秒，两份完整JSON逐字节一致，11项manifest前后通过。有限模型不认证原FULL几何或一般预算。

## 精简归档

[research.zip](research.zip)：12成员、76346字节，SHA256 `706197a368e9c03d9820d7018b7765134fc1dfdfc19cc5aeeaca216d54727981`。仓库新增2文件，更新项目README。项目根目录821项，Git对象包366.94 MiB。

展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round223_20261016_final_checkpoint`；独立ZIP：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round223_20261016_research.zip`。工作目录、展开备份、ZIP逐成员相同后才移除本轮工作目录。

复现：`python3 reproduce.py --verify-only`；完整运行 `python3 reproduce.py --output-dir <新的包外目录>`，需要Python及MATLAB Symbolic Math Toolbox。`--python-only`仅部分复现。最终回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round223_20261016_driver_replay/receipt.json`。

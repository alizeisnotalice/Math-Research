# 第273轮：终端高来源恢复与全近域根阶中心转移误差

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。

返回仍保原 terminal-high 的正母式，恢复真实软来源权重，以精确三因子加载 ell 同时处理共同硬窗口和 halo。全部 near 中心转移的加性误差按 O(sqrt(n)) W_nu 支付，不需局部来源质量帽或大维度吸收。主项 Q_ell 保 high、near、终端载体与连续增长折扣，仍未取得一般根阶预算；未声称 halo 主交通或一般端点闭合。

[研究与复现包](research.zip)：39,079字节，10成员；SHA256 `29fdb60b878526d434b4d948804815db70059864fcc316672c04b6dfcbfc44a0`。

- proof.tex.txt、report.md、screen.md、audit.md：详细LaTeX、来源账、独立先筛与终审PASS；不编译PDF。
- numeric_checks.py/json、loading_checks.py/json、reproduce.py、manifest.json：最小复现材料。

本机16项精确接口例子、8个维度菜单的50/80位标量价格（166字段最大归一化差约3.45e-49），以及2/4/8源加载检查通过。加载检查含12项因子等式、14项上界、3项展开恒等式。完整reproduce.py实际重跑，两JSON逐字节一致；同一36参数族的两组审计有理断言通过。有限接口算术不认证FULL或一般阶数。

解压后用 Python 3（依赖mpmath）运行 `python3 reproduce.py`；`--check-only`核哈希与独审算术。展开备份与ZIP全部成员逐字节核对。本轮仅新增本索引和ZIP、更新总览。

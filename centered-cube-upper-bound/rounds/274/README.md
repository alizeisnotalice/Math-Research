# 第274轮：原坐标抵消与全当前正交通的双因子改进

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合；本轮没有新增统一根阶主项付款。

保原终端good–bad两槽的坐标差，在同一完整后验下同时改进连续质量因子和目标损失，得到全当前正交通的较小上界 Q_cross<=Q_ell。没有裁出新小交通份额。另保时钟幂族与零幂公式，核准实际near硬见证对共同衰减因子的限制。原R273已付中心转移误差沿用；Q_cross及其它原分支的全局来源预算仍未证明。

[研究与复现包](research.zip)：66,362字节，11成员；SHA256 `a0e107d56177275c864999ebcee71562a800ce5483324df8aaa276d61f4af7a9`。

- proof.tex.txt、report.md、screen.md、audit.md、numeric_audit.md：详细LaTeX、费用、独立筛查与终审；主稿PASS，不编译PDF。
- numeric_checks.py/json、cut_checks.py/json、reproduce.py、manifest.json：最小复现材料。

实际本机三阶段9族时钟+6族三群cut检查；40/70位最大差分别约7.30e-36和4.59e-41。完整reproduce.py两JSON逐字节一致；12组坐标律/45点、8组零幂及180项端点/乘积核验通过。有限模型未认证真实P核、FIRST/high、历史或FULL；大相对clock gap同时保H²价格，不当成原交通反例。

解压后用 Python 3 和 mpmath 运行 `python3 reproduce.py`；`--check-only`核清单与精确/存储数字身份。本轮仅新增本索引和ZIP、更新总览；展开备份与ZIP逐字节核对。

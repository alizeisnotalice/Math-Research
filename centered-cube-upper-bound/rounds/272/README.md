# 第272轮：全尺度目标包络与共同硬窗口的中心转移

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。

新目标核包络覆盖两物理尺度及全部soft，真实dx列价≤24√n exp(6||y−z||1/b)a^(−n)，无需尺度菜单。精确连续逆degree增长将时钟常数4改为1并保数据折扣。在共同硬窗口，把纯损失中心移到真实硬来源，全部加性误差以维度无关的一次来源费支付，不需高补源局部质量帽；原硬带帽进一步消去主项中的外层补源。

同一ν的二来源主项Q和硬窗口外halo仍未付，原所有来源/选择器/signed账保持。全adaptive核列不等于一般来源预算。主稿先独筛再终审PASS。

[研究与复现包](research.zip)：61,604字节，10成员；SHA256 `2c66b197c1ab7ef3133ad6586c65c071bbf31711e116101eb4ace6f47c350526`。

- `proof.tex.txt`、`report.md`、`screen.md`、`audit.md`：详细LaTeX、费用、独审与缺口；不编译PDF。
- `numeric_checks.py/json`、`anchor_checks.py/json`、`reproduce.py`、`manifest.json`：必要复现材料。

本机58个arcsine系数菜单、27组双Abel、252项幅度比、81项profile检查及18个真实C中心转移菜单通过；50/80位最大差分别约2.394e-44、5.318e-52。完整reproduce.py重跑两个JSON逐字节一致，另5+44项独审有理算术通过。有限工具验证不认证FULL或一般阶数。

解压后 `python3 reproduce.py`，依赖mpmath；`--check-only`核哈希与独审精确算术。仓库外展开备份已逐字节核验，本轮仅新增本索引和ZIP并更新总览。

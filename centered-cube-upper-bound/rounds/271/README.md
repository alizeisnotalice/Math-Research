# 第271轮：共同存活来源的双时钟压缩与近位移几何

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合，本轮没有将新空间积分当作已付预算。

原near hard来源在终端载体共同存活，使后验平方保留一次真实硬来源并产生degree跨时钟衰减。Schur估计给常数4的时钟压缩，保删除空档和零度数处理。连续损失D_Z对已删除目标也满足纯权差界，进一步获得目标距离因子，near目标有1/√n损失。固定参数target-only列保原clamp，真实l1近位移只付exp(3/(8√n))，不产生指数维度费用。

原FIRST、两后验、真实dx、唯一j及共同叶费用完成回代；仍需证明实际空间联合积分Aint的根阶来源预算。固定列工具不能自动处理接收点依赖的选择器。详细主稿已独立先审及终审PASS。

[研究与复现包](research.zip)：34,996字节，9成员；SHA256 `4a7a0b12dabf5a99c1849b442d34ba7bb1426b24857e935b32fb7ef36869958a`。

- `proof.tex.txt`、`report.md`、`screen.md`、`audit.md`：详细论证、原资格、费用及缺口，不编译PDF。
- `numeric_checks.py/json`、`kernel_checks.py`、`reproduce.py`、`manifest.json`：必要复现材料。

本机12组区间并集三阶段验证，40/70位每精度36行独立求积，132字段最大差约3.71525e-40；另36项Fraction接口算术及6项真实目标核检查通过。主代理重跑主JSON逐字节一致并复核目标核。有限数值不是一般证明或FULL/区间证书。

解压后 `python3 reproduce.py`，依赖mpmath；`--check-only`只核验哈希。仓库外展开备份已逐字节核验，本轮仅新增本索引和ZIP并更新总览。

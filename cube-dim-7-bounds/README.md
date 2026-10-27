# 中心立方体极大算子：一般高维覆盖改进

本项目起于七维，目前按用户最新授权集中一般高维覆盖。固定Lebesgue背景、中心轴平行立方体、所有坐标共同半径、一般L1输入及全部正尺度。七维专项已停止，下界暂停；三四维来源及其他项目不改。

当前调度见[STATE](STATE.md)，严格的一般维数结果见[覆盖账本](results/general_cover_ledger.json)，历史固定维数的准确证书见[原账本](results/bounds_ledger.json)。所有内部结果均不称公开发表或全球最优。

第035轮已证明并获完整独立审核：对共同充分大的n，

\[
C_n\le n\log n-n\log\log n-
\frac{n}{131072}\sqrt{\frac{\log\log n}{\log n}}.
\]

完整机制是反对角多壳匹配和环境行/列度数的共享缺额支付，见[报告](upper/round035/REPORT.md)、[证明](upper/round035/energy_theory/PROOF.md)、[独立全范围审核](audit/round035/energy_theory/proof_review.md)。首项系数仍为1。[第036轮](upper/round036/REPORT.md)已独立认证显式起点loglogn≥128；该门槛不声称最小。

第034轮还有较低起点的准确参数族，见[报告及证书](upper/round034/REPORT.md)。新旧界只在各自适用范围比较，始终取适用界最小值；新渐近阶不意味着在所有旧门槛处均更小。严禁从展示小数重造证书。

## 最小复核与证据

解析证明承担一般输入、全尺度、全部充分大维数的量词；有限数值样本不能替代它。实际计算和源码hash随轮次保存。部分历史生产程序会覆盖自己的输出，复核应在可丢弃副本中进行。

- [R035准确比较](upper/round035/energy_theory/check_comparison.py)：标准库Fraction；比较完整有理参数。
- [R035核代数](upper/round035/typical_distribution/verify.py)：SymPy；区分单次cut障碍与新能量机制。
- [联合支付独立小例](audit/round035/joint_cover/verify.py)：标准库Fraction；验证module改善与必要边界罚项。
- [R034十二点独审](audit/round034/constant_family/check_matrix.py)：标准库Fraction；保留实际主机来源路径，不是重新搜索。
- [来源清单](upper/round035/MANIFEST.json)：锁定证明、代码、原实际运行输出及源身份。

本机生产计算由Terra负责，Sol推导/独审，主代理负责证明整合和共享MATLAB调度。不用空转或填充内存制造负载；只将真实运行作为计算证据。

## 历史与论文

各轮完整报告位于[rounds](rounds/)。旧下界构造及完整并集审计保留于[第006轮](lower/round006/proof.md)，本轮未启动下界计算。原七维上界点保存在准确账本，停止专项后只封存已完成的输出。

[现有论文PDF](paper/main.pdf)及[LaTeX](paper/main.tex)是已审核的13页旧稿，仍有效但尚未吸收R030之后的覆盖改进。没有声称重新编译或将新结果公开发表。旧轮次“暂停”记录属于历史；当前持续研究授权以STATE为准。

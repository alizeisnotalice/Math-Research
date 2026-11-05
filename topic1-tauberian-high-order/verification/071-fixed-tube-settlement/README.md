# R071 最小标量核

完整全输入证明在[固定整管推导](../../proofs/071-fixed-tube-settlement.md)，[独立数学重建](../../reviews/071-fixed-tube-settlement.md)逐项核集合包含、目标壳、同一预算及新N9。这里不以算核证明几何。

作者实际运行：`wolframscript -file bootstrap_scalar.wl`，Wolfram15 ARM PID79446。核零参数Hessian、移动最大点一阶梯度及峰值二阶余项；文件另保留早期保守bootstrap指数核，它不承担最终直接吸收证明。

独审实际运行：`python3 scalar_check.py > scalar_check.json`，仅标准库Fraction。核一般n峰值二次多项式、预算份额、q3…24指数样本及全m传递阶数；一般q恒等式在独审正文中解析证明。两个程序都是短小本机符号/有理核，不需要大规模重复计算。`manifest.json`绑定实际输出及冻结证明字节。

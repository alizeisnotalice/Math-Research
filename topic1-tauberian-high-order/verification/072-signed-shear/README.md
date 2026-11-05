# R072：正负剪切的真实比较与局部全图核

[完整证明](../../proofs/072-signed-shear-comparison.md)和[独立重构](../../reviews/072-signed-shear-comparison.md)给每个固定n≥3的两类真实有限输入九阶相对改进。构造先证明原半径、共同设计质量和全算子安全，再用行列式差；不是从det单独推halo，也不证明C_n的绝对九阶系数。

实际运行：

1. `wolframscript -file determinant.wl`：一般n的三阶det差残差0；固定同孔质量的n3归一值4只是代数效应。
2. `wolframscript -file neighbors_batch.wl`：依赖同目录neighbors_common；n2..4、两符号共228个非零小整数向量，e∈(0,1/[4(n−1)])整个区间精确Reduce。主PID80823，实际15子核、零mismatch。全整数及任意n分类由解析界证明。
3. `wolframscript -file full_charts_batch.wl`：依赖full_charts_common，三个真实global-min图×三个核项共9个独立精确积分，实际9子核，保留pair外最浅轴。再运行`full_charts_scaled.wl`、`full_charts_samples.wl`。
4. `python3 full_charts_fraction.py > 新输出.json`：标准库Fraction独立多项式积分，两组有理参数×9项全部等于已有WL输出；不依赖SymPy。输入是本包full_charts_samples.json。

全图精确式的7/4首项仅为一个相容方向的isolated single-pair线性模型，不能代替joint源层。完整证明通过weighted堆叠/jump条/clip账本保留a u⁸交叉退款，再由两轴边forest的局部反射证明正负同价，余量O(a u⁹)。三边triangle不能这样同时反射，它的密度先在av⁹区分。

最早plane-only核遗漏offmin及d₂≥s限制，不能作完整3D cross价格；scope更正原文保留，旧plane原始脚本/日志留在本地不作为本包成功证据。完整single-F积分已替代它，仅用于检查图域/方向；一般几何不由有限计算承担。原计算元数据中的WL内存是分配峰值，不冒称OS RSS或持续满CPU。所有文件SHA与冻结源绑定于manifest。

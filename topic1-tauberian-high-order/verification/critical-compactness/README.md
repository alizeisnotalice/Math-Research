# 第17轮：临界紧性的局部符号与横截面核验

`local_bridge.wl` 在本机 Mathematica/Wolfram Language 15.0.0 ARM 实跑四项恒等式：两向CDF面差、对数侧长Hessian、保留独立正预算的五层望远镜、标量盒体积对质量盈余的准确一阶导数，全部True、FAILURES=0。还核验一般n的各向异性周期胞体积差、归一化负差系数−κ²及稀有相混合贡献趋零。`local_bridge-result.txt`为实际输出。主代理在归档副本又实际运行同脚本一次，结论相同，无消息；单内核，无时间限。

`transverse_partition.py` 使用标准库Fraction，精确检查三维错相竖列铺砌的一张角片：总权重1，两中间胞各1/2，没有单一胞覆盖整张片。主代理实际从归档副本运行，输出为`transverse_partition.log`。它只核有理区间，不证明一般铺砌计数。

复现：从本目录运行 `wolframscript -file local_bridge.wl > local_bridge-result.txt` 与 `python3 transverse_partition.py > transverse_partition.log`。脚本没有私有绝对路径，标准Python无需额外库。

全部输入的正预算、自由面几何、有限输入反例及条件性LT计数的证明见`proofs/017-critical-compactness.md`及独审`reviews/017-compactness.md`。有限/符号计算不代替这些证明，不证明原输入能进入LT，更不证明真实Cn的精确六阶系数。独审报告中的Decimal有限样本只作数值检错，最终定理不用这些小样本作为严格证书。

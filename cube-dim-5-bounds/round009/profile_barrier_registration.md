# 主新增：全rho的rational Young形状障碍

实际运行前登记。将030的两个固定rho排除推进为完整连续参数族：alpha1,beta1/2；Z=a0[1−rho(W_x+2)/12]^(−3/2)，rho<1,a0≥0，允许空间自适应选择。只排除此明确族，不排除一般正lambda或真实Green。

新问题取x0=(1,1,1,0,0)，lambda有4个1−5rho/12和12个1−rho/12。利用030已推且主独立重推的U²≤8pi detD/C7及V=(128/35)xadjD x/C5，按lambda共同比例归一后令t=sqrt((12−rho)/(12−5rho))，范围0<t<9/7。用16个实际方向构D，SymPy准确核det/adj、Mellin与最终预算平方的有理恒等式；不枚举网格、不重跑两个旧rho。

一次有限CPU批：1个5x5符号矩阵、2个新Mellin指数、4个Fraction比较，目标证明点上族预算<1/10000，再用明确相对lambda/Loewner稳定性扩到正四维measure面patch，排除a.e.绕过单点的误读。eps=1/1000，patch factor=((1009/991)^9)(1000/999)^4，须准确<6/5，故patch预算<3/25000<当前I所需15198/1e8。所有常数来源/代码hash、PID/runtime保存。只涉及证明预算，不声称G或C5上界/下界；无MATLAB并行、无窗口，最多120秒检查点。

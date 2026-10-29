# 第38轮：条件矩小例与压缩接口

本阶段只做确定的小例及一般条件矩证明，不追加旧来源大样本。

来源沿用37轮的精确 quantile-cell：n=3,5，K4，R2,4；按 Bin(n,1/2) 的 lower-CDF quantile 分段，空段完整保留。表θ_b=(b+.5)/R，plus/minus cell offset 为θ_b/1−θ_b；全宽δ=1/(64KR)真实box，每源格质量K^{-n}，总质量1。

冻结长度10/9、7/6、11/9。每坐标用两个等概率有理query位置：围绕某表源峰进入或退出cube的中心转折点，偏移±δ/8，因此两个点落在宽δ/4的同一圆周细bin。坐标 i 使用表 i mod R、cell i mod K、边方向交替±，不依据实验结果改选。独立枚举全部2^n隐藏query。

每个query直接K^n源格求和，并以正PB段系数公式对照；全部用Fraction核对（含完整真实box重叠），另以NumPy浮点结果对照。再比较直接条件均值与替换各class核为两点均值后的PB多项式，计算真实方差与prod_i(1+max_{b,±}relativeVar(z_{i,b,±}))−1。

预计最大case仅32个query×1024源格，内存远小于2GB、总预算几分钟。执行一次确定计算，无随机seed、无筛选。先记录实际结果，再写Tonelli独立均值与中心化monomial的L2 Minkowski证明。Fraction证明覆盖这些有限输入，不是通用浮点区间认证。

随后给高维bin-count接口设计：common cell law、实际pχ及Bin(n,pχ)先验、核一二阶矩、条件Cantelli。高维PB近似必须等待独立Berry–Esseen明确误差；本阶段不以Gaussian proxy代替真实平均，也不立即启动大n实验。

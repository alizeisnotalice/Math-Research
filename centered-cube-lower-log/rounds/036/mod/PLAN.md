# 第36轮：模类别joint-cell来源的有限筛选计划

K固定12，R=2、3、4、6，n=3、7、15。每joint粗格质量K^−n，classifier为 \(b=\sum_i(j_i\bmod R)\bmod R\)。source格均匀，单坐标residue均匀，故实际先验精确为1/R。

标签只抽一次后固定：seed202610290360得到permutation [11,6,4,5,1,7,2,3,9,0,8,10]，用s_j=permutation_j mod R。第b张表 \(u_b(j)=((s_j+b)\bmod R+1/2)/R\)，点a=(j+u)/K，box宽δ=1/(64KR)。box留在格内，每格质量及总L1质量严格正确。

真正joint source抽样顺序：先uniform抽j_vec，再由该source格算b，再抽该点的box。prior对照则先独立uniform抽b，再抽j_vec。比较joint、全部R张固定表、实际先验prior；不得将两种抽样混淆。

真实共同cube L平均按residue循环卷积DP：每表b分别卷积n个坐标的R类别权重，取最终residue b项，再对b求和。fixed表取总类别权重乘积，prior按1/R平均。

先落盘最小实现并运行n3/K12/R2的全12^3格直接有限求和与DP；再扩R及小有理控制。之后才冻结N、seed与query grid执行MC。无Gaussian。

L预定包含短尺度δ/16至1/(2K)的13点几何段、至1.25的41点线性段及1/K、2/K、3/K、4/K、6/K、1、10/9、7/6、11/9。grid冻结后不随收益更改。L1是weak≥1确定锚。

proposal=.5 uniform+.5 uniform选全部来源与全部L，按真实source抽T，再X=T−U_cube mod1；实际RN由全部来源全部L的DP平均给出。权重≤2，避免稀有层仅均匀抽样。best weak报告全样本/ESS支持扫描，并以分半选择阈值后验证，不能只报低阈值CDF。

每R/n固定一个seed，case顺序n3/7/15、每n内R2/3/4/6，seed202610290361起。最坏n15/R6、batch128全L定标后按成本选择：优先N4096；若保守估计>180秒则N2048，仍>180秒就停止并报告。RAM2GB，预计实际远小于此限。禁止因阳性增加N、seed、R或改标签。

额外描述诊断固定64个查询中心，记录10/9、7/6、11/9下joint/prior−1及residue Fourier乘积残差；先数值再解释混合。一般严格估计由主线处理，本轮不重复宣称未证障碍。

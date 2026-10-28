# 第37轮：plus计数分位数joint-cell来源

先小例数值，再估计。只冻结一个quantile convention，不按收益改分类。

对S~Binomial(n,.5)，q_b取CDF首次达到b/R的整数：R*sum_(s≤q_b) C(n,s)≥b*2^n；末cut为n，初始cut为−1。第b段包含前cut+1至本cut。重复cut的空段保留；π_b用精确整数Binomial质量/2^n记录，绝不改成1/R或默删表。

K偶，χ(j)=1_{j<K/2}，source格uniform，每格质量K^-n。b由本source格的S确定。θ_b=(b+.5)/R，plus格u=θ_b，minus格u=1−θ_b；a=(j+u)/K。box宽δ=1/(64KR)，留在格内，质量1。

prior先独立按实际π抽表，与joint的source-cell关联不同。实现prior分类时用一组独立n个fair bits取S再classify，精确生成π；不是用近似Gaussian或把π=.1等浮点CDF当定义。

真实共同L：每表PB多项式取该表quantile段系数和，再跨表相加，无额外π。fixed表是对应坐标总权重乘积，prior才按π平均。先n3/K4/R2、4全64格+Fraction核对，再扩大。

大动态范围计算采用log(product单坐标总均值) × 正归一PB概率，段概率直接正系数求和，不做两个接近1的CDF相减。极小下溢的误差预算单列，不能冒充正态近似。

冻结候选n15/31/63，R2/4/8，K固定12，不看收益择K。L含δ/16短段、至1.25的共同线性段、cell尺度、L1及10/9,7/6,11/9。proposal=.5 uniform+.5全部来源和全部目标L的真实source抽样后卷积cube核；实际RN同模型同grid。

seed顺序202610290371起（每n内R2/4/8）。按最坏n63/R8 batch64全L定标，只依成本在N4096/2048/1024中选预计≤180秒者；若都超则停止。RAM2GB，预计远低。一次/case，不加seed追阳性，不做大n正态。

记录所有fixed表、实际prior、joint的best-weak支持/未过滤扫描，分半选阈值/heldout验证；不以两个下界排序声称weak改善，不要求昂贵fullweak上CI。

# 第38轮高维来源表与BE区间独立审查

2026-10-29。只读审查 `../round038_highdim/` 的 `source_tables.py`、`build_cache.py`、`probe_highdim.py`、`SOURCE_TABLES.md`、两版source验证及脚本、HIGH_DIM_PLAN/REGISTRATION/RESULTS。未重建1.46GB cache，未抽新MC，未修改冻结文件。moment组合、完整RN及统计族由另一独审负责，本报告主要核查来源核与BE参数转移。

**结论：未发现必须修正的来源定义、full-count端点或BE数学误差传播错误。数学包络在所声明的class均值误差前提下成立。当前FP证据不足以将实际保存结果称为严格有限真概率证书。**

## 1. 保存元数据与小核对

实际保存n=6704800、K255、R4、v8、w4、M3、Q48960、J233、N128；没有扩大样本。来源表构建12.39963秒；不是对各J另抽query。所有三个源码SHA256与PLAN匹配，RESULTS的REGISTRATION SHA256匹配，标准化reference对原版验证的hash也匹配。精确Fraction尺度集去重后233点，包含10/9、11/9；固定cuts为3338016、3339311、3340606，彼此不同且严格在计数支持内。

自有 `highdim_small_checks.py/json` 已实际运行，不调用build_cache、不用RNG：

- 288个真实GF来源的峰/bin边界位置，以Fraction trapezoid质量独立求和比较算法，最大观察count误差1.8119×10^-12，未出现差一个完整峰的端点错误。
- 54个Fraction class比例及Bernoulli方差扰动例通过；162个CDF参数矩形角点/内点例通过。
- BE的1/√n常数检查、η在[.49,.51]导出p在[.3,.7]内的有理检查通过。

这些是有限小核对与源码审查，不是完整FP认证。

## 2. 来源、class质量与周期分解

GF(2^8),poly0x11d的指数表、B4 Sidon/六defect/短lag及非邻接检查均在构造时执行；这只验证该有限来源，不补足v≥16logn或B≥v等原匹配概率参数。保存元数据明确两者不满足。

每表将原四标签在前127格作+b、其余作−b的mod4旋转。峰中心为(j+(label+.5)/4)/K，δ盒留在原格。共同背景为1−p+(p/K)D_M，每格都有相同峰质量p/K，所以各表全周期class质量相同：

`η=a+(p/(2K))Σ c_m sin(2πma)`，a=127/255。

η≈.498047800888639，代码未将其改为.5。确定系数界 `πm c_m≤20/(3m)`、`Dsup≤(20/3)H_M`、`Dprime≤(880/21)M` 给背景下界约.4760348584>0，方向正确。c_m的实际三角函数/平方根浮点表示与理想系数的误差尚未区间认证。

L=1+ℓ时，一个完整周期加中心x+.5、宽ℓ的余窗口，是原对称cube窗口的精确周期分解。`evaluate`的class质量只加一次；L=1单独返回全周期class质量和零range/variance，符合恒等式。

## 3. peak计数与端点

δ盒与宽h query-bin卷积的精确平均是两个中心均匀变量差的trapezoid CDF；差与和在此有相同分布。支撑半宽s=(δ+h)/2。完整峰区间 `[lo+s,hi−s]` 用左端包含、右端包含计数。

左partial候选取 `[lo−s,lo+s)`，排除已计入full的内端点；右partial取 `(hi−s,hi+s]`，同样排除full内端点。两个外端点的贡献恰为零，可保留。候选transition计数允许把内端点再记作“可能变化”，这仅使range帽更宽，不会把均值再加一次。ℓ>δ+h保证两端候选带分离；每端最多一个atom由峰间距界保证。源码的searchsorted side选择与这些开闭约定一致。

288个小核对刻意让atom落在左右窗口端点及±s支撑边界，而非只抽一般位置。观察到的极小浮点误差不能代替精确端点舍入证明，但足以排查漏/重计完整峰这类单位级错误。

## 4. 背景中点误差与range

class/window边界先分片，常数背景overlap分片线性，中点规则在实数算术下精确。D积分在每片内二阶导绝对值≤2Dprime；除以bin宽h后，中点总误差≤Dprime h²/12。乘p/(KL)即代码analytic_mean_error_abs，最大约7.71246×10^-12。边界恰落bin端点时不需要新增正长度片；重复分片计算只会覆盖重算同一bin，不会累加其均值。

range帽包含常数背景、两端transition峰及D-mask贡献。mask在class边界跳跃，但窗口积分连续且a.e.导数是两端masked density之差，因此D部分2Dsup的Lipschitz界合法。真实均值下端为nominal−e>0时，range²/[4(nominal−e)²]可作class相对方差帽。

上述解析中点误差与经验FP guard=10^-13分别保存。guard不是由区间舍入推导，不能把两者之和的数值标签当成完整计算证书。

## 5. 从class误差到PB三矩

关键常数 `p_sum_error=n e/minS` 正确，不缺额外2倍或分母膨胀。设nominal class均值为a,b且a,b>e，真实值分别在[a−e,a+e]、[b−e,b+e]内。比例a/(a+b)对a递增、对b递减，最坏角点为(a+e,b−e)和(a−e,b+e)，两者总和仍为a+b。因此

`|p_true−p_nom|≤e/(a+b)≤e/minS`。

函数p(1−p)在[0,1]上Lipschitz常数1，所以PB总均值和总方差各有误差≤n e/minS。此结论不依赖误差的独立性。总kernel乘积的log误差则由每坐标total±2e得到，统一使用 `n[-log(1−2e/minS)]`，同时包住上下方向。保存最大PB参数误差约5.25260×10^-5，最大logmean膨胀约1.05052×10^-4。

这里e必须真包含各class误差；该传播证明不会替代对e所含FP guard的认证，也不会自动包住随后BLAS点积的舍入。

## 6. BE、CDF角点、固定cuts及prior

正周期源且η∈[.49,.51]、L≥1时，严格宏观界给p∈[49/151,102/151]，包含于[.3,.7]。条件bin均值也满足同一界。源码的rational check

`(.5606)^2(.58)^2/.21<1`

确实保证一般非同分布BE的CDF绝对误差≤1/√n≈.00038619539。它不使用iid常数。

真实μ、variance位于nominal各±p_sum_error构成的矩形内。normal-CDF对μ单调；标准差方向取决于cut+.5−μ的符号。`cdf_enclosure`分别对两种标准差取min/max，恰是这个矩形上标准化参数的全角点包；不要求μ与variance误差能独立实现。之后外加±BE并截到[0,1]正确。

整数cuts固定后，即使其定义来自nominal μ±sd，它们仍是合法固定source分类，不能称精确分位数。四段为 `[0,c0]`、`[c0+1,c1]`、`[c1+1,c2]`、`[c2+1,n]`，无重叠/遗漏；两段使用c1和表0/2。首尾CDF用精确0/1，中间段下包为upperCDF的下端减lowerCDF的上端，上包反向相减，BE端点预算未遗漏。

prior用的是**实际类别质量η**的Bin(n,η)，含n×10^-13的参数不确定guard；未用Bin(n,.5)。各表η相同，四段（或两段）划分计数支持，故理想joint来源质量1；prior权重也按这些固定cuts计算。保存prior区间可以总宽度不归一为1，这不会破坏正加权上下包；不能把区间端点再强行归一而改变其意义。

positive表权乘以PB段上下包、随后求和保留真实条件均值的方向。没有把normal近似裸值替换成真概率后遗漏BE费用。首尾段仅一个有限CDF端点；内段付两端。若未来采用σ=0、小μ近边界等其他参数，这个driver的variance>error断言及当前p范围需重新核验，不能直接泛化。

## 7. 验证轨迹及剩余证据资格

主一维验证为 `SOURCE_VALIDATION_normalized_reference.json`：4表各52个预定bins、共64个class/bin有理peak核对，最大class均值与Gauss reference差约1.416×10^-12，Fraction peak折算class误差约4.19×10^-15。原版reference权重未归一导致常数积分的浮点消减偏差；原文件保留，新版对**同一批确定位置**重新算reference，source hash未变。两次都没有高维抽样，不能称第二次独立确认。

normalized reference修正常数权重是合理数值修复，但Gauss/三角函数和浮点bin端点仍未区间认证；Fraction证据只覆盖所测有理peak部分。当前已知解析误差、经验FP guard、BE误差三者应继续分列。

尚未认证：理想系数及η的完整FP包围；各kernel mean/range在所有Q×R×J上的向外舍入；高维BLAS点积、log/exp及erfc；实际proposal归一化/RN的逐点FP误差。这些不能仅由观察误差小、C_Z接近1或source验证有限点数补齐。因此当前保存统计可称“真实源条件矩/BE包络算法的普通浮点pilot”，不可称严格有限真概率证书或已证明新弱收益。

另外，`best_lower_CI`内置的理论锚1来自允许L=1并取λ↑1；L=1不在233点列表里。解释为完整cube算子（或233点再加L=1）下界合法；若只定义严格233点的弱范数，应明确这项额外锚点。它不扩展已抽样h0至10的阈值曲线为全弱曲线上包。

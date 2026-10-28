# R088：自由选择阈值降低覆盖上界

FULL ACCEPT（内部严格结果）。令T=nlogn−nloglogn。对每个固定s∈[4,5]，设

    mu(s)=(s−1)+(s+1)log((s+1)/(2s))，
    Cfree(s)=log[s/(1+mu(s))]。

每个c<Cfree(s)都有一个固定参数组与有限门槛，此后Cn≤Htau<T−cn。
固定s=9/2的明确有理元组给出Cn<T−.9n（eventual）。此外取s=5得

    limsup (Cn−T)/n ≤ −log[5/(5−6log(5/3))]，
    log[5/(5−6log(5/3))]>.9493。

所有新门槛仍未量化。对象保持Lebesgue背景、中心轴平行立方体、共同正半径、
一般L1输入及全部正尺度。R089的L=2^27结果此时仍为待审候选。

## 机制

原始有限选择定理允许任意tau>0，之前平方根选择并非算法条件。
固定x=(1−sigma)logn、t=n/x，先选chi>mu、chi(1−sigma)<d<sigma，
在真实历史开始前固定tau=n exp(−dlogn)。它与原几何核配合时，

    tau/k=rhi (logn)^2 exp[−(d−chi(1−sigma))logn]→0。

保留全部成本Htau=(1+A/tau)(tau+B−gamma n)，令Q=B−gamma n>0，
则Htau=Q+tau+A+AQ/tau，后三项除以n被
exp(−dL)+exp(−sigma L)+L exp(−(sigma−d)L)控制并趋零。
新条件sigma>chi/(1+chi)放宽旧平方根条件sigma>2chi/(1+2chi)，
由此降低准确价格−log(1−sigma)。绝不从Cn≤Htau倒推Cn≤F。

局部轮廓、中央删除、五事件、全匹配、真实源面方差及weightedcap使用已核且
与tau无关的R085几何；有限history、行列/缺失质量和原算子桥对新tau逐项重核。
先固定所有严格参数，再统一放大n；没有端点参数或n依赖配比替代。

## 完整比较

R087的旧完整成本满足Phi87>T−[c87+200000/(37073L)]n，
c87=619/400+log(37073/100000)<.556。L≥65536时误差<.0001。
本轮固定H88<T−.9n，故真实共同有效域上Phi87−H88>.3439n。
共同域须同时满足旧L≥2^80及新未知门槛；65536仅比较辅助条件。
旧公式下包不是Cn下界。原所有较低维数有效界保留按域取小。

公开比较延续此前已核文献口径：前导nlogn系数仍为1，本轮改善内部覆盖线性项。
未新检索文献，未称公开发表、全球最佳或所有覆盖最优。

## 实际运行与审查

作者预登记一次本机Fraction/log24运行PID32025，13项PASS；独审直接读
原R027及v3审核，独立log20运行PID32499，5组PASS；主另注册完整新旧比较
PID32481，4项PASS。三者原始receipt和来源hash分开保留，无旧基线重跑。
最小复核在new_tau/check.py、primary_comparison/check.py及audit/round088/full/check.py。
全文全域证明与有限算术的职责分别记录；来源hash见MANIFEST.json。

R086诊断已完成，无新界。R089显式候选独审进行中；R090任意固定s推广与R091
本机并行模型尚未接受。论文未改，下界暂停，七维专项停止，覆盖研究继续。

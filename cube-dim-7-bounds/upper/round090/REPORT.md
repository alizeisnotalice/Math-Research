# R090：任意固定配比与新的渐近线性系数

FULL ACCEPT（内部严格结果）。设T=nlogn−nloglogn，

    c∞=−log(1−log2)=1.181387061856…。

对每个epsilon>0，存在有限门槛，使此后所有整数维n满足

    Cn<T−(c∞−epsilon)n，
    limsup (Cn−T)/n≤log(1−log2)。

本轮没有给出显式维数门槛，未断言有限域Cn<T−c∞n。一般L1输入、Lebesgue背景、
中心轴平行共同半径立方体及全部正尺度保持。旧显式.9的log n≥2^27界继续保留。

## 构造与证明

对每个先固定的有限s≥4，令mu(s)=(s−1)+(s+1)log((s+1)/(2s))。
自由tau构造给每个c<Cfree(s)=log[s/(1+mu(s))]的固定参数上界。
s>5的局部估计重新推导：ell≥1/2，Q≤log2，中央删除完整计入J，
|J−1|≤2s/rho、|dlogJ|≤(4s²+6s+1)/rho、J1'≤5s/rho。
新五事件矩预算、source真variance≤2s²x、e依s的count floor及target/source
尺寸共同给出实际history上的缺失质量。新tau在任何history前固定，最终成本
始终为Htau=Q+tau+A+AQ/tau，不优化回F。

(1+mu(s))/s=1−(1+1/s)log2+(1+1/s)log(1+1/s)→1−log2。
每个epsilon先选一个有限s，再选全部严格参数，最后令n超过该tuple门槛；
这条顺序从实际有限构造推出limsup，不是把无限s形式代入局部定理。

## 准确新旧比较

旧R088族在s∈[4,5]的系数≤Cfree(5)<.95：
Cfree'(s)=[1+log((s+1)/(2s))]/[s(1+mu(s))]>0。
本轮c∞>1.18，故存在一个固定有限s和tuple，使完整Hnew<T−1.18n。
旧固定H88≥B88−1.5n>T−[c88+40/(11L)]n，c88=1.5+log(.55)<.903。
L≥4096时准确比较给H88−Hnew>.275n，须同时满足两函数真实有效门槛。
这不是在4096直接激活新界，也不把旧成本下包当Cn下界。各低域界按域取小。

公开比较沿用先前已核原始来源：前导nlogn系数仍为1，本轮改进内部覆盖线性项。
未新增文献检索、未称公开发表或全球最佳；c∞不是所有覆盖方法的最优性定理。

## 实际证据与复核

作者一次注册本机Fraction/局部积分辅助PID32534：新常数准确记账、log24包围，
s=10,q=.6的四个局部配置含g/R>.8，积分诊断误差最大2.22e−16。
独审另一次exact polynomial/log16运行PID33534，5组PASS，未重跑作者s10诊断。
主注册一次log24 Fraction完整比较，4项PASS。数值诊断不替代全s与全域证明。
来源hash、最小代码、原receipt、独立全文审查及主接受均收入本轮MANIFEST。
复核程序为general_s/check.py、primary_comparison/check.py、audit/round090/full/check.py；
在副本运行以保留PID原始回执。

R091的1.1/2^33准确候选正在单独审核，R092单侧源帽机制仅诊断。
论文未改，下界暂停，七维专项停止；研究继续。

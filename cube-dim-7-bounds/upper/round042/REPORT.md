# 第042轮：二阶矩列截断与新的显式高维上界

本轮严格得到：对每个整数n满足L=logn≥2^108，

Cn≤F(A_n(t),B_n(t)−2^-31 nL^-1/3)
  <nlogn−nloglogn−2^-32 nL^-1/3，

t=n/[L−L^(1/3)]，F(A,B)=(sqrt A+sqrt B)²。

这是完整原算子上界，适用Lebesgue、中心轴平行立方体、所有坐标共同半径、一般L1及全部正尺度；不是finite model cost或样本外推。主/两个独立Sol分别重建operator链、一维矩与显式门槛，整链及比较附录均已ACCEPT。leading nlogn系数仍1，不声称阈值/参数最优或公开发表。论文R041与本轮并行，round编号不是完成顺序；各自独立提交。

## 完整构造与两份预算

设s=L^(1/3)、d=1/s、x=L−s、t=n/x、β=1+1/t，h=d²/100、ζ=n^-3、κ=β−ζ、p=exp(−s/4)、M0=nhp。保持原largest-eligible selector及冻结真实increments，原corner未捕获测度μa≤λa；限制到D∞={||z||<t:存在zi<1−βmax(1,||z||)}后λa恰为α^-nλ，α=a/(a+1)。不是借由无控制包络极限。

后序radius r=1+u、source radiusρ=tq=t(1−v)，u,v∈[d/4,3d/4]，仅连接|u+v−d|<h及||z−w||<κρ。Target good set A用真实固定面统计Sw=Σ(|wi|−(r−d))+、准确均值d+(n−1)d²/(2r)、容忍t/100及Kw=#{|wi|>r−2d}≤12nd。固定面条件分布的Chebyshev/Markov及原cap扣除给nd/16≤λA≤nd；完整row deleted products在[e^-2p,e²p]，原source cap总量可单独支付，得到raw row≥e^-2M0/4、≤36M0，Eraw≥2^-10 ndM0。

新尺度中source的普通坏质量不能再相对M0小，不能照搬R035/040的全局坏source扣法。本轮改为精确column上包与二阶矩的energy扣除。固定source z，q=ρ/t，r_-=q+d−h。每个匹配target r≥r_-，归一化交长ell_i(r)=min(1,(r+κρ−|zi|)/(2r))随r递减且≥1/2；interval长度≤2h，所有面与deleted product共同给

Craw(z)≤Y(z)=(4nh/r_-)∏ell_i(r_-)。

这是逐点几何upper，不把无条件平均代替extreme或good-set bound。Source每张固定face有一个|zi|=ρ的pin坐标，其ell_pin=1−g/(2r_-)≤1，g=d−h+ζρ；其余n−1自由坐标均匀分布。准确积分给

E ell²=1−g²/(2r_-ρ)+g³/(12r_-²ρ)。

因为Ld³=1，正三次项和rq扰动有常数级补偿，不能省作o1。统一moment审计对全部L≥64给∫Y²dν_source≤48ndM0²，主完整证明保留更宽C2=144；未把C设置成1或用(Eell)²替代Eell²。

取H=2^20，P=Praw∩{Y≤HM0}。被删除的不是source总质量，而是其column energy：Eraw−Egood≤∫_{Y>HM0}Y≤(HM0)^-1∫Y²≤(144/H)ndM0≤(2^-10/4)ndM0。于是Egood≥(3·2^-10/4)ndM0、row/column均≤HM0。

## 同一真实history的缺额及一般输入

同一Borel核K的实际rowμ≤τ由某个真实target j的一次eligibility支付，不能拼接不同history。对共同缺额D=λ−μ，Eλλ−Eμμ=E(D_A,λ_P)+E(μ_A,D_P)≤M[D(A)+D(P)]≤M D(total)，A/P不交。故总缺额≥(Egood−τλA)/M。

选τ=sqrt(A0(B0−δ))，δ=2^-31nd、A0=(1+1/t)^n、B0=原cap预算。新x使A0≤n e^-s、B0≤nL，故τ/M0≤100s^(7/2)e^-s/4≤2^-12；配合λA≤nd、M=HM0得到缺额≥2^-31nd。全部正性由B0>8n、δ<n保证。

对每个实a≥8n^4，target实际半径<3、source实际半径>4；κ+2/(a+1)<αβ把增量点捕获传给同一个target centre，源确为它的前序。λa两侧精确共乘α^-n，energy乘其平方、degrees与target mass乘一次，有限a缺额不减。加回τ和自身a^-n、接原finite coverfactor，再固定n,t,τ仅令已证算子数值界的a→∞，得到F(A0,B0−δ)。任意L1和所有正radius用真实witness有限覆盖紧水平集及内正则性传递，不设support或尺度上下限。本路线不需要Green核、n≥5奇异子空间或四维线外连续论证。

## 显式门槛和scalar成本

显式s≥2^36逐项量化所有eventual条件，完整证明explicit_threshold/THRESHOLD_PROOF.md。目标面均值误差<1/3、good-stat和shift误差<1/5、平方log/deleted误差<1/10，所以row products保e±2。删除坐标upper严格使用logΠ_except≤−Tw+vj，不使用错误的−log(1−vj)≤vj。

Source cap/M0≤819200/s³<1/36；τ/M0由Taylor8化为<2^38/s^4<2^-12。矩审计L≥64及H≥4·144·1024、H≥36均为明确充分条件。最终B0≤T+2nd²+1，而A0+2sqrt(A0B0)可用两个Taylor8支付；除以nd的正余项为

R(s)=2/s+2/s^5+8!/s^7+2^9·8!/s^5<2^-32。

所有尾部是递减power，21作者和17独立Fraction/integer基点检查均实际PASS；没有有限网格替代全s≥2^36的解析尾。n≥s6/2、n≥6、h+ζt<d/2、ρ>4及finitea guard等隐含条件已全部补齐。因此最终严格负修正有明确域logn≥2^108，不声称最小起点。

## 与旧完整成本的准确比较

旧族用t_old=n/(L−3logL)，新族用t_new=n/(L−L^(1/3))，所以不能只比扣减δ或展示小数。比较附录读取旧12行原certificate.q准确分数及已证B_old≥T−4n/L，得每个自身有效的列明old完整root cost

F_old≥T−4n/L−2n sqrt(logL/L)。

在L≥2^240，2sqrt(logL)/L^(1/6)严格小于2^-35、4/L^(2/3)≤2^-158，两项和<2^-32，全尾单调。因此new full F_new<U42=T−2^-32 nL^-1/3严格低于old完整cost的共同lower，不是只把一个宽展示上界取紧。R034 q/Y门槛和R040参数依赖未知N一律只取自身valid域交集；跨全部t的unit-core族比较仅用已证asymptotic Uold−2>T给eventual结论，未把2^240强加给未知O项起点。

当前报告界取所有适用旧新界的最小值。公开对照仍限实际核过原文的范围，不宣称公开发表、首次、全球最好或穷尽。负修正nL^-1/3的阶强于此前n sqrt(logL/L)，不等于真实Cn的渐近等式，也不降低首项系数1。

## 证据、复核及后续

本轮大结论由完整解析证明和独立全链审查负责；本机只执行必要SymPy/Fraction一维矩、17/21门槛及准确比较，未执行n7、旧网格、下界或生产LP。最小复核代码分别在second_moment_audit/verify.py、explicit_threshold/verify_threshold.py、comparison/verify.py与独审目录。代码可能更新自己的receipt，应在可丢弃副本运行；源身份部分保留当时主机路径。原件hash在MANIFEST，冻结full/moment/threshold及比较分别锁定，不混报局部矩审核与operator审核。

R041论文已经把本轮完整proof和显式门槛写入，按单独24页稿提交；主仍需正式整合编译/渲染，未在此报告提前声称main PDF替换。R043继续检验分数矩/entropy的具体模板限制，不把它当Cn不可能定理；其他new-selector/commonmass方向保持。goal active，七维专项及下界停止，其他项目未改。网络失败时独立commit保留待推送，不强推。

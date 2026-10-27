# R006 四矩Radau独立Sol审计

2026-10-27。独立clone起点 `73646ed3d9818eb1d87a157aac09b8a709718ad7`。仅写 `audit/round006/radau4/`；旧轮次、作者目录及其他工作者文件只读。全文读取Sol `radau4_proof.md`/registration及Terra生成器、assembly、checker、certificate、ledger、report、receipt和manifest。Sol自己的code_review明确属于作者自审，本审未将其接受意见当作前提。

**结论：接受四矩Radau理论、解析严格改善及Terra端点到完整I的接口。** 本审未发现实质修复项。独立准确账本与768-bit标量检查证明

`I < 7723627498344/10^18 = 0.000007723627498344`，

且该端严格低于R005真实矩定义的确定预算表达式的有效下端。所有根、权重和修正方向有独立证明/区间检查，未依赖打印小数。没有全域g证书，没有新C7上界。没有重算4096个科学面板；数值系数生产仍由Terra负责，主集成人保留最终采用和证明责任。

## 1. 原模型、辅助变量和四个矩

原七维64个weight-one未单位化body方向、一个共同Cauchy时间及八个不交Hadamard块保持不变。T=detΣ_block^(1/7)独立同分布，U=T1+T2，q(r)=EU^-r；旧q均是完整Mellin积分，q0=1来自概率归一，不能在r0调用有Γ(r)极点的公式。

沿用a=1/8倾斜，令Δ=logU1−logU2，W=cosh(Δ/32)−1。四块目标仍是

`q_four(1/4)=2^-1/4 q(1/8)^2 E_a g(W)`，`g(w)=T16(1+w)^-1/4`。

因为16·Δ/32=Δ/2，T16(cosh(Δ/32))=cosh(Δ/2)，没有改变目标或漏常数。八块determinant Minkowski与两个四块和AM–GM最后给 `M≤q(1/8)^4(E_ag)^2/2`。

Rl=E cosh(lΔ/32)=q(a−l/32)q(a+l/32)/q(a)^2，l0..4。逐项展开(cosh−1)^j，四矩正是作者式(2)，尤其m3系数1/4,−3/2,15/4,−5/2，m4系数1/8,−1,7/2,−7,35/8。本审准确符号核验旧变量Y=2W(W+2)的 `EY=4m1+2m2=R2−1`、`EY²=16m2+16m3+4m4=R4/2−2R2+3/2`。

新q仅1/32、3/32、5/32、7/32；q1/8来自R004，q1/16、3/16、1/4来自R005，q0准确1。涉及的最小指数0，所有正指数负矩有限，所以只需四个矩，没有暗需第五矩或complex q。

条件于块内部角度，T是正系数/Gamma4半径，具有(0,∞)完整连续支撑；U、正tilt、log差及cosh变换使W在(0,∞)每个开区间有正质量、在0无原子。这足以支持下述strict moment/weight证明，不依赖经验矩正定性。

## 2. Weighted orthogonality、两个正根

令D=m1m3−m2²，它是正测度w dP(W)中1,w的strict Gram determinant，故D>0。给定

`alpha=(m2m3−m1m4)/D`，`beta=(m2m4−m3²)/D`，

p2(w)=w²+alpha w+beta满足E[W p2(W)]=E[W²p2(W)]=0。右端分别仅用至m3/m4，方程和参数分母正确。

作者自足root论证成立：若p2在(0,∞)至多有一个sign-changing zero，以所有这些根的线性因子作degree≤1乘子，则p2乘该乘子在正半线固定非零符号。乘子把所有奇重正根变成偶重，其他根在正半线不改变符号。完整支撑使其对w dP积分严格非零，与两条正交式矛盾。因此p2有两个不同正根0<d1<d2，alpha<0、beta>0、discriminant>0。

该证明没有把区间中点的根当作真节点，也无须第五数值矩。beta正也可由m2m4−m3²的strict Cauchy–Schwarz核对。数值稳定formula d2=(−alpha+sqrtDisc)/2、d1=beta/d2正确，已用准确区间算术检查其分母和顺序。

## 3. 三个正权、degree4 exactness

节点0,d1,d2的degree2 cardinal Lagrange polynomials记Li，定义λi=E Li(W)。对degree≤4多项式P，除以w p2，得到P=w p2·ell+r，degree ell≤1、degree r≤2。正交式消去首项期望，首项在三节点也为0，所以quadrature exact through degree4。

显式λ1=(m1d2−m2)/[d1(d2−d1)]、λ2=(m2−m1d1)/[d2(d2−d1)]，λ0=1−λ1−λ2均正确。对degree4的Li²应用exactness，得λi=E Li²>0，因完整支撑不集中在有限多项式零点。常数P1给Σλi=1。

可选λ0=detH2/(m2m4−m3²)同样成立：L0=p2/beta，E p2·D等于detH2；H2是1,W,W²的strict Gram matrix。Terra本轮没有用此公式修补权重，也没有clip负端。source中权正、mass及0..4矩residual均通过；这些是实现诊断，exactness本身由上述证明承担。

## 4. CM T16和五阶Hermite全域上界

T16的16个根cos((2j−1)π/32)都严格介于−1,1，leading2^15。因此g是常数2^-15/4乘16个 `(w+c_j)^-1/4`，所有c_j=1−root>0。每个因子完全单调，Leibniz使乘积各阶导交替严格符号，特别g^(5)<0、g0=1、g正严格递减。

P4匹配g在0,d1,d2的值及在d1,d2的一阶导，共五个独立条件。对每个w≥0，重复节点Hermite余项

`g(w)−P4(w)=g^(5)(ξ) w(w−d1)²(w−d2)²/5!≤0`

在整个非负半线成立，包括w>d2；ξ仅需在相关有限区间，g在那里C5。节点处等号。P4只用degree4，故其期望可积，不需要积分某个未知第五矩来支付余项。

Exactness给 `E_ag≤B6=λ0+λ1g(d1)+λ2g(d2)`，0<B6<1。区间计算的职责是包住真矩定义的节点/权重，并非证明某个任意区间或中点插值多项式本身全域majorizes g。Terra实现遵守这条分界。

## 5. 相对R005复合majorant的严格改善

旧Y=cosh(Δ/16)−1=2W(W+2)，T8(T2(x))=T16(x)，所以旧f(Y)=g(W)。本审准确符号核了Chebyshev复合恒等式。

R005的真矩quadratic Hermite P_old全域支配f；由于f'''严格负，其接触点恰0和d_old>0。Q(w)=P_old(2w(w+2))的degree≤4，支配g，正半线接触点恰0及唯一w*>0。新的两个正节点不同，因此至少一个不等于w*，其权严格正。于是

`B6=Σλi g(di)<Σλi Q(di)=E Q(W)=E P_old(Y)=B5`。

这严格比较两个expectation budget表达式，exactness在最后两步起作用；并不需要另证明新P4在每个w点都小于Q。倾斜law和四块目标两轮相同，故比较合法。B5是旧真实矩定义的Radau表达式，不能把R005某个rounded上预算的lower误标为同一对象。

## 6. δ6下贡献与R005真实表达式lower方向

corr(d)=1−T16(1+d)^-1/4对d>0严格增加。若root和weight区间分别包含真值且权正，则

`δ6,L ≤ λ1_L corr(d1_L)+λ2_L corr(d2_L) ≤1−B6`。

下节点、下权必须组合为正乘积，方向正确；取整个sum的lower后得到B6≤1−δ6,L。Terra的correction(d,4)由四次doubling计算T16，差v每步用2v(t+1)，随后四次根及正分母 `(T16−1)/[R(R³+R²+R+1)]`。本审独立重算声明下输入的表达式，严格证明存储δ下端低于其768-bit下端。

R005比较方向独立审查如下：真实旧μ=R2−1、ν=R4/2−2R2+3/2；旧δ(μ,ν)随μ增、随ν减。用μ.upper与ν.lower给δ旧upper；验证0<δupper<1后，1−δupper是B5 lower。再用q1/8.lower及正单调幂、正F，得到旧真实确定预算I5的lower。作者代码这条方向正确。本审以准确Fraction极端矩端、768bit重新构造更紧有效I5safe.lower，并证明交付声明的旧lower≤独立有效lower≤真实I5。

因此最终严格链是

`I6_768upper < Terra dyadic upper <7723627498344/10^18 < declared I5 lower ≤ independently valid I5 lower`。

最后两项不是actual I的lower。若只与旧rounded upper比较，本来不足以给真预算改善；本轮确实完成更强的真预算表达式比较。

## 7. 四组q参数、完整误差与旧证书复用

两个旧certificate复制件均与原文件逐字节相同，哈希同时符合复制manifest和原manifest的certificate项。reused_q_1_8及reused_q_R005字段逐项相同。没有旧q面板重算、没有q0 Mellin调用。

|r|a=r/2|rho=3a/4|halfwidth|面板数|eta|
|---|---|---|---|---:|---|
|1/32|1/64|3/256|1/128|2048|2/3|
|3/32|3/64|9/256|1/64|1024|4/9|
|5/32|5/64|15/256|1/32|512|8/15|
|7/32|7/64|21/256|1/32|512|8/21|

所有网格覆盖[0,32]、rho<a，准确中心(2j+1)halfwidth。每组Gamma normalization πΓ(r)、strip bound Γ(a−rho)Γ(a+rho)m(a−rho)m(a+rho)、eta^129/(1−eta)及Γ(a+16)tail全部分别支付，不能复制R004的2/3误差给后三组。

复域包络来自Gamma(v)m(v) log-convex对称端点乘积；Cauchy圆在Re参数>0域内。degree128遗漏129，raw总误差32Mstar eta^129/(1−eta)，32是总长度。N16 recurrence尾为 `m(a)^2Γ(a+16)^2/[πΓ(r)31·32^31]`，涵盖全部未积分t≥32。正tail仅加q上端，Taylor误差双侧支付。

源码ctx.prec512、cap129、signed even coefficients准确积分、coefficientwise conjugate实现analytic paired series均与规格一致。16 pilot终态被计入4096且剩余从index4开始；ledger的4096个(group,index)唯一完整、全部positive completed。本审未重算任何系列系数，系数区间依赖Terra producer及certified FLINT/Arb库；完整分析证明、账本与标量复核是本审接受范围。

## 8. 独立有限运行、可携带入口及失败记录

本审先登记后执行，入口minchecker.py只读、stdout-only，默认合并后项目相对目录，也可用--numeric-dir/--theory-dir/--r004-dir/--r005-dir指定独立源。没有--write或隐形result写入。shell仅将stdout保存到本审results.json。

精确Fraction汇总四组dyadic积分上下端，768-bit log-Gamma重组q的完整余尾；均被Terra交付区间包含。进一步以准确有理区间算术重建R、四矩、D、alpha/beta/discriminant和weights；sqrt分别从区间两端用Arb外包，根区间保持有理向外。不采用midpoint替代。独立根/权区间逐项被交付区间包含，全部正且root1.upper<root2.lower；mass和degree0..4residual包含理论值。lower节点/权的合法性因此先于δ验证建立。

两次最初执行都在所有数学断言已通过后的JSON组装阶段失败：Python对超过4300digits的大Fraction字符串施加限制。先改根/权摘要为768-bit向外dyadic，再修helper避免对已有巨大Fraction重复转str，第三次完成。该修复只改变序列化路径，没有改变准确内部区间、参数、科学面板或验收不等式。数学断言失败0、实现失败2/修复2。成功脚本0.2850806670030579秒，Python3.15.0、SymPy1.14.0、python-flint0.9.0。共享基础库，不称独立算术内核。

所有Sol/Terra/旧证书输入在运行前后hash相同，Terra manifest逐项验证，proof/registration source锁一致。results保存完整源hash、exact moment接口、向外dyadic根权摘要、δ及新旧预算比较。Terra receipt本轮未报告科学或checker失败；本审不凭缺少历史日志补造成功过程，只验证交付final文件与实际本审运行。

## 9. 完整I transfer和停止边界

原full-boundary Gaussian/Minkowski界 `I≤F M^(7/6)`仍属同一个body64共同时间过程，F=`1/[105·2^(25/6)π^(13/12)]`。本审从原C^(7/6)H/896路径重组F，先全boundary后原law的14face×64orthant对称；没有把单shape正face假称对称，也没有遗漏short time或singular neighbourhood。

因此接受新的完整I上预算。一般L1、同一个D对全部共同r>0、正部核全空间有限源的既有接口保持适用；仍须独立global essential g0证书才能用 `C7≤7^(7/6)I0/g0^(7/6)`。本轮不提供g，也不声称新C7界或model exclusion。

本批结束并停止，不开启新策略、额外指数、mesh或MC。无MATLAB、子代理、commit/push；未修改旧轮次或其他工作者文件。

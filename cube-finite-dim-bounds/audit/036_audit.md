# 036 三维完整连续上界候选最终验收

2026-10-25。接受固定body+face模型、n=3、τ=1/3给出的严格结论

\[
 C_3\le\sqrt{\frac{1547906585017734375000000000}{24707828464688747938727861}}
 <\frac{1979}{250}=7.916<8.
\]

这里C3是恒定Lebesgue测度下中心立方体极大算子的弱(1,1)常数，所有坐标用同一半径。没有假设真实minimum等于有限样本minimum，没有以原先目标g=.00524替代033的实际完整下端。账本由主代理负责，本审计不改账本。

## 锁源和本次实际执行

`036_source_snapshot.json`锁46项主目录输入，包括最终manifest SHA256 `19665768679a50ecb1ae2efed8ec5bb05ceb4a4762abfc08883f971975b578a4`及其25个payload；032、033、034的既有最终证据、031迁移、wholeface短time证明以及准确成本文件均锁定。最终报告的历史pending字段只表示生成时尚未验收，不能作为当前数学断言；本次将其所需接口逐项连接。

新`036_check.py`只用stdlib Fraction/JSON/hash及python-flint Arb，不导入、不执行作者checker、selfcheck、aggregate、final_verify、Green evaluator或point函数。独立512bit标量路线首次实际运行通过，无实现修补，内部约0.343秒。它复算全部8256条凸弦、索引/面积、准确总和、短时合成及成本，并读取外部原raw12份文件核全部SHA和保留字段。新增Green节点评估0；未重跑8990784节点、未变更精度/网格、未重搜参数、未读写节点cache。

共享输入明确为032已审真实96锥/ray表达、033已独立fresh的59976节点下界、034已审连续插值与短时理论，以及原作者384bit全批执行产生的范围/mean证据。共享算术库为python-flint/Arb；本次8256个标量包围在512bit独立计算。本次不将作者全批称为另一套独立节点实现；也没有仅因作者或主selfcheck PASS自动接受结论。

## 完整覆盖、原始记录和逐块计算

固定配置是ε=9/5000、d=32、grid128、rho=5/2、node_bits128、precision384、τ=1/3、cone SHA23258e3d…。冻结checker SHA `4f89a377fbdf263061c2d6525d2101a773f7a6ba4127bb367bd3b7328ed57eb5`的实际源码就是034所审接口；本批逐份核config和该hash。

独立重建所有0≤j≤i<128的8256个indices；序号i(i+1)/2+j mod12分成12个两两不交的688块集合。每part的实际集合与预定集合完全一致，全集无缺漏、重复或越界。j<i取权2、i=j取整方块权1，权和16384；每方块面积1/128²，所以总面积准确为1。对角线只是零测边界，diagonal整方块没有再减半。

每条compact记录的L、U、mean上端m和contribution均是零radius准确dyadic。独立确认0≤L≤m≤U、U>0、contribution的dyadic与其Fraction文本相等。令l=√L、u=√U，以512bit正根重算

\[
 \operatorname{mean}H^{3/2}
 \le l^3+\frac{u^2+ul+l^2}{u+l}(m-L).
\]

全部8256个独立贡献上端均不大于已保存的准确contribution。该式由凸上弦和非负斜率成立；U=L时m=L仍正常，不除U−L。L=0也正常。range和mean来自已审的连续Bernstein变换及准确mean，而非仅离散节点min/max。

本次还实际读取scratch内原raw档12份，逐份SHA与最终provenance相同；每条compact保留字段与同序raw逐项相等。原whole_square有理端点、weight、nodes=1089、chord_mean_upper与乘面积后的contribution均另核一致。原四个复用pilot确为(127,57)、(127,127)、(127,0)、(0,0)，对应4356节点；其余8986428节点是原全批fresh，总计8990784。这里只核原记录与执行来源，不把这些原节点误报为本次新算。主原12processes均exit0的实际执行来源在锁provenance中记录；raw的外部路径仅用于本次可选增强审阅，便携默认checker不依赖该外部档。

所有contribution用Fraction准确求和，恰等于staging/all_result的IH_upper_exact及provenance准确总和。其展示值为

\[
 B=0.0005620627584450029912931962438471980490\ldots.
\]

报告小数不作严格端点，准确分子分母在`036_results.json`。压缩并未改变B或最终packedI；未使用浮点求和。

## 连续误差、全部短时及最终I

034已经独立检查并接受本次d32参数的tube、Cauchy、roots alias、DCT、Bernstein和准确mean接口，本次重新读冻结实际代码确认仍是同一链。Hε=∫ε∞ptdt，完整48对反锥有符号Hermite实DD合并，零/重节点均依032正常处理。eps正则化不是删除某个signed项或忽略Γ项。

ψ≥3(|ξy|+|ξz|)给polyellipse tube，h=1/256、rho5/2的高度21/5120、δ=277/640000>0。Epoly=176M(2/5)^32 <7.706294e−8。128bit有理roots替代由Cauchy导数给2Dh2⁻128，并在node ball内先付后传播；Arb mag半径向上舍入已包含。准确DCT和degree32 Bernstein给整个闭方块range，不是只给33×33节点。mean(P)用μ0=1、奇数0、偶数−1/(k²−1)，range与mean同时加Epoly；m取U与mean(P)上端+Epoly的min，并要求m≥L。全批没有不合格记录被容差修补。

这一连续代码/误差证明已由034接受，本次证据验收核其冻结代码、执行配置、每块保留结果及准确聚合。没有出现需额外节点校准的疑点，因此按预登记不重跑约899万节点。该范围说明不意味着离散节点本身足以证明连续正则化函数的range。

对Tε=∫0εptdt，034的wholeplane density事件界及真实marginal给正象限norm误差

\[
 E=2^{-4/3}\frac34
 \left(\frac{256}{3\pi^3}\right)^{1/3}
 \left(\frac{16}{3\pi}\right)^{2/3}
 \epsilon^{4/3}.
\]

该式已包含方向线附近、整个plane及空间无穷尾；没有另遗漏奇性帽。独立512bit以正三次根重新计算E，而非调用作者exp/log路线，E≈.00012996912954638279224。Minkowski给

\[
 I=\int_{[0,1]^2}G(1,y,z)^{3/2}dy\,dz
 \le (B^{2/3}+E)^{3/2}.
\]

独立该表达式的上端不超过冻结packedI；packedI与manifest逐字一致，是准确dyadic，展示约.00057822817482354580122732。因此可向外取

\[
 I\le I_0=\frac{23129127}{40000000000}=.000578228175.
\]

此处E是未知实际短时norm的上界，未将其下端当成实际误差下界；组合不取不利端点。

## 完整g下端及准确成本

033已fresh独立计算136闭块59976节点，完整min下端为.00524264797871796797512876…。本次重新从033全部136条H_lower_exact准确取min，等于其保存全局min，故

\[
 G\ge g_0=\frac{2621323989}{500000000000}=.005242647978
\]

在完整cube成立：033认证的H1/50≤G，在正face的闭块覆盖和cube群/齐次向内延拓均已审。这里只取已有准确下端的向外有理数，没有新minimum搜索，也没有倒退到旧预设.00524。

cube对称和face射线Jacobian给J=24I。设c=g0/3，并取

\[
 K=\frac{(G-c)_+}{8(g_0-c)}.
\]

齐次层积分给∫K≤√27 I0/g0^(3/2)。独立Fraction得到cost²为本报告首式，且

\[
 (1979/250)^2-27I_0^2/g_0^3
 =\frac{22590856465666228399282662701}
 {1544239279043046746170491312500}>0.
\]

所有量正，平方比较严格保序，cost≈7.9150759240454368739853，故严格低于7.916。另以有限交错atan准确核Machin上界π<355/113；角度分支由tan(4atan(1/5)−atan(1/239))=1及0<角<.8<π/2确定。该π上界不是本n3成本必要条件，也不提前接受后续n4锐化。

## 合法核、全尺度和一般L1闭环

本次重新核031的各实际前提仍是同一固定τ=1/3过程：正权未单位化方向、正保守自伴卷积半群、ψ~|ξ|、form域H1/2/operator域H1、所有小大cutoff及weighted jump率。032完整真实96锥只用于该同一pt/G，不能用另一谱分数核替代。

由真实pt的短时质量和长时Ct⁻3，G weak-L3/2、locallyL1、a.e.有限；其正齐次及正性保证g>0和J有限。半群Jensen的全空间分布证明给H=(G−c)+∈L1、LH=μ−δ0，其中μ是正概率、μ({0})=0。该全分布链同时处理body/face奇线、非光滑level和无穷尾，不能只在线外检查点源号。新的g0/I0只改阈值和成本，不改该源证明；c>0任意本来就在031许可范围。

因此K≥0、偶、L1、K≥1cube/8且LK≥0于原点外分布。031/026的核比较对一般K的Lévy源允许无限：CND符号、ρ≤Mψ、有限ν截断的L2极限和接触集正性，不作∞−∞。同一个障碍u、σ和接收集D控制所有共同半径r；ψ的一齐次给ρr=r⁻1ρ(rξ)≤Mψ。先有理r的可数坏集，再由L1核缩放连续扩到全部r，非负有界紧支撑输入最后单调逼近一般|f|∈L1。没有把各坐标独立半径、多参数强极大或非恒定背景混入结论。

故实际g/I证据、合法majorant、质量成本与完整一般L1/all commonr链已闭合，可接受本报告首式。此前031/034的“targets only”并不是当前仍有未付前提，而是历史有限任务的正确范围；本次正是完整候选验收。

## 复现、MATLAB及终止

主整合后最小复核：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/036_check.py
```

审计clone实际运行加`--main /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds --raw-archive /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/round011o-run/bridge-raw-archive`。默认只读仓内compact和既有锁依赖，不需要原raw或cache。`--save`可指定自己所属输出位置。

MATLAB最小脚本仅核area、secant恒等式、Machin算术与严格cost²差，不替代8256准确包围或连续理论；实际回执另见036_matlab_receipt.json。有限036任务完成后锁manifest停止。本报告不接受11Q sharper n4结论，不启动新研究轮；后续完整第二审计及两篇TeX/PDF按主另行安排。

主串行MATLAB已实际exit0（session41508），同SHA762a728f3b68f9fa450e5155cabf1cff38b952e4df4139290f6f7e1655e1448a；审计实际读取log，area1/chord/Machin及精确cost²和严格差均PASS。vpa展示7.9150759240454368739853216081571569262533555181128，log SHA6da935c9e251ba56df803dc985cecc6c09d44bae5368cff8f4f25f9a32160369。收尾全部锁源未变。

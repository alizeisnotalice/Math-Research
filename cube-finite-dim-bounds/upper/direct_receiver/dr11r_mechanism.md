# 11R：直接接收权对偶、连续最大值及两盒准确证书

2026-10-25。固定Lebesgue R^n，中心轴平行方体Q(x,r)=x+[-r,r]^n，每个接收点只用一个共同radius。此文件不涉及Green生成元。六个已登记配置限n3/4、r∈{1,2}、h每receiver box一个常数；没有新的Cn上界。

## 1. 完整量词及一般L1条件接口

给f∈L1，F=|f|≥0，λ>0，有限measure可测E及可测r:E→(0,∞)，满足

(2r(x))^-n ∫_(Q(x,r(x)))F(y)dy>λ，a.e.x∈E。

任意非负可测h:E→[0,∞]，定义

Th(y)=∫_E h(x)1_(Q(x,r(x)))(y)/(2r(x))^n dx。

kernel joint measurable。若Th≤1 a.e.dy，则非负Tonelli给

λ∫_E h dx≤∫_E h(x)avg_Q F(x)dx=∫F(y)Th(y)dy≤||f||1。

特别∫h自动有限。h可以依赖f、λ、E、r，不能省略实际容量Th≤1。

若能对每个这样的见证数据、每个有限E统一构造h，使∫_E h≥|E|/U且Th≤1，则得到整个中心cube maximal的弱型常数≤U。一般完整bad set先取与[-N,N]^n交再MCT；不先假定bad set有限measure。

所有正尺度的measurable选择无额外问题：对f∈L1，固定x的avg_Q随r连续（边界Lebesgue-null与DCT）；在bad set上可按预列正有理r的第一个严格见证选r(x)，得到countable-valued measurable选择。**未证明的是该任意E/任意选定r的一致mass构造**，不是选择器本身。

恒定Lebesgue实际用于(2r)^n正规化、receiver volume、dxdy Tonelli及输入质量。换w(x)dx后的分母与adjoint不同，不能直接移植。这里未排除原子dual pricing；但本文一般算子结论的输入按f∈L1表述。

## 2. 两receiver boxes的真正连续最大值

E=B1∪B2，interiors不交，边界零measure任意归属；r_i恒于B_i，h=h_i≥0。设B_i=Π_j[l_ij,u_ij]，则

a_i(y)=|B_i∩Q(y,r_i)|/(2r_i)^n=Π_j ell_ij(y_j)/(2r_i)^n，

ell_ij(t)=max(0,min(u_ij,t+r_i)−max(l_ij,t−r_i))。

Th=h1a1+h2a2。ell为连续分段affine，其全部断点包含l_ij±r_i,u_ij±r_i。对每坐标取两盒断点之并，在每个网格cell中Th为multi-affine，等于其2^n顶点值的tensor线性插值（权重非负、和1）。所以cell内最大值在顶点；外部越过最外断点，所有相应因子为0。于是

sup_(y∈R^n)Th(y)=max_(y属于有限断点笛卡尔积)Th(y)。

这不是receiver中心采样或任意均匀grid近似。闭盒r_i等于半边长时重复断点只去重，不漏端点。所有配置均rational，负载向量v(y)=(a1,a2)准确为有理。

## 3. 常数-per-box的准确LP与弱对偶最优证书

maximize m=|B1|h1+|B2|h2，subject to v(y)·h≤1于全部断点、h≥0。只保Pareto最大负载向量合法，因为h≥0。对任何有限非负input prices alpha_k，若

Σ alpha_k a_i(y_k)≥|B_i|（i=1,2），

则每个可行h均有m≤Σalpha_k。这只是有限weak dual，未调用无限LP strong duality。给出可行h与相同目标的prices便准确证明该constant-per-box类最优；不证明对任意可变h(x)也最优。

## 4. 六个配置的全部连续约束与结果

每维n=3,4，B2=[0,4]^n、r2=2。三配置分别为：

- volume_critical：B1=[−16,0]×[0,4]^(n−1)，r1=1。
- face_adjacent：B1=[−4,0]×[0,4]^(n−1)，r1=1。
- vertex_adjacent：B1=[−4,0]^n，r1=1。

两盒interiors不交，都是正volume。记V=4^n。前两配置所有最强连续容量约束恰为

h1+(1/4)h2≤1，h2≤1。

第一个在y_A=(−1,2,...,2)，第二个在y_B=(2,...,2)。所有横向因子≤1且在y_j=2同时取1；沿分离轴，半径分别1与2的断点负载Pareto仅(1,1/4)、(0,1)，故这是完整约束而非只必要条件。

vertex配置同理得到h1+rho h2≤1、h2≤1，rho=4^-n，见证y_A=(−1,...,−1)、y_B=(2,...,2)。各坐标共同断点上左盒normalized overlap只能0或1；若左盒product为1，右盒product至多(1/4)^n；否则左product0、右至多1。multi-affine定位补足所有连续点。

统一取h=(1−rho,1)。设v1=|B1|、v2=|B2|。prices为alpha_A=v1、alpha_B=v2−rho v1≥0，满足全部列等式及sum prices=m，所以准确最优。

|配置|n|h1,h2|最优m|receiver measure|局部receiver比值|
|---|---:|---|---|---|---|
|volume_critical|3|3/4,1|256|320|5/4|
|face_adjacent|3|3/4,1|112|128|8/7|
|vertex_adjacent|3|63/64,1|127|128|128/127|
|volume_critical|4|3/4,1|1024|1280|5/4|
|face_adjacent|4|3/4,1|448|512|8/7|
|vertex_adjacent|4|255/256,1|511|512|512/511|

表中比值仅这些E/r及constant-per-box类的准确接收mass结算，不是Cn。volume_critical中alpha_B=0，所有可行constant h都不能取得比m更多质量；这是此限制类的损失例，不能排除更细h(x)。

**中心-only反例。** 六例在h=(1,1)时，两个receiver centers的负载都准确等于1；真正global max却是5/4（前两）或1+4^-n（vertex），在上述断点y_A取得。只查receiver centers会错误批准不可行h。优化h的global max准确为1。

## 5. 哪些是已覆盖基准，哪些只在有限类改进

仅两固定半径时h≡1/2的U≤2是平凡finite-scale计数。更精细的uniform envelope也已有普通解释：对任意disjoint scale labels E_j及r1<...<rm，

Σ_j |E_j∩Q(y,r_j)|/(2r_j)^n≤∫max_j 1_(Q(y,r_j))(x)/(2r_j)^n dx

=1+Σ_(j≥2)[1−(r_(j−1)/r_j)^n]。

因此固定r=1,2给2−2^-n。它等价于有限averaging kernels的pointwise max envelope，不是新机制突破；其常数随尺度数/跨度增大可发散，不能处理全尺度。

对“最多两个disjoint axis boxes”的更强uniform baseline是5/4：两盒必有一坐标分离，沿该坐标用left/right halfline上包围。r1=1左侧、r2=2右侧时，normalized lengths和在y=−1最大5/4，其他坐标因子≤1；相反次序由反射。两radius相同时和≤1。所以h≡4/5适用于整个两盒限制类。前表equal-face和vertex的非uniform优化严格改善该baseline，volume_critical则在constant-per-box类准确达到5/4。

该改善不稳定于任意拆盒/并盒：分离轴因两盒而存在，一般scale标签集可交错。有限kernel envelope的worst load可由nested cube shells实现，而shell不是单receiver box。给每对盒单独可行h再相加会重复花同一input容量。因此本批没有任意E/任意r的一致mass lemma。

## 6. 有限空间逼近的正确预算及准确障碍

若坚持保留任意给定的可测E/半径标签，不能默认它可由内含full boxes在measure上逼近。例fat Cantor S⊂[0,1]由第k阶段删除2^(k−1)段各长4^-k构造，总删除measure1/2；S闭、empty interior、measure1/2。E=S×[0,1]^(n−1)有正measure却不含positive-volume axis box。该为解析量词反例，不是第七个数值配置。

但这**不是实际严格bad set的障碍**：每个固定r的平均随x连续，故U_r={avg_rF>λ}开放。对有界open A_N={MF>λ}∩(−N,N)^n，取closure包含于某个U_q∩A_N的admissible dyadic boxes（q正有理）。小盒总存在；最大admissible dyadic boxes内部互不交并覆盖A_N除dyadic边界零集。每盒选一个对应q作共同radius。有限prefix的volume增加到|A_N|，不用近似原有measurable selector。

因此一个足够且更精确的缺口是：对**任意有限数量**的interior-disjoint receiver boxes和任意assigned positive radii，能否统一找constant-per-box h，使容量≤1及mass≥union volume/U，U不依赖box数量/尺度？若有，前述有限prefix、N→∞即给全L1/all-scale上界。本批只解决box数2和r=1,2，远未证明这个uniform finite-family mass引理。它是充分条件，不宣称与最优Cn等价。

一个正确、可检验的有限转移接口如下。令E_j为实际disjoint可测receiver标签，固定r_j>0，surrogate F_j为可计算box unions（允许surrogate重叠，但按labels求和），weights a_j≥0，epsilon_j≥|E_j△F_j|。若

T_F≤1−delta，Σ_j a_j epsilon_j/(2r_j)^n≤delta，

则T_E≤1；而实际接收mass至少Σ_j a_j|F_j|−Σ_j a_j epsilon_j。两个结论分别由kernel L∞界和triangle inequality准确给出。

该预算有限或countable可用，前提是两类误差和真实收敛；不要求E本身包含盒。全尺度r_j↓0时普通measure误差不够，须支付带(2r_j)^−n的容量误差。**还缺的核心几何引理**是对任意E/可测r统一产生满足这些容量和mass预算的权重，并有固定U；本批没有证明这一存在性，也没有无限LP strong-duality结论。

## 7. 计算与范围

实际Python Fraction以6进程分别构造六例完整断点集，总8232顶点；准确Pareto/2变量LP、可行h及finite prices全通过。没有uniform数值grid、随机点、额外配置或Green调参。MATLAB由主串行同样断点/已知h与prices交叉，来源和回执在report/results。Python输出是准确有限证书，Matlab浮点只作独立实现交叉。

本批停止，不追加实验。一般L1/所有共同radius的定理只是第1节条件接口；初次锁定输入时main为U3=8、U4严格14753/1000；交付前再读主账本已由独立036/037更新为严格U3<1979/250、U4<14371/1000。这些主成果不来自本receiver局部比值，本批未改变账本。

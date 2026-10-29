# 连续uniform坐标的sharp Lipschitz MGF候选

先运行 `uniform_probe.py`，再给如下自含证明。它不修改已冻结MCD证明、upper_math、receipt或任何CI；此处为独立候选，尚未应用于tail数据。

## 1. 1-Lipschitz函数的单调重排仍1-Lipschitz

设U~Unif[0,1]，g:[0,1]→R为1-Lipschitz，q为g(U)的非降quantile。g连续，其range为一个闭区间。任取range中a<b，沿连续路径从值a到值b，取首次到b前最后一次到a的两时刻。中间整段取值都在(a,b)，而Lipschitz性使该段长度至少b−a；路径方向反过来时同理。因此

`Leb{x:a<g(x)<b}≥b−a`。                       (1)

对0<u<v<1，若a=q(u)<q(v)=b，则quantile定义给F(a)≥u、F(b−)≤v，故中间质量≤v−u。结合(1)，

`0≤q(v)−q(u)≤v−u`。                         (2)

原子对应q的平坦段，不使该论证失效。q连续延拓至0,1仍为1-Lipschitz。由quantile定义，q(U)与g(U)同分布，不是假设g本身单调。

## 2. 中心化后的convex order

令m=∫q，a(u)=q(u)−m，b(u)=u−1/2。二者均非降、均值0。由(2)，q(u)−u非增，所以D(u)=a(u)−b(u)非增且∫D=0。非增函数的prefix平均至少为总平均，故

`∫_0^v a(u)du≥∫_0^v b(u)du`，所有v∈[0,1]。

等价地，各suffix积分a不超过b。对任何实t，非降性给

`∫_0^1(a(u)−t)_+du = sup_v ∫_v^1(a(u)−t)du`。

逐v用suffix比较，得到stoploss支配

`E(g(U)−Eg(U)−t)_+ ≤ E(U−1/2−t)_+`。         (3)

同均值加(3)即可逐个正铰链组合证明所有piecewise-linear convex函数的期望比较，再由紧区间上的凸函数逼近推广至任意有限convex函数。因此

`g(U)−Eg(U) ≼cx U−1/2`。                    (4)

没有用一般浓度定理作为黑箱。由q的Lipschitz性，m−q(0)≤1/2、q(1)−m≤1/2，中心化变量也确落在[−1/2,1/2]内。

## 3. Exact uniform MGF与1/24

对任意实t，(4)应用于exp(tx)给

`E exp[t(g(U)−Eg(U))] ≤ sinh(t/2)/(t/2) ≤ exp(t²/24)`。 (5)

t0时商定义为1。后一不等式可直接逐系数证明：令z=t/2，

`sinh(z)/z=Σ_{k≥0} z^{2k}/(2k+1)!`，

`exp(z²/6)=Σ_{k≥0} z^{2k}/(6^k k!)`。

而(2k+1)!≥6^k k!，由相邻比(2k+3)/3≥1归纳即得。因此无需无限乘积或未核验特殊函数定理。g(u)=u的variance=1/12，说明中心MGF的普适二次系数1/24不能再减。

一般宽h的uniform坐标X及Lipschitz常数L的函数f，只需单位化，(5)变为

`E exp[t(f(X)−Ef(X))] ≤ exp[t²(Lh)²/24]`。   (6)

## 4. Conditional Doob与factor6 tail

隐藏X_i独立、各自在自己的bin上连续uniform。设H具有坐标Lipschitz常数L_i，κ_i=L_i h_i，S=Σκ_i²。条件于prefix，未来坐标期望仍是当前X_i的L_i-Lipschitz函数；独立性保证当前X_i的条件law仍uniform。对Doob增量逐步使用(6)：

`E exp[t(H−EH)]≤exp(t²S/24)`。

Markov优化t=12a/S得到

`P(H−EH>a)≤exp(−6a²/S)`，a>0、S>0。         (7)

S0时H恒定，另行处理strict边界。这里改进的是连续uniform加Lipschitz条件下的条件MGF，不能从仅bounded oscillation的旧McDiarmid结论直接换常数。

对非空positive gated项F_b，H_b=log F_b。若已有真实均值上包u_b，Jensen仍给EH_b≤log u_b，故λ>u_b时

`P(F_b>λ)≤exp[−6log(λ/u_b)²/S_b]`。          (8)

空段、zero weights、u0、S0与strict阈值规则照旧MCD证明处理。

## 5. 当前source为何能用同一C log product surrogate

当前h=δ/6，旧source的range cap v_{i,s}确为一个整bin Lipschitz上界L_{i,s}乘h。endpoint候选覆盖所有可能的非零slope，N0仍含完整background项。令l_{i,s}为真实均值下端，q_i=max_s v_{i,s}/l_{i,s}=2sqrt(r_old,i)，q_i≤q0。

z_s最小值≥l_s−v_s≥(1−q0)l_s，因此

`Lip(log z_s)·h ≤v_s/[(1−q0)l_s]≤q_i/(1−q0)`。

固定其他坐标后F_b=a_0z_0+a_1z_1，a_s≥0，故log F_b的Lipschitz常数不超过两类log核的最大值。于是

`Σ_i (Lip_i(log F_b)h_i)² ≤[4(1+r0)/(1−q0)²] Σ_i log(1+r_old,i)`。

原预算r0=10^−6、C=81/20=4.05已经覆盖该系数，所以保存的

`S_b=C log(1+Rold,b)`

仍可用。这里相同的是这个确定surrogate；不能把旧的实际Σosc_i²误当成Lipschitz cost，它可能更小。若source换成variance/12预算，仍须保留旧quarter r或调整q系数。

总和threshold分配照旧：U=Σu_b、λ>U，对t_b=λu_b/U逐项用(8)，得到

`upper≤min(1,Σ_b exp[−6log(λ/U)²/S_b])≤min(1,R exp[−6log(λ/U)²/Smax])`。

固定产品用R1、Sfixed。再逐共同L组合；要保留旧importance cap仍须与既有Markov/old upper取min。真实prior权重和非空正poly条件不变，不需重新求PB或source核。

## 6. 输入资格和数据边界

探路的72函数包含affine/reverse、clipped、asymmetric plateau、tent、oscillatory以及固定seed选出的64个rational piecewise-linear函数；3528个stoploss比较用Fraction精确成立，792个MGF数值检查通过。该有限验证不替代上面的全函数证明。

若只保留振幅，结论错误：uniform域的step函数取±1/2，或affine函数作用于两端点law，中心MGF为cosh(t/2)，严格大于uniform MGF；stoploss在0分别为1/4与1/8。它们缺少Lipschitz或连续uniform前提。之前旧MCD小例使用两点law，不能验证本次factor6；当前真实query-bin条件law是连续uniform，才能迁移。

本稿不修改冻结upper_math的factor2，也未用factor6计算任何保存tail/CI。root若采纳，可从已保存logUppermean与Smax/Sfixed做独立离线分支，无须source重建；但仍需独审、冻结后处理与统计family，携带原mean/BE/r0和普通FP未认证前提。不得把该候选直接称为新weak或C_n结论。

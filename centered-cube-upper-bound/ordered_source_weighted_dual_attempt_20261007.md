# 原 Bernoulli 来源加权对偶：有界凸缺陷的一次质量预算

2026-10-07。固定同一任意非负 L1 密度 f、W=||f||1、物理尺度和原 G_c；不改主账，不重做共同 Green/帽反例。本轮得到一般 n 的 **真实卷积边正缺陷预算**，及含 receiver r_*(x) 的明确交换余项。它不是原 weak endpoint 的闭合；正缺陷的时空预算不能直接支付选时图。

## 1. 先使用正确的 capped 对偶

令 nu_A=G_Af，p_A(r)=r^|A|(1−r)^(n−|A|)，
F_r=sum_A p_A(r)nu_A，F_*=max_r F_r。对 lambda>0，
E={F_*>lambda}，取可测最小 maximizing r_*(x)。有限 mask 系数在共同零集外均有限，F_r 为连续多项式。置
\[
 \theta(x)={\bf1}_E(x)\lambda/F_*(x),\qquad
 q_A(x)=\theta(x)p_A(r_*(x)).
\]
E 外 theta=0。原对称性和正 Tonelli 给 **精确式**
\[
 \lambda|E|=\sum_A\int q_A\nu_A
           =\int f(y)\sum_A G_Aq_A(y)\,dy.       \tag{1}
\]
这收的是同一 f；没有逐 y 界，也没有把 r_*(x) 换成 r_*(y)。

若去掉 theta 的 lambda/F_*，质量加权平均变为 ∫_E F_*，是 strong tail，不是所需 weak 交通。本文保留 capped 式(1)，不把该更强量作为待证目标。E 为空时(1)为零；有限系数总和的 L1 可积性保证 E 有限测度，或者先作有限空间截断再用单调极限。

## 2. 一般 Markov 交换图上的有界凸势

以下定理不需要 f 的 product 结构或有限 entropy，也不需要 kernel 具有 Lebesgue 密度。

充分假设：在 sigma-finite measure space (X,m) 上，G_1,...,G_n 是 commuting、positive、normal conservative Markov operators：G_i1=1，保非负函数的积分，并可对非负函数作 monotone limits/Jensen（例如 doubly stochastic probability kernels）。它们是 L1 contraction，且 A 的顺序无关。原对称轴向卷积 G_c^(i) 满足这些假设；本节预算本身不要求对称，对称性仅在§3交换公式使用。

固定 lambda>0，定义
\[
 \Phi_\lambda(t)=t-\lambda\log(1+t/\lambda),\quad
 \Phi'_\lambda(t)=t/(\lambda+t),\quad
 \Phi''_\lambda(t)=\lambda/(\lambda+t)^2.
\]
因此 0≤Phi≤t、Phi 单调凸。即使 ∫f log f 为无穷，Phi(f) 仍可积且积分≤W。每个 nu_A 的质量 W，故
\[
 D_A=\int\Phi_\lambda(\nu_A)\in[0,W].
\]
对 i∉A，v=nu_(A∪i)=G_i nu_A，定义真实正 Jensen defect
\[
 e_{A,i}(x)=G_i\Phi_\lambda(\nu_A)(x)
                  -\Phi_\lambda(\nu_{A\cup i})(x)\ge0,\quad
 \delta_{A,i}=\int e_{A,i}=D_A-D_{A\cup i}\ge0. \tag{2}
\]
对不有界输入，先截 nu_A∧M；其及 Phi 单调增至原函数，normality 和 Tonelli 延续 Jensen。Phi≤nu_A 保证所有积分有限；(2)不是“无穷减无穷”。

### 一次来源定理

令 j=|A|，
\[
 \beta_A=\frac1{{n\choose j}(n-j)}.
\]
对均匀 random permutation of coordinates，从空 mask 到全部 mask，每条边 (A,i) 被走到的概率恰 beta_A。每一条路径均严格 telescopes，平均后
\[
 \boxed{\sum_{A,i\notin A}\beta_A\delta_{A,i}
       =D_\varnothing-D_{[n]}\le D_\varnothing\le W.} \tag{3}
\]
不是每条路径另领 W 后相加；(3)是同一来源的路径平均。

等价地，nonnegative beta integral
\[
 \int_0^1 r^j(1-r)^{n-j-1}dr=\beta_A
\]
和 Tonelli 给
\[
 \boxed{\int_0^1\!\int
   \sum_{A,i\notin A}r^{|A|}(1-r)^{n-|A|-1}e_{A,i}(x)\,dm(x)\,dr
       \le W.}                                 \tag{4}
\]
这里 mask 求和有限，函数非负；端点形式按积分的极限理解。所有来源相关性留在 nu_A，随机排列是算子坐标的访问顺序，不是输入 source labels 的独立性假设。

还可令 D(r)=sum_A p_A(r)D_A。求导得
−D'(r)=sum_edges r^j(1−r)^(n−j−1)delta_(A,i)≥0，
故(4)是一个真实正流，而不是 signed derivative 的绝对值预算。

## 3. 一个可付的边交换项，以及保留选时的 signed 缺陷

先令 u=nu_A、v=G_i u，h_lambda(v)=lambda/(lambda+v)=1−Phi'(v)。
逐点凸性及两者同质量给
\[
 \int h_\lambda(v)(v-u)
 \le\int[\Phi_\lambda(u)-\Phi_\lambda(v)]
 =\delta_{A,i}.                                \tag{5}
\]
因 h≤1，左边绝对可积。这是一条真实 source-weighted 边收费，按 beta_A 平均后总费≤W。

对实际 child B=A∪{i} 的 capped test q_B（§1定义），由于
F_*(x)≥p_B(r_*(x))v(x)，
\[
 q_B\le\min\{1,\lambda/v\}\le2h_\lambda(v).
\]
零 v 时用直接极限。设 t=q_B/h_lambda(v)，则 0≤t≤2。对原 G_i 对称卷积，凸性加自伴性得到
\[
 \boxed{\int q_B(v-u)
 \le \int t\,e_{A,i}
 +\lambda\int(G_it-t)\log(1+u/\lambda).}         \tag{6}
\]
证明将(5)的逐点凸性式乘 t，再写
Phi(u)−Phi(v)=[Phi(u)−G_iPhi(u)]+e；
两组交换项合成
(G_it−t)[u−Phi(u)]=lambda(G_it−t)log(1+u/lambda)。
全部项有限，因为 t bounded，lambda log(1+u/lambda)≤u。

第一项在任何这一类实际 test family 下，按 beta_A 平均总量≤2W，**确实是一次来源的正费用**。第二项是一个明确的 signed nonlocal commutator，含 E、F_* 与 r_*(x)，没有把它删掉。它的 trivial 绝对值界只是每条边≤2W，沿路径给 O(n)W；不是 polylog。驻点均值在同一个 receiver 成立，不自动给这个空间交换项正确符号。我们未证明其正部可由(3)累计支付。

一维 resolvent 的单边 superlevel exchange 确有正确符号，但这里只把它视为局部正确接口；一般 E={max_rF_r>lambda} 不是任一 v 的 superlevel。将该一维符号应用到(6)的任意 E 或 t 会丢掉 commutator。本轮不以一维情形作为主目标。

## 4. 为什么 (3) 尚未支付真正的 graph

对原完整 family，逐点 fundamental theorem 的真实 edge flow 为
\[
 F_{r_*}-f
 =\sum_{A,i\notin A}
       \left[\int_0^{r_*(x)}s^{|A|}(1-s)^{n-|A|-1}ds\right]
        [\nu_{A\cup i}(x)-\nu_A(x)].
\]
故(1)也给
\[
 \lambda|E|=\int\theta f+
       \sum_{A,i\notin A}\int h_A(x)[\nu_{A\cup i}-\nu_A],\quad
 h_A(x)=\theta(x)\int_0^{r_*(x)}s^{|A|}(1-s)^{n-|A|-1}ds
       \le\beta_A.                             \tag{7}
\]
第一项≤W，保持了 receiver 选时。

但是 h_A **不是** child test q_(A∪i)，也没有证明
h_A/h_lambda(v)≤2beta_A。只知道 h_A≤beta_A 并不能将(7)免费放进(6)的可付第一项：h_lambda(v) 在 large v 区域可能很小。相反，(6)的 q_B 有 child prior p_B(r_*) 保证 relative cap，而(7)含 incomplete beta integral，二者不能静默互换。

因此当前两项真实缺口具体为：

* receiver-dependent edge test 的 signed commutator（6）；
* 从真正(7)的 h_A 到(6)/(3)可支付 tests 的无损转换，或直接控制其 graph-weighted defect。

这是保持 r_*(x) 的精确交换缺陷。本文没有把(4)的 dr dm 预算取到选时图，也没有要求逐 y column bound；没有用 conditional weak averaging。新的 bounded entropy 至少避免了不受控的 ∫f log f 和每层 W，但它本身尚不回答这两项。

## 5. 原 G_c 两相位、相关 spike/proper-coordinate 来源守卫

预登记后新 guard 一次运行。n=8,32,128，Fejér degree d=3,5,7；phase grids 64/128/256；c=1,3/4,1/2；lambda=1/8,1/2,1,2,8。

在 full n-torus 的坐标 x_i∈[0,2pi] 上，theta=v·x、eta=w·x。第一两列 (v_i,w_i)=(1,0),(0,1)，其余按
(1,1),(1,−1),(2,1),(−1,2),(0,1),(1,0)循环。因此 Haar 的相位 pushforward 为均匀二维 torus，输入坐标可高度相关。

F_d(theta)=|sum_(j=0)^d exp(ijtheta)|²/(d+1) 非负、均值1。两份固定质量1来源：
F_d(theta)F_d(eta)，及
.4 F_d(theta)F_d(eta)+.4 F_d(theta+pi/3)F_d(eta+pi/2)+.2 F_d(theta)。
前者有相关 spike；后者为同一 mixed 来源，包含对 v_i=0 的坐标不活动的 proper-coordinate component。没有据此约束一般源复杂度。

真实 G_i 对 phase mode (k,l) 的 multiplier **恰为**
1/[1+cB((kv_i+lw_i)²)]。来源有限 Fourier，没有 source 频谱 truncation；zero mode恰1。没有 reset 或代理生成元。natural/reverse/两个 seeded permutations 均从同一来源生成全部 ν_A，未对每条路径另收来源质量。

三轮记录72个 path NPZ；168个 Fraction beta identity checks，27108个 numerical predicates，终态 PASS_EXACT_AND_NUMERICAL_SCREEN，实际10.00秒exit0，session57197已结束。min integrated edge defect −2.23e−16、selected pointwise Jensen −2.28e−15，均在注册1e−7 screen 内；telescoping discrepancy≤3.34e−16。

|n|64→128/128→256的最大差|仅128→256最大差|最大path fee显示|
|---|---:|---:|---:|
|8|6.081e−7|4.869e−11|0.31531|
|32|2.790e−5|3.569e−8|0.42825|
|128|2.128e−4|7.996e−6|0.49716|

非线性 Phi 的 trapezoid/FFT Gi action 仍是浮点，不是 interval；特别 n128 的 grid 差真实保留，没有用1e−7的 positivity screen 把它宣称为积分误差上界。点态 Jensen 的连续真值由正 Markov 理论负责，FFT 的 nonlinear aliasing 仅是守卫。这些是原核模型上 (2)(3) 的针对检验，不覆盖全 n维任意输入，不认证 actual FIRST/CPGP/hardband 或 Rn 截断源。数值不拟合次数，也不验证(6) remainder或(7) graph payment。

## 6. 终态

一般 Markov theorem (3)(4) 与原对称边交换 (5)(6) 已证；它们都收同一来源 W，允许任意非负 L1、奇异轴向 jump measures。原 kernel 三轮 guard终态，脚本/注册及72个NPZ SHA256保存。旧 stopping 稿 §3严格号已按 root要求改成≤，未改旧数值或主账。

原 A_ord 的 polylog 弱端点仍未闭合。新的正费用是 graph 转移前的真实 convex edge 预算；两处 receiver-dependent 缺陷明确保留，不能把“≤W”写作原输出已经付清。

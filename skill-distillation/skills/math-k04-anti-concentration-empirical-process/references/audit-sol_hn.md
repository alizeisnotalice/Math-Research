# K04 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-7c5071d32603e0f2 Thm2 lines169–202允许可分中心高斯绝对上确界 Z finitea.s.且可退化，Qbounds依VarZ；非中心上确界扩展另需边缘方差的下确界>0。 CCK P-fbc49cdca0a78f82 Thm3 lines350–406与densityproof633–798：finite centered Gaussian，逐分量σj>0但covariance可奇异，a_p=Emax(Xj/σj)，Lévy半宽ε。Thm2 lines180–223用逐元素covariance errorΔ比较CDF，Y正方差而X可退化，打印a_p记号variance/SD不一致须明确定义标准化b_p。 Tu–Boczar P-6d836a3f30fbd5fb Thm1.2 lines19–47、Lemmas2.1–2.4 lines53–208：标准Gaussian输入全局非负且非零的二次多项式，PSD增广Q≠0、Ef=trQ>0，P(f≤εEf)≤sqrt(2eε)。它是相对均值的原点小球，不是任意偏移区间界。 一维logconcave P-9dec1113f304933f：continuous moment/Orlicz Thm3.8/Cor3.9 lines411–486、discrete variance/fourth Thm4.5/4.6 lines558–632；Lebesgue density M 与整数pmf M不同单位。 Cohen–Conze P-a150ab8a0390f583 Thm1.9 lines567–665：空间iid且F支持[0,1]，确定采样序列占用N_n、M_n=maxN_n、V_n=ΣN_n²，严格M_n²/V_n→0；Y_n归一化√V_n，弱极限W°∘F。 CCK P-52fde473d714f0b2 Thm2.1 lines257–331、proof1325–1390：iid、centered、pointwise measurable、包络Lq且q≥3、紧pre-Gaussian、n≥3、κ³≥E||E_n|f|³||；每个ε∈(0,1],γ∈(0,1)给标量sup耦合。

## 证明骨架与已检查推导

先定statistic/version，再用Gaussian浓度function bound；经验过程近似及熵另行控制后才合分布函数误差。 对max|X|用(X,−X)双倍向量保ε；光滑indicator+softmax interpolation再以anti-concentration恢复CDF。p17远尾不能用EmaxX≤σmax a_p，改用{maxX≥a}⊂{max(Xj/σj)≥a/σmax}与标准化Borell。 q11≤trQ22时用Laplace/trace branch及Ef≤2trQ22；q11>trQ22时直接Schur补f≥(sqrtq11+q12ᵀX/sqrtq11)²，分q12=0确定情形与正方差scalarGaussian。该已检查修补绕开原γunion严格<与≤的端点问题。 单交叉用扣去常数，双交叉用扣去割线并以质量/均值矩消除affine积分；连续固定M由uniform/centeredExp比较得M²Var范围[1/12,1]，离散asym度量Laplace给M²Var+M≤1。 将重复值合成独立空间点上的确定权重，max权重/√V_n→0给Lindeberg，协方差是F(s∧t)−F(s)F(t)；tightness与generalized inverse引用Billingsley，不能升级为定量耦合。 Δ_n含φ_n+γ^{−1/q}ε||F||2+n^{−1/2}γ^{−1/q}||M||q+n^{−1/2}γ^{−2/q}||M||2+n^{−1/4}γ^{−1/2}(E||G_n||_{F·F})^{1/2}H_n^{1/2}+n^{−1/6}γ^{−1/3}κH_n^{2/3}。截断δ_n只在失败概率γ(1+δ_n)+C logn/n，视觉p7核阈值只有K(q)Δ_n。网格最大值与两侧局部振荡分开再标量Strassen；有耦合(a,η)时直接事件包含给CDF夹逼，并加半宽a的Gaussian反集中。

本次实际检查：Z=|G|,ε.1，Var1−2/π=.36338，Q=.0796557，lower.04783upper.57400；σ→0Q→1。Var0Q1，原lower=1、upper=12(再截1)。 CCK例X=(G,G)：a_p0，L(maxX,.1)=.0796557≤.4。绝对max的双倍向量a′=sqrt(2/π)，L(|G|,.1)=2Φ(.2)−1=.1585194≤4·.1(1+a′)=.719154；这里区间宽.2，不能抄另一来源Q宽.1数值。X=(G,2G)给EmaxX=.398942而标准化a_p0，反证该中间式。 ε.01的f=X²，真概率erf(sqrt(.01/2))=.0796557≤sqrt(.02e)=.233164；Q0反例概率1，须排零。正constant f1在ε<1事件空，说明trQ22=0不能用除trace优化。 实算continuous uniform[-.5,.5]的M1,Var1/12与Exp1的M1,Var1；geo度量 pmf2^{−k−1}的M.5,Var2，使M²Var+M1。Section5 n2,f=x⁴−x²,points±.5,P=−3/16，在0 remainder product−3/64；n1右Unif[1,2]与左Unif[-2,-1]满足printed sign而均值3/2>−3/2。另实际读取lines603–632，从几何幂和原始矩逐系数核中心四次矩有理式和导数；cubic精确因子为((1+M)q+M−1)[(M²+12M+18)q²+(22M²−36)q+M²−12M+18]，二次判别式96M²(5M²−12)<0，0<M<1且q<(1−M)/(1+M)时导数严格负；M1点质量单列。脚本和系数证书在SOL_HN/check-logconcave-fourth-algebra.py与logconcave-fourth-algebra.json。 实际枚举序列(0,0,1,2)、三个独立Bernoulli(1/2)指标：N=(2,1,1)、M2、V6，原中心和方差3/2，除√V后方差1/4，错误除√n后为3/8。所有采样同一点时M_n²/V_n=1，Y_n仍为一个中心Bernoulli，方差1/4但四次矩1/16而同方差Gaussian为3/16，不能声称FCLT。 单函数Rademacher f=X、包络1、q3：A1/A2/A3均满足，κ1、φ_n0、H_nlogn、M1且G_n(f²)=0，可具体核误差输入而不冒称K(q)=1。已知耦合Y~N01、Z=Y+.1时η0、d_K=2Φ(.05)−1≈.0398776≤L(Y,.1)≈.0796557，半宽/全宽保持区分。

## 正反例与失败处理

正例：单索引 X_u=G~N(0,1)，Z̃=|G|，ε=0.1；Var(Z̃)=1−2/π≈0.36338，精确 Q=2Φ(0.1)−1≈0.07966。Theorem 2 界给 0.1/(√12√(0.36338+0.1²/12))≈0.0478 ≤Q≤0.1√12/√(0.36338+0.1²/12)≈0.574。

条件缺失例：取单索引 X=σG 且σ→0，对固定ε>0有Q(|σG|,ε)=2Φ(ε/σ)−1→1；不能宣称存在与σ无关的C使Q≤Cε。中心绝对上确界Theorem2仍成立但退化时无小区间衰减，不能误称定理不适用。

## 中心立方体迁移推断

中心立方体迁移接口：将目标 M_c 统计量的有限网格高斯极限写为 X_Q，并证明可分性及 Var(sup_Q|X_Q|) 的上下界；按 Theorem 2 得到 Gaussian concentration-function 界后，再独立控制非高斯到高斯的近似误差与重叠 cube 依赖。本文不替 cube 指数类给出方差或近似误差。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

退化允许不等于一致的小区间衰减；立方体重叠与高斯近似不由反集中本身证明。 dimension-free表述仅通过a_p，非n一致；σmin趋0不可uniform，Gaussian-to-empirical耦合及Stein/Borell外引仍待核。Comment5漏ε为已检查修补，不认证源全部proof。 原文一般Carbery–Wright degree/logconcave结果是外引未核；本二次版本的修补不证明一般degree或一般偏移区间，也不自动覆盖退化协方差使f在支持上恒零的输入。 literalSection5函数类应A_{n−2}及奇数阶product方向需修，禁止照抄。四次矩局部代数现已检查，比较归约与等号分类仍引用源陈述，不能借局部因子认证全族比较或全篇。无对数凹性的双区间混合能使M²Var>1。 该文仅定性GC/FCLT，不供应sup的反集中或Gaussian近似速率；序列/field依赖或random sampling另需合法条件化，不能将空间iid静默改成stationary依赖。 δ_n尾截断、κ和乘积类波动不得删掉；scalar supremum coupling不提供过程sup-norm coupling。Strassen、Borell、Rademacher矩/熵等外引未独立核完；目标cube函数类仍需证明这些具体输入。

版本决定非退化前提：Theorem 2 的中心绝对上确界允许零方差，此时 Var(Z)=0、Q(Z,ε)=1，下界为1而上界为12（可再截为1）；其界不会随ε趋0。非中心sup版本的正方差下界不可省略，经验过程高斯近似也须另证。 CCK p17中EmaxX≤σ_max a_p的中间式不能普遍使用，例如X=(G,2G)的a_p=0而EmaxX=1/√(2π)>0；可按原文证明目标用标准化max事件和Borell另行修补，不能据此否定主定理。Borell/Stein/非高斯CLT外部证明未独立核完。 对数凹来源的核心一维界不能继承到任意Gaussian上确界或cube统计量；必须先证明目标分布的对数凹性。Section5已检查两个字面反例，不能认证该节全部推广。离散原子不满足连续密度界。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。

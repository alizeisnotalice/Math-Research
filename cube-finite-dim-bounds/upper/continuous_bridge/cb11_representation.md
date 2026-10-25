# cb11：n3 body+face Cauchy 的完整96锥/Green接口

2026-10-25。固定真实生成元

    L_tau=−Σ_body |v·D|−tau Σ_face |w·D|，tau≥0，
    V_body={(1,±1,±1)}, V_face={ei±ej:i<j}。

方向未单位化，所有Cauchy项共享时间t，face项的真实尺度为tau。Lebesgue及Fourier inverse (2pi)^−3不变。本分支没有做成本、g或J搜索，四维不在本接口内。

## 1. 全空间96锥的准确覆盖

令a≥b≥c≥0是频率绝对值的排序。body符号为max(4a,2(a+b+c))；face的六项是

    Σ_(i<j)(|xi+ xj|+|xi−xj|)=2(2a+b)，

所以psi_tau=max(4a,2(a+b+c))+2tau(2a+b)。设A=(1,0,0)、B=(1,1,0)、C=(2,1,1)、D=(1,1,1)。排序正chamber拆为ABC和BCD：

    R_ABC^−1(a,b,c)=(a−b−c,b−c,c)，
    R_BCD^−1(a,b,c)=(b−c,a−b,b+c−a)。

两组系数非负分别恰为a≥b+c和a≤b+c；它们覆盖整个排序chamber，公共BC仅为边界。两个det绝对值均1。48个signed permutation chambers覆盖R³且内部不重叠，故得到96个simplicial cones，边界重复只在Lebesgue零测集。准确列表、det和每ray的body/face rate整数系数在cb11_cones.json；这不是仅凭体积或画图断言覆盖。

四ray rates为

    kA=4+4tau，kB=4+6tau，kC=8+10tau，kD=6+6tau。

每cone上psi_tau(Σu_jr_j)=Σk_j u_j，u_j≥0。它们对tau≥0均正。所有96锥由SymPy准确核det、两个inverse及谱体积；没有预先相信design草案。

## 2. 真半群、复锥密度与精确大时间尾

取独立标准Cauchy变量，令X_t=Σ_body v C_t^v+Σ_face tau w C_t^w，其特征函数正是exp(−tpsi_tau)。body的三条basis已张成R³，故与其余概率测度卷积给严格正、连续、保守的Lebesgue密度p_t。tau=0就是body-only，face退为零向量。

对96锥分别积分Fourier inverse，得到真实全空间密度

    p_t(x)=1/(8pi³)Σ_cones |detR| Π_j(t k_j−i x·r_j)^−1.       (1)

各项可能有负实部；完整复项保留。相反锥恰共轭，数值实现仅以2Re合并成48对，并未删项或取abs。该配对在zero/repeated knots处也精确。

谱单位球体积及峰值准确为

    Vol{psi_tau≤1}=5/[6(1+tau)(2+3tau)(4+5tau)]，
    p_t(0)=5/[8pi³(1+tau)(2+3tau)(4+5tau)t³]。

第一式由48个ABC tetrahedra及48个BCD tetrahedra求和det/(6prodrates)，第二式由齐次layer-cake或直接(1)得到。Fourier三角界和概率正性给0<p_t(x)≤p_t(0)。因此对任何T>0，

    0≤∫_T^∞p_t(x)dt≤5/[16pi³(1+tau)(2+3tau)(4+5tau)T²].    (2)

没有遗漏time tail。tau0给5/(64pi³t³)，与旧body-only正常数一致。本接口不重复主对face奇性邻域的角向可积幂/cutout预算。

## 3. 从截断时间积分取真实Green，不把zero knot直接塞log

每锥定义b_j=(x·r_j)/k_j，W=|detR|/prodk_j。对epsilon>0，所有cone time积分都有限，部分分式给

    ∫_epsilon^∞Π_j(t−ib_j)^−1dt
      =DD2(f_epsilon;b1,b2,b3)， f_epsilon(b)=log(epsilon−ib)。

同值节点按Hermite导数解释；f_epsilon对所有实b光滑，principal log在右半平面没有branch歧义。无穷端项在DD中准确抵消。因此

    G_epsilon(x)=∫_epsilon^∞p_t(x)dt
      =(8pi³)^−1Σ W Re DD2(f_epsilon;b)。                 (3)

由于真实p_t≥0，epsilon↓0时G_epsilon单调趋向真实扩展值G，不允许逐锥丢掉发散项后猜极限。

令L_b=log|b|。对非零x，每锥的三条ray独立，不能三个b同时0。各锥的实部具有准确形式

    Re DD2(f_epsilon;b)=gamma log epsilon + finite + o(1)。

完整规则如下（不含本锥权W）：

|nodes|gamma|finite|
|---|---|---|
|a,b,c全非零互异|0|Σ L_bj/Π_(k≠j)(bj−bk)|
|a,a,b，均非零、a≠b|0|[L_b−L_a−(b−a)/a]/(b−a)²|
|a,a,a，a≠0|0|−1/(2a²)|
|0,b,c，bc≠0、b≠c|1/(bc)|L_b/[b(b−c)]+L_c/[c(c−b)]|
|0,b,b，b≠0|1/b²|(1−L_b)/b²|
|0,0,b，b≠0|−1/b²|L_b/b²|

例如f′_epsilon(0)=−i/epsilon的发散是纯虚，恰在反锥共轭配对中消去；不能把其复项当作单独可积项。表内gamma及finite在配对后完整保留。cb11_limits_check.py实际SymPy核zero/repeated形式的极限和系数和；这是形式身份验证，没有新增空间样点。

将全部权重准确相加，得到

    G_epsilon(x)=(8pi³)^−1[Gamma(x,tau)log epsilon
                              +C(x,tau)+Σ_q A_q(x,tau)log q+o(1)]。

q为非零knots的绝对值；只在log参数处合并±q，之前分母/导数的符号全部保留。若Gamma<0，真实G=+infinity；若Gamma=0，真实G就是有限log表达。Gamma>0会迫使G_epsilon趋负无穷，与p_t≥0矛盾，因此精确表达不可能出现；代码将其当错误报警，不取abs或截为0。x=0则p_t(0)~t^−3，直接G=+infinity。

这是zero/repeated knots的完整处理，不要求每锥积分在0单独有限，也不要求预先猜哪条physical line奇异。每非零点的实部渐近只有log项与常数，因此Gamma=0时不存在漏掉的更慢无穷项。全局Gamma及finite使用Fraction准确聚合，避免用“小epsilon数值”判有限或无穷。

## 4. 可调用数值接口与稳定性范围

模块cb11_frame_green.py，SHA256

    ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849

接口：

    green(x,tau)                     # auto，float输出
    green(x,tau,method='mp',dps=80)  # 固定高精度
    green_expression(x,tau)         # exact Gamma、constant、log系数
    green_face(y,z,tau)
    green_array(X,tau)
    heat_kernel(t,x,tau)
    peak_coefficient(tau), green_time_tail_upper(T,tau)

内部先max|x|齐次归一，最后付scale^−2。rates、nodes及所有聚合系数为准确Fraction；float输入解释为其十进制拼写，若要准确tau=1/3应传Fraction(1,3)或字符串"1/3"。Gamma准确判号；有限时ΣA_q=0，可以将所有log替换为log(q/q_ref)减少共同正常数的数值取消。

auto先作fsum和取消条件检查，条件差时只在固定80位精度评价同一完整log表达；未用正负项截断。method='float'是原始双精度调试路径，极端取消时不能当严格值；mp也只是高精度数值而非interval证书。近奇性、coalescence的连续box评价尚需另行严格接口，本模块的pointwise exact分类不等于完整face g/J包围。

## 5. 预登记36参数点实际校准

固定tau=1/10,1/3,1,3，各同一8个有理点，tau0只前4点校准，共36参数点；点坐标与全部结果在cb11_results.json。第二点原登记3/8在tau3造成A/B knots相等，数值前纯有理预检发现，已透明勘误为3/11；没有对旧点作G/p或成本实验。全部36个最终点满足无zero/重复knots的预定非退化条件。

实际一次数值批：双精度解析值、完整96复锥heat与48对压缩heat、独立直接time求积；随后仅一组80位重算。同tau0四点再对旧cb9正有理积分校准，不重跑旧证书。结果：

- float解析G对80位最大相对差1.02×10^−15；
- 直接time积分对80位最大相对差2.27×10^−14；
- 无time quadrature warnings，tau0四点校准通过；
- 全批约0.14秒；没有成本/min/g/J搜索。

这些固定非退化点不是全域认证，zero/coalescence的合法性来自第3节准确极限规则而非偷偷补空间网格。SymPy/Fraction/mpmath/Scipy及定义共享如实披露，自检不代独审。主已实际串行运行便携cb11_matlab.m，源码SHA `a9d9fecf349c66498e35bb2ee4bfed56cd02d71c5ed9290d200abdb69a10d86a` 未变，exit0。96个det1锥、同36个参数点及peak校准PASS；Green最大相对差1.78×10^−13、heat最大相对差1.16×10^−15。实际执行来源和stdout/hash见cb11_matlab_receipt.json/cb11_matlab_stdout.txt。只是共享公式的浮点cross-check，未冒作全域interval或独立完整公式审计，没有重跑本批数学实验。

## 6. 证明边界与下一接口建议

严格解析结论：96锥覆盖、真实p及G截断积分、所有重复/zero knots的完整合并规则、准确谱体积和time大尾。代码给pointwise真实G的数值接口及精确分类，未给数值舍入的严格区间。

来源hash见cb11_sources.json，包含主登记与不可变输入；设计草案仅起猜想作用。固定Lebesgue用于Fourier inverse、概率密度/卷积、锥Jacobians及time积分；原026一般L1/全共同尺度仍是比较框架依赖，不从固定样点重新证明。四维、主成本pilot、全cube minimum与J认证均不由本文件自动推出。

为后续连续interval evaluator，不能把generic DD有理系数直接在跨coalescence的box上除小差。可采用epsilon>0的真G_epsilon：每锥用固定下限time积分，分母epsilon−ib实部>0，zero/repeated knot不再有奇点；再用已证off-jump-line小time界支付G−G_epsilon，近真实jump-line帽另用主cutout预算。这一接口可绕过无害的coalescence/zero planes，但还需验证其区间收敛/box误差，当前没有试跑或声称可完成J界。另一选择是按zero/coalescent分层作confluent log算法，跨层需Taylor/积分余项，不能仅generic grid。

本首批已按有限预算停止；任何后续J连续认证须另登记。不宣称新U3/U4或固定frame族障碍。

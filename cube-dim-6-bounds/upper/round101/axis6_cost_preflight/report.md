# 101：六轴真实 Green 正部族的严格成本预检

## 结论

**同一等权 axes6 模型的理想 Green 空间正部截断成本严格满足66<U<12000，
所以这一族不能改进64。** 下端来自真实正面平均值与 Jensen，不是一个粗
上预算大于门槛。因此即使把该模型的真实 g、真实 I 精确计算出来，也无法
用其最优空间截断 (G-c)_+ 改进当前上界。

精确对象是六条单位坐标轴 I6、共同 Cauchy 时间、alpha1，
p_t(x)=prod_i t/[pi(t^2+x_i^2)]。G=int_0^infinity p_tdt，q=6/5，
I=int_[0,1]^5 G(1,y)^qdy，B=I^(5/6)。新证明全 cube 最小值在角点，
已有角点 Beta 公式由此成为准确全域 floor：

    g=3/(512pi^5),   U=(6B/g)^(6/5)=6^(6/5)I/g^(6/5),
    1/100000<I<1/300,   66<U<12000。

I 和 U 的区间有意很宽；严格下端已越过64，按预登记立即停止，不再扩大
窗口、积分、参数或精度。本轮只做上界模型预检，没有任何 C6 下界推进。
它不是 C6>=66，也不排除该半群衍生的任意其他核、修正、不同轴权或时间
截断方案；结论针对准确写明的等权 Green 空间正部族。全局上界仍为64。

## 严格去重和真正新增的对象

已读084 registry/schema、006/008/067/072 指定轮次与模型资料，并补查14个
六维/通用源目录；再以 Markdown、Python、JSON 搜索补查关键 n6/theory/
numerics/scratch 根。lookup_scan.json、lookup_all_formats.json、LOOKUP.md
记录范围、命中及 hash。没有因084漏列一个路径就把模型认作新颖。

补查确实发现一个容易漏掉的旧件：

|旧件|SHA-256与准确状态|
|---|---|
|theory clone 的 axis_frame_cost_probe.py|f845cf942707c27412c9d6ad04eb59f89fe3db1cd95f956cf4b26b576f3b4f28|
|相邻 axis_frame_cost_probe.json|27aadd616636150c761711a570affbaf9a7a30054aee189b93fef101cc256d51|

它是同一 I6/unit-Cauchy 模型，已有角点3/(512pi^5)，65536个 Sobol 点、
96个时间 Gauss 节点，形式成本750.7412233。旧文件明确标记 numerical_candidate_only，
并声明未证明角点是全 cube 最小、积分估计没有认证误差。这个历史值只作
去重证据，未参与任何本轮不等式，也未重新运行旧点或网格。

072 已证明同一轴模型的局部半凹及一条真正 product 路径；052/069 是 tiny
真实核参考，不是完整成本排除。006/008/067 的 alpha3/2 结果属于 body32，
012/018/081 是 body32+axes，通用 finite-dim round009 的 axis anchor 是
三维 body-diagonal 点值，均不能替代本项成本。U08 提供 n3/4 的坐标 Cauchy
接口参考，没有六维轴模型成本证书。

在上述明确范围内没有找到同一模型的完整严格实际成本排除，所以没有触发
“提取旧证书、0science停止”的条件。本次新增的是旧数值 pilot 缺失的解析
严格成本证书；不是新 frame，也不主张范围外或公开文献的穷尽新颖性。

## 真实模型、角点 floor 与全空间接口

完整推导在 theorem.md。每个 c_t(s) 随 |s| 单调减，故在整个[-1,1]^6
p_t(x)>=p_t(ones)，正积分给 G(x)>=G(ones)，包括扩展值奇异点。
角点附近的小时间 O(t^6)、大时间 O(t^-6) 支配给连续性，cube 内正体积邻域
可逼近角点，所以 essential infimum 也正好等于该角点值，而非只证边界 inf。
Beta(7/2,5/2)/2=3pi/512 给上述准确 g。角点常数是旧公式的来源核对和
解析解释，没有作新的旧参数点 kernel 评价。

坐标生成元 L=-sum|D_i| 是真实对称保守半群。psi 与欧氏范数满足
|xi|<=psi<=sqrt(6)|xi|，L2 域 H1；热核上界为 pi^-6 t^-6，不能混入
body32 的1/(160pi^6)。瞬时 t^-6、维数6>1、同时间保正/守恒，满足001
完整有限 spanning-Cauchy 障碍证明所需的每个维数假设。

G 局部 L1、R_c=(G-c)_+全局 L1；后者也由全 face I 有限与 cube-radial
计算直接得到。正部源直接在全空间证明，未删除支撑<=2的真实奇异集合：
G=P_tG+U_t、正部凸性给 P_tR_c>=R_c-U_t，从测试极限得到
mu_c=LR_c+delta0>=0。大 cutoff 用 R^-1 缩放和 R_c 的 L1 回收总质量1；
写 min(G,c) 并用 epsilon^5 的小 cutoff 缩放，去除原点正原子。

一般 L1、全部半径接口重新核查了同一个障碍坏集 D：有界紧支撑 f 的
f=sigma-Lu 分解与正概率源卷积，在同一个 D 补集控制所有 K_(c,r)。
先合并正有理尺度的零集，再由 L1 dilation 连续性扩到全部正实尺度，
最后单调截断 |f| 扩到任意 L1。K_(c,r)>= (2r)^-6 1_[-r,r]^6，r 是六轴
共用的单一半径。product semigroup 不等于六个独立尺度的 strong maximal
operator，也没有使用一维最大算子的 composition。

符号/面归一化亦精确匹配：六个最大坐标 faces、64个反射给 J=384I，
radial Jacobian r^5，正积分得 int R_c=320c^-1/5 I。
K_c=R_c/[64(g-c)] 的质量为5I c^-1/5/(g-c)，唯一最小点 c=g/6，最小值
恰为U。这说明真实成本下端排除的是所有 c，不是只试一个截断值。

## 决定性的实际成本下端

因正面体积1，B>=A=int_[0,1]^5G(1,y)dy。真实 product 与 Tonelli 给

    A=int_0^infinity t/[pi(1+t^2)] * [atan(1/t)/pi]^5 dt。

在固定 [1/4,1/2] 上 atan(1/t)>=atan2>pi/3；窗口外贡献非负，所以

    A>log(20/17)/(486pi)。

log(20/17)>93/578 来自 log(1+x)>x-x^2/2、x=3/17，方向严格。
复用只读E098已证有理 pi 区间，取粗下端 pi>157/50，得到

    6B/g > R_L=(512/243)(157/50)^4(93/578)。

本机 Fraction 实际核验 R_L^6>66^5，故 U>66>64。这下包的是实际 A、B、I
及理想模型成本；没有从真实核的上包函数反推实际积分下界。其证明不使用
旧750.74数值，也不依赖 body32 尚未闭合的全域 g 或098/100局部补项。

## 完整全时间上端与为何不再细算

为给出完整成本区间，也必须支付全 face、全部时间与所有奇异邻域。
对0<t<=1，用真实半线 Cauchy 质量1/2与峰值插值，得

    ||p_t(1,.)||q <=2^-25/6 pi^-11/6 t^1/6。

对 t>=1，pointwise p_t<=pi^-6t^-6，face体积1。空间 Minkowski 与两段
准确时间积分于是给

    B <= (6/7)2^-25/6 pi^-11/6+1/(5pi^6)
      <3/448+1/3645<1/128。

这里既未漏0时端，也未漏无穷尾。pi>157/50 的十一幂大于8^6、2^-25/6<1/16，
各 fractional-power 比较均转为正有理数整数幂检查。pi<22/7 与准确 g 给
6B/g<8(22/7)^5；实际核验其六次幂小于12000^5，故 U<12000。
同样有理幂检查给整个真实 I 的粗区间(1/100000,1/300)。它证明成本有限，
不是竞争性上预算；真实成本下端已失败，不应为缩小宽上端而盲目扩算。

## 真实运行与入口修复账目

去重读取在新科学前完成，registration.json 先登记3项/上限4项，CPU总上限
30秒，每程序含保存/receipt/stdout/exit/cleanup墙钟上限180秒。

首个启动确有一个可避免的入口错误，原件全部保留：guard 在 Popen 前把
CPU hard 设为2，child继承后不能设为28，故在数学 try 之前的 setrlimit
抛出 ValueError。该次数学配置0、断言0；不是超时或数学反例。child已回收，
guard内部到最终receipt前墙钟0.043824秒，177秒SIG_DFL hard fuse当时仍在。
原工具 raw wait 字段约1e-5秒不是整个程序计时；该次未记录精确完整tool
roundtrip，不伪称它等于guard内部计时。原 traceback、run_receipt.json、
tool_execution_receipt.json、执行绑定及两份源代码字节冻结。

主线程明确授权一次另名修复，并在执行前全文静读 run_guard_v2.py 与原 child，
核v2 hash及限定diff、22源pins和9执行pins后正式放行。修复只涉及CPU继承：
spawn前guard soft2/hard30，spawn后parent降hard2；child原数学代码和参数未改，
新收据独占run_receipt_v2.json。177/150秒默认硬fuse、160秒wait、5秒own-group
kill/reap策略均不变，直到最后receipt/stdout/os._exit都没有cancel。

放行的唯一 continuation 一次通过：3/4配置、77断言，数学CPU0.002046秒；
child OS CPU0.026661秒、guard CPU0.029646秒；guard内部到最终receipt写前
墙钟0.032274秒。Date.now在完整await exec前后测得0.091秒tool roundtrip，
包含guard启动、数学child、guard最终receipt/stdout/exit/reap及工具传输开销。
该完整观察值与guard内部值分字段保存，不能把两种scope混作同一时间。
CPU限额child28+guard2<=30，hard fuse直到退出保持生效；无实际超时。

所有77项断言通过后即停止数学。22个只读源在数学前后hash一致，未运行
旧corner kernel、旧Sobol/time grid、任何body native、MATLAB或J100，也未
推进075/099或其他下界工作。历史axis数值仅提取原值与hash。未递归分派、
stage、commit或push，所有写入仅在本人新101目录。

本地实际工具是标准库 Fraction，精确整数幂与已有有理pi证书；数值决定
不依赖浮点、MC、QMC或未认证求积。完整解析证明、最小运行代码、新旧
入口/修复diff、全部执行绑定和真实收据均列在manifest。首失败后的修复是
主明确放行的一次独立continuation，不是自动重试；其后无失败或再运行。

本轮为这个实际上界构造族的严格失败证书，不改变原C6<=64。按登记封存后
停止；最终理论接受、账本和是否归档推送均由主线程处理。

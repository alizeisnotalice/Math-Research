# 原完整 Bernoulli family：非产品、多尺度输入的弱水平集探路

2026-10-07，tensor 代理。本轮只研究固定物理尺度、原完整
\[
 T_r^{c}=[(1-r)I+rG_c]^{\otimes n},\qquad0\le r\le1,\quad
 \widehat G_c(\xi)=\frac1{1+cB(\xi^2)},\quad
 B(v)=\frac v{\log(1+v)}-1.                    \tag{1}
\]
本轮不重复强 D_n、first-event 强下界、L2 或投影压缩。采用 L02 的显式可行输入/有限搜索/回建及误差方向工作流；技能不提供弱型结论。物理尺度归一到 1。目标是判断
\[
 \|\sup_rT_r^cf\|_{1,\infty}\le A_{\rm ord}(n)\|f\|_1 \tag{2}
\]
的 polylog 版本是否被原核输入直接排除，而不是拿 reset 代理或一维核替代 (1)。

**结果：本轮没有排除 polylog，也没有证明 (2)。** 三轮 n=8,32,128，c=1,1/2 的四类非产品/多尺度输入，原核实测弱比值约 0.78–1.03。有限输入的高置信筛选上界最高约 1.61，但浮点核实现不是区间算术；更不能把这些上界当任意输入、任意维数的上界。新增实质是原 kernel 的完整连续 r、空间水平集和有限 L1 回建接口均明确可审，而不是再报告一个行条件。

## 1. 完整来源输入：同一周期空间上的计数约束

取共同周期空间 [0,4)^n，以归一化 Lebesgue measure。对 P=1,2,4，置
\[
 b_{P,i}(x)={\bf1}_{x_i\bmod P\in[0,P/2)},\quad
 K_P(x)=\sum_i b_{P,i}(x),\quad
 \pi_j=2^{-n}{n\choose j}.
\]
各模型平均质量一，完整源固定后不依赖 receiver、r 或 c。n 为偶数，m=n/2，d=floor(sqrt(n)/2)。预设四个来源：

1. f0=1_{K4=m}/pi_m；
2. f1=(1_{K4=m−d}+1_{K4=m+d})/(pi_(m−d)+pi_(m+d))，等价于两个对称 shell 各自归一后等权混合；
3. f2=1_{K4∈{m−1,m,m+1}}/(pi_(m−1)+pi_m+pi_(m+1))；
4. f3=(period1 的 f2 规则 + period2 的 f1 规则 + period4 的 f0)/3。

计数约束使输入坐标相关；f3 的三种分辨率由同一源 mixture 决定，**不能先对三个分量分别取最大再当原输出**。这些是有界密度，不是原子。输入计数只是具体压力族定义；没有据此假定任意一般源有有限复杂度或独立坐标，也不按来源数量收费。

对固定 receiver，令 g_(c,P)=G_c^{per P}*1_[0,P/2)，则
\[
 p_i(r)=(1-r)b_{P,i}(x)+r\,g_{c,P}(x_i).
\]
原 kernel 对不同坐标的乘积作用给 exact conditional count polynomial
\[
 \prod_i[(1-p_i(r))+p_i(r)z].
\]
其第 j 系数是 T_r^c 1_{K_P=j}(x)。这不是输入来源独立性；我们对相关函数 1_{sum bits=j} 应用一个张量卷积。holding 的 r=0、全连续 r=1、全部中间 faces 均保留。C++ 正系数递推计算这些概率；f3 在每一个同一 r 合成三分量，然后才取最大。

## 2. 连续 r 有纯有理外包，没有未核峰值间隙

Root 已独立完成并由 gated 审读
[positive_bernstein_rational_mesh_20261007.md](positive_bernstein_rational_mesh_20261007.md)。
本轮提出的核心驻点身份是：对 F(r)=sum a_j r^j(1−r)^(n−j)，a_j≥0，在内部最大 p，代数后验满足 E j=np，因此
\[
 F(q)\ge F(p)e^{-nD(p\Vert q)}.
\]
Root 的有理节点 q_j=j²/[j²+(M−j)²]、χ²覆盖半径 4/M²（M 偶）把它变成完全有理网格证书。三轮 (n,M)=(8,128),(32,256),(128,512) 都有
\[
 G(x):=\max_jT_{q_j}f(x)\le\sup_rT_rf(x)\le(512/511)G(x). \tag{3}
\]
这适用于任意非负源，不要求我们的计数压力族。它控制参数幅度误差，**不是空间弱界**；不能从 O(sqrt n) 个节点免费推 polylog 弱费。本轮使用 root 已证结论，没有重复其 4939 项精确实验。

## 3. 原 G_c 的 Fourier 计算和显式误差

没有替代原 G_c 的正概率核。period P 的响应为
\[
 g_{c,P}(x)=\frac12+\frac2\pi
 \sum_{\substack{k\ge1\\k\ {\rm odd}}}
 \frac{\widehat G_c(2\pi k/P)}{k}\sin(2\pi kx/P). \tag{4}
\]
原 B 在所有使用的 Fourier modes 上直接用 log1p 公式计算。只截到 K=16384，FFT table 有 262144 点，并作周期线性插值。

写 a=2pi/P。对遗漏 modes v≥4，有 B(v)≥v/[2log(1+v)]；由 log(1+a²k²)≤2log k+log(1+a²) 及递减正级数积分，
\[
 \epsilon_{\rm series}\le
 \frac4{\pi c a^2K^2}
 [\log K+\tfrac12+\tfrac12\log(1+a^2)].          \tag{5}
\]
这些都是原 symbol 的高频项，不是用另一个核近似原低频谱。

当 c≥1/2，原 Exp 从属表示给 G_c≤2G_1=2w。periodic peak≤2[1+2/(e^P−1)]。半周期 indicator 卷积的导数是两个平移 periodic G 的差，故 Lipschitz 常数至多该 peak。因此 table 插值误差
\[
 \epsilon_{\rm interp}\le
 [2+4/(e^P-1)]\,P/(2\cdot262144).              \tag{6}
\]
series 与 interpolation 误差相加约为 8.26e−6、1.003e−5、1.587e−5（P=1,2,4；以 c=1/2 的较大项显示）。

若各 p_i 的误差≤epsilon，product Bernoulli probability measures 的 TV coupling 给任意 event 概率差≤n epsilon，故输入高度 H 的响应误差≤H n epsilon。对 f3 按固定三个来源权重相加，仍先合成同 r。结合 (3)，计算 profile Ghat 对真正 continuum 输出的方向为
\[
 \max(0,\widehat G-\Delta)\le T_*f
 \le (512/511)(\widehat G+\Delta).             \tag{7}
\]

**浮点限制：** FFT、DP、表格插值和 log1p 都为 double；另加 1e−10 的一维浮点 screen allowance，但它不是已证舍入外包。故 (5)(6) 是解析 truncation/discretization 界，实际保存数据的“置信上下界”还条件于浮点 screen，不能称全区间认证。没有隐藏这个区别。

## 4. 弱水平集不是平方能量或 rare-event 盲采

每轮独立 uniform receiver batches 为四份；每份样本数 n8 为 1024，n32/n128 为 512。c=1 和 c=1/2、四模型在同一 batch 共用 receiver，四 batch 彼此独立。因此不同 c 的结果不是额外独立重复。

阈值预设为 lambda=2^(j/8)，j 从 −8 到 ceil(8log2 H)+8，仅由输入高度 H 决定，不事后选择有利阈值。保存全部阈值行，同时显示它们的最大 weak ratio。对聚合 profile CDF 用 DKW，24 个 dimension/c/model 的联合失败概率≤.01：
\[
 \epsilon_{\rm DKW}
 =\sqrt{\log(2\cdot24/.01)/(2N_{\rm total})}.
\]
分别为 .03216696、.04549095、.04549095。虽然同轮模型相关，union bound 不要求它们独立；每个 profile 的 receiver samples 独立即可。

对每个阈值，lower probability 使用 Ghat−Delta、upper 使用 (512/511)(Ghat+Delta)，再加减 DKW。网格之间的完整 weak supremum upper 另乘 2^(1/8)。lambda<1/2 的 ratio≤1/2；lambda≥H 的水平集由 Markov kernel 的精确输入高度帽为空。这是选定有界输入的有限概率处理，**没有**从少量样本声称检查了任意指数高输入高度的极稀事件。

主要显示结果如下。括号为本次解析误差+统计 screen 的 lower/upper；浮点未区间认证的限制如上。

|n|模型|c=1 的 empirical weak ratio|lower/upper screen|
|---|---|---:|---:|
|8|central shell|0.9148|0.8066 / 1.1156|
|8|shifted double shell|0.9771|0.9069 / 1.1420|
|8|central band|0.9337|0.8920 / 1.0637|
|8|shared multiscale mixture|0.8563|0.8268 / 0.9687|
|32|central shell|0.8935|0.6848 / 1.3080|
|32|shifted double shell|1.0224|0.8239 / 1.3313|
|32|central band|0.9209|0.8127 / 1.1223|
|32|shared multiscale mixture|0.8686|0.8247 / 0.9976|
|128|central shell|0.8672|0.6721 / 1.6131|
|128|shifted double shell|0.9372|0.8368 / 1.5367|
|128|central band|1.0034|0.7870 / 1.3302|
|128|shared multiscale mixture|0.8763|0.8220 / 1.0000|

c=1/2 的前三类最高 ratio 相同，主要因为最高阈值由 r=0 的输入层驱动；它不是 fullfuture 对两 c 的独立支持证据。multiscale mixture 为 .7844、.8019、.8207。完整 profiles 及 winner node 已保存，future r 没有从实验中删除。

这一输入族没有找到迫使 (2) 超过 polylog 的反例，且其高阈值主要受原输入层影响，说明这轮族对 maximal propagation 的压力有限。既没有拟合阶数，也没有用这些上界宣称 A_ord 有界。任意非产品来源、多数其它周期和中间 c 仍未覆盖；c 的两个端点不能免费扩成共同 sup_c。

## 5. 回到 Rn 的有限 L1 输入：只回代下界方向

Periodic f 不是 Rn 的 L1 源，故不能直接当反例。将它截到 Q=[0,4Mbox]^n，Mbox=65536 n² 为整数，nu=f1_Q。其质量恰为 |Q|，输入高度 H≤2^n。对距 Q 边界至少 R=8n 的 receiver，任意 r∈[0,1]、c∈[1/2,1]：
\[
 0\le T_rf-T_r\nu
 \le2Hn e^{-R}\le n\,2^{1-7n}.                \tag{8}
\]
因为 holding 不移动、每个 visited 原 G_c 位移的 tail≤2e^(−R)，按坐标 union bound，仍是原 kernel。该 lifting 误差全 r 一致。删去边界 strip 的相对 Lebesgue 体积
\[
 1-(1-R/(2Mbox))^n\le nR/(2Mbox)=1/16384.      \tag{9}
\]
先在整个大盒按整数周期复制 torus level set，再扣 (9)，所以一个 torus lower weak level 可提升为
\[
 \lambda|\{T_*\nu>\lambda\}|/\|\nu\|_1
 \ge\lambda[\Pr_{\rm torus}\{T_*f>\lambda+
                 n2^{1-7n}\}-1/16384]_+.       \tag{10}
\]
没有条件化省略大盒边界，也不把 torus 局部概率直接写作有限来源质量。**本轮没有把 torus upper 回代为 Rn upper**；大盒外部 receiver 和有限截断后的全空间输出仍须另付。所选源也未核 hardwinner、hardband、原 FIRST、CP/GP、SC-F/R 或完整历史门；(2) 本身是任意初始密度接口，不能称本实验为 actual geom 余项。

## 6. 可重审收据和当前判断

先保存 bernoulli_weak_endpoint_registration_20261007.json，再运行新 guard 一次，实际 8.78 秒 exit0；session 96850 已结束。没有改旧文件、旧实验或主账。主结果 JSON 和 24 个 NPZ 保存完整 profile、winner、receiver、seed、r nodes、输入高度与文件 SHA256。

主脚本 exact_checks=902 的含义是 896 个有理相邻网格区间值加六个 max-radius/constant 等式，并非 902 个独立 Boolean asserts；另外 48 个有限性/高度诊断。Root 的4939项有理网格守卫独立负责完整标量定理。

再只读已保存数据作独立接口核验：n8 的三点、两个 c、全部129 r，用直接 Fourier series（不用 FFT 插值）和显式 2^8 个 source bit events 枚举，复核 DP/来源 mixture/interface。2348 项 PASS_SAVED_INTERFACE_REVIEW；max profile 差 3.992e−11 / 4.557e−11。24 个 NPZ hash 全匹配。这不是重跑主 Monte Carlo，更不是一般弱型证明。

**终态判断：** 完整原 family 的弱端点 polylog 仍是可行但未证的策略；本轮未发现原核反例，不能判断它过强。新增连续网格和原核有限 L1 回建使该接口可直接压力，而数值所选输入尚未显示空间占用难点。Root 的 hazard/source transfer 可以合法条件于 A_ord；它与本实验都不支付原 cube 的 R_dagger。主目标继续 active。

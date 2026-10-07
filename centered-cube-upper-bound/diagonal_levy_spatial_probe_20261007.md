# 原 c=1 的斜向 Lévy 测度增量：解析充分链与三轮空间诊断

日期 2026-10-07；模型 gpt-6.1-sol high。使用 J01，已读其 method、provenance、cube-interface；原输入接口承接 [full_rank_source_probe](full_rank_source_probe_20261007.md) §5、§7。不使用 Fourier 点态序冒充正 measure 序，不重跑旧 fixed K 或两相位空批。

本轮得到可证明的原空间合同：**c=1、L 递增时，dr/dlog L>=2r(1-r) 保证原 P 的 Lévy measure 增量为正。** 等号是 odds(r)=C L²。它建立一条真正原核 Markov 斜向链，但不能免费合并所有 C，也尚不支付 actual FIRST/fullfuture 的完整 receiver/source 预算。

## 1. 原 G_a 的 cut 密度与 pole 审计

使用 Fourier 约定 hat G_a(xi)=f_a(xi²)，其中

\[
 B(z)=\frac z{\log(1+z)}-1
 =\int_0^1[(1+z)^t-1]dt,
 \quad f_a(z)=\frac1{1+aB(z)}
 =\frac{\log(1+z)}{az+(1-a)\log(1+z)},\quad0<a\le1. \tag{1}
\]

主支 cut 为 (-infinity,-1]。上半平面的 B 严格有正虚部：arg(1+z) 在 (0,pi)，每个 0<t<1 的 (1+z)^t 有正虚部，因此 1+aB 不可能为零。下半平面由共轭相同。在 (-1,0)，z/log(1+z) 在 (0,1)，于是 B 在 (-1,0)，1+aB>0（a=1 也在开区间严格正）；在非负轴更不为零。z=0 是可去点，f_a(0)=1。因此没有遗漏 (-1,0) pole 或 off-cut pole。z=-1 对 a<1 时 f_a 趋 1/(1-a)，对 a=1 时只有 log branch singularity，均不给孤立原子 pole。

在 cut 上 z=-v+i0、v>1，ell=log(v-1)，

\[
 -\pi^{-1}\operatorname{Im}f_a(-v+i0)
 =\rho_a(v):=\frac{av}{[av-(1-a)\log(v-1)]^2+\pi^2(1-a)^2}>0. \tag{2}
\]

直接 keyhole Cauchy 积分给 f_a(z)=integral_1^infinity rho_a(v)/(z+v)dv。外圆因 f_a(z)=O(log|z|/|z|) 消失，branch endpoint 小圆因至多 log 奇异消失，且已排除其它 poles。rho_a(v)=O(1/(av)) 在无穷远，积分收敛。逆 Fourier 使用 (xi²+s²)^(-1) 的密度 e^{-s|x|}/(2s)，再令 v=s²，得到原候选公式确实为

\[
 g_a(x)=\int_1^\infty
 \underbrace{\frac{a s^2}{[a s^2-(1-a)\log(s^2-1)]^2+\pi^2(1-a)^2}}_{R_a(s)}
 e^{-s|x|}\,ds. \tag{3}
\]

其质量由 f_a(0)=1 给出；非负且偶。a=1 还原 R_1(s)=s^{-2}，所以 g_1=w 为旧原密度。这里 a 是原 c 参数在 c=1 的积分内变量，不是物理窗口左端点。

## 2. 先对 a 精确积分，得到原 Lévy 的 phase 表示

对 0<r<1，原 P_r=K_r^(1/2) 的单坐标 Lévy density 为

\[
 \Lambda_r(x)=\frac12\int_{1-r}^1g_a(x)\frac{da}{a}
 =\int_1^\infty F_r(s)e^{-s|x|}ds. \tag{4}
\]

Tonelli 合法。记 ell=log(s²-1)，代换 u=(a/(1-a))s²-ell，du=s² da/(1-a)²，得

\[
 F_r(s)=\frac12\int_{1-r}^1\frac{R_a(s)}a da
 =\frac1{2\pi}\operatorname{atan2}\!\left(\pi,
 \frac{1-r}{r}s^2-\log(s^2-1)\right). \tag{5}
\]

atan2 在 (0,pi)，避免 arctan 选支；等价 F=(pi/2-arctan(u_r/pi))/(2pi)。s->1+ 和 s->infinity 都有 F->0，后者 F_r(s)~r/[2(1-r)s²]。前者约 1/(2|log(s²-1)|)，因此 integration-by-parts 下界项也为零。原 measure 总质量

\[
 |\Lambda_r|=-\frac12\log(1-r),\qquad
 P_r=\exp\{\Lambda_r-|\Lambda_r|\delta_0\}. \tag{6}
\]

式 (6) 是真实原 compound Poisson 核，不是新的 reset 或 Gaussian 代理。物理尺度 L 的 density 是 Lambda_{r,L}(x)=L^{-1}Lambda_r(x/L)。n 维原 P 是各坐标 tensor，相同单坐标 measure 结论可独立叠加，维数 n 不进入斜率条件。

## 3. 真正正空间 measure 的斜向充分条件

令 t=log L，q=dr/dt>=0。在 L=1 的 scaled coordinate x，真实 measure 生成元密度是

\[
 H_{r,q}(x)=q\partial_r\Lambda_r(x)-\Lambda_r(x)-x\Lambda_r'(x). \tag{7}
\]

它不是冻结系数的 P_r generator；第二、三项正是物理 dilation 的变化。对 x>0 用式 (4) 分部积分（上述两端项为零），

\[
 H_{r,q}(x)=\int_1^\infty
 [q\partial_rF_r(s)+s\partial_sF_r(s)]e^{-s x}ds. \tag{8}
\]

令 u_r(s)=((1-r)/r)s²-log(s²-1)，直接求导为

\[
 q\partial_rF_r+s\partial_sF_r
 =\frac{s^2}{2(u_r^2+\pi^2)}
 \left[\frac q{r^2}-\frac{2(1-r)}r+\frac2{s^2-1}\right]. \tag{9}
\]

因此 q>=2r(1-r) 时括号严格正，H>=0。这是**正谱混合系数推出正空间 measure**，不是 Fourier multiplier 的点态比较。q=2r(1-r) 时剩下 s²/[(u_r²+pi²)(s²-1)]，在 s=1 附近约 1/[(s-1)log²(s-1)] 可积；在无穷远为 O(s^{-4}) 可积。式 (7) 在 x=0 也由极限成立。q 更大时增加正的 partial_rLambda。

特别 root 提议 q=2(1-r) 也充分，但不是本轮得到的较小充分斜率；本轮没有证明 2r(1-r) 是空间 measure 层面的必要或最优斜率。次临界斜率的负谱系数本身不证明空间 density 有负部。

另有更直接的有限步证明，完全不依赖数值或局部 differentiation。令 odds(r)=r/(1-r)=C L²，则

\[
 \Lambda_{r,L}(x)=\int_{v>1/L}F_r(Lv)e^{-v|x|}dv,\qquad
 u_r(Lv)=v^2/C-\log(L^2v^2-1). \tag{10}
\]

L 增大时已有 support 上的 u 严格下降，atan2(pi,u) 严格增加；谱 support (1/L,infinity) 扩张。因此任意两个有限尺度 L2>L1 有 Lambda_{r2,L2}>=Lambda_{r1,L1} 作为**原空间 density**，不仅在 Fourier 意义上。一般只要 L2>=L1 且 odds(r2)/L2²>=odds(r1)/L1²，同样逐点谱 domination 成立。

等号链在 L 翻倍时

\[
 r_2=\frac{4r_1}{1+3r_1}. \tag{11}
\]

它保持 0<r<1。root 的较大 q=2(1-r) 链则 r2=1-(1-r1)/4。有限步原概率商是 compound Poisson with jump measure Lambda_{r2,L2}-Lambda_{r1,L1}；它们的商与 cocycle 由 exponent 差分直接给出。无需把同 r 已否定的比例商重跑。

## 4. 真生成元、质量与适用范围

在紧 r 区间内、L 有界的 C1 斜向路径满足上述充分条件，则物理 density

\[
 H_t(z)=L(t)^{-1}H_{r(t),q(t)}(z/L(t))
\]

非负、偶、有限质量，质量为 q/[2(1-r)]。它有有限一阶矩：谱率最低 1/L，空间尾指数衰减；原 Ga 及 compact r 的 mixture 在原点有限。对 bounded C1 测试函数可用

\[
 (\mathcal A_t\varphi)(x)=\sum_{i=1}^n
 \int_{\mathbb R}[\varphi(x+z e_i)-\varphi(x)]H_t(z)dz. \tag{12}
\]

无 killing、无 boundary cemetery；这是时间非齐次、空间平移齐次的原增量过程。A_t1=0，偶性给合法坐标测试的 drift 为零；单坐标总 rate=q/[2(1-r)]，logistic 等号链为 r，n 坐标总 rate=nr。有限区间积分有限，不爆炸；演化核是原 P 参数之间的 compound-Poisson 商。r=1 的无限总 Lévy 质量端点不由本有限强度构造处理。

logistic 等号链的 normalized 单次 jump 是 rates 至少 1/L 的 Laplace probability mixture：式 (9) 的正系数 Gamma(s) 将 H(x) 写成 integral Gamma(s)e^{-s|x|}ds，归一化时混合权为 2Gamma(s)/(rate*s)。因此任意整数 k>=1 有 E|Z|^k<=k!L^k。原 P 的小频率系数由 B(v)=v/2+O(v²) 给出，marginal variance=rL²/2；沿 t=logL 求导为 r(2-r)L²。除以单坐标 rate r，精确 normalized jump variance 为 (2-r)L²。这些是原位置 jump 的矩，不能直接变成任意非负 f、捕获后的 mask、receiver selector 或 actual 门的 Bernstein 增量界。

若 X0 用固定完整来源 nu/W，先 P_{r1,L1} 后沿这条原链，§4 的 receiver sensor count 表示可以使用它；但**表示不等于 overlap fee**，source-once 停止 mass<=W 仍留投影 D。固定来源可以沿核演化，不能将它误称每个新 L 都对应同一饱和 obstacle 的 nu_b、mu_b identity。

## 5. 三轮已注册原空间数值：完整保留负部与 refinement

专属 registration 在执行前保存。三轮均测试 r=.01/.1/.5/.9/.99，各有 5 个预设 q：.25、.5、1、2 倍 logistic slope 及 root 的 2(1-r)，每 round 共 25 records。物理有限步 logL=.01/.05/.1；logistic倍数使用对应 ODE 精确 endpoint，root 候选使用其精确 endpoint。未根据结果优化 slope。

Gauss spectral nodes 为 128/256/512，a-integral 对照 nodes 为 64/128/256，x nodes 为 129/257/513。s=1+t/(1-t) 将无穷远映到有限区间；没有截去尾部，也无 clipping。x 包括 0 和 [1e-6,1] 对数格。式 (7) 以 q R_{1-r}/[2(1-r)] 与 F 的直接积分求 H，避免式 (9) 的可积 endpoint singularity。generator 对 |x|>=1 解析非负（所有 s>=1，式 (7) 中 (sx-1)F>=0），因此负部只可能在 [-1,1]，保存偶对称 trapezoid 诊断。finite-step 差的数值负部只检查该区间，不把有限网当全空间认证。

三轮共 75 generator records、225 finite-step profiles、120 公式交叉验证，运行约 0.1205 秒：

|指标|128 nodes|256 nodes|512 nodes|
|---|---:|---:|---:|
|充分正斜率所有网点的最小 H|0.0018384111|0.0018384095|0.0018384091|
|充分正斜率 finite-step 最小差|1.8567216e-5|1.8567199e-5|1.8567196e-5|
|原 Ga Fourier 符号最大残差|5.421e-6|3.383e-7|7.061e-8|
|闭式 F 与独立 a 积分最大残差|6.262e-11|8.299e-15|2.526e-15|
|出现负 generator 的 records|10/25|10/25|10/25|

每轮负 generator 均为全部 5 个 r 的两个次临界倍数，最小值均在 x=0。最细轮 .25/.5 倍候选在 r=.99 时 H(0)约 -8.2133/-4.5860，负部网积分约 .37788/.17104；r=.01 的 .5 倍则只有 H(0)约 -1.034e-5、网积分约 2.586e-9。这些是浮点诊断，尤其小负部不能当 interval-certified 反例。正候选的普遍成立来自 §3 的解析证明，不来自三轮数值正。refinement 差没有认证 quadrature error；无随机 seed、无统计 CI、无独立样本声明。

保存全部原 x 网、s 网、weights、各 H 和每个 finite-step profile，均有 SHA-256。公式交叉验证使用原 Ga 的 spectral integral against 2s/(s²+xi²)，xi=0,.5,2,8，对照原 log(1+xi²)/(a xi²+(1-a)log(1+xi²))；另以独立 a 求积对照 s=1.0001,1.1,2,8 的闭式 F。登记与脚本的 hashes 写在结果 JSON；receipt 覆盖本文与全部 profiles。

收到 root 的原位置 moment 补充后，另登记并执行一个有理 Taylor 守卫，没有重跑上述 quadrature。三轮 n4/16/64、L=1/2,1,2，各用全部五个原有理 r，15 个 exact cases、75 项 identity checks 全通过。核单坐标 rate=r、原 P variance=rL²/2、其物理导数、normalized variance=(2-r)L²<=2L²；全部 Fraction，无随机样本。它核原 marginal coefficient 的代数和位置 moment，不把位置增量浓度推广到原输入输出交通。

## 6. 实际全参数覆盖仍缺什么

一条 odds=C L² 的 curve 在 [a,2a] 覆盖的是它自己的 r 窗口 [(r1),(4r1/(1+3r1))]。所有 C>0 的 curves 可以**集合意义上** foliating 全部 (r,L) 参数，但不是一条可共享一次 stopping-source 的序列。沿每条曲线复制一份 W、混合 conditional weak 或对 C 免费取 sup 都不合法。

原 receiver 还会选择不同 sigma,L、hard winner 和所有 fullfuture，原首跳/continuation 的物理 L 同变。本文只控制原辅助 P_r=K_r^(1/2) 的 c=1 斜向关系；没有证明这些真实选点满足式 (10) 的 partial order，亦没有将 fullfuture 压成一条 diagonal。hard sensor 可接到这条链上，但 full FIRST/history 与 CP/GP/source pair 资格不随表示自动恢复。

因此本轮是**原 Markov 跨尺度接口的一项真实进展**，不是实际 geom sqrt(n) 闭合。c<1 时 Lévy measure 是两个 phase 的差，c=1 的正谱证明不可免费延用；本轮只登记并核 c=1，没有以负谱系数冒充 c<1 空间反例。关于 c<1 或多条 diagonal 的 source-weighted overlap，需要另一个明确全覆盖合同，当前不运行旧 fixed-family weak 测试来填补它。

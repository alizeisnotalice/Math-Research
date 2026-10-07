# 物理尺度的有界熵正变差：一般线性上界与同源平方根必要下界

2026-10-07。只新增本稿与同前缀数值/证书，不改主账或旧数据。结果：一般任意正 L¹ 输入的总正变差有 source-once n log2·W 粗上界；存在完整正 L¹ 输入使它至少为 c√n·W，所以不能用维数无关或 polylog 的熵变差免费连接不同物理尺度。**尚未证明**所需 O(√n polylog n) 一般上界；下界不排除该目标，也不认证任何 actual FIRST/geom 样本。

## 1. 查重、文献与量词

本地搜索了 bounded entropy、正变差、entropy variation 等，已有 `ordered_continuous_bounded_entropy_20261007.md` 与独审只给**固定物理尺度**的 ordered 时间耗散；它们不包含本题物理尺度的标量总正变差。2026-10-07 外部检索了 “cube averages entropy variation”“convex norm-variation cubic averages”等组合，未定位到可直接引用的同一结果；这不是不存在性断言。

已读取 [Bourgain–Mirek–Stein–Wróbel, On dimension-free variational inequalities for averaging operators in R^d](https://arxiv.org/pdf/1708.04639) 的定义与 Theorems 1.1–1.3（PDF前两页）。其 cube 结论是 p>1、r>2 的 Lp(r-variation) 估计；本题是对任意 f≥0,L¹ 的非线性积分熵作 scalar V¹ 正变差，费用要求 W=∫f。它不直接给本题端点。没有把该文结果改名为新熵定理，也未宣称读审其全部证明。

也核对了 [Norm-variation of cubic ergodic averages](https://arxiv.org/abs/1903.04370) 的原文摘要与出版方 introduction；它研究多个 commuting transformations 的 cubic/multilinear ergodic averages 之 norm convergence/variation，不是本题同一正 L¹ cube 卷积的 bounded-entropy scalar variation。此项只作查重，不引用其技术结论。

应用 [L03](/Users/zhengzhihao/.codex/skills/math-l03-rational-and-interval-certificates/SKILL.md) 与已读 provenance/method/cube-interface 的证据守卫：数值 profile 与有理区间认证分开，完整输入与源质量先重构。Skill 不提供下述新数学结论。

固定 λ>0、任意 f≥0,f∈L¹(R^n)、W=∫f，u_L=h_L^⊗n*f，h_L 为全边長 L 的归一化中心 cube。定义

\[
 \Phi_\lambda(t)=t-\lambda\log(1+t/\lambda),\quad
 E_\lambda(L)=\int\Phi_\lambda(u_L)dx,
 \quad V^+_{[a,2a]}E=\sup_{a=L_0<\cdots<L_m=2a}
 \sum_{j=1}^m[E(L_j)-E(L_{j-1})]_+.
 \tag{1}
\]

0≤Φ≤t、0≤Φ'≤1、Φ''=λ/(λ+t)²。故 E 有限、≤W，不要求 L log L 或有界输入。不同 L 始终使用同一个 f；不按 signed Fourier 层分别领取 W。

## 2. 一般 n log2·W 粗费与真实端点下降

令 s=log L。归一化 cube 的尺度导数为

\[
 \partial_s u_L=n(v_L-u_L),
 \tag{2}
\]

v_L 是 2n 个面的等权概率平均与 f 的卷积，v_L≥0、∫v_L=W。对光滑 f 由 boundary differentiation 得到；一般 L¹ 输入用正卷积和 L¹ 近似，导数范数≤2nW，将 (2) 延拓为 L¹ 绝对连续身份。Φ 为1-Lipschitz，链式微分同样成立。因此

\[
 [\partial_sE_\lambda(L)]_+
 =\left[n\int\Phi_\lambda'(u_L)(v_L-u_L)dx\right]_+
 \le nW,
 \quad
 \boxed{V^+E_\lambda\le n\log2\,W.}
 \tag{3}
\]

source-once 但仍为线性阶。利用质量守恒还可准确写为

\[
 \partial_sE_\lambda(L)=n\lambda\int\frac{u_L-v_L}{\lambda+u_L}dx.
 \tag{4}
\]

这定位了真正待证的非线性 face 敏感度预算，不把 signed 正负部取消省掉。

端点则有严格的正平均关系：

\[
 h_{2a}=2^{-n}\sum_{\varepsilon\in\{-1,1\}^n}
 h_a(\,\cdot-a\varepsilon/2),
 \quad E_\lambda(2a)\le E_\lambda(a).
 \tag{5}
\]

第二式为 Jensen 与平移不变性，适用于全部完整 f。它不推出区间内单调；下面同源输入有多次正回升。

## 3. 完整正 L¹ 输入：√n 个共同频率侧瓣

先取 a=λ=1，n=M²，整数 M≥4，ε=1/8、D=10⁶n²。令

\[
 f_n(x)=1_{[-D/2,D/2]^n}(x)
 \prod_{i=1}^n[1+\varepsilon\cos(2\pi Mx_i)].
 \tag{6}
\]

这是完整同源输入：每因子≥7/8，非负且 L¹；MD 为整数，故每一维 integral=D，完整 W=D^n。没有 restricted source 的再归一化、signed 独立层、输出 retuning 或分尺度改变输入。

在距离输入盒边界至少 L/2 的接收点，真实卷积为

\[
 u_L(x)=\prod_i[1+\alpha_L\cos(2\pi Mx_i)],
 \quad\alpha_L=\varepsilon\frac{\sin(\pi ML)}{\pi ML}.
 \tag{7}
\]

在完整周期盒上按 normalized Lebesgue 积分，相位独立均匀，记 U=∏(1+α cos Θ_i)，EU=1。
谷点 L_k=k/M、k=M,…,2M 时 α=0，周期熵恰为 Φ_1(1)。峰测试点 L_k^+=(k+1/2)/M、k=M,…,2M−1 时

\[
 s=\alpha^2=\frac1{64\pi^2(k+1/2)^2},\quad
 \frac1{256\pi^2}\le ns\le\frac1{64\pi^2}<\frac1{576}.
 \tag{8}
\]

不需要取真实最大侧瓣峰；测试点给正变差下界即可。

## 4. 有界熵回升的严格原输入矩下界

原同源 U 的完整矩为

\[
 EU^2=(1+s/2)^n,\quad EU^3=(1+3s/2)^n,
 \quad EU^4=(1+3s+3s^2/8)^n.
 \tag{9}
\]

令 V=E(U−1)²=(1+s/2)^n−1，

\[
 F_4=E(U-1)^4=(1+3s+3s^2/8)^n
 -4(1+3s/2)^n+6(1+s/2)^n-3.
 \tag{10}
\]

定义 G(U)=Φ_1(U)−Φ_1(1)−Φ_1'(1)(U−1)。因 EU=1，周期熵回升恰为 EG。G≥0；在 |U−1|≤1 的区间，Φ''≥1/9，故 G≥(U−1)²/18。其补集 U>2 满足 (U−1)²≤(U−1)^4，因此

\[
 EG\ge(V-F_4)/18.
 \tag{11}
\]

下面给统一下界，而不是根据三轮拟合。置 q=ns≤1/(64π²)<1/576。作为 s 的函数，F_4(0)=F_4'(0)=0。对 0≤t≤s，用 (9) 逐项二次微分、(3+3t/4)²≤10 及各幂≤e^{4q}，得

\[
 F_4''(t)\le22n^2e^{4q}<44n^2,
 \quad F_4(s)\le22q^2.
 \tag{12}
\]

这里二阶系数可保守求和为10+3/4+9+3/2<22；e^{4q}<2 是初等严格界，无浮点决定。又

\[
 V\ge ns/2\ge1/(512\pi^2),\quad
 F_4/V\le11/(4\pi^2)<11/36<1/2.
 \tag{13}
\]

所以每个测试侧瓣的周期熵回升均满足

\[
 EG\ge V/36\ge1/(18432\pi^2)>1/294912.
 \tag{14}
\]

这些是同一原 f 的真实矩，未把 Φ 的 signed 级数分量当独立正来源。

## 5. 从周期计算严格还原到 R^n 的原 L¹ 输入

记周期响应在 [-D/2,D/2]^n 上的熵为 E_per(L)。接收 core 为 [-D/2+L/2,D/2−L/2]^n；其中原完整 f 的响应与 (7) 相同。

周期响应的一维因子≤2，core 外但原盒内的周期响应质量≤2nL·W/D。对于原响应，若来源距各输入边界至少 L，则加一个长度 L 的均匀 cube 后必落在 core；原来源的一维密度因子≤2，所以原响应在 core 外的质量≤4nL·W/D。0≤Φ≤u 因而给

\[
 |E_1(L)-E_{per}(L)|\le4nL W/D\le8nW/D.
 \tag{15}
\]

这个 full-space boundary 误差是原完整来源概率的 union bound，没有模拟边界或把 torus 当 L¹ 输入。每一 valley→test peak 的实际熵增量由 (14)–(15) 至少为

\[
 W[1/294912-16n/D]\ge W/589824\qquad(n\ge16).
 \tag{16}
\]

把 M 次 valley→peak 放入一个递增分割（随后 valley 的负增量不会扣除正变差），得到

\[
 \boxed{V^+_{[1,2]}E_1(f_n)\ge\sqrt n\,W/589824,
 \quad n=M^2,\ M\ge4.}
 \tag{17}
\]

任意固定 λ>0、a>0 可取 f(x)=λ f_n(x/a)，则 E 与 W 都乘 λa^n；同一 normalized 下界保持。故一般 o(√n) 预算不可能，尤其 dimension-free/polylog 预算不可能。它**不排除 O(√n polylog n)**，也不是一般 cube weak/geom 反例。

## 6. 三轮 L03 新证书与数值探路

执行前保存同前缀 `_registration.json`：λ=1、ε=1/8、D=10⁶n²，n/M=16/4、64/8、256/16；每维两份独立 seed，样本8192、16384、32768，nested nodes129、257、513。各尺度采用共同相位，完整 f 全程冻结。保存 `_round1/2/3_power_sums.npz`、输入/数组/script hashes，足以重构 profile。

严格证书：Machin π=16 atan(1/5)−4 atan(1/239)，24项 alternating rational 余项区间；所有 s、原2/3/4阶矩、(11)、boundary误差、实际 peak rise 以 Fraction 运算。保存方向明确的十进制有理外包（整数 floor/ceil），不以浮点残差判正。28个峰测试点全部 actual L¹ rise lower>0，另独立从保存 π/原输入参数重算全部28项，不调用主脚本。

|n|认证实际正增量和/W 下界|两份经验 grid V+/W|fine−coarse（诊断）|
|---:|---:|---:|---:|
|16|0.0000830732|0.000202323、0.000197076|4.39e−7、4.28e−7|
|64|0.000173306|0.000400192、0.000391394|5.42e−20、浮点0|
|256|0.000350203|0.000791052、0.000784682|两份浮点0|

数值采用真实乘积 U 的 degree10 log Taylor 加速并保留解析 remainder，64个原乘积样本在多尺度重算；该 profile 的浮点、trigonometric、Monte Carlo 误差未做区间认证。保存 iid sample standard errors，未冒称 uniform confidence interval；正变差从经验均值取正部，有选择/非线性偏差，coarse/fine差不是连续变差误差证书。严谨下界只来自独立的有理矩证书和 (15)，不靠这些经验 profile。

原输入验证另用 rational clipping limits 的一维真实有限盒 antiderivative，在 interior/boundary接收点核对 (7) 并保留完整 source；trigonometric数值检查仅是浮点实现复核。正性、MD整数、W=D^n 与 boundary预算是解析/有理证书。总运行22.82秒，无旧数据重跑。

## 7. 统一 Fourier/sinc 侧瓣：能给什么，不能给什么

对任意 f∈L²，二次能量 Q(L)=∥h_L*f∥²_2 有一个一般 dimension-free **L²-input** 正变差界：

\[
 \boxed{V^+_{[a,2a]}Q\le(2\log2/\pi)\|f\|_2^2.}
 \tag{18}
\]

证明使用 angular Fourier convention，m_L(ξ)=∏sinc(Lξ_i/2)。固定 ξ，将 i 按 |aξ_i/2|<π/2 与≥π/2分组。前一组在全部窗口 t=L|ξ_i|/2≤π 上 sinc²(t) 非增，不能贡献正导数。后一组有 b 个坐标，各因子≤ρ=4/π²<1/2，且

\[
 [\partial_{\log L}\operatorname{sinc}^2(t)]_+
 \le1/t\le2/\pi.
\]

乘积的正导数≤(2/π)bρ^{b−1}≤2/π，因为 b≤2^{b−1}。b=0 的正导数为0。由 Plancherel、正部放大、log窗口长度log2得 (18)。对固定 n，multiplier导数有 uniform O(n) 界，足以对任意 L² f 用 dominated convergence/绝对连续积分；不是只检查单 mode 后外插。

但 (18) 不给所求 W 预算：任意正 L¹ 输入可以不在 L²，即使同时在 L²，也没有 ∥f∥²_2/λ≲W 的 uniform关系。本轮完整 source 恰有 ∥f_n∥²_2/W=(1+ε²/2)^n，随 n 指数增大。更根本地，∫Φ_λ(u_L) 的导数是 (4)，不是一个 Fourier diagonal quadratic form；不能由 Φ≤u²/(2λ) 推断它的总正变差被 Q 的总正变差支配。点态/标量大小比较不控制正变差。用 Φ 的 signed power expansion 再逐层计费也会丢同一原来源预算。

目前没有从 (18) 推出 bounded nonlinear entropy 的 O(√n polylog n) source-once上界。统一 sinc 侧瓣可以完成 quadratic L² 接口，但连接到 (4) 的任意 f、固定 λ 的非线性行权重仍是具体缺项。

## 8. 对一般 geom 连接的结论

物理尺度 entropy route 不能期待 polylog variation 免费转接；平方根阶本身已经是一般接口的必要成本。本轮 upper/lower 范围为 c√nW 的必要下界与 nlog2W 粗上界，夹在中间的 O(√n polylog n) 上界未证。即使后者成立，还需把原 actual FIRST、不同 hard/soft winners、门控差分负部/commutator与该总变差准确相接，不能把本稿的完整未门控 f 变差自动算作 R_dagger 已付。

本稿不改当前联合主账，不把 endpoint熵下降误写成连续单调，也不把产品构造推广为所有输入的上界。

# 完整 all-visited 核的空间 Hörmander 界与固定 family 弱型接口

2026-10-07。本轮不重做已付 holding/proper-face 分支，不修改总稿或实际主账。新结果是原冻结轴向 family 的 **all-visited 高频 Banach-valued Hörmander 界**。证明保留整个时间乘积及跨跳层取消；对全部 \(0<c\le1\)、\(\ell>0\) 一致。在预设 \(\Lambda=1024n^2\) 下，结合先前低频 Abel/obstacle 工具，单个冻结 semigroup 的 maximal 具有 \(O(\log(n+2))\) 弱型常数。实际 cube 的时变参数、共同来源与实际门的 transfer 仍未证明，本稿不向 cube 主账计费。

## 1. 对象、阅读范围与来源

仍使用 F02 的正性/张量工作流；它不提供 endpoint 定理。此前逐层 triangle 失败、reset 链及 proper-face 收据不重跑。唯一外部定理输入为 [Gomilko–Tomilov, arXiv:1410.1505, Corollary 3.5, (3.11)](https://arxiv.org/pdf/1410.1505)：本轮实际核读该处、其 \(c_0=0\) 说明及相关假设，不声称全文验证。对 \(L^1\) 的一维热 semigroup，\(M(A)=1\)、Yosida \(c_0=0,c_1\le1\)；下文自行核 BF 和这些热核常数，再代入该定理。空间估计、矩形传播和轮廓移位是本稿推导，不是该文现成结论。

取空间 Fourier 约定 \(\widehat f(\zeta)=\int e^{-ix\cdot\zeta}f(x)dx\)。先设 \(\ell=1\)，最后用统一 dilation 回代。令
\[
B(\lambda)=\frac{\lambda}{\log(1+\lambda)}-1,
\qquad \psi_c(\lambda)=\frac{B(\lambda)}{1+cB(\lambda)},
\qquad \widehat G_c(\zeta)=\frac1{1+cB(\zeta^2)}.
\tag{1}
\]
固定生成元与 semigroup 为
\[
L=\sum_{i=1}^n\psi_c(-\partial_i^2)
=c^{-1}\sum_i(I-P_i),\qquad H_t=e^{-tL}.
\tag{2}
\]
总跳率 \(n/c\)，内部 Poisson 坐标均匀标签与输入来源坐标独立性无关。任意非产品、相关、XOR 或 Cantor 来源作用同一线性算子；all-visited 不是来源复杂度限制。

一维 \(H_t^{(i)}\) 的 holding 质量 \(e^{-t/c}\)，continuous 密度记 \(p_t\)。其 Gaussian-clock 正混合表示给出 \(p_t\ge0\)、偶且递减，质量 \(1-e^{-t/c}\)。完整 all-visited 密度为
\[
P_t^{\rm full}(x)=\prod_{i=1}^np_t(x_i),\qquad
H_t^{\rm full}=\bigotimes_{i=1}^n(H_t^{(i)}-e^{-t/c}I).
\tag{3}
\]
这是每个跳跃坐标至少访问一次的全部标签；等于前稿 \(C_k\) 汇总后的完整时间核。证明从此保留 (3)，不先对各 \(k\) 取绝对值。

## 2. 原符号的全尺度 \(2/3\) scaling

由
\[
B(\lambda)=\int_0^1[(1+\lambda)^s-1]ds
\]
可见 \(B\) 是 BF：每个 \(0\le s\le1\) 因子的导数完全单调。\(x/(1+cx)\) 也是 BF，复合得 \(\psi_c\) BF，且 \(\psi_c(0)=0\)。因此正 Gaussian subordination 及下面时间解析定理合法。

进一步，对全部 \(\lambda>0\)，
\[
\boxed{\quad \frac23\le\frac{\lambda B'(\lambda)}{B(\lambda)}\le1.\quad}
\tag{4}
\]
上界为 BF 凹性。下界的具体推导如下。令 \(w(s)=(1+\lambda)^s-1\)；因 \(w\) 凸且 \(w(0)=0\)，\(w(s)/s\) 递增。在概率底 \(2s\,ds\) 下对 \(s\) 与 \(w(s)/s\) 用正协方差，得到
\(\int sw(s)ds/\int w(s)ds\ge\int 2s^2ds=2/3\)。逐 \(s\) 又有
\[
\lambda s(1+\lambda)^{s-1}\ge s[(1+\lambda)^s-1]
\]
（两边相减为 \(s[1-(1+\lambda)^{s-1}]\ge0\)）。积分即为 (4)。故对 \(r\ge1\)，
\[
r^{2/3}\le B(r\lambda)/B(\lambda)\le r.
\tag{5}
\]
不能把此正指数 global scaling 原封不动用于饱和的 \(\psi_c\)。证明始终在 \(B\) 的实际阈值处分区。

## 3. 真实 continuous 密度的统一峰值与空间差

记 \(\beta_t=B^{-1}(1/t)\)、\(\beta_c=B^{-1}(1/c)\)。由 (5)，
\(\int_0^\infty\widehat G_c(\zeta)d\zeta\le4\sqrt{\beta_c}\)，故 \(G_c(0)\le(4/\pi)\sqrt{\beta_c}<\infty\)。每个 convolution power 的峰值不超过 \(G_c(0)\)，因此 \(p_t(0)\le(1-e^{-t/c})G_c(0)\)。

更关键的统一大时间界为
\[
\boxed{\ p_t(0)\le5\sqrt{\beta_t},\qquad t\ge c.\ }
\tag{6}
\]
Fourier 反演中必须保留 \(1/\pi\)：
\[
p_t(0)=\frac1\pi\int_0^\infty
[e^{-t\psi_c(\zeta^2)}-e^{-t/c}]d\zeta.
\]
在 \(\zeta\le\sqrt{\beta_c}\)，有 \(\psi_c\ge B/2\)。再以 \(\sqrt{\beta_t}\) 分开，前段至多 \(\sqrt{\beta_t}\)，后段由 (5) 至多
\(\sqrt{\beta_t}\int_1^\infty e^{-u^{4/3}/2}du\le2\sqrt{\beta_t}\)。
在 \(\zeta\ge\sqrt{\beta_c}\)，令 \(a=t/c\ge1\)、\(g=(1+cB)^{-1}\le1/2\)：
\[
e^{-a}(e^{ag}-1)\le ag e^{-a/2},\qquad
\int_{\sqrt{\beta_c}}^\infty g\,d\zeta\le3\sqrt{\beta_c}.
\]
又 \(\sqrt{\beta_c/\beta_t}\le a^{3/4}\)，所以该段至多
\(3a^{7/4}e^{-a/2}\sqrt{\beta_t}\le12\sqrt{\beta_t}\)；这里使用
\(a^{7/4}\le a^2\)、\(\sup a^2e^{-a/2}=16e^{-2}<4\)。合计除以 \(\pi>3\) 便得 (6)。所有常数独立于 \(c\)。

固定非零平移 \(h\)，令 \(r=\|h\|_\infty\)、\(A=\psi_c(r^{-2})\)。一维整核的真实尾有
\[
\mathbb P(|X_t|>r)\le8tA.
\tag{7}
\]
证明为在 \(0\le\zeta\le2/r\) 平均 \(1-\cos(\zeta X_t)\)：当 \(|X_t|>r\) 时其平均至少 \(1/2\)，而 Fourier 期望至多 \(t\psi_c(4r^{-2})\le4tA\)。holding 原子并未被漏掉。

定义真实外部空间差
\[
D_t(h)=\int_{\|x\|_\infty>2r}
|P_t^{\rm full}(x-h)-P_t^{\rm full}(x)|dx.
\]
小时间用坐标 union、(7) 及 \(\|x-h\|_\infty>r\)，得 \(D_t\le16ntA\)。大时间 \(tA\ge1\) 蕴含 \(t\ge c\)、\(\beta_t\le r^{-2}\)。BF 凹性给
\[
r\sqrt{\beta_t}\le[tB(r^{-2})]^{-1/2}\le(tA)^{-1/2}.
\]
偶递减密度满足 \(\|p_t(\cdot-h_i)-p_t\|_1\le2|h_i|p_t(0)\)。按坐标逐次平移 tensor，其他因子的质量均 \(\le1\)，从 (6) 得全空间差 \(\le10n(tA)^{-1/2}\)。因此
\[
\boxed{\ D_t(h)\le16n\min\{tA,(tA)^{-1/2}\}.\ }
\tag{8}
\]
这是真正的空间差，非时间解析界的重命名。原生成元内部的坐标独立性用于核表示，未用于来源模型。

## 4. 时间解析如何合法转成复时间空间衰减

一维热核直接给
\(\|tAe^{-tA}\|_{1\to1}=\tfrac12\mathbb E|Z^2-1|\le1\)，其中 \(Z\) 为标准 Gaussian；正热 semigroup 的 \(M(A)=1\)。由 §1 已核的 Corollary 3.5，代入 BF \(\psi_c\) 得
\[
\|t\psi_c(A)e^{-t\psi_c(A)}\|_{1\to1}\le4.
\tag{9}
\]
这里用的是一维常数，不用 \(4n\) 免费推出空间差。
分成 \(k\) 个可交换 semigroup 因子得
\(\|L_1^kH_t^{(1)}\|\le(4k/t)^k\)。Taylor、\(k!\ge(k/e)^k\)、\(e<3\) 给当 \(|z-t|/t\le1/64\) 时
\(\|H_z^{(1)}\|\le\sum_k(12|z-t|/t)^k<2\)。
特别，令
\[
\theta=1/128,
\qquad |\operatorname{Im}v|\le\theta.
\]
因为 \(|e^{i\operatorname{Im}v}-1|\le\theta\)，此 strip 上每个 continuous 因子的 convolution-measure 范数 \(\le2+|e^{-e^v/c}|\le3\)。整个 (3) 的 \(L^1\) kernel 范数 \(\le3^n\)。该结论来自含取消的解析 semigroup；直接对复 Poisson 系数取绝对值产生的指数增长不能替代它。

令 \(F_h(v)\) 为 (3) 的平移差限制于 \(\|x\|_\infty>2r\)，视为 \(L^1_x\)-值解析函数。strip 内
\(\|F_h(v)\|_1\le M=2\cdot3^n\)，实轴上有 (8)。在以 \(v\) 为横向中心、宽 \(2\theta\)、高 \(\theta\) 的上矩形中，其底边最大范数至多
\[
m=16e^\theta n\min\{s,s^{-1/2}\},\qquad s=e^vA.
\]
从点 \((v,\theta/2)\) 的底边 harmonic measure 至少 \(1/4\)：内接正方形 \([v-\theta/2,v+\theta/2]\times[0,\theta]\) 的中心对四边退出概率各 \(1/4\)，其底边属于外矩形底边。这同时适用于 \(0\le\operatorname{Im}v\le\theta/2\)：或用相应矩形 harmonic function 的底边下界。严格的后一结论可由显式 subsolution
\[
\cos\!\frac{\pi x}{2\theta}\,
\frac{\sinh[\pi(\theta-y)/(2\theta)]}{\sinh(\pi/2)}
\]
在 \(x=0,0\le y\le\theta/2\) 给出下界
\(1/[2\cosh(\pi/4)]>1/4\)。

对每个范数 \(\le1\) 的 dual functional 用标量 two-constants/最大原理，再取范数；若 \(m>M\) 则直接用 \(M\)。上下 strip 对称，得到
\[
\|F_h(v+iy)\|_1
\le M^{3/4}(32n)^{1/4}
\min\{s^{1/4},s^{-1/8}\},\quad |y|\le\theta/2.
\tag{10}
\]
因 \(e^\theta<2\)。于是完整复时间空间差可积，且
\[
\int_{\mathbb R}\|F_h(v\pm i\theta/2)\|_1dv
\le J_n:=48\,3^{3n/4}n^{1/4}.
\tag{11}
\]
这里小/大时间积分分别为 \(4\)、\(8\)；常数从 \(12(2\cdot3^n)^{3/4}(32n)^{1/4}\) 得出。端点的衰减也在全部中间高度一致，后续轮廓左右边真正趋零。

## 5. 全核 Mellin 匹配、局部 Banach 核及 Hörmander

使用 logtime Fourier \(\widehat f(\xi)=\int e^{-2\pi i\xi v}f(v)dv\)，令 \(\omega=2\pi\xi\)。完整 all-visited 的准确空间核为
\[
\mathcal K_\omega^{\rm full}(x)
=\frac{i\omega}{1+i\omega}
\int_{\mathbb R}e^{-i\omega v}P_{e^v}^{\rm full}(x)dv.
\tag{12}
\]
它不是把 \(H_{e^v}\sigma\) 的初始常量当普通 Fourier。对 \(\sigma=Lu\)，原 derivative Mellin 定义仍沿用前稿。

首先核 (12) 的定义与 operator 匹配。固定 \(c\) 时，小 \(t\) 的 full 质量 \(\le(t/c)^n\)；大 \(t\) 时 (6) 及 \(B(\lambda)\ge\lambda/4\) 对 \(\lambda\le1\) 给 \(p_t(0)\le10t^{-1/2}\)（\(t\ge4\) 足够）。这里凹性给 \(B(\lambda)\ge\lambda B(1)\)，而 \(\log2<4/5\) 蕴含 \(B(1)>1/4\)，故 \(\beta_t\le4/t\)；没有在低频使用未经核实的渐近。因而对 \(f\in L^1\cap L^2\)，
\[
\|H_t^{\rm full}f\|_2\le(t/c)^n\|f\|_2\quad(t\le c),
\qquad
\|H_t^{\rm full}f\|_2\le10^{n/2}t^{-n/4}\|f\|_1\quad(t\ge4).
\]
故其 logtime Mellin 是真正的 \(L^2\) Bochner integral。空间 Fourier 上对每个非零 \(\zeta\)，(3) 乘积的小 \(t\) 消去和大 \(t\) 衰减允许积分；有限 inclusion-exclusion 或正则化 binomial 给 (12) 等于前稿
\[
\widehat R_{\rm full}(\xi)
=\widehat R_\sigma(\xi)-\widehat R_{\rm miss}(\xi),
\qquad
\widehat R_\sigma=-\frac{\Gamma(1-i\omega)}{1+i\omega}L^{i\omega}\sigma.
\tag{13}
\]
原 \(C_k\) 的全和仍是 Abel-limit，未宣称全 \(k\) 的 \(L^1\) 绝对收敛。\(\omega=0\) 时 full 为零；空间零频是 Lebesgue 零集，固定欧氏 family 没有 \(L^2\) 常数模。有限链的零模不在本定理范围。

其次，(12) 给合法的局部 Banach kernel。对任意有界空间集 \(E\)，
\(\|\mathbf1_E P_t^{\rm full}\|_1\le(t/c)^n\) 于小 \(t\)，
\(\le|E|10^nt^{-n/2}\) 于大 \(t\)。同 §4 的局部 \(L^1(E)\) two-constants 和轮廓移位，得
\(\int_{|\xi|>\Lambda}\|\mathbf1_E\mathcal K_{2\pi\xi}^{\rm full}\|_1d\xi<\infty\)。局部常数允许依赖 \(c,E\)；这只是定义核，不把它当统一源费用。卷积在测试输入上与 (13) 一致，再由 \(L^2\) 延拓匹配。因此离对角 kernel representation 不再是假设。不能省略这一点，只凭 multiplier 估计宣称 CZ。

最后，对平移差应用 (10)–(11)，将 logtime 轮廓向
\(-\operatorname{sign}(\omega)i/256\) 移位。左右边在 \(L^1\) 趋零，故
\[
\int_{\|x\|_\infty>2\|h\|_\infty}
|\mathcal K_\omega^{\rm full}(x-h)-\mathcal K_\omega^{\rm full}(x)|dx
\le J_n e^{-|\omega|/256}.
\tag{14}
\]
频率连续积分由正 Tonelli 合法，得到原待证 (15) 的具体值
\[
\boxed{
H_\Lambda\le\frac{256}{\pi}J_n e^{-\pi\Lambda/128}
<4096\,3^n n e^{-3\Lambda/128}.}
\tag{15}
\]
所有空间差、核表示和 cutoff 费对全部 (c\in\(0,1]\)、\(\ell>0\) 一致。恢复 \(\ell\) 时，取 \(x/\ell,h/\ell\)；Hörmander 的空间积分 Jacobian 正好抵消 kernel dilation，时间相位也不改变模。

## 6. 明确阈值的弱型高尾，而非强 \(L^1\) 尾

由旧 face 稿已证的强 \(L^2_x(L^1_\xi)\) bound，
\[
B_\Lambda\le8e^{-\Lambda}
+\frac{20}{\pi}n^{5/2}e^{-\pi\Lambda/(5\sqrt n)}.
\tag{16}
\]
新 (15) 和已核 kernel representation 使其 CZ 合同现在可以调用；具体仍在 height
\(h_0=\tau(3/2)^{n/2}/(2B_\Lambda)\) 分解：三倍 bad cubes 体积费 \(3^n\|f\|_1/h_0\)，good 的 Chebyshev 费 \(4B_\Lambda^22^nh_0\|f\|_1/\tau^2\)，bad 的外部 \(L^1\) Banach 费 \(2H_\Lambda\|f\|_1\)。因此
\[
\tau|\{\|\mathcal T_\Lambda f\|_{L^1_\xi}>\tau\}|
\le4(6^{n/2}B_\Lambda+H_\Lambda)\|f\|_1.
\tag{17}
\]
这处理连续频率的强 Banach norm，不使用 \(L^{1,\infty}\) Minkowski 或未经核验的 tensor weak endpoint。

固定输入无关的多项式 cutoff
\[
\boxed{\ \Lambda_n=1024n^2.\ }
\tag{18}
\]
由 \(\log2<1,\log3<2,\log n\le n\)，(15) 右侧严格小于
\(e^{12+3n-24n^2}\le e^{-9}<1/512\)。又 \(6^{n/2}<e^n\)：(16) 的第一项乘上该因子 \(<8e^{n-1024n^2}<1/1024\)，第二项 \(<e^{2+(7/2)n-(3072/5)n^{3/2}}<1/1024\)。所以
\[
H_{\Lambda_n}<1/512,\quad
6^{n/2}B_{\Lambda_n}<1/512,
\qquad
\tau|\{Q^{\rm full}_{\Lambda_n}>\tau\}|
<\|f\|_1/64.
\tag{19}
\]
取原 \(f=\sigma=\nu_{\rm b}-\mu_{\rm b}\)，\(\|\sigma\|_1\le2W_{\rm b}\)，full 高频费用 \(<W_{\rm b}/32\) 为弱型。旧 miss 在更大 cutoff 保留强 \(L^1\) 费用 \(<W_{\rm b}/8\)。显式以阈值 \(\tau/2\) 分配，得
\[
\boxed{\quad
\tau|\{Q_{\Lambda_n}>\tau\}|\le5W_{\rm b}/16.
\quad}
\tag{20}
\]
未证明 full \(Q\) 的强 \(L^1\) 预算，也未证明普适 \(L^1\) imaginary powers；这两个更强结论不用于 (20)。

## 7. 固定 family maximal 的回代与真正仍缺的接口

沿用已审 [固定 jump 障碍接口](fixed_jump_obstacle_interface_20261007.md)，对非负 \(\nu\in L^1\cap L^2\) 及任意 cap \(\kappa\) 构造
\(\nu=\nu_{\rm b}+\nu_{\rm g}\)、\(\mu_{\rm total}=\mu_{\rm b}+\nu_{\rm g}\le\kappa\)、\(\sigma=Lu=\nu_{\rm b}-\mu_{\rm b}\)、\(\kappa|\Omega|\le W_{\rm b}\le W\)。只使用同一原来源；外部初始 good 源不进入 \(Ku=\nu_{\rm b}\)。在 \(\Omega^c\)，
\[
H_t\nu=H_t\mu_{\rm total}
+(u-H_tu)/t+R_t\sigma\le\kappa+R_t\sigma.
\tag{21}
\]
前稿 fixed Abel/derivative-Mellin 仅在原域外给 \(\|T_\Lambda\sigma\|_{L^2(\Omega^c)}^2\le C_0\log^2(e\Lambda)\kappa W_{\rm b}\)，其中 \(T\) 支配全部时间的低频 inverse；高频 inverse 最大值由 \(Q\) 支配。这里 \(C_0\) 沿用已审谱常数，不从本轮数值拟合。

对输出阈值 \(\alpha>0\)，预设
\(\kappa=\alpha/[2\log(e\Lambda_n)]\le\alpha/2\)。域内体积费至多 \(2\log(e\Lambda_n)W\)。域外以 \(T>\alpha/4\)、\(Q>\alpha/4\) 覆盖，由 Chebyshev 和 (20) 分别付
\(8C_0\log(e\Lambda_n)W\)、\(5W/4\)。因此
\[
\boxed{\quad
\alpha|\{\sup_{t>0}H_t\nu>\alpha\}|
\le[(2+8C_0)\log(e\Lambda_n)+5/4]W
=O(\log(n+2))W.
\quad}
\tag{22}
\]
初始仅 \(L^1\) 的非负密度可用递增的有界紧支撑截断，先对有理时间单调极限，再用 bounded-\(L\) 幂级数的 a.e 紧时间一致连续代表回代全时间；不增加源质量或 smoothing 费。有限原子输入的奇异 holding/轴面输出不在此密度 maximal 声明内，不把原子 mollification 默认为免费等价。

这是一个可复用的**单个冻结生成元 family 定理**，常数与 \(c,\ell\) 无关。它不自动估计
\(\sup_{c,\ell,t}H^{c,\ell}_t\nu\)，更不估计原有序时变 cube 的 actual FIRST 核。不同冻结算子各自的 \(u,\Omega,\mu_{\rm b}\) 不能并为共同一份 \(W\)；原 cube 响应与此 maximal 的正点态支配、同参数来源分配、原实际门 transfer 均仍需独立证明。故此轮不更改 actual cube 的未付空间余项或总账。

## 8. 三轮原符号注册守卫与终态

先保存 `frozen_jump_allvisited_hormander_registration_20261007.json`，后执行同前缀新 guard。三轮 \(n=8,32,128\)，GL 配对阶数分别 \(16/32,32/64,64/128\)；每轮预设 \(c=1,1/16,1/4096\)。核原 \(B\) 的全尺度 elasticity/scaling、正 continuous 峰值在 \(t/c=1,2,4,16,64\) 的方向，以及 inverse-clock 空间幅度。为避免小参数消去，数值 \(B\) 在 \(\lambda<10^{-3}\) 用固定四阶 Taylor，多余项未作区间证书；这与浮点 inverse/GL 一并属于实现诊断。峰值 Fourier 积分保留显式正尾上界；GL 部分非 interval。

另用原 \(\psi_c\) 的完整 (3) multiplier，取 balanced、单坐标大频率、五尺度三种 Fourier 配置及 \(\omega=1/2,1,2,4\)，检查实轴积分与向下 \(i/256\) 的**完整时间积分**一致。108 条 contour 记录均保留低端 \(e^{-40}/n\) 和高端 \(2^ne^{-T\cos(1/256)\sum\psi}/[T\cos(1/256)\sum\psi]\) 尾，不用 reset 或逐 \(k\) 截断代替原核。GL 粗细差和 contour 差为浮点诊断，不是区间认证；全无穷频率、全 \(c\) 以及 all-dimensional 常数由解析证明负责。

三轮一次运行 `PASS_EXACT_AND_NUMERIC`，708 个谓词，保存脚本/注册 SHA256。精确 Fraction 只核 cutoff 代数、Taylor 余量和弱型预算；符号、Fourier 密度与 contour 为数值守卫，不拟合 dimension 阶数，不认证实际 FIRST/历史门。旧实验和 reset 均未重跑。

保存结果只读核对：两个 SHA256 均匹配；45 条峰值记录最大 GL 配对差 \(1.918\times10^{-8}\)，108 条 contour 最大实轴/移位差 \(2.048\times10^{-15}\)，最大 contour GL 配对差 \(1.007\times10^{-15}\)。这些是数值实现诊断，不能解释为原无穷空间积分的区间误差。

**终态：原 all-visited Hörmander (15) 已证；预设多项式 cutoff 下完整高尾弱型 (20) 与固定 family maximal (22) 闭合。未证明 full 高尾强 \(L^1\)，未证明参数 supremum、原时变实际 cube transfer；cube 主账仍不计本工具费。**

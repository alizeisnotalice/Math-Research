# 原有序噪声向固定 jump 工具的局部转移与时变 forcing 缺口

2026-10-07。只新增本稿及同前缀注册/守卫/结果，不修改总稿、主账或他人文件。可证的新接口是指定窄窗口中的正核支配；全 early 的冻结误差有一份不刷新来源的 \(O(\sqrt n)W\) **forcing 时空总变差**预算。后者尚非可选输出时刻的已付误差，不能用固定 maximal 的条件弱界作平均。本稿没有闭合原时变实际 cube 余项。

## 1. 已读对象、技能和禁止偷换

本轮实际读取 J01、J05 的 SKILL、provenance、method、cube-interface，按其流程核真实生成元、冻结误差与源预算；技能没有提供转移定理。读取原 \(\texttt{概率接口阶段证明.tex}\) 4878–5130 的 NB、有序退出及原平均反例，8353 后的共同帽障碍；再读新 fixed all-visited endpoint 及物理尺度转移稿。旧反例/守卫不重跑。

固定物理尺度 \(\ell\)，原有序噪声为
\[
N_{s,\ell}=\big[(1-s)\delta_0+s w_\ell\big]^{\otimes n},
\qquad
\widehat w(\zeta)=m(\zeta^2)=\frac{\log(1+\zeta^2)}{\zeta^2}.
\tag{1}
\]
以 \(B(v)=1/m(v)-1\)、\(c_s=1-s\) 记
\[
\widehat G_{c_s,\ell}(\zeta)=\frac1{1+c_sB(\ell^2\zeta^2)},
\qquad
L_{s,\ell}=c_s^{-1}\sum_i(I-G_{c_s,\ell}^{(i)}).
\tag{2}
\]
在每个 \(s\le S<1\) 的窗口，\(\partial_sN_s=-L_sN_s\) 在 \(L^1\) 合法，\(L_s1=0\)。全空间卷积生成元彼此可交换，但不是同一个固定生成元；killed 域和实际标签不能据此免费交换。

当前 actual early 首跳入口满足 \(t\le1/2\)，所以对应冻结参数只需
\[
\boxed{c_t=1-t\in[1/2,1].}
\tag{3}
\]
不把后续 continuation 的终点 \(\sigma\) 误写成仍 \(\le1/2\)：入口 \(t\le1/2\) 与剩余相对窗口长度是不同条件。以下全 early forcing 预算只覆盖演化时段 \(0\le s\le1/2\)，不覆盖 arbitrary continuation 的全部 \([t,1]\)。

fixed theorem 对每个指定 \((c,\ell)\) 给 \(O(\log(n+2))\) 弱 maximal，uniform \(c,\ell\) 不等于这些参数的 supremum 定理。已核物理尺度 Gaussian 网给固定 \(c\)、同一密度输入、\(\sup_{\ell\in[a,b],\tau>0}\) 的 \(O(\sqrt n\log(n+2))\) 弱界；本稿不额外重收该费用，也不把其初始输入换为逐尺度 \(h_\ell^{\otimes n}*\mu\)。

## 2. 原 cocycle 的精确 Bernoulli 表示

给定入口 \(t_0\)、\(c_0=1-t_0\)，令
\[
r=(S-t_0)/c_0,\qquad 0\le r<1.
\]
单坐标原 multiplier 直接约去：
\[
\frac{c_S+(1-c_S)m(v)}{c_0+(1-c_0)m(v)}
=\frac{1+c_SB(v)}{1+c_0B(v)}
=(1-r)+r\widehat G_{c_0}(v).
\]
所以原完整转移为
\[
\boxed{
K_{S,t_0,\ell}
=\big[(1-r)I+rG_{c_0,\ell}\big]^{\otimes n}.}
\tag{4}
\]
保持 holding、所有部分坐标面、混合来源和每个 coordinate convolution；不是重新 FIRST 的对象。这不是 frozen compound-Poisson semigroup，也不是一份共同时间的普通半群平均。

尤其在 \(t_0=0\)，\(G_1=w\)，(4) 正是 \(N_S\)。在任意 fixed \(t_0,\ell\)，可对同一 \(\nu_0=N_{t_0,\ell}*\nu\) 使用 (4)，但 \(\|\nu_0\|_1=W\) 并不授权跨入口重新领取 \(W\)。

## 3. 一个指定窄窗的真实正支配

令 frozen semigroup
\[
H_\tau^{c_0,\ell}=\exp[-\tau c_0^{-1}\sum_i(I-G_{c_0,\ell}^{(i)})],
\quad a=\frac{r}{1-r},\quad \tau=c_0a.
\]
其 Poisson 展开保留每个坐标恰零跳或一跳的序列。对任意给定掩码 \(A\)、\(j=|A|\)，(4) 的原权为 \((1-r)^{n-j}r^j\)，上述 Poisson 子项权为 \(e^{-na}a^j\)，两者比例恰为
\[
C_n(r)=(1-r)^ne^{nr/(1-r)}.
\]
全部其他 Poisson 序列为正，因此
\[
\boxed{
K_{S,t_0,\ell}\le C_n(r)H^{c_0,\ell}_{c_0r/(1-r)}
}
\tag{5}
\]
为正卷积 measure 不等式。没有源坐标独立假设，没有丢掉重复跳的负项；重复跳在右边只提供额外正质量。

常数有明确连续预算：
\[
\log C_n(r)
=n\int_0^r\frac{u}{(1-u)^2}\,du
\le\frac{nr^2}{2(1-r)^2}.
\tag{6}
\]
预设 \(m_n=\lceil\sqrt n\rceil\)、\(r_n=1/(2m_n)\)。对 \(0\le r\le r_n\)，\(nr^2\le1/4\)、\(r\le1/2\)，故
\[
C_n(r)\le e^{1/2}<2.
\tag{7}
\]
这是参数和输入无关的前向覆盖；不是拟合一个有利时钟。冻结时钟 \(\tau=c_0r/(1-r)\) 与真实 elapsed \(S-t_0=c_0r\) 不同，不能在 Duhamel 中混用。

若一个指定 \(t_0,c_0,\ell\) 窗口确用共同非负 \(L^1\) 密度 \(\nu_0\)，(5) 给
\[
\alpha|\{\sup_{t_0\le S\le t_0+c_0r_n}
K_{S,t_0,\ell}*\nu_0>\alpha\}|
\le2C_{\rm fr}(n)\|\nu_0\|_1.
\tag{8}
\]
若同一 \(\nu_0\) 还独立于物理尺度 \(\ell\)，可用已付尺度网给单窗 \(\sup_\ell\) 的 \(O(\sqrt n\log(n+2))\|\nu_0\|_1\) 费用。入口 \(t_0=0\) 的 noise-only 小软度分支有这种共同原输入；\(t_0>0\) 时 \(N_{t_0,\ell}*\nu\) 随 \(\ell\) 变，不自动满足它。完整原响应 \(h_\ell^{\otimes n}*K\) 的硬平均缺口同样仍在。

实际 early 入口的 \(c_0\in[1/2,1]\) 不保证其 future remaining fraction \(r=(\sigma-t_0)/c_0\le r_n\)。故 (8) 是原短 continuation 的合法局部接口，不是对全 early-source continuation 的覆盖。超过 \(r_n\) 后，(6) 会产生维数损失，不能把其系数仍写常数。

## 4. 冻结 Duhamel 的准确误差与源预算

这里先固定一个物理尺度和同一指定入口。若 \(U_s=K_{s,t_0}*\nu_0\) 为原 free 演化，(2)、(4) 的精确消去给
\[
(L_s-L_{t_0})K_{s,t_0}
=\frac{s-t_0}{c_0^2}
\sum_i(I-G_{c_0}^{(i)})^2
\bigotimes_{j\ne i}[(1-r)I+rG_{c_0}^{(j)}].
\tag{9}
\]
此处右边是带符号算子，不是正核。因 \(\|I-G\|_{1\to1}\le2\)，
\[
\int_{t_0}^{t_0+c_0r}
\|(L_s-L_{t_0})U_s\|_1\,ds
\le2nr^2\|\nu_0\|_1.
\tag{10}
\]
这比对变化率直接作粗估计好，但仅因 free product 的精确消去成立；不能从 \(\eta_s\le U_s\) 推出带符号项同样支配。

对原正 killed/live 演化 \(\eta_s\)，完整全空间式为
\[
\partial_s\eta_s+L_s\eta_s=-\rho_s,
\quad \rho_s=V_s\eta_s+\text{整个外域 jump 落点},\quad
\rho_s\ge0.
\tag{11}
\]
同一初始 source/domain/物理尺度的 NB 合同给
\(\|\eta_s\|_1\le W_0\)、\(\int\|\rho_s\|_1ds\le W_0\)；外部初始 good 部分须按原合同另拆，不能遗漏。
真正通用的差算子因式分解为
\[
L_s-L_{t_0}
=\frac{s-t_0}{c_sc_0}
\sum_i(I-G_{c_s}^{(i)})(I-G_{c_0}^{(i)}).
\]
所以在一个 specified 窗口内，forcing 总变差满足
\[
\begin{aligned}
\int\|(L_s-L_{t_0})\eta_s\|_1ds
&\le4n[-r-\log(1-r)]W_0\\
&\le\frac{2nr^2}{1-r}W_0.
\end{aligned}
\tag{12}
\]
窄窗 \(r\le r_n\) 给上界 \(\le W_0\)。这些是实际 rate/kernel 的带符号 forcing 预算，尚不是其传播后 maximal 的输出费用。

## 5. 全 early 可有一份 forcing 预算，但不是已付输出误差

不在每个窗口重启 source。对同一 \(\eta_s\) 在 \(0\le s\le1/2\)，预先按 \(c_j\) 几何划分：\(c_0=1\)，每步 \(c_{j+1}=\max\{1/2,(1-r_n)c_j\}\)，入口 \(t_j=1-c_j\)。定义 piecewise frozen \( \widetilde L_s=L_{t_j}\)。此划分覆盖 entire early，参数与输入无关。

在每格，
\[
\|L_s-\widetilde L_s\|_{1\to1}
\le\frac{4n(s-t_j)}{c_sc_j}
\le8nr_n.
\]
故在同一原 live occupation 上
\[
\boxed{
\int_0^{1/2}\|(L_s-\widetilde L_s)\eta_s\|_1ds
\le8nr_n\int_0^{1/2}\|\eta_s\|_1ds
\le4nr_nW
=\frac{2n}{m_n}W\le2\sqrt n\,W.}
\tag{13}
\]
**(13) 仅是 forcing 的时空总变差。** 它只用一次原初始质量和整个时间积分，未将 \(N_{t_j}\nu\) 在每格当作新源收费。可与 \(\int\|\rho_s\|_1ds\le W\) 并列记录，但不能命名为已付 output error。

令 \(\widetilde U(s,v)\) 为这些 piecewise frozen free generators 的 time-ordered 正传播（在共同物理尺度时它们虽可交换，仍不是单个 generator）。精确 Duhamel 是
\[
\eta_s=\widetilde U(s,0)\nu
-\int_0^s\widetilde U(s,v)
\{\rho_v+(L_v-\widetilde L_v)\eta_v\}\,dv.
\tag{14}
\]
新增 forcing 为带符号、派生于原来源的空间时间 measure；其源标签不能在新空间点重新解释 FIRST/父组资格。

固定 endpoint 下 \(L^1\) contraction 可检查 (14)，但 receiver 选择 \(s=s(x)\) 时需要控制完整右侧的 propagation/maximal。fixed theorem 只控制 \(\sup_\tau H_\tau^{c_j}f\)，不控制 \(\widetilde U\) 的全乘积。更不能把
\(\int \sup_\tau H_\tau^{c_j(v)}|F_v|\,dv\)
的条件弱常数积分后，当成其空间弱常数；弱 \(L^1\) 没有该 Minkowski。时空 forcing mass \(\le2\sqrt nW\) 不解除这一障碍。

就算在每格使用 (8)，输入 \(\eta_{t_j}\) 的质量会在多格反复存活，\(\sum_j\|\eta_{t_j}\|_1\) 不能由 \(W\) 控制；并且要再使用 \(O(\sqrt n)\) 个物理尺度网，不能把两个独立网各自费用隐藏为一次 \(O(\sqrt n)\) 付款。这正是本轮仍未付的 receiver 参数选择接口。

## 6. 查重后的最小未付接口

旧 5072 附近已有严格 capped-obstacle 反例：
\(A_0\int_0^1N_sds\) 的全维空间核在原点附近可为正，从而域外平均可为正。故把 \(N_s\) 当 \(e^{-sA_0}\)，或免费继承 fixed obstacle 的 ergodic sign，都不合法。本轮不复测该例。8353 后共同尺度帽反例否定的是特定共同饱和分解；不把它泛化为所有 c-window 方法不可能。根端另审不同 c 的矩约束时，也须保留其矩交换假设，本文不依赖尚未纳入的全一般共同 c 帽否证。

| 接口 | 已核范围 | 尚未授权的推广 |
|---|---|---|
| (5) 正支配 | 一个指定 c0、同一初始密度、相对长度 ≤r_n | 把整个 continuation 或所有入口合并成一次 fixed fee |
| (8) 局部弱费 | 同一入口/源；含尺度网时源须独立于尺度 | 逐窗或逐条件 c 的弱界作平均 |
| (10)/(12) | free/killed 的指定窗 forcing 时空总变差 | 将带符号 forcing 当正原 source 或已付选时输出误差 |
| (13) | 同一 early live occupation、一次原 W；forcing TV ≤2√nW | 当前 receiver-chosen 输出 maximal 或 actual traffic 的已付误差 |
| (14) | 正完整 time-ordered propagation 的带符号 Duhamel | 单个 fixed semigroup maximal theorem 自动适用 |

真正可回代的下一条需证明：对原共同来源、原完整门，(14) 的 piecewise propagation（或其实际 source-labelled 子项）在 receiver 参数选择后有可支付空间界；费用须按原 absorption/forcing 的一份时空 measure，而非每窗另领 \(\|\eta_{t_j}\|_1\)。若只得到弱水平集，还须用原完整行上界和阈值比较接回 actual traffic。物理尺度、硬平均 \(h_{L_s}\)、参数 \(c_t\in[1/2,1]\) 和原历史门均仍保留，不能在新落点重判。

本轮有用的精确增益是 (5) 的合法局部 bridge 与 (13) 的一次源 forcing TV；当前一般空间费用仍缺后续传播/选择，而不是再缺一个 frozen family endpoint。没有新 cube 主账费用。

## 7. 三轮原模型注册守卫

先保存同前缀 registration，再执行新 guard。\(n=8,32,128\)，\(c_0=1,3/4,1/2\)，相对窗口 \(r=0,r_n/2,r_n\)。每个非零 r 核全部掩码大小 \(j=0,\ldots,n\) 的原 Bernoulli/冻结 no-repeat Poisson 权比例，保持 holding 与全部重复跳的正剩余。原 Fourier symbol 在 \(\lambda=1/16,1/4,1,4,16,256\) 直接以 80 位 Decimal 的 \(\log(1+\lambda)\) 核原 cocycle、\(G_cN_s=w\) 和 free forcing (9) 的精确代数形式。

Fraction 核相对窗口、log 常数的解析上预算、free/killed forcing 质量、覆盖全 early \(c\in[1/2,1]\) 的几何网与 (13) 的一次 occupation 预算。Decimal 是高精度实现诊断，非 interval；正核及无穷 Poisson 剩余的证明由 (5)–(6) 负责。没有 reset 替代、随机来源、阶数拟合或实际 FIRST/CPGP 认证。

终态一次运行 **PASS_EXACT_AND_DECIMAL**，1947 个谓词；注册/脚本 SHA256 保存于结果。只新增本轮实验，旧固定/尺度/共同帽数据未重跑。可证局部支配与 forcing TV 保持；完整原过程转移尚未付，不改 actual 总账。

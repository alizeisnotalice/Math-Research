# 原 ordered winner 的定向响应差收费与小 contrast 接口

2026-10-07。独立新稿；不再修改已审 [capped test-energy 稿](ordered_capped_test_energy_20261007.md) 及其脚本/hash。固定物理尺度、固定入口 \(c>0\)，任意非负 \(f\in L^1(\mathbb R^n)\)，\(W=\int f\)，有限 \(Z\)。本文仅研究完整原 ordered 最大化的接口；不自动支付 actual geom/FIRST/历史门，不进 cube 主账。可付结论是已证有界熵预算的应用；一般小 contrast 收费仍未证明。

## 1. 保留 receiver winner 的原量

令
\[
v_z=\prod_i[e^{-z}I+(1-e^{-z})G_{c,i}]f,\qquad
J_z=\sum_i(G_{ce^{-z},i}-I),\quad v'_z=J_zv_z.
\]
各 \(G_d\) 是原正对称保质量轴核，乘子
\[
\widehat G_d(\xi)=\frac1{1+d[\xi^2/\log(1+\xi^2)-1]}.
\]
取 \(M=\max_{0\le z\le Z}v_z\)、\(E=\{M>\lambda\}\)、可测原最大点 \(\tau(x)\)。任意可测 \(S\subset E\) 上，
\[
\theta=\mathbf1_S\lambda/M,\quad
H_z=\theta\mathbf1_{\{z<\tau\}},\quad
t_z=H_z(1+v_z/\lambda),\qquad0\le t_z<2.
\tag{1}
\]
原 flow、熵与 commutator 是
\[
q_S:=\lambda|S|,\qquad
F_S=q_S-\int\theta f=\int_0^Z\!\int H_zJ_zv_z,
\]
\[
D=\frac\lambda2\int d\mathfrak m\,
\frac{(\Delta v)^2}{(\lambda+v(x+h))(\lambda+v(x))}
=\int\Phi_\lambda(f)-\int\Phi_\lambda(v_Z)\le W,
\tag{2}
\]
\[
C_S=-\frac\lambda2\int d\mathfrak m\,\Delta t\,
\Delta\log(1+v/\lambda),\qquad
q_S\le\int\theta f+2D+C_S.
\tag{3}
\]
这里 \(d\mathfrak m=\sum_i dz\,dx\,G_{ce^{-z},i}(dh)\)，\(\Delta\) 为 \(x+h\) 减 \(x\)，\(\Phi_\lambda(v)=v-\lambda\log(1+v/\lambda)\)。它不是 raw \(f\) 能量或抽象 independent test；所有 gate、winner 仍来自同一 \(f\)。有限时间绝对可积由 bounded tests、\(\lambda\log(1+v/\lambda)\le v\)、\(\|J_zv_z\|_1\le2nW\) 保证。

## 2. 大 regularized contrast 同边子支已经可付

给定 \(\delta>0\)，设
\[
A_\delta=\left\{(z,x,h,i):
\frac{\lambda+\max(v_z(x),v_z(x+h))}
{\lambda+\min(v_z(x),v_z(x+h))}\ge1+\delta\right\},
\qquad N_\delta=A_\delta^c.
\tag{4}
\]
这是对称的真实同边事件，包括 equality 的大 contrast 支。定义
\[
D_A=\frac\lambda2\int_{A_\delta}d\mathfrak m\,
\frac{(\Delta v)^2}{(\lambda+v(x+h))(\lambda+v(x))}
\le D.
\]
把大支 commutator 的正费用定向到 \(u=v(x+h)>v=v(x)\)、\(t(x)>t(x+h)\)。记
\[
a=\frac{\lambda+u}{\lambda+v},\qquad
Q=\frac{(u-v)^2}{(\lambda+u)(\lambda+v)}
=\frac{(a-1)^2}{a}.
\]
该支的**未归一化**正积分是
\[
C_A^+=\lambda\int d\mathfrak m\,
\mathbf1_{\{a\ge1+\delta,\ u>v,\ t(x)>t(x+h)\}}
[t(x)-t(x+h)]\log a.
\tag{5}
\]
式(5)的系数是 \(\lambda\)，不是 \(\lambda/2\)：反向边给相同正费用。对 \(a\ge1+\delta\)，
\[
\log a\le a-1
\le\frac{1+\delta}{\delta}\frac{(a-1)^2}{a}.
\]
又 \(t(x)-t(x+h)\le2\)。对称性使 \(D_A=\lambda\int_{\{u>v\}\cap A_\delta}Q\,d\mathfrak m\)，因此
\[
\boxed{C_A\le C_A^+
\le\frac{2(1+\delta)}{\delta}D_A
\le\frac{2(1+\delta)}{\delta}W.}
\tag{6}
\]
其中 \(C_A\) 保留其 signed 定义。负方向可以使它更小；没有把每条边/每个时间重新领取一份 \(W\)。例如 \(\delta=1\) 的系数是4。

## 3. 剩余准确 signed 与正定向接口

定义真实小支
\[
C_N=-\frac\lambda2\int_{N_\delta}d\mathfrak m\,
\Delta t\,\Delta\log(1+v/\lambda).
\tag{7}
\]
它仍可正可负。对应的正定向上界是
\[
\boxed{C_N^+=\lambda\int d\mathfrak m\,
\mathbf1_{\{\,v(x+h)>v(x),\ 1<a<1+\delta,\ t(x)>t(x+h)\,\}}
[t(x)-t(x+h)]\log a.}
\tag{8}
\]
精确无遗漏地 \(C_S=C_A+C_N\)，且 \(C_N\le C_N^+\)。\(a=1\) 没有费用；\(a=1+\delta\) 放在已付大支；strict 端点不造成额外集合。

故 full \(S=E\) 可保留更准确的账
\[
q_E\le W+2D+\frac{2(1+\delta)}{\delta}D_A+C_N
\le(5+2/\delta)W+C_N.
\tag{9}
\]
gap \(S=\{f\le\lambda/2,M>\lambda\}\) 则
\[
q_S\le(8+4/\delta)W+2C_N.
\tag{10}
\]
若先付初始高集 \(\{f>\lambda/2\}\)，其费用为 \(2W\)；不得再把它遗漏或无说明重复领取。

可选 \(\delta=1/\log(n+2)\)，式(6)及(9)的大支费用是 \(O(\log(n+2))W\)。这只是费用分配，既不证明 optimal 阶数，也不证明 \(C_N^+\) 或 \(C_N\) 已付。由 \(\log a<\delta\) 直接估 jump 总活动会保留 \(nZ\) 因子，不能称为吸收。当前最小未付接口可以写为
\[
(C_N)_+\le B_n W
\quad\text{或更强 } C_N^+\le B_nW,
\tag{11}
\]
其中 \(B_n=\mathrm{polylog}(n)\) 是待证目标，所有边、gate、\(M,\tau\) 必须仍是(1)的原对象。只需 signed 版本，不必先证明正定向强版本。

## 4. 原 winner 的同边符号分类

准确乘积差分是
\[
\Delta t=\frac{H(x)+H(x+h)}{2\lambda}\Delta v
+\left(1+\frac{v(x)+v(x+h)}{2\lambda}\right)\Delta H.
\tag{12}
\]
第一项对 commutator 的贡献非正。如果 \(\Delta H=0\)，整条边是免费负支。

两端都活跃时，\(t=(\lambda+v)/M\)。在 \(u=v(x+h)>v(x)\) 的方向，adverse test 顺序恰为
\[
t(x)>t(x+h)
\quad\Longleftrightarrow\quad
\frac{M(x+h)}{M(x)}>
\frac{\lambda+u}{\lambda+v}.
\tag{13}
\]
若两端活跃且 \(M\) 相等，则不是正费用。若高响应端的 \(t\) 为0，则低响应端必须仍属于 \(S\)、\(z<\tau(x)\)，而另端可能已过自己的 winner 或从未属于 \(S\)。这些情况不能统称为同一“过去退出”。

式(12)保留一个重要取消：**不能先删第一项，再沿用 \(\Delta t\le2\) 对余项收费**。余项单独可很大；(6)是对完整 \(\Delta t\) 的上界。未来 winner 只限制 \(v_z\le M\)，没有保证 \(\Delta H\) 与 \(\Delta v\) 的符号，不能凭(13)断言 source-once 空间费用。

小支还可量化“取消以后回到原 flow”这一现象。令 \(a_0=\lambda+u,b_0=\lambda+v\)，
\[
R(a_0,b_0)=(a_0+b_0)\log(a_0/b_0)-2(a_0-b_0).
\]
定义
\[
F_N=-\frac12\int_{N_\delta}d\mathfrak m\,\Delta H\Delta v,\quad
P_N=\frac14\int_{N_\delta}d\mathfrak m\,
[H(x)+H(x+h)]\Delta v\,\Delta\log(1+v/\lambda)\ge0,
\]
\[
R_N=\frac14\int_{N_\delta}d\mathfrak m\,\Delta H\,R(a_0,b_0).
\]
则准确
\[
C_N=F_N-P_N-R_N.
\tag{14}
\]
取 \(s=(a_0-b_0)/(a_0+b_0)\)。由
\(\operatorname{artanh}s-s=\sum_{k\ge1}s^{2k+1}/(2k+1)\)，
\[
|R(a_0,b_0)|\le\frac{|a_0-b_0|^3}{6a_0b_0}.
\]
另一方面 \(H\le\min(1,\lambda/v)\le2\lambda/(\lambda+v)\)。在 \(N_\delta\)，设 \(r=\max(a_0,b_0)/\min(a_0,b_0)<1+\delta\)，则 \(r-1/r\le2\delta\)。这些给
\[
|R_N|\le\frac{\delta}{3}D_N,\qquad
D_N=D-D_A.
\tag{15}
\]
所以小支的非线性余项可由熵支付；主项 \(F_N\) 却仍是同一 receiver 门下的原小 contrast flow，没有来源一次支付。这是准确 sign bookkeeping，不是新的空间闭合，也不能把 \(F_N\) 当成独立可收费来源。

## 5. \(H\) 时间单调不提供免费 Markov 停时

对每个固定 receiver，\(H_z=\theta\mathbf1_{z<\tau}\) 随 \(z\) 单调下降。但 \(\tau\) 依赖完整未来确定性响应 \(v_\cdot(x)\)；沿原 Markov 路径改变 receiver 后，它不是自动可预测，也不是已给超鞅。

以下 bounded primitive 核对说明单调时间分部积分到底保留了什么。令
\[
e_z=J_z\Phi_\lambda(v_z)-\Phi_\lambda'(v_z)J_zv_z\ge0,
\qquad\int_0^Z\int e_z=D\le W.
\]
逐 receiver 积分到其原 \(\tau\)，再用对称性，准确得到
\[
\boxed{\int_0^Z\!\int\Phi_\lambda(v_z)J_zH_z
=\int\theta[\Phi_\lambda(M)-\Phi_\lambda(f)]
+\int_0^Z\!\int H_ze_z.}
\tag{16}
\]
所有项可积，\(\Phi_\lambda(v)\le v\)、\(H\le1\) 即足够；不需 \(\Phi\) 无界 entropy 的 \(L\log L\) 假设。式(16)亦可直接由 \(\int H_z\partial_z\Phi(v_z)\) 的基本微积分证明，避免把单点时间 gate 强行微分成未说明的随机补偿。

终端项满足
\[
(1-\log2)q_S\le\int\theta\Phi_\lambda(M)\le q_S.
\tag{17}
\]
因为 \(M/\lambda>1\) 且 \(1-\log(1+r)/r\) 对 \(r>0\) 递增；下界在 \(|S|>0\) 时严格，空集仍保留非严格式。初始项 \(\int\theta\Phi(f)\le W\)，正缺陷项 \(\int H e\le D\le W\)。因此简单 adjoint/time 分部积分并没有消灭 receiver 水平集费用：未知的 \(q_S\) 以正终端项留下。不能把左侧 \(\int\Phi(v)JH\) 因 \(H\) 的时间单调而称为负或已有来源费用。

式(16)不排除进一步利用真实核及(13)的空间取消；它只排除这一未经证明的符号捷径。一般(11)仍待证，完整 future \(Z\to\infty\) 的统一预算也未建立。

## 6. 排除全空间时间配对，保留同边问题

另代理的 [空间错配稿](paired_entropy_spatial_mismatch_20261007.md) 已只读核其范围：同一原 \(c=1\) 输入由远移的两个有限大盒组成，一处保留真实 high-frequency winner test、另一处保留单轴低频熵。正边积分与 winner 稳定的欧氏提升给全空间时间配对 \(\int\sqrt{d(z)a(z)}dz\ge\sqrt nW/2048\)，初始 gap test 版本亦有 \(\sqrt nW/4096\)。这里是先全空间积分后 Cauchy 的错配，不能解释成 signed commutator 的下界，也不是 actual geom 反例。

因此本稿不继续一般全空间时间配对的 polylog 猜想。未付方向限定为原同边 signed \(C_N\)，或真正保留空间/坐标/jump 的配对；不把 fast test 与远处 slow 熵人工交叉。该他稿的900项有理收据不计入本文新守卫，也未由本作者重跑。

## 7. 专属纯有理常数收据

先写 [注册](ordered_oriented_small_contrast_registration_20261007.json)，再执行 [脚本](ordered_oriented_small_contrast_guard_20261007.py) 一次，输出 [结果](ordered_oriented_small_contrast_results_20261007.json)。三轮细化参数 \(N=8,32,128\)，\(\delta=1,1/2,1/8\)，比值网同时覆盖 equality、大差及趋近1；再核 capped node 的差分、同 \(H\) 符号、两端 alive 的最大值比值分类及 near remainder 的有理 envelope。所有数是 Fraction，不运行旧原核实验、不生成新模型。

终态 exit0：三轮17040、19848、31080项，共67968项 exact 检查通过。注册 SHA256 为 \(\texttt{686e5e085dfe5ca7a701dbe1f596e3b8714f113270615bc88784ce5eb2a75688}\)，脚本 SHA256 为 \(\texttt{babe2e2336a6ba835767d6b5ade34526934c297da20e243476879a59196226cb}\)。

有理守卫认证的是 \(\log r\le r-1\) 之后的比值常数及精确代数；log/atanh 不等式由上面的解析证明负责，没有把浮点 log 当 interval 证书。守卫中的 capped 双节点不宣称能被同一实际输入同时实现；它们只校核必要 scalar envelope。原输入近并列及 gap 压力已在前稿保留，不重跑。大 contrast 子支是合法熵应用，小支及原 cube 主目标仍未付。

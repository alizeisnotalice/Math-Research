# R033 · 典型外壳点、变量前序区与增长扣减

本报告独立检查并闭合主协调者提出的一般高维构造。结论是**通过**，但后序集必须明确使用 `1≤||w||∞≤s`，而不是自动成立的 `1≤m(w)`。证明采用显式 κ 与有限 a guard，不依赖连续选择器或未量化的紧性。没有七维特异性搜索，没有 Cn 下界研究。

## 1. 定理与依赖

记所有对数为自然对数，

\[
L=\log n\ge4096,\quad y=\sqrt L,\quad l=\log L,
\quad x=L-3l,\quad t=n/x,\quad\beta=1+1/t,
\]

\[
A_n(t)=\beta^n,\qquad
B_n(t)=1-\left(\frac\beta2\right)^n+n\log t
-\frac n2\int_1^t\ell_t(\rho)^{n-1}\frac{d\rho}{\rho},
\quad\ell_t(\rho)=1+\frac1{2t}-\frac1{2\rho}.
\]

**定理。** 令 δ=n/(128√L)，则

\[
C_n\le\left(\sqrt{A_n(t)}+
\sqrt{B_n(t)-\delta}\right)^2. \tag{1}
\]

它严格改进同参数旧公式超过 δ。结合 R030 已独审的上界

\[
F_n(t)\le W_n
=T+(3l+1)n/L+(9l^2+1)n/L^2,
\quad T=nL-nl,
\]

得到

\[
C_n<W_n-\frac{n}{128\sqrt L}. \tag{2}
\]

特别地，对 L≥2^30，

\[
C_n\le n\log n-n\log\log n
-\frac{n}{256\sqrt{\log n}}. \tag{3}
\]

依赖只包括：已经证明的原有限覆盖/负载架构、其角点包络，以及 R030 的显式 Wn 上界。下面重新给出所需的新几何和所有有限 a 连接；不使用单一 unit-core cut 的结论。前导系数仍为 1，(3) 改进的是其后的负修正项。

## 2. 原有限实例与变量前序区引理

固定 n,t，令 a>2t、α=a/(a+1)、c=(aβ−1)/(a+1)。在任意负载支撑点处作原归一化：支撑点为 0，最后一个最小正贡献立方体 Qk 的半径为 1。其中心归一化为 u=yk/(a+1)，反射坐标使 0≤ui<α。

相对末项未捕获的全部前序实际增量给出测度

\[
d\mu(z)=\alpha^{-n}2^{-n}\sum_jr_j^{-n}
\mathbf1_{I_j/(a+1)}(z)\,dz,
\qquad m(z)=\max(1,\|z\|_\infty).
\]

总质量 U 恰为这些项的原权重和。真实增量不交与位置半径耦合给出 μ≤λa，其中

\[
d\lambda_a(z)=\alpha^{-n}2^{-n}m(z)^{-n}
\mathbf1_{H_a}(z)\,dz,
\]

\[
H_a=\{\|z\|_\infty<t:\exists i,\ z_i\le\alpha-cm(z)\},
\qquad\lambda_a(H_a)=\alpha^{-n}J_a(t).
\]

实际末中心的包络包含于 Ha，是因为 zi−ui≤m(z)<cm(z)，所以未捕获必须发生在负方向。

令极限角点域为

\[
D_\infty=\{\|z\|_\infty<t:\exists i,\ z_i<1-\beta m(z)\},
\qquad d\lambda_\infty=2^{-n}m(z)^{-n}\mathbf1_{D_\infty}(z)\,dz.
\]

在所有 m(z)≥1 上，

\[
\alpha-cm(z)-[1-\beta m(z)]
=\frac{(\beta+1)m(z)-1}{a+1}>0.
\]

故 D∞⊂Ha。对任意 Borel 集 E⊂D∞，有准确等式

\[
\lambda_a(E)=\alpha^{-n}\lambda_\infty(E). \tag{4}
\]

无需 TV 收敛或集合选择器连续性。

**变量前序区引理。** 设一个后序集 A⊂D∞ 满足 λa(A)≥δ。对每个 w∈A，有一个 Borel 前序区 P(w)⊂D∞，满足 λa(P(w))−τ≥p≥δ。若任意在 A 中有正实际增量质量的项 j，选择其任意 w∈Ij/(a+1)∩A 后，所有在 P(w) 中有实际质量的项均半径严格大于 rj、并由 bQi 捕获 yj，则

\[
U\le\alpha^{-n}J_a(t)-\delta. \tag{5}
\]

**证明。** 设总缺额 Δ=λa(Ha)−U<δ。则 μ(A)≥λa(A)−Δ>0，有限原族中至少有一个真实项 j 与 A 有正增量交集。固定该 j 和一个实际增量点 w。由前序与捕获假设，同一个 j 的一次 eligibility 给 μ(P(w))≤τ。但

\[
\mu(P(w))\ge\lambda_a(P(w))-\Delta>\tau,
\]

矛盾。没有对 w 积分，没有假设所有 w 共用同一 P，也没有使用可测选择器。证毕。

## 3. 典型外壳后序集的质量

定义

\[
s=1+\frac1{16y},\quad
\eta=1-\frac1{2y},\quad
\eta'=1-\frac1{4y},\quad
d=s-\eta=\frac9{16y},
\]

\[
S(w)=\sum_{i=1}^n(|w_i|-\eta)_+,
\quad K(w)=\#\{i:|w_i|>\eta\},
\]

\[
\mathcal A=\{w:1\le\|w\|_\infty\le s,
\ \exists i\ w_i=-\|w\|_\infty,
\ S(w)\le4nd^2,\ K(w)\le8nd\}.
\]

该集是闭有界集。K 是下半连续整数值函数，所以它的上述 sublevel 是闭集；负最大面条件是有限个闭集的并。由于 s<t，每个负最大坐标满足 1+||w||∞>β||w||∞，故 A⊂D∞。

固定壳半径 r∈[1,s] 和一个负面 wi=−r，其余坐标均匀分布于 [−r,r]。线性期望给

\[
\mathbb E S\le d+\frac n2d^2\le nd^2,
\qquad
\mathbb E K\le1+nd\le2nd.
\]

这里 nd≥2 即够，而下面的粗阈值保证它远大于 2。Markov 给 S、K 两个坏事件各≤1/4，故每个负面至少一半面积是好点。负面本身占全部 2n 个面的半数，因此 A 至少占所有壳面积的 1/4。排除最大坐标并列的零测集合后，

\[
\lambda_\infty(\mathcal A)
\ge\frac n4\log s
\ge\frac{n}{128y}=\delta. \tag{6}
\]

最后使用 log(1+v)≥v/2（0≤v≤1）。由 (4)，有限 a 还有 λa(A)=α^-nλ∞(A)≥δ。

## 4. 每个典型点接受的大前序区

对 w∈A，定义

\[
P_\kappa(w)=D_\infty\cap
\{\eta t<\rho=\|z\|_\infty<\eta't,
\ \|z-w\|_\infty<\kappa\rho\}.
\]

记 Φκ(w)=λ∞(Pκ(w))。本节先在 κ=β 下下界，随后给显式 κ<β 的误差。

### 4.1 局部捕获的完整面系数

在固定 ρ 壳面上，ρ≫s，所以局部捕获盒 w+(−βρ,βρ)^n 在第 i 坐标与 [−ρ,ρ] 的归一化交长为

\[
\ell_i^\beta
=1-\frac{(|w_i|-\rho/t)_+}{2\rho}.
\]

因为 ρ/t>η，

\[
\sum_i(1-\ell_i^\beta)
\le\frac{S(w)}{2\eta t}
\le\frac{2nd^2}{\eta t}
=\frac{81}{128\eta}\frac xL
\le\frac{81}{127}<\frac23. \tag{7}
\]

所以任何删去一个坐标的交长乘积都大于 1/3。每个坐标至少有一个面被捕获，而丢失第二个面只可能发生在 |wi|≥ρ/t>η。因此捕获面数至少为

\[
2n-K(w)\ge2n-8nd\ge15n/8.
\]

壳面积归一化系数为每面 1/2 乘其余交长乘积，故完整立方体的局部捕获系数至少 5n/16。

### 4.2 扣掉原角点捕获区

原角点捕获区在该壳上只有 n 个正面，其系数为

\[
\frac n2\ell_t(\rho)^{n-1}.
\]

对 ρ≤η't，

\[
\frac{n-1}{2}\left(\frac1\rho-\frac1t\right)
\ge\frac{(1-1/n)x}{8y\eta'}
\ge\frac{9y}{128}\ge\frac92.
\]

故 ℓ_t(ρ)^(n−1)≤exp(−9/2)<1/50，原系数<n/100。原角点捕获区与局部捕获区的交不大于整个原区，因此 D∞ 中的局部捕获系数至少

\[
5n/16-n/100=121n/400>n/4.
\]

积分得所有 w∈A 的统一下界

\[
\Phi_\beta(w)
\ge\frac n4\log(\eta'/\eta)
\ge\frac{n}{16y}=8\delta. \tag{8}
\]

### 4.3 显式 κ 与捕获损失

取

\[
\zeta=n^{-3},\qquad\kappa=\beta-\zeta.
\]

有 ζ≤(β−1)/2，因此 1<κ<β。完整立方体的第 i 个归一化坐标交长为

\[
\ell_i^\kappa=
\min\left(1,\frac{1+\kappa}{2}-\frac{|w_i|}{2\rho}\right),
\]

对 κ 的 Lipschitz 常数≤1/2，始终处于 (0,1]。从 β 缩小到 κ 后，交长乘积最多损失 (n−1)ζ/2。

每个坐标至多丢一个面。丢面发生的 ρ 集合包含于

\[
\frac{|w_i|}{\beta-1}\le\rho\le
\frac{|w_i|}{\kappa-1},
\]

当 wi=0 时不丢面；否则其 logρ 长度≤log[(β−1)/(κ−1)]≤2tζ。每个面的归一系数≤1/2，所以全部丢面的积分损失≤ntζ。存留面的产品变化积分损失≤n(n−1)ζ/2·log(η'/η)≤n²ζ/2。由于 t<n，

\[
0\le\Phi_\beta(w)-\Phi_\kappa(w)
\le nt\zeta+n^2\zeta/2
\le\tfrac32n^2\zeta=\frac3{2n}
\le\delta. \tag{9}
\]

这里对完整立方体捕获区计损失。再与固定 D∞ 相交只会减小损失，不会放大；故 (9) 对 Φ 也成立。得到 Φκ(w)≥7δ，对所有 w∈A 统一有效。无需 compact-uniform continuity。

## 5. 显式有限 a 捕获与原算子不等式

选择

\[
a\ge8n^4.
\]

则 a>2t，β<2，且

\[
\alpha\beta-\kappa-2/(a+1)
=\zeta-\frac{\beta+2}{a+1}>0. \tag{10}
\]

若某真实未捕获项 j 的增量含有 w∈A，则

\[
cr_j\le\|u-w\|_\infty\le\alpha+s,
\qquad r_j<\alpha+s<3,
\]

因为 c>1、s<2。若某实际未捕获项 i 的增量含有 z∈Pκ(w)，则 ri≥||z||∞>ηt>4，故 ri>rj，半径优先规则强制 i<j。中心 vi=yi/(a+1)、vj=yj/(a+1) 满足

\[
\|v_i-v_j\|_\infty
<\|z-w\|_\infty+(r_i+r_j)/(a+1)
<\bigl(\kappa+2/(a+1)\bigr)r_i
<\alpha\beta r_i.
\]

所以同一个 j 接受全部 Pκ(w) 实际项，eligibility 给 μ(Pκ(w))≤τ。

取 τ=√(An(t)[Bn(t)−δ])。下面的粗界保证 Bn>δ，且

\[
\tau\le n/L\le n/(64y)=2\delta.
\]

用 (4)、(8)–(9)，对每个 w∈A，

\[
\lambda_a(P_\kappa(w))-\tau
=\alpha^{-n}\Phi_\kappa(w)-\tau\ge5\delta.
\]

变量前序区引理遂给，对所有 a≥8n^4 的全部原有限实例，

\[
U\le\alpha^{-n}J_a(t)-\delta.
\]

再加末中心捕获预算 τ 和自身项 a^-n，接入原覆盖公式，

\[
C_n\le
\left(1+\frac{(\beta+1/a)^n}{\tau}\right)
\left[\tau+a^{-n}+\alpha^{-n}J_a(t)-\delta\right].
\]

固定 n,t,τ，对这一组对全部原实例成立的数值上界令 a→∞，得到 (1)。不要求选择历史收敛，也未将 κ=β 直接用于一个无 guard 的有限算法。

## 6. 阈值、预算正性与最终渐近修正

L≥4096 时，x≥3L/4>1：函数 L/4−3logL 此后递增，在 L=4096 时用 log4096=12log2<12 即得正值。又 x<L、n=e^L≥max(2L,L²/2)，所以 t≥n/L≥L/2、t<n、s<t、ηt>4、nd≥9y/8>2、n≥4、n²≥192y。还有 η≥127/128。这些阈值正好支撑 (7)–(10)。

旧积分上界 I_n(t)≤2t/(n−1)，非负核心及 x∈[L/2,L] 给

\[
B_n(t)\ge T-4n/L>8n.
\]

而 Bn≤1+nlogt≤nL，An≤e^x=n/L³。因此 Bn>δ 且上述 τ≤n/L 成立。

展开两个根式的准确差，

\[
F_n(t)-\left(\sqrt{A_n(t)}+\sqrt{B_n(t)-\delta}\right)^2
=\delta+2\sqrt{A_n(t)}
\left[\sqrt{B_n(t)}-\sqrt{B_n(t)-\delta}\right]>\delta.
\]

接入已审 Wn，得到 (2)。这里没有把新公式族的下界当成真实 Cn 的下界。

最后令

\[
R(L)=\frac{3\log L+1}{L}+
\frac{9(\log L)^2+1}{L^2}.
\]

√L·R(L) 的两个加项在 L≥2^30 上分别严格递减，导数符号来自 5−3logL<0 和 18logL−(3/2)[9(logL)^2+1]<0。在 L0=2^30，有 logL0=30log2<21，所以

\[
\sqrt{L_0}R(L_0)
<\frac{64}{2^{15}}+\frac{3970}{2^{45}}
<\frac1{256}.
\]

由此 Wn<T+n/(256√L)，与 (2) 合并得 (3)。log2<7/10 可由 exp(7/10)>1+7/10+(7/10)^2/2+(7/10)^3/6>2 直接证明。

## 7. 实际检查与范围

`check_thresholds.py` 使用本机 Python 标准库 Fraction 实际检查关键有理阈值：81/127<2/3、面系数差>1/4、exp(9/2) 的正项截断>50、log2 上界所需的截断，以及 L0 处最终余项包围。结果保存在 `threshold_checks.json`。这些有限计算只核代数和阈值基点；全维数、全真实历史的结论由本报告解析证明给出。

仅写 `upper/round033/typical_shell/`。R030–R032 保持冻结；没有参数优化、共享 MATLAB、七维细调、commit/push/fetch/rebase。新的增长扣减来自真实外壳后序点与随点变化的前序区，明确超出单 unit-core cut 的小于 2 障碍。


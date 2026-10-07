# 原 resolvent 的 Stieltjes 密度与斜向 Lévy 正增量

2026-10-07。独立解析审计；只新增本稿，不修改实际主账。使用 J01 的真实时变生成元合同，已读 SKILL、method、cube-interface；F02 的谱正性/空间正测度区分，已读 SKILL、method。本稿自含支切积分证明，不借未核自分解定理。radial agent 另负责新有限守卫，本稿不运行或重跑其数值。

主结论：对原 \(c=1\) soft square root，曲线
\[
 \frac{r(L)}{1-r(L)}=C L^2,\qquad C>0
\tag{1}
\]
确实给**原 Lévy 测度序**，所以有限步商是真正 Markov convolution。每轴沿 \(\log L\) 的跳率为 \(r\le1\)。这严格区别于 fixed-r 跨 L 的不可能商。但条件首跳后的 \(c_t<1\)、任意接收 selector、原 FIRST/history 与 count 尾预算尚未由此解决。

## 1. 约定与没有漏掉的 pole

Fourier 约定 \(\widehat f(\xi)=\int e^{-ix\xi}f(x)dx\)。令
\[
 f(z)=\frac{\log(1+z)}z=\int_0^1\frac{dv}{1+vz},\quad
 G_a(z)=\frac{f(z)}{a+(1-a)f(z)}
 =\frac{\log(1+z)}{az+(1-a)\log(1+z)},\quad 0<a\le1.
\tag{2}
\]
取主支 log，域 \(\mathbb C\setminus(-\infty,-1]\)。z=0 的值为1，是可去点。上半平面 \(\operatorname{Im}f<0\)，故 \(a+(1-a)f\) 不为0；当 a=1 它恒为1。下半平面由共轭；实轴 z>−1 的 f>0，也无零点。因此 (2) 没有域内 pole。特别 \((-1,0)\) 上不能凭未消去的分母在0取零而添加原子。

在 z→−1 的 branch 端，a<1 时 \(G_a\to1/(1-a)\)，a=1 时至多 log 发散；均不给 pole。无 cut 外负实 pole，也无额外 exponential 项。固定 a>0 时大圆上 \(G_a=O_a(\log|z|/|z|)\)；这是固定 a 的估计，不声称 a→0 一致。

## 2. cut jump、空间密度与归一化

写 \(t>1,\ \ell=\log(t-1)\)，上侧 cut 有
\[
 G_a(-t+i0)=\frac{\ell+i\pi}{-at+(1-a)\ell+i(1-a)\pi},
 \qquad
 -\frac1\pi\operatorname{Im}G_a(-t+i0)
 =\rho_a(t):=\frac{at}{[at-(1-a)\ell]^2+[(1-a)\pi]^2}>0.
\tag{3}
\]
绕 cut 的 Cauchy 公式给
\[
 G_a(q)=\int_1^\infty\frac{\rho_a(t)}{q+t}\,dt,\qquad q\ge0.
\tag{4}
\]
具体收敛检查：无 pole；大圆贡献趋0；−1 小圆贡献在 a<1 时 \(O_a(\varepsilon)\)、a=1 时 \(O(\varepsilon|\log\varepsilon|)\)。cut 密度在无穷为 \(O_a(1/t)\)，所以 (4) 可积；在 t↓1，a<1 时为 \(O_a(|\log(t-1)|^{-2})\)，a=1 时为1/t。先固定 q>0 作 contour，再用正单调极限到 q=0。由 \(G_a(0)=1\)，严格
\[
 \int_1^\infty\rho_a(t)\frac{dt}t=1.
\tag{5}
\]
反 Fourier \(1/(\xi^2+s^2)\leftrightarrow e^{-s|x|}/(2s)\)，变量 t=s² 给
\[
 \boxed{g_a(x)=\int_1^\infty
 \frac{as^2 e^{-s|x|}}
 {[as^2-(1-a)\log(s^2-1)]^2+[(1-a)\pi]^2}\,ds.}
\tag{6}
\]
正 Tonelli 及 (5) 证明 \(\int g_a=1\)，\(\widehat g_a(\xi)=G_a(\xi^2)\)。固定 a>0 时 g_a(0)<∞，连续、偶、严格正，并是 exponential 正混合。候选 (6) 的分子 s² 与全部常数均正确。

a=1 给 \(g_1(x)=\int_1^\infty e^{-s|x|}s^{-2}ds\)，正是原 w。a→1 可在 a≥1/2 下以 \(Cs^{-2}\) 控制 (6)，从而 pointwise 与 L1 收敛。a→0 则 \(G_a(q)\to1\)，概率测度 \(g_a(x)dx\Rightarrow\delta_0\)；不能把 a=0 直接代入 (6) 的零分子后声称零 kernel。弱收敛可不用引用特征函数定理：对任意 T，(5) 与 \(1-G_a(T^2)\to0\) 控制 \(t\le T^2\) 的质量，再由 exponential 混合的尾概率 \(\exp(-s\varepsilon)\) 得 \(\Pr(|X|>\varepsilon)\to0\)。这里不额外声称 a→0 的全空间 density 一致极限。

## 3. 精确 phase，原 Lévy 测度

取 \(0<r<1,\ A=1-r\)，原每轴 P 的有限 Lévy density 为
\[
 \Lambda_r(x)=\frac12\int_A^1g_a(x)\frac{da}a,\qquad
 |\Lambda_r|=-\frac12\log A.
\tag{7}
\]
这是 \(P_{r,L}\) 的 compound-Poisson 表示：\(\widehat P_{r,L}
=\exp\{\widehat\Lambda_{r,L}-|\Lambda_r|\}\)，\(\Lambda_{r,L}(x)=L^{-1}\Lambda_r(x/L)\)。

固定 s>1，令 \(t=s^2,\ell=\log(t-1)\)、
\[
 u_r(s)=\frac{1-r}{r}s^2-\ell,\qquad
 F_r(s)=\frac1{2\pi}\operatorname{atan2}(\pi,u_r(s))\in(0,1/2).
\tag{8}
\]
由 \(u_a=as^2/(1-a)-\ell\) 与 \(du_a/da=s^2/(1-a)^2\)，正积分直接给
\[
 \frac12\int_A^1
 \frac{s^2\,da}{[as^2-(1-a)\ell]^2+[(1-a)\pi]^2}
 =\frac1{2\pi}\left[\frac\pi2-\arctan\frac{u_r(s)}\pi\right]
 =F_r(s).
\]
故无近似地
\[
 \boxed{\Lambda_r(x)=\int_1^\infty F_r(s)e^{-s|x|}ds.}
\tag{9}
\]
F_r(s) 在 s↓1 趋0，s→∞ 为 \(O_r(s^{-2})\)，所以 \(\Lambda_r(0)<∞\)。r↓0 时 Λ 的总质量趋0。r↑1 时 phase 有点态极限，且 \(\Lambda_1(x)\le e^{-|x|}/(2|x|)\) 在 x≠0；总质量无限，但 \(\int(1\wedge x^2)\Lambda_1(x)dx<∞\)。本文有限步正商证明只需 r<1；不为有限时间强行到达 r=1。

## 4. 斜向充分条件：保留原 measure，不仅 symbol

令 \(\theta=\log L,\ q=dr/d\theta\)。对固定 x 与 y=x/L，原生成元的候选 density 是
\[
 \partial_\theta\Lambda_{r,L}(x)
 =L^{-1}\left\{\frac q{2(1-r)}g_{1-r}(y)
                   -\Lambda_r(y)-y\Lambda_r'(y)\right\}.
\tag{10}
\]
q 在此只是斜率，不是原 soft response q。为避免 x=0/积分分部混淆，可先在 (9) 中换 v=s/L：
\[
 \Lambda_{r,L}(x)=\int_{1/L}^\infty F_r(Lv)e^{-v|x|}dv.
\tag{11}
\]
下端 phase 为0。逐点在 s>1，
\[
 \boxed{q\,\partial_r F_r+s\,\partial_sF_r
 =\frac{s^2}{2(u_r(s)^2+\pi^2)}
 \left[\frac q{r^2}-\frac{2(1-r)}r+\frac2{s^2-1}\right].}
\tag{12}
\]
因此 \(q\ge2r(1-r)\) 是真实 density (10) 非负的充分条件。特别 logistic 等号时 (12) 留下严格正项；更大的斜率也保持非负。原所询 q=2(1-r) 是更强但同样充分的条件，较小的 \(2r(1-r)\) 已足够。

导数合法性可在固定紧参数 \(r\in[\varepsilon,1-\varepsilon],L\in[\ell_0,\ell_1]\) 核：高 s 为 \(O(s^{-2})\)；s↓1 的可积项是 \(O(1/[(s-1)\log^2(s-1)])\)。L¹密度与总质量可微，(10) 在 x>0 成立，x=0 由同一正 spectral integral 给规范代表。不存在遗漏的 support 边界原子，因为 \(F_r(1+)=0\)。

还可完全不依赖微分：在 (1) 上，
\[
 u_{r(L)}(Lv)=\frac{v^2}{C}-\log(L^2v^2-1).
\tag{13}
\]
对原来的 v≥1/L 区间，L增大使 u 减小，F=atan2(π,u)/(2π) 增大；同时 spectral support 向下扩张，新增部分非负。因此任意有限 \(L_2\ge L_1>0\) 都有
\[
 \Lambda_{r(L_2),L_2}-\Lambda_{r(L_1),L_1}\ge0
\tag{14}
\]
作为原空间正测度。商 \(P_{r(L_2),L_2}/P_{r(L_1),L_1}\) 是该有限测度差的 compound-Poisson 概率 kernel；不是仅因 symbol 小于1便称 Markov。

沿等号曲线每轴总率
\[
 \partial_\theta|\Lambda|=\frac{r'}{2(1-r)}=r\le1.
\tag{15}
\]
所以在有限 θ区间真正非齐次生成元为 \(\mathcal J_\theta\phi(x)=\int[\phi(x+z)-\phi(x)]\,\partial_\theta\Lambda_{r,L}(dz)\)，轴向 tensor 时取轴和，总率≤n，保常数、无 killing、无 drift。新状态路径是这条物理尺度曲线的 Markov 演化；不能称它就是原 fixed-L ON 路径。

这还给原 position increment 的有限矩工具。(11)/(12) 的正 exponential mixture 每个 rate v≥1/L，所以 rate-r 的单轴瞬时跳核满足
\[
 \int |z|^k\,\partial_\theta\Lambda(dz)\le r\,k!L^k,\qquad k\ge1.
\tag{15a}
\]
偶性给零漂移。由 \(B(q)=q/2+O(q^2)\)，原 P 的方差为 \(rL^2/2\)；沿 logistic 求导，瞬时 variance 为 \(r(2-r)L^2\)，故 rate-r 条件下每次跳的 variance恰为 \((2-r)L^2\)。

对有限区间 \([\theta_0,\theta_1]\)、固定向量 \(\mathbf a\)，n轴独立非齐次跳的线性增量 \(Y=\mathbf a\cdot(X_{\theta_1}-X_{\theta_0})\) 满足
\[
 \log\mathbb E e^{sY}
 \le\frac{s^2V}{2(1-|s|b)},\quad
 V=2\|\mathbf a\|_2^2\int_{\theta_0}^{\theta_1}rL^2d\theta,\quad
 b=L_{\max}\|\mathbf a\|_\infty,\quad |s|b<1.
\tag{15b}
\]
直接证明是原 Poisson cumulant中逐轴展开，使用零一阶矩与 (15a)，再以 \(\sum_i|a_i|^k\le\|\mathbf a\|_2^2\|\mathbf a\|_\infty^{k-2}\) 求几何级数；没有把不同来源归一化后视为新输入。这只控制真实 position 的线性增量；任意 \(f(X)\)、receiver winner、历史门不是可免费代入的 coordinate-additive/Lipschitz测试，(15b) 不支付原几何交通。

## 5. 首跳后 c<1 不能免费继承

原条件 continuation 的 \(c=1-t,\ p=(\sigma-t)/(1-t)\) 使用
\[
 \widehat P^{c}_{p,L}
 =\sqrt{\frac{1+c(1-p)B(L^2\xi^2)}{1+cB(L^2\xi^2)}}
 =\frac{\widehat P^1_{\,1-c(1-p),L}}
        {\widehat P^1_{\,1-c,L}},
\quad
 \Lambda^c_{p,L}=\frac12\int_{c(1-p)}^c g_a(\cdot/L)\frac{da}{aL}.
\tag{16}
\]
两条 c1 正链相减，不能因此断言其 Lévy 差也递增。这一点与保持原首标签条件化直接相关。

具体看 \(0<c<1,0<p<1\)，令 \(A=c(1-p)/(1-c(1-p)), B=c/(1-c)\)，t=s²、\(\ell=\log(t-1),v=1/(t-1)\)，\(D_A=(At-\ell)^2+\pi^2,D_B=(Bt-\ell)^2+\pi^2\)。条件 phase 是 \(F^{c}_p=F^1_{1-c(1-p)}-F^1_{1-c}\)。若仍取 \(p'=2p(1-p)\)，其 spectral 导数除以 t 精确为
\[
 \frac{v-A^2/B}{D_A}+\frac{B-v}{D_B}.
\tag{17}
\]
通分后的分子为
\[
 (B-A)\left\{
 v[(A+B)t^2-2t\ell]-2At\ell+(1+A/B)(\ell^2+\pi^2)
 \right\}.
\tag{18}
\]
固定 c,p 的 t→∞ 时，括号主项 \(-2At\log t\) 为负。因此 c1 的**逐 spectral coefficient 正性证明**确实失效。这尚不是条件空间 Lévy derivative 负的反例：signed exponential mixture 可以仍产生非负空间函数，不能拿 (18) 的高频负权冒充 measure-order 否定。更大的条件斜率或其它路径仍待审。

## 6. 能接什么，尚不能接什么

完整 pre-first 原 soft \(K_{\sigma,L}\) 是 c1；(1) 上其原 P square root 有共同源/正增量，任意移动 outer hard sensor 可以进入已证 acceptance-count 表示。此处“可以进入”限于同一固定 C 的 family，所用响应/winner必须来自它本身。

原 actual \(L_s(x),\sigma(x)\) 不是预先固定 C 的一条曲线，C=\(\sigma/[(1-\sigma)L_s^2]\) 可随 receiver 变化。不能按每个 C 重领同一来源，或免费并无穷 family。沿较大 L 的 logistic 增大 σ 落在原 fullfuture 测试域中，故 fullfuture cap仍可测试；但 cap不说明原 selected kernel被未来核点态支配，也不支付跨 C 的联合 selector。

更关键的是原 accepted traffic条件在 \(t,C_{\rm source},i,\rho,y,w,z,\theta\) 与原 ON/FIRST 历史下，continuation 是 (16)。从整份 c1曲线 path 抽到某一 first事件，并未证明其条件 law等于该固定物理 L 下原首跳 t条件 law；t改变时基点 \(P^1_{t,L}\) 本身改变。没有同一 canonical joint RN coupling，就不能凭全 kernel可卷积商，继承原 LCA、seed、失败权与所有实际门。

本轮解析新增的是**真实斜向 Lévy 正链**，以及 c<1迁移的具体未证接口，不是 actual空间付款，不改变 \(R_{\rm angle}\) 或主账，也未给 Palm小 count/stop后尾预算。

## 7. 数值与证据范围

本稿没有运行数值。radial agent负责独立注册/执行 c1 phase 与 finite-step guard，本稿只独审解析，完成后可引用其专属终态收据，不能冒称本稿跑过或其有限样点证明了全参数序。(14) 已由原 positive density 的全参数解析给出。

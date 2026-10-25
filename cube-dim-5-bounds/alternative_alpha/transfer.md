# 五维 16-body-direction α=3/2 的有限源传递

本证明固定 Lebesgue 测度、中心共同半径方体和
\[
 V=\{(1,\epsilon_2,\epsilon_3,\epsilon_4,\epsilon_5):\epsilon_j\in\{-1,1\}\},
 \quad L=-\sum_{v\in V}|v\cdot D|^{3/2},\quad P_t=e^{tL}.
\]
它证明的是条件性上界接口：最后仍需要独立的 essential 下界 (g) 及面
积分 (I)。本文件不使用、也不需要任何小时间奇异层的点态渐近。

Fourier 约定为 \(\widehat f(\xi)=\int e^{-ix\cdot\xi}f(x)dx\)，逆变换带
\((2\pi)^{-5}\)。以下每个等式都是全空间等式；不存在把 α=1 的 domain
论证照搬到分数阶生成元的步骤。

## 1. 半群、分数阶域和 cutoff

令 \(Z_t\) 是标准对称 (3/2)-stable 变量，
\(\mathbb E e^{isZ_t}=e^{-t|s|^{3/2}}\)，并令
\(X_t=\sum_{v\in V}vZ_t^v\)，其中各变量独立。其特征函数为
\[
 e^{-t\psi(\xi)},\qquad \psi(\xi)=\sum_{v\in V}|v\cdot\xi|^{3/2}.
\]
有限方向张成 \(\mathbb R^5\)，所以 ψ 在单位球上的连续正函数，从而存在
\(0<a\le b<\infty\)，使
\[
 a|\xi|^{3/2}\le\psi(\xi)\le b|\xi|^{3/2}. \tag{1}
\]

从 16 个方向中取任一五方向基 (B)（已有整数几何给 \(|\det B|=48\)）。
五项和的密度是 \(|\det B|^{-1}\prod_{j=1}^5s_t((B^{-1}x)_j)\)，其中
\(s_t(z)=t^{-2/3}s_1(t^{-2/3}z)\) 是一维 stable 密度。其余十一项只是概率
卷积。因此 (p_t) 严格正、连续、保质量、自伴，且
\[
 p_t(x)=t^{-10/3}p_1(t^{-2/3}x),\qquad
 \|p_t\|_\infty\le C t^{-10/3}. \tag{2}
\]
这里 (C=48^{-1}\|s_1\|_\infty^5) 即可；不需要它的数值。卷积近似恒等性给
\(P_tf\to f\) 于 (L^1,L^2\)，并对 (C_0\) 函数一致收敛。

在 (L^2\) 中，(L) 是自伴非正乘子 \(-\psi\)，故
\[
 \operatorname{Dom}_{L^2}(L)=H^{3/2}(\mathbb R^5),\qquad
 \operatorname{Dom}((-L)^{1/2})=H^{3/4}(\mathbb R^5), \tag{3}
\]
其中范数均加上 (L^2) 项；这是 (1) 和 Plancherel 的直接等价。

对 \(0<\alpha<2\) 的一维稳定符号，存在正的常数 \(c_\alpha\)，使每个
\(\phi\in C_c^\infty\) 满足
\[
 |v\cdot D|^\alpha\phi(x)=c_\alpha\int_\mathbb R
 \frac{2\phi(x)-\phi(x+sv)-\phi(x-sv)}{|s|^{1+\alpha}}\,ds. \tag{4}
\]
本处取 \(\alpha=3/2\)。二阶 Taylor 公式和三角不等式分别给其 (L^1) 及
\(L^\infty\) 控制：被积二差分不超过
\[
 \min\{s^2|v|^2\|D^2\phi\|_{L^p},\;4\|\phi\|_{L^p}\},\qquad p=1,\infty. \tag{5}
\]
将 (5) 对 \(|s|^{-1-3/2}ds\) 积分，在 (0) 端使用 (s^2)，在无穷端使用
常数，说明
\[
 L\phi\in L^1\cap C_0\cap L^\infty. \tag{6}
\]
连续性由 (4) 的平移连续性和同一可积支配给出。

取光滑径向、径向非增的 \(0\le\chi\le1\)，它在原点邻域为一、紧支撑。
于是 \(R\mapsto\chi_R(x)=\chi(x/R)\) 对每个 \(x\) 单调增加到一。从 (4) 的
变量替换得到**精确**尺度式
\[
 \|L\chi_R\|_\infty=R^{-3/2}\|L\chi\|_\infty,\qquad
 \|L\chi_\varepsilon\|_1=\varepsilon^{7/2}\|L\chi\|_1. \tag{7}
\]

## 2. 有界紧支撑源的 obstacle 解与质量

先令 \(f\ge0\) 有界且紧支撑，记 \(F=\|f\|_\infty\)、\(m=\|f\|_1\)。
由 (2)
\[
 \|P_sf\|_\infty\le\min\{F,Cms^{-10/3}\}. \tag{8}
\]
对 (0<t\le1) 定义
\[
 V_t=t\sum_{j\ge0}P_{jt}f,\qquad
 u_{t,0}=0,\qquad u_{t,k+1}=(P_tu_{t,k}+t(f-\kappa))_+ . \tag{9}
\]
由单调性，右 Riemann 和被 (8) 的积分控制：
\[
 \|V_t\|_\infty\le tF+\int_0^\infty\min(F,Cms^{-10/3})ds
 \le A:=F+\frac{10}{7}(Cm)^{3/10}F^{7/10}. \tag{10}
\]
这里第二项由在 \(s=(Cm/F)^{3/10}\) 处分裂积分直接得到。级数在除首项外
于 (C_0) 一致收敛，且
\[
 P_tV_t=V_t-tf. \tag{11}
\]
归纳使用 (u_{t,k}\le V_t) 和 (11) 得
\(u_{t,k+1}\le(V_t-\kappa t)_+\)。该上界有紧支撑：写
\(V_t-\kappa t=tf+(V_t-tf)-\kappa t\)。在一个包含 \(\operatorname{supp}f\)
的紧集外，首项为零，而 \(V_t-tf\to0\) 于无穷；故其和最终至多
\(-\kappa t/2<0\)。故极限 \(u_t\) 属于 \(L^1\cap L^2\)，
并满足 (9) 的固定点式。

令 \(L_t=(P_t-I)/t\)、\(\sigma_t=f+L_tu_t\)。固定点式的正、零两支分别给
\[
 0\le\sigma_t\le\kappa,\qquad \sigma_t=\kappa\ \text{a.e. on }\{u_t>0\}. \tag{12}
\]
因为 (u_t\in L^1)，保质量可合法用于 (9)，从而 \(\int\sigma_t=m\)。又由
\(\int u_t\sigma_t=\kappa\|u_t\|_1\) 和 \(\sigma_t=f+L_tu_t\)，
\[
 \langle u_t,-L_tu_t\rangle+\kappa\|u_t\|_1=\int fu_t. \tag{13}
\]
结合 \(0\le u_t\le A\) 产生统一界
\[
 \|u_t\|_1\le mA/\kappa,\quad
 \|u_t\|_2^2\le mA^2/\kappa,\quad
 \|\sigma_t\|_2^2\le\kappa m,\quad
 \|L_tu_t\|_2\le\|f\|_2+\sqrt{\kappa m}. \tag{14}
\]

写 \(\psi_t=(1-e^{-t\psi})/t\)。基本不等式
\(1-e^{-z}\ge(1-e^{-1})\min(z,1)\) 给：若 \(t\le Q^{-1}\) 且
\(\psi(\xi)\ge Q\)，则
\[
 \psi_t(\xi)\ge(1-e^{-1})Q. \tag{15}
\]
由 (14) 和 Plancherel，\(\widehat u_t\) 在 \(\{\psi\ge Q\}\) 的 (L^2)
尾至多是常数乘 (Q^{-1})。另一方面，在固定 \(\{\psi\le2Q\}\) 上，
频率截断的 (u_t) 因 (14) 而在 (H^1) 有界（紧频率支撑把 (L^2) 控制提升到
任意固定整数 Sobolev 阶）。Rellich 在每个有界空间集上给紧性。先固定 (Q)，
再令 (Q\to\infty) 的对角抽取与 (15) 的统一尾控制给某个 \(t_j\downarrow0\)：
\[
 u_{t_j}\to u\quad\hbox{于局部 }L^2\hbox{ 强收敛且 a.e.，}
 \qquad u_{t_j}\rightharpoonup u,\quad\sigma_{t_j}\rightharpoonup\sigma
 \quad\hbox{于 }L^2. \tag{16}
\]
Fatou 与 (14) 保留 (u\ge0\)、(u\in L^1\cap L^2)。

对 \(\phi\in\mathcal S\)，
\(0\le\psi_t\le\psi\) 和 dominated convergence 给
\(L_t\phi\to L\phi\) 于 (L^2)。在
\(\langle u_t,L_t\phi\rangle=\langle\sigma_t-f,\phi\rangle\) 中使用 (16)，得到
\[
 Lu=\sigma-f\quad\text{于分布及 }L^2,\qquad u\in H^{3/2}. \tag{17}
\]
若 \(h\in L^2\) 有界支撑且 \(h=0\) a.e. 于 \(D^c\)，其中 \(D=\{u>0\}\)，
则 (16) 保证在 \(D\) 上最终 \(u_{t_j}>0\)，故 (12) 和 dominated convergence
给 \(\int h\sigma=\kappa\int h\)。于是
\[
 0\le\sigma\le\kappa,\qquad\sigma=\kappa\quad\text{a.e. on }D. \tag{18}
\]

对固定 (R)，以 \(\chi_R\in L^2\) 测试 \(\sigma_t\rightharpoonup\sigma\)，
再令 (R\uparrow\infty)，得 \(\sigma\in L^1\)、\(\int\sigma\le m\)。以
\(\chi_R\) 测试 (17)，并用 (7)，
\[
 \left|\langle u,L\chi_R\rangle\right|
 \le\|u\|_1R^{-3/2}\|L\chi\|_\infty\longrightarrow0.
\]
因此 \(\int\sigma=m\)，而 (18) 给 \(\kappa|D|\le m\)。至此没有使用全空间
强紧性、统一支撑或局部性断言。

## 3. Green 势的正部是全空间有限源

令 \(G=\int_0^\infty p_tdt\)，允许取值 (+\infty\)。设
\(q=5/(5-3/2)=10/7\)。由 (2)
\[
 \Big\|\int_0^T p_tdt\Big\|_1=T,\qquad
 \Big\|\int_T^\infty p_tdt\Big\|_\infty\le C' T^{-7/3}. \tag{19}
\]
以阈值 \(\lambda/2\) 分解并取 \(T\asymp\lambda^{-3/7}\)，得到
\(G\in L^{10/7,\infty}\)。对任意有限测度集 (E)，Tonelli 和 (2) 给
\[
 \int_EG\le\int_0^\infty\min\{1,C|E|t^{-10/3}\}dt<\infty. \tag{20}
\]
所以 (G) 局部 (L^1)、a.e. 有限；对于每个 (c>0)，
\(H=(G-c)_+\in L^1)，因为 weak-(L^{10/7}) 的分布尾从 (c) 到无穷可积。

由 (6)、(2)，对 \(\phi\in C_c^\infty\)
\[
 \int G|L\phi|
 \le\int_0^\infty\min\{\|L\phi\|_\infty,
 C t^{-10/3}\|L\phi\|_1\}dt<\infty. \tag{21}
\]
故 Fubini 合法。半群导数
\(\partial_tP_t\phi=P_tL\phi\) 以及
\(\|P_T\phi\|_\infty\le CT^{-10/3}\|\phi\|_1\to0\) 给
\[
 \langle G,L\phi\rangle=\int_0^\infty P_tL\phi(0)dt=-\phi(0),
 \quad\text{即 }LG=-\delta_0. \tag{22}
\]

Tonelli 和半群律（右端由 (19) 有限）给
\(P_tG=\int_t^\infty p_sds\)，因此 a.e. \(G=P_tG+U_t\)，其中
\(U_t=\int_0^tp_sds\)。保正和 Jensen 的正部不等式给
\[
 P_tH\ge(P_tG-c)_+\ge H-U_t. \tag{23}
\]
对非负 \(\phi\in C_c^\infty\) 配对、除以 (t)。因为
\[
 \frac{P_t\phi-\phi}{t}-L\phi
 =\frac1t\int_0^t(P_sL\phi-L\phi)ds\longrightarrow0
 \quad\text{于 }L^\infty
\]
（使用 (6) 和 (C_0) 上的近似恒等性），而
\(t^{-1}\langle U_t,\phi\rangle\to\phi(0)\)，令 (t\downarrow0\) 得
\[
 \mu:=LH+\delta_0\ge0 \tag{24}
\]
作为整个空间的正 Radon 测度。没有在任一奇异线或无穷远删去测试函数。

由 (7) 与 (H\in L^1)，
\[
 \int\chi_Rd\mu=1+\int H L\chi_R\longrightarrow1.
\]
单调收敛于是给 \(\mu(\mathbb R^5)=1\)。再令 \(h=\min(G,c)\)。由 (22) 和
\(H=G-h\) 有 \(\mu=-Lh\)。因 (0\le h\le c)、\(\mu\ge0\)，(7) 给
\[
 0\le\int\chi_\varepsilon d\mu
 =-\int hL\chi_\varepsilon
 \le c\|L\chi_\varepsilon\|_1
 =c\varepsilon^{7/2}\|L\chi\|_1\longrightarrow0.
\]
故 \(\mu(\{0\})=0\)。这是正部势的全空间有限源结论。

## 4. 同一接收集、全部尺度及一般 (L^1)

假设现在有一个经认证的 essential 下界 (G\ge g>0) a.e. 于
\(Q=[-1,1]^5\)，并取 (0<c<g)。令
\[
 N=32(g-c),\qquad K=H/N,\qquad M=\|K\|_1.
\]
则 (K\ge2^{-5}1_Q) a.e.。由 (24) 的 Fourier 变换，
\[
 -\psi\widehat K=(\widehat\mu-1)/N. \tag{25}
\]
先令 (u\in\mathcal S)，把 (25) 乘 \(\widehat u\) 并作逆 Fourier 变换；
所有项都在 (L^2)，且 \(K\in L^1\)、\(\mu\) 有限。再用
\(\mathcal S\) 在 (H^{3/2}) 稠密、Young 不等式延拓，得到
\[
 K*Lu=(\mu*u-u)/N\quad\text{于 }L^2,\quad u\in H^{3/2}. \tag{26}
\]
对于 \(K_r(x)=r^{-5}K(x/r)\) 和 \(\mu_r=(x\mapsto rx)_\#\mu\)，
\(\psi(r\xi)=r^{3/2}\psi(\xi)\) 使 (25) 在 (r\xi) 处化为
\[
 K_r*Lu=\frac{\mu_r*u-u}{N r^{3/2}}. \tag{27}
\]
这就是分数阶版本的缩放因子；它不能误写为 (r^{-1})。

对第 2 节 obstacle 取 \(\kappa=\lambda/M\)。由 \(f=\sigma-Lu\)、(27)，
在同一个集合 (D^c\) 上，对每个固定 (r>0) 有
\[
 K_r*f=K_r*\sigma+\frac{u-\mu_r*u}{Nr^{3/2}}
 \le K_r*\sigma\le\kappa M=\lambda. \tag{28}
\]
这里 \(u=0) a.e. 于 (D^c)，且 \(\mu_r*u\ge0) 由 μ 的正性及
Tonelli 给出。先对正有理 (r) 合并零集。对有界紧支撑 (f)，
\[
 \|K_r*f-K_s*f\|_\infty\le\|f\|_\infty\|K_r-K_s\|_1\to0\quad(s\to r), \tag{29}
\]
其中最后一步是 (L^1) dilation 的强连续性。因此 (28) 同时对全部 (r>0)
成立，且
\[
 \lambda\,|\{\sup_{r>0}K_r*f>\lambda\}|
 \le\lambda|D|\le M\|f\|_1. \tag{30}
\]

对一般 (f\in L^1)，取
\(f_j=\min(|f|,j)1_{[-j,j]^5}\uparrow|f|\)。非负卷积的单调收敛给
\[
 \sup_rK_r|f|=\sup_r\sup_jK_rf_j=\sup_j\sup_rK_rf_j.
\]
右端是可测函数的递增极限；对严格超水平集取单调极限便将 (30) 延至全部
\(L^1\) 输入和全部正尺度。又 (K_r\) 支配相同尺度的中心方体平均，故
\(C_5\le M\)。

## 5. weak-​(L^q) 自动给 a.e. 面积分；条件成本

令 \(d=5-3/2=7/2\)、\(q=5/d=10/7\)。缩放给
\(G(rx)=r^{-d}G(x)\) a.e.。记 \(V(\lambda)=|\{G>\lambda\}|\)。第 3 节的
weak-​(L^q) 界给 (V(1)<\infty)，而齐次性给精确式
\[
 V(\lambda)=\lambda^{-q}V(1). \tag{31}
\]
将空间分成 (5\cdot2^5) 个以坐标面为底的锥，并以
\(x=r(1,y)\)、Jacobian \(r^4\) 参数化。Tonelli（允许先取值 (+\infty)）给
\[
 J:=\sum_{\text{all faces }F}\int_FG^q=5V(1)<\infty. \tag{32}
\]
因此不需要连续性或点态排除奇异线：weak-​(L^q) 已使每个面 a.e. 的
\(q\) 次积分有限。body 对称使每个面、每个正正交象限相同；若
\[
 I=\int_{[0,1]^4}G(1,y)^{10/7}dy,
\]
则 \(J=5\cdot2^5I\)，故 \(V(1)=32I\)。由层蛋糕和 (31)，
\[
 \|H\|_1=\int_c^\infty V(s)ds
 =\frac{7}{3}\,32I\,c^{-3/7}. \tag{33}
\]
于是
\[
 M=\|K\|_1=\frac{7}{3}\frac{I c^{-3/7}}{g-c}. \tag{34}
\]
令 (c=3g/10)，即一般驻点 \(c=\alpha g/5\)，得到
\[
 \boxed{\displaystyle C_5\le\left(\frac{10}{3}\right)^{10/7}
 \frac{I}{g^{10/7}}.} \tag{35}
\]

这里 (g) 是 essential cube 下界。若某面上的正测集满足 \(G<a<g\)，取
\(r<1\) 且接近一使 \(r^{-7/2}a<g\)，其射线扫出 cube 内正体积，矛盾；
故 essential cube/face 下界接口匹配，而没有把 a.e. 下界错误升级为逐点边界下界。

结论：在尚未认证 (g,I) 前，(35) 只是严格、可承接的传递接口，并不构成
任何新的数值上界。

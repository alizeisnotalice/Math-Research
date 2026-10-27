# 八维有限方向 Cauchy Green 势：奇异分层、概率源与全尺度接口

2026-10-27；独立范围 `cube-dim-8-bounds/upper/round001/green/`。

本文是本轮新写的证明，不把三维、四维审计的接受标签当作八维定理。允许的四份来源的 SHA-256 和阅读范围见 `sources.json`。证明用到的标准分析定理是 Fourier 反演与 Plancherel、Tonelli/Fatou/支配收敛、局部 Rellich 紧性、正分布为 Radon 测度以及有限测度卷积的 Young 不等式。有限代数检查见 `check.py` 和实际运行的 `checks.json`；它只核常数、指数和模型的有限代数，不代替下述分析证明。

## 1. 对象、归一和结论范围

设 \(V=\{v_1,\ldots,v_m\}\subset\mathbb R^8\setminus\{0\}\) 是有限 spanning frame，\(c_j>0\)，不把 \(v_j\) 单位化。向量可以平行或重复。令
\[
 q_a(s)=\frac{a}{\pi(a^2+s^2)},\qquad
 X_t=\sum_{j=1}^m v_j Z_{j,t},\quad Z_{j,t}\sim q_{c_jt}
\]
各项独立。取
\[
 \psi(\xi)=\sum_jc_j|v_j\cdot\xi|,\qquad
 L=-\psi(D),\quad P_tf=p_t*f,
 \quad G(x)=\int_0^\infty p_t(x)\,dt.
\]
Fourier 变换为 \(\widehat f(\xi)=\int e^{-ix\cdot\xi}f(x)dx\)，逆变换因子 \((2\pi)^{-8}\)。所有输入、输出及卷积均用同一 Lebesgue 测度。

令
\[
 S_3=\bigcup_{A\subseteq V,\ \operatorname{rank}A\le3}\operatorname{span}(A).
\]
包括空集张成的 \(\{0\}\)。这是有限个闭线性子空间的并，故闭且零八维测度。本文证明：

1. \(G=+\infty\) 恰好在 \(S_3\) 上；在 \(\mathbb R^8\setminus S_3\) 上 \(G\) 有限且局部连续。
2. \(G\in L^{8/7,\infty}\cap L^1_{\rm loc}\)，\(H_c=(G-c)_+\in L^1\)；全空间中
   \[LH_c=\mu_c-\delta_0,\quad \mu_c\ge0,\quad\mu_c(\mathbb R^8)=1,\quad\mu_c(\{0\})=0.\]
3. essential cube 下界 \(0<g_0\le\operatorname{essinf}_{[-1,1]^8}G\) 给非负 \(L^1\) 主核，进而给一般 \(L^1\) 输入和所有共同半径的弱 \((1,1)\) 上界。
4. 在面 \(x_1=1\) 上，短时 Green 部分的全平面 \(L^{8/7}\) 误差是显式常数乘 \(\varepsilon^{9/8}\)。长时统一上确界是 \(C/(7T^7)\)，其简单 \(L^{8/7}\) 推论只适用于有限面域。

第 2、4—6 节的半群/源/比较机制实际在所有整数 \(n>1\) 成立；第 3 节的精确奇异层数和第 8 节的指数专属于 \(n=8\)。本文不提供八维的连续 \(g_0,I_0\) 数值证书，不导出新的数值常数。

## 2. 真实密度、测试生成元及缩放

以下先写任意 \(n>1\)。取方向基 \(B=(v_{i_1},\ldots,v_{i_n})\)，则这 \(n\) 项之和有密度
\[
 p_t^B(x)=\frac1{|\det B|}\prod_{k=1}^n
 \frac{c_{i_k}t}{\pi((c_{i_k}t)^2+(B^{-1}x)_k^2)}.
\]
其余项和的分布记为 \(\rho_t\)，故 \(p_t=p_t^B*\rho_t\)。因此 \(p_t>0\)、\(p_t\in C_0\)、\(\int p_t=1\)，且
\[
 \|p_t\|_\infty\le Ct^{-n},\qquad
 C=\frac1{\pi^n|\det B|\prod_{k=1}^nc_{i_k}}.
\]
各 Cauchy 项的共同缩放给 \(p_t(x)=t^{-n}p_1(x/t)\)，从而 \(p_t(x)\) 在 \(t>0\) 联合连续。概率构造给自伴、保正、保常数和质量的卷积半群。\(p_t\) 是近似恒等核，且 \(t\mapsto p_t\) 在正时间 \(L^1\) 连续。

单位球上连续函数 \(\psi\) 严格正：若 \(\psi(\xi)=0\)，所有 \(v_j\cdot\xi=0\)，spanning 迫使 \(\xi=0\)。于是存在 \(0<a\le b<\infty\) 使
\[
 a|\xi|\le\psi(\xi)\le b|\xi|.
\]
故 \(L\) 的 \(L^2\) 算子域为 \(H^1\)，形式域为 \(H^{1/2}\)；均指带 \(L^2\) 项后的等价范数。没有使用高维 \(H^1\) 到连续函数的错误嵌入。

对 \(\phi\in C_c^\infty\)，真实生成元是
\[
 L\phi(x)=\frac1{2\pi}\sum_jc_j\int_\mathbb R
 \frac{\phi(x+sv_j)+\phi(x-sv_j)-2\phi(x)}{s^2}\,ds.
\]
用二差分的上确界估计
\[
 \|\Delta_{sv_j}^2\phi\|_\infty
 \le\min\{s^2|v_j|^2\|D^2\phi\|_\infty,4\|\phi\|_\infty\},
\]
及同样的 \(L^1\) 估计（方向二阶导的 \(L^1\) 范数）积分，得 \(L\phi\in L^1\cap L^\infty\cap C_0\)。系数是 \(c_j\)，不是 \(c_j^2\)，未单位化向量的长度不能删去。一阶齐次性给
\[
 L[\chi(\cdot/R)](x)=R^{-1}(L\chi)(x/R),
\quad \|L\chi_R\|_\infty=R^{-1}\|L\chi\|_\infty,
\quad \|L\chi_\varepsilon\|_1=\varepsilon^{n-1}\|L\chi\|_1.
\]
对 Schwartz \(\phi\)，由 \(\psi|\widehat\phi|\in L^1\) 及 Fourier 支配收敛，
\[
 (P_t\phi-\phi)/t\longrightarrow L\phi\quad\text{于 }L^\infty,
 \qquad \partial_tP_t\phi=P_tL\phi.
\]

## 3. 八维奇异分层的完整证明

### 3.1 补集：大系数集合必须有秩至少四

设 \(E\Subset\mathbb R^8\setminus S_3\)，\(d=\operatorname{dist}(E,S_3)>0\)，取
\[
 D_V=\sum_{j=1}^m|v_j|>0,\qquad \delta=\frac d{4D_V}.
\]
考虑开邻域 \(O=\{x:\operatorname{dist}(x,E)<d/4\}\)。对一组系数 \(s=(s_j)\)，置
\[
 J(s)=\{j:|s_j|>\delta\}.
\]
若 \(\operatorname{rank}\{v_j:j\in J(s)\}\le3\)，大项和属于 \(S_3\)，小项和范数至多
\[
 \Big|\sum_{j\notin J(s)}s_jv_j\Big|\le\delta\sum_j|v_j|=d/4.
\]
故 \(\operatorname{dist}(\sum_js_jv_j,S_3)\le d/4\)，不可能落入 \(O\)，因为 \(O\) 中点到 \(S_3\) 的距离大于 \(3d/4\)。因此，落入 \(O\) 时，**整个大系数方向集合的秩至少为四**。这之后才可取其中四个独立方向。仅仅证明“有四个大系数”不够，因为这四个方向可能相关。

令 \(\mathcal I_4\) 为全部四元独立方向指标集，逐个 \(I\in\mathcal I_4\) 选含 \(I\) 的方向基 \(B_I\)。在 \(I\) 的四个系数上限制 \(|Z_{j,t}|>\delta\)，保留**未归一化**的子概率律；其密度界为
\[
 q_{c_jt}(s)1_{|s|>\delta}\le\frac{c_jt}{\pi\delta^2}.
\]
基中另四项密度界是 \(1/(\pi c_jt)\)。基变换及余项概率卷积给该受限随机和密度上界
\[
 C_I=\frac{\prod_{j\in I}c_j}
 {\pi^8|\det B_I|\,\delta^8\prod_{j\in B_I\setminus I}c_j}.
\]
这里幂为 \(t^4t^{-4}=1\)，并未除以受限事件的概率。由事件并集，对每个 Borel \(F\subset O\)，
\[
 \mathbb P(X_t\in F)\le\sum_{I\in\mathcal I_4}
 \mathbb P(X_t\in F,\ |Z_{j,t}|>\delta\ (j\in I))
 \le C_E|F|,\quad C_E=\sum_IC_I.
\]
所以 \(p_t\le C_E\) 在 \(O\) 几乎处处；真实 \(p_t\) 连续，将该界升级为 \(O\) 上逐点界。对 \(0<t\le1\) 用 \(C_E\)，对 \(t\ge1\) 用 \(Ct^{-8}\)。这给 \(t\mapsto p_t(x)\) 的局部统一可积控制，支配收敛证明 \(G\) 在 \(E\) 上有限连续。因 \(E\) 任意，补集结论成立。

### 3.2 低秩空间：下界的小时间幂不可积

若 \(x\in S_3\)，可选独立方向 \(v_{i_1},\ldots,v_{i_r}\)，\(0\le r\le3\)，使
\[
 x=\sum_{k=1}^ra_kv_{i_k}.
\]
将其补成八方向基 \(B\)，在基坐标中 \(B^{-1}x=(a_1,\ldots,a_r,0,\ldots,0)\)。余项和 \(Y_t\overset d=tY_1\)。取 \(M<\infty\) 使事件
\[
 E_t=\{\|B^{-1}Y_t\|_\infty\le Mt\}
\]
概率 \(\theta>0\)，与 \(t\) 无关；若没有余项则取 \(\theta=1\)。这样的 \(M\) 存在，因为 \(B^{-1}Y_1\) 几乎处处有限。

在 \(E_t\) 上，\(0<t\le1\)，前 \(r\) 个基密度均有下界 \(k_xt\)：分母由 \((c_j)^2+(|a_k|+M)^2\) 统一控制；后 \(8-r\) 个零坐标的基密度均有下界 \(k/t\)。取余项期望得
\[
 p_t(x)=\mathbb E[p_t^B(x-Y_t)]\ge k'_x t^{2r-8}.
\]
\(2r-8\le-2\)，故 \(\int_0^1p_t(x)dt=\infty\)。零点取 \(r=0\) 同样成立。这证明 **\(G^{-1}(\{+\infty\})=S_3\)**，包含所有相关方向的低秩张成空间，而不只方向线。

## 4. 全空间弱 \(L^{8/7}\)、可积正部及有限概率源

这一节对所有 \(n>1\) 成立，最后取 \(n=8\)。设 \(q=n/(n-1)\)，写
\[
 U_T=\int_0^Tp_tdt,\qquad W_T=\int_T^\infty p_tdt,
\quad\|U_T\|_1=T,\quad0\le W_T\le\frac C{n-1}T^{1-n}.
\]
取 \(T=(2C/((n-1)\lambda))^{1/(n-1)}\)，则 \(W_T\le\lambda/2\)，Chebyshev 给
\[
 |\{G>\lambda\}|\le A_n\lambda^{-q},
 \quad A_n=2\left(\frac{2C}{n-1}\right)^{1/(n-1)}.
\]
有限体积 \(E\) 满足
\[
 \int_EG\le\int_0^\infty\min\{1,C|E|t^{-n}\}\,dt<\infty.
\]
故 \(G\in L^1_{\rm loc}\)、几乎处处有限。\(H_c=(G-c)_+\) 满足
\[
 \int H_c=\int_c^\infty|\{G>s\}|ds
 \le (n-1)A_nc^{-1/(n-1)}<\infty.
\]
这一步处理整个空间，包括 \(S_3\) 邻域和无穷远；不是由齐次性独自推出的。

对 \(\phi\in C_c^\infty\)，
\[
 \int G|L\phi|\le\int_0^\infty
 \min\{\|L\phi\|_\infty,Ct^{-n}\|L\phi\|_1\}\,dt<\infty.
\]
所以半群导数可合法积分。由自伴性及 \(P_T\phi(0)\to0\)，
\[
 \langle G,L\phi\rangle=\int_0^\infty P_tL\phi(0)dt=-\phi(0),
 \qquad LG=-\delta_0.
\]
正 Tonelli 和半群律给
\[
 P_tG=\int_t^\infty p_sds,\qquad G=P_tG+U_t.
\]
后者在几乎处处有限点使用，无须在 \(S_3\) 上相减无穷量。Jensen、保常数和正部的 1-Lipschitz 性质给
\[
 P_tH_c\ge(P_tG-c)_+\ge H_c-U_t.
\]
对非负紧光滑测试配对，除以 \(t\) 并令 \(t\downarrow0\)。左侧使用 \(H_c\in L^1\) 与测试生成元的 \(L^\infty\) 收敛，右侧使用 \(U_t/t\rightharpoonup\delta_0\)，得
\[
 LH_c+\delta_0\ge0
\]
于**全空间分布**。正分布定理给正 Radon 测度 \(\mu_c=LH_c+\delta_0\)。没有排除任何奇异子空间或 level 界面，也不要求 \(H_c\) 紧支撑。

取径向单调的 \(0\le\chi_R\uparrow1\)，且原点附近等于 1。由大截断缩放及 \(H_c\in L^1\)，
\[
 \int\chi_Rd\mu_c=1+\int H_cL\chi_R\longrightarrow1.
\]
正测度单调收敛给 \(\mu_c(\mathbb R^n)=1\)。置 \(h=\min(G,c)\)，几乎处处有 \(G=H_c+h\)，合法配对给 \(\mu_c=-Lh\)。小截断则给
\[
 0\le\int\chi(x/\varepsilon)d\mu_c(x)
 =-\int hL\chi_\varepsilon
 \le c\varepsilon^{n-1}\|L\chi\|_1\longrightarrow0.
\]
因此 \(\mu_c(\{0\})=0\)。真正使用 \(n>1\) 的地方包括长时可积、\(q>1\) 的正部层积分和小截断的消失。

八维常数明确为
\[
 q=8/7,\quad A_8=2(2C/7)^{1/7},\quad
 \|W_T\|_\infty\le C/(7T^7),\quad
 \|L\chi_\varepsilon\|_1=\varepsilon^7\|L\chi\|_1.
\]

## 5. 同一接收集的障碍分解：逐条核查 \(n>1\)

设非负 \(f\in L^\infty\) 紧支撑，\(F=\|f\|_\infty\)、\(m_f=\|f\|_1>0\)，\(\kappa>0\)。取 \(0<t\le1\)，定义
\[
 V_t=t\sum_{j\ge0}P_{jt}f,\quad
 u_{t,0}=0,\quad u_{t,k+1}=(P_tu_{t,k}+t(f-\kappa))_+.
\]
由 \(\|P_sf\|_\infty\le\min(F,Cm_fs^{-n})\)，正时间尾在 \(L^\infty\) 收敛，\(V_t-tf\in C_0\)。递减函数的右 Riemann 和给
\[
 \|V_t\|_\infty\le tF+
 \frac n{n-1}(Cm_f)^{1/n}F^{1-1/n}\le C_f.
\]
Tonelli 给 \(P_tV_t=V_t-tf\)，迭代单调，受 \(V_t\) 控制，极限 \(u_t\) 满足固定点方程。更强地
\[
 0\le u_t\le(V_t-\kappa t)_+.
\]
固定 \(t\) 时右侧有 essential 紧支撑，因此 \(u_t\in L^1\cap L^2\)。这里不要求 \(V_t\in L^1\)，也没有声称各 \(u_t\) 支撑一致。

令 \(L_t=(P_t-I)/t\)、\(\sigma_t=f+L_tu_t\)。固定点的正分支和零分支分别给
\[
 0\le\sigma_t\le\kappa,\quad
 \sigma_t=\kappa\text{ 于 }\{u_t>0\},\quad
 \int\sigma_t=m_f.
\]
最后的等式是在 \(u_t\in L^1\) 已确立后由质量守恒得到。与 \(u_t\) 配对给
\[
 \langle u_t,-L_tu_t\rangle+\kappa\|u_t\|_1=\int fu_t.
\]
自伴收缩性使第一项非负，从而
\[
 \|u_t\|_1\le m_fC_f/\kappa,\quad
 \|u_t\|_2^2\le m_fC_f^2/\kappa,\quad
 \|\sigma_t\|_2^2\le\kappa m_f,\quad
 \|L_tu_t\|_2\le\|f\|_2+\sqrt{\kappa m_f}=C_2.
\]

写 \(\psi_t=(1-e^{-t\psi})/t\)。若 \(t\le1/Q\)，在 \(\psi\ge Q\) 上
\[
 \psi_t\ge(1-e^{-1})Q,\qquad
 \int_{\psi\ge Q}|\widehat u_t|^2
 \le\frac{C_2^2}{(1-e^{-1})^2Q^2}.
\]
这里可用 unitary Fourier，因此没有 Plancherel 归一因子。固定低频投影的 \(H^1\) 范数有界；局部 Rellich 加高频小尾，给 \(t_j\downarrow0\) 的子列使
\[
 u_{t_j}\to u\text{ 于 }L^2_{\rm loc}\text{ 及 a.e.},\quad
 u_{t_j}\rightharpoonup u,\quad\sigma_{t_j}\rightharpoonup\sigma\text{ 于 }L^2.
\]
Fatou 给 \(u\ge0\)、\(u\in L^1\cap L^\infty\)。测试上的 \(L_t\phi\to L\phi\) 于 \(L^2\) 允许传递等式，得 \(Lu=\sigma-f\in L^2\)；Fourier 域和 \(\psi\asymp|\xi|\) 得 \(u\in H^1\)。弱极限仍在 \([0,\kappa]\)。在 \(D=\{u>0\}\) 的几乎处处点，\(u_{t_j}\) 最终为正，故 \(\sigma_{t_j}\to\kappa\)；在 \(D\) 的有限测度截集上用有界支配收敛与弱极限，得 \(\sigma=\kappa\) 于 \(D\)。

先以正截断测试弱极限得 \(\sigma\in L^1\)、\(\int\sigma\le m_f\)。再测试 \(f=\sigma-Lu\)：
\[
 |\langle Lu,\chi_R\rangle|
 =|\langle u,L\chi_R\rangle|
 \le R^{-1}\|u\|_1\|L\chi\|_\infty\to0.
\]
故 \(\int\sigma=m_f\)，并有 \(\kappa|D|\le m_f\)。这明确排除质量逃逸。

因此所有 \(n>1\) 的所需障碍结论是
\[
 f=\sigma-Lu,\quad u\ge0,\quad u\in L^1\cap L^\infty\cap H^1,
 \quad0\le\sigma\le\kappa,\quad\sigma=\kappa\text{ 于 }D,
 \quad\int\sigma=\int f.
\]

## 6. 有限源 Fourier 结算、全部尺度与一般 \(L^1\)

设非负 \(K\in L^1\)，\(M_K=\int K>0\)，满足全空间等式
\[
 LK=a(\mu-\delta_0),\qquad a>0,\quad \mu\text{ 为正概率测度}.
\]
对 Schwartz 测试，\(K\in L^1\) 和 \(\psi|\widehat\phi|\in L^1\) 允许直接 Fourier 逆积分与 Fubini，因而
\[
 -\psi(\xi)\widehat K(\xi)=a(\widehat\mu(\xi)-1)
\]
是连续函数恒等式。具体地，原来的紧支撑测试等式可延伸到 Schwartz 测试：令 \(\phi_R=\chi_R\phi\)，则 \(\phi_R\to\phi\) 于 Schwartz 拓扑，二差分估计给 \(L\phi_R\to L\phi\) 于上确界范数；\(K\in L^1\) 和有限 \(\mu\) 允许取极限。故没有不合法地将一般分布与非紧支撑核卷积。

乘任意 \(u\in H^1\) 的 unitary Fourier 变换，并使用 \(L^1*L^2\to L^2\) 及有限测度的 \(L^2\) Young 界，得
\[
 K*Lu=a(\mu*u-u)\quad\text{于 }L^2.
\]
令 \(K_r(x)=r^{-n}K(x/r)\)、\(\mu_r=(z\mapsto rz)_\#\mu\)。精确缩放为
\[
 LK_r=\frac ar(\mu_r-\delta_0),\qquad
 K_r*Lu=\frac ar(\mu_r*u-u).
\]
因子 \(1/r\) 不能漏掉。

先取非负、有界、紧支撑 \(f\)，以 \(\kappa=\lambda/M_K\) 作上一节的障碍分解。对每个固定 \(r>0\)，在同一个 \(D^c=\{u=0\}\) 上，\(\mu_r*u-u\ge0\) 几乎处处。即使 \(\mu_r\) 低维奇异，\(dx\,d\mu_r\)-Fubini 仍保证卷积的合法性、对代表的无关性以及非负性。于是
\[
 K_r*f=K_r*\sigma-K_r*Lu\le K_r*\sigma\le\kappa M_K=\lambda
 \quad\text{a.e. 于 }D^c.
\]
先对正有理 \(r\) 合并零集，再由
\[
 \|K_r*f-K_s*f\|_\infty
 \le\|f\|_\infty\|K_r-K_s\|_1\to0
\]
及 \(L^1\) dilation 连续性，扩展到全部正实 \(r\)。严格坏集模零集包含于同一个 \(D\)，从而
\[
 \lambda|\{\sup_{r>0}K_r*f>\lambda\}|
 \le\lambda|D|\le M_K\|f\|_1.
\]

一般 \(f\in L^1\) 取
\[
 f_j=\min(|f|,j)1_{[-j,j]^n}\uparrow|f|.
\]
非负卷积的单调收敛及 \(\sup_r\sup_j=\sup_j\sup_r\) 给
\[
 \sup_rK_r*f_j\uparrow\sup_rK_r*|f|.
\]
每级的上确界由有理尺度表达，因而可测，极限同样可测。严格超水平集和质量单调收敛，得到一般 \(L^1\)、全部共同尺度的同一弱型界。没有为一般 \(L^1\) 输入假定未证明的尺度一致连续，也没有为不同 \(r\) 重选 \(D\)。不同逼近输入 \(f_j\) 可以有不同的障碍集。

取 \(Q=[-1,1]^n\)，若 \(0<g_0\le\operatorname{essinf}_QG\)，\(0<c<g_0\)，置
\[
 N=2^n(g_0-c),\qquad K=H_c/N.
\]
第 4 节给 \(a=1/N\) 及所需有限正源，而
\[
 K\ge 2^{-n}1_Q\quad\text{a.e.}
\]
给逐尺度 cube 平均支配。非负 Lebesgue 卷积不受核在零测度集上重定义影响，因此 essential 下界已足够。于是中心轴平行、各坐标使用同一半径的 cube 极大算子满足 \(C_n\le M_K\)。固定 Lebesgue 背景用于自伴卷积、质量守恒、Fourier、测度卷积 Fubini、\(K*\sigma\le\kappa\int K\) 及 \(\kappa|D|\le\int\sigma\)；没有转移到任意背景或独立矩形尺度。

## 7. 八维质量公式和 essential 面接口

真实缩放给 \(G(rx)=r^{-7}G(x)\)，包括无穷值。由弱 \(L^{8/7}\)，\(V_1=|\{G>1\}|<\infty\)。置
\[
 J=\sum_{F\subset\partial Q}\int_FG(y)^{8/7}\,dS(y).
\]
沿 cube 射线 \(x=ry\) 的 Jacobian 为 \(r^7drdS(y)\)，不同开面的锥只有零体积交界。Tonelli 给
\[
 V_1=\frac18J,\quad
 |\{G>s\}|=V_1s^{-8/7},\quad
 \int H_c=7V_1c^{-1/7}=\frac78Jc^{-1/7}.
\]
特别 \(J<\infty\)，并未先假定面迹可积。

正时间 \(p_t\) 在紧 cube 上严格正且联合连续，故 \(\min_{x\in Q}\int_1^2p_t(x)dt>0\)。局部可积性又使 essential cube 下确界有限。其与 essential boundary 下确界相等：若正面测度的一组 \(y\in\partial Q\) 满足 \(G(y)<a<g_{\rm cube}\)，则取 \(r<1\) 充分接近 1 并用齐次性，会在 cube 内生成正体积且 \(G<g_{\rm cube}\) 的点集，矛盾。严格说先把这组点按 \(G(y)\le a\) 及一个固定的接近 1 的 \(r\) 区间截取。反向由 \(r\le1\) 的齐次性直接得出。下半连续性自身不足以把内部 a.e. 下界升级为边界逐点下界。

在本节的有限 frame 下，还可用第 3 节升级为逐点 cube 下界：补集的真实连续性排除低于 essential 下界的点，\(S_3\) 上为无穷。

固定 \(g_0\) 后，
\[
 \int K=\frac{7J}{8\cdot256}\frac{c^{-1/7}}{g_0-c}.
\]
对数导数是 \((8c-g_0)/(7c(g_0-c))\)，唯一标量最优点为 \(c=g_0/8\)，于是
\[
 C_8\le\frac{8^{1/7}J}{256g_0^{8/7}}.
\]
若模型在坐标反射下不变且其对称群传递地作用于全部面，令
\[
 I=\int_{[0,1]^7}G(1,y)^{8/7}dy.
\]
则 \(J=16\cdot128I=2048I\)，从而
\[
 C_8\le8^{8/7}I_0/g_0^{8/7},\qquad
 C_8^7\le8^8I_0^7/g_0^8
\]
只要有同一真实模型的连续证书 \(I\le I_0\)。一般有限 frame 没有这些对称性时，必须保留 \(J\)，不能擅自换成单面 \(I\)。

## 8. 固定面短时预算与长时预算

### 8.1 全平面峰值与固定切片质量

令 \(J_1=\{j:(v_j)_1\ne0\}\)、\(D_1=\sum_{j\in J_1}|(v_j)_1|>0\)。对每个 \(j\in J_1\)，选含 \(v_j\) 的方向基 \(B_j\)，记
\[
 b_j=\frac{c_j}{\pi^8|\det B_j|\prod_{i\in B_j\setminus\{j\}}c_i},
 \qquad B=D_1^2\sum_{j\in J_1}b_j.
\]
对任意 \(0<\beta<1\)，若所有 \(j\in J_1\) 的系数满足 \(|s_j|\le\beta/D_1\)，则
\[
 \left|\sum_js_j(v_j)_1\right|\le\beta.
\]
故在开区域 \(|x_1|>\beta\) 中，至少一个法向非零系数绝对值大于 \(\beta/D_1\)。对该项用受限密度上界 \(c_jtD_1^2/(\pi\beta^2)\)，基中其余七项用 \(1/(\pi c_it)\)，余项概率卷积，事件并集给
\[
 p_t(x)\le\beta^{-2}B t^{-6}\quad\text{a.e. 于 }\{|x_1|>\beta\}.
\]
先以真实密度连续性升级为逐点界，再于 \(x_1=1\) 令 \(\beta\uparrow1\)，得
\[
 \boxed{\sup_{y\in\mathbb R^7}p_t(1,y)\le Bt^{-6}.}
\]
全过程都在空间密度层面进行，没有将零概率平面截面当成正概率事件。

第一坐标的真实一维边际是尺度 \(At\) 的 Cauchy，其中
\[
 A=\sum_jc_j|(v_j)_1|>0.
\]
Tonelli 及边际分布识别只直接给几乎处处 \(s\) 的恒等式
\[
 \int_{\mathbb R^7}p_t(s,y)dy=\frac{At}{\pi(s^2+A^2t^2)}.
\]
**不能仅凭该 a.e. 恒等式代入 \(s=1\)。** 固定 \(t>0\)，选属于恒等式满测度集的 \(s_k\to1\)。\(p_t\) 连续，Fatou 给
\[
 \boxed{\int_{\mathbb R^7}p_t(1,y)dy
 \le\liminf_k\int p_t(s_k,y)dy
 =\frac{At}{\pi(1+A^2t^2)}\le\frac A\pi t.}
\]
该上界已足够，并不声称仅凭 Fatou 得到等号。

### 8.2 \(L^{8/7}\) 短时范数及 \(8/9\) 系数

对非负 \(h=p_t(1,\cdot)\)，由 \(h^{8/7}\le\|h\|_\infty^{1/7}h\)，
\[
 \|p_t(1,\cdot)\|_{8/7}
 \le\|p_t(1,\cdot)\|_\infty^{1/8}
       \|p_t(1,\cdot)\|_1^{7/8}
 \le B^{1/8}(A/\pi)^{7/8}t^{1/8}.
\]
关键时间幂为 \((-6)/8+7/8=1/8\)，不是 \(0\) 或负数。设
\[
 T_\varepsilon(1,y)=\int_0^\varepsilon p_t(1,y)dt.
\]
先对有限时间与有限空间区域用 Minkowski，再单调极限，得到全平面预算
\[
 \boxed{\|T_\varepsilon(1,\cdot)\|_{L^{8/7}(\mathbb R^7)}
 \le\frac89 B^{1/8}(A/\pi)^{7/8}\varepsilon^{9/8}.}
\]
这包括全部 \(S_3\) 的面交集邻域和空间无穷尾。该奇异集在面内零七维测度，其无穷点不会破坏这个有限范数。

若另有逐自由坐标反射对称性，则正象限的 \(8/7\) 次积分是全平面的 \(2^{-7}\)，范数因子为 \(2^{-49/8}\)。对单位正面可使用该改进，但这需模型对称性，不能对一般 frame 自动加入。

### 8.3 长时：有限面域与全平面应分开

任意 \(x\in\mathbb R^8\)，
\[
 R_T(x)=\int_T^\infty p_t(x)dt\le\frac C{7T^7}.
\]
因此对有限七维测度的面域 \(E\subset\mathbb R^7\)，
\[
 \boxed{\|R_T(1,\cdot)\|_{L^{8/7}(E)}
 \le |E|^{7/8}\frac C{7T^7}.}
\]
特别单位正面 \(|E|=1\) 可直接使用；整个面 cube 的面积为 \(2^7\)。**不能把统一点态常数直接在全平面上积分而声称同一个有限预算**，因为全平面测度无穷。

如确实需要全平面长时范数，可另用固定切片质量
\[
 \|p_t(1,\cdot)\|_1\le\frac1{\pi At}
\]
和全空间峰值 \(Ct^{-8}\)，得到
\[
 \|p_t(1,\cdot)\|_{8/7}\le C^{1/8}(\pi A)^{-7/8}t^{-15/8},
\quad
 \|R_T(1,\cdot)\|_{L^{8/7}(\mathbb R^7)}
 \le\frac87C^{1/8}(\pi A)^{-7/8}T^{-7/8}.
\]
这个独立全平面预算的衰减率与紧面域的 \(T^{-7}\) 不同。固定切片质量在这里同样只用已证的 Fatou 上界。

## 9. 十六方向 \(4+4\) 直和模型的严格规格

令
\[
 W=\{(1,\epsilon_2,\epsilon_3,\epsilon_4):\epsilon_i\in\{-1,1\}\},
\quad V=(W\times\{0\}^4)\cup(\{0\}^4\times W),\quad c_v=1.
\]
这是十六个长度 \(2\) 的未单位化方向，两个四维块各八条 body 方向。以 \(x=(x',x'')\) 写
\[
 \psi_8(\xi',\xi'')=\psi_4(\xi')+\psi_4(\xi''),\quad
 p_t^{(8)}(x',x'')=p_t^{(4)}(x')p_t^{(4)}(x''),
\]
\[
 \boxed{G_8(x',x'')=\int_0^\infty
 p_t^{(4)}(x')p_t^{(4)}(x'')dt.}
\]
两块共享同一个 \(t\)。这个表达式不是 \(G_4(x')G_4(x'')\)，后者有两个独立时间积分。不能据此将四维弱型常数平方，不能把四维的 \(g_4,I_4\) 直接代作八维证书。

取四维 Hadamard 基
\[
 H_4=\begin{pmatrix}
 1&1&1&1\\1&-1&1&-1\\1&1&-1&-1\\1&-1&-1&1
 \end{pmatrix},\qquad |\det H_4|=16.
\]
八维块对角基行列式绝对值 \(256\)，因此可用
\[
 C=\frac1{256\pi^8},\quad
 \|p_t^{(8)}\|_\infty\le\frac1{256\pi^8t^8},\quad
 R_T\le\frac1{1792\pi^8T^7}.
\]
每块 Gram 为 \(8I_4\)，全 frame Gram 为 \(8I_8\)。范数界例如 \(\sqrt8|\xi|\le\psi_8(\xi)\le8\sqrt2|\xi|\) 成立。该粗谱界不承担任何数值 \(g/I\) 认证。

第一面上恰有八条方向的首坐标为 \(1\)，所以 \(D_1=A=8\)。每条可由块内坐标反射映射到 \(H_4\) 的一列，再带同一变换的其余列补成 \(|\det|=16\) 的四维基；另一块任选 Hadamard 基。因此每个 \(B_j\) 可取 \(|\det B_j|=256\)。第 8 节给
\[
 B=8^2\frac8{256\pi^8}=\frac2{\pi^8},\quad
 \|T_\varepsilon(1,\cdot)\|_{8/7,\mathbb R^7}
 \le\frac89\,2^{1/8}8^{7/8}\pi^{-15/8}\varepsilon^{9/8}.
\]
自由坐标反射允许单位正面再乘 \(2^{-49/8}\)。两块内置换、反射和交换两块使全部十六个面等价，故本模型可以用第 7 节的 \(J=2048I\) 和 \(C_8^7\le8^8I_0^7/g_0^8\)。这只给正确的证书接口；本轮未产生 \(g_0,I_0\) 的数值证书。

## 10. 来源、计算和未完成部分

来源只读自冻结包的 `paper_finite_source.md`、`audit.md`、`papers/n3/main.tex`、`papers/n4/main.tex`，以及 clone 和来源根的 `AGENTS.md`。两篇正文的相关分析和短时证明已重读。八维的 \(S_3\) 分层、四独立大系数论证及 \(\varepsilon^{9/8}\) 预算由本文重新证明；不是把三四维的方向线结论按名称迁移。

专用 clone 是 `git clone --no-hardlinks --no-checkout` 从本地七维 clone 建立；HEAD 为 `79f8acb20c418e60c5196a3c830103253e34c900`。稀疏检出报告一项与本轮无关的缺失历史 blob：`48d436008eaeecca94834fc826350c6e94e26490`，路径 `centered-cube-upper-bound/adversarial_shell_pressure_20261007/adversarial_shell_search.py`。`AGENTS.md` 已正常检出；本轮仅创建新 `green/` 目录，未读取或修改该无关路径，未提交或推送。

初次本机 `python3` 尝试导入 SymPy 得到 `ModuleNotFoundError`。随后使用标准库 `Fraction` 的准确消元和有理指数运算，不安装依赖。`check.py` 的实际执行核查十六方向 Gram、Hadamard/块基行列式、每条法向方向的补基、源缩放幂、短时/长时积分幂、质量公式和成本幂；结果在 `checks.json`。未调用 MATLAB，未计算 Green 节点、未开展优化或新的大批实验。有限检查的通过仅支持这些明示有限断言；\(S_3\) 分类、源正性及全尺度接口的严格依据是上文证明。

剩余工作明确为：针对所选八维真实模型构造连续 \(g_0\) 下界与 \(I_0\) 上界证书，以及由主代理完成整体证明的最终独立审阅。本文件不宣称该数值阶段已完成。

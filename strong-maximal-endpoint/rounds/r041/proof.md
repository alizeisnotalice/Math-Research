# 第41轮证明：导数 Mills、中心化桥与正中偏差单步比较

日期：2026-11-06。以下证明使用实际有限 N 终端；最后一节只闭合正侧中间带，不声称全空间上界已经完成。

## 1. 记号及精确终端权

\(\mu=1+e^{-\mu}\)，\(p_c=(1+e^c)^{-1}\)，\(Z_c=\mathrm{Exp}(1)+c\mathrm{Bern}(p_c)\)。记密度 q_c、矩母函数 M_c、\(\Lambda_c=\log M_c\)、\(K_ch(r)=Eh(r-Z_c)\)。直接卷积给
\[
 q_c(z)=\frac{e^{-z}}{1+e^{-c}}(1+\mathbf1_{z\ge c}),\quad z\ge0,
 \qquad q_c(z)\le2e^{-z}.                                      \tag{1}
\]
令 a 在 mu 的固定正紧邻域，\(y=\Lambda_a'(\eta)\)、\(\nu=\Lambda_a''(\eta)\)、\(\kappa_3=\Lambda_a'''(\eta)\)，\(0<\eta\le\delta_0<1/2\)。\(f_{t,a}\) 为参考和密度，\(v_{t,N}=K_a^t\phi_N\)。

定义 m(c)=E Z_c=1+c p_c，且 phi_N(u)=exp(−u)（u>0）、phi_N(u)=1−u/N（u≤0）。实际 payoff 的精确倾斜总质量为
\[
 h_N(\eta)=\int\phi_N(u)e^{\eta u}du
 =\frac1{\eta(1-\eta)}+\frac1{N\eta^2}.                         \tag{2}
\]
归一化终端权 \(d\rho_{N,\eta}=\phi_N(u)e^{\eta u}du/h_N\) 是正侧 Exp(1−eta)、负侧 Exp(eta)、负侧 Gamma(2,eta) 的正混合。因此对每个固定 j≥0，
\[
 E_\rho|U|^j\le(j+1)!\eta^{-j} \quad(\eta\le1/2),              \tag{3}
\]
对全部 N≥1 一致；无需假设 N eta 很大。eta 微分归一化权一次，只增加一次对应的加权矩和 O(eta^(−1)) 系数。

## 2. 带参数与时间导数的统一 Mills 引理

设 g_t 为参考增量经 eta 倾斜、中心化后的和密度。对充分大 t、\(t\eta^2\ge A^2\log t\)，定义
\[
 \Delta_N(t,\eta,a)=\frac{v_{t,N}(ty)}{f_{t,a}(ty)h_N(\eta)}-1.
\]
则
\[
 |\partial_t\Delta_N|\le C/(t^2\eta^2),\qquad
 |\partial_\eta^j\partial_a^q\Delta_N|
 \le C_{j,q}/(t\eta^{2+j})\quad(j+q\le3).                     \tag{4}
\]
常数对 N≥1 和 a 的该紧邻域一致。时间导数是下面实时间 Fourier 延拓的导数，N、a、eta 固定；整数时间复原真实参考值。

**Fourier 证明。** 中心特征函数为
\[
 \psi(u)=e^{-iuy}\frac{1-\eta}{1-\eta-iu}
 \frac{(1-p_a)+p_ae^{a(\eta+iu)}}{(1-p_a)+p_ae^{a\eta}}.
\]
由于 \(1-p_a-p_ae^{a\eta}>0\)，Bernoulli 圆盘不包围零；指数因子也无零。取共轭对称连续对数 L，用 exp(tL) 定义实 t 延拓。它不必处处非负，但中心展开给 \(g_t(0)\asymp t^{-1/2}>0\)。近零有 |psi|≤exp(−c u²)，外频 |psi|≤rho<1，且指数因子给统一 L² 可积界。参数导数在中心化后满足 \(\partial_\eta L,\partial_a L=O(u^2)\)，所以固定阶参数微分保持同阶 Gaussian Fourier 界；外频微分仅增加固定次多项式权，取固定充分大 t 即可积分。

由 Fourier 积分及中心奇次抵消得到
\[
 \|g_t''\|_\infty\le Ct^{-3/2},\quad
 \|\partial_tg_t''\|_\infty\le Ct^{-5/2},\quad
 g_t'(0)=O(t^{-3/2}),\quad\partial_tg_t'(0)=O(t^{-5/2}).
\]
中心导数的后两项必须保留 cubic cumulant：Gaussian 主项是奇函数，积分精确为0。仅取绝对值会损失半阶。

令 \(R_t(x)=g_t(x)/g_t(0)\)。商式微分先给 \(\|\partial_tR_t''\|_\infty=O(t^{-2})\)，以及中心 \(\partial_tR_t'(0)=O(t^{-2})\)。因 R_t(0)=1，Taylor 在全部实 x 上给
\[
 |R_t(x)-1|\le C(|x|+x^2)/t,\quad
 |\partial_tR_t(x)|\le C(|x|+x^2)/t^2,                         \tag{5}
\]
固定阶 eta、a 归一化参数导数也有相应 t^(−1) 界。该论证不截断 x 在 log t 内。

精确地 \(\Delta_N=E_{\rho_{N,\eta}}[R_t(-U)-1]\)。以 (3) 积分 (5)，再微分归一化终端权，得到 (4)。这包括有限 call 项，避免将固定正 eta 的尾估计错误地直接用于 eta趋零。

## 3. 精确中心化的密度桥

置 \(\widehat q_a(z)=q_a(z)e^{\eta z}/M_a(\eta)\)。倾斜恒等式给
\[
 \frac{q_a(z)f_{t-1,a}(ty-z)}{f_{t,a}(ty)}
 =\widehat q_a(z)\frac{g_{t-1}(y-z)}{g_t(0)}.                  \tag{6}
\]
中心 Fourier 展开与四阶导数界一致给
\[
 \frac{g_n'(0)}{g_n(0)}=-\frac{\kappa_3}{2\nu^2n}+O(n^{-2}),\quad
 \frac{g_n''(0)}{g_n(0)}=-\frac1{\nu n}+O(n^{-2}),\quad
 \frac{g_n'''(0)}{g_n(0)}=O(n^{-2}),\quad
 \frac{\|g_n^{(4)}\|_\infty}{g_n(0)}=O(n^{-2}).
\]
其中物理第三导数的 Gaussian 奇项同样消失。再用 \(g_{t-1}(0)/g_t(0)=1+1/(2t)+O(t^{-2})\)，全轴 Taylor 给
\[
 \frac{g_{t-1}(y-z)}{g_t(0)}=1+B(z)/t+R(z),\quad
 B=\frac12-\frac{\kappa_3(y-z)}{2\nu^2}-\frac{(y-z)^2}{2\nu},
 \quad |R(z)|\le C(1+z^4)/t^2.                              \tag{7}
\]
这是全局加性展开，不声称有不可能的全局相对误差。\(E_{\widehat q_a}(y-Z)=0\)、二阶矩 nu，故 EB=0；精确 \(K_af_{t-1}=f_t\) 又强制 ER=0。裸的 1/t 质量误差被精确消去。

## 4. 实际有限 payoff 的中心化桥

令 \(H_n(x)=E_{\rho_{N,\eta}}g_n(x-U)\)、\(H(x)=H_{t-1}(x)/H_t(0)\)。则
\[
 \frac{q_a(z)v_{t-1,N}(ty-z)}{v_{t,N}(ty)}=\widehat q_a(z)H(y-z).
                                                                    \tag{8}
\]
(4) 给 H_t(0)≥c/sqrt(t)。由 (3) 及上述物理密度导数界，
\[
 A_1=H'(0)=O((t\eta)^{-1}),\quad A_2=H''(0)=O(t^{-1}),\quad
 A_3=H'''(0)=O((t^2\eta)^{-1}),\quad\|H^{(4)}\|_\infty\le Ct^{-2}.
\]
写 \(X=y-Z\)，精确 \(K_av_{t-1}=v_t\) 把 Taylor 常数项消去，得到
\[
 H(x)=1+A_1x+\frac{A_2}2(x^2-\nu)+\frac{A_3}6(x^3+\kappa_3)
       +[R_4(x)-E_aR_4(X)],\quad |R_4(x)|\le Cx^4/t^2.       \tag{9}
\]

以下动作比较需 \(|a-\mu|\le K_0\eta\)。若 |c−a|≤K eta，则在 \(\widehat q_c=q_ce^{\eta z}/M_c\) 与 \(\widehat q_a\) 之间，X 的一阶矩差为 O(eta²)，二、三阶矩差为 O(eta)，四次多项式加权 TV 为 O(eta)。这些阶数来自显式 (1)：均值的 c 导数在 mu、eta=0 为0；其余矩在紧动作集光滑，而移动跳点区间长度是 |c−a|。质量比 M_c/M_a 不被丢弃，其与1之差为 O(eta³)。于是
\[
 \frac{K_cv_{t-1}}{v_t}=\frac{M_c(\eta)}{M_a(\eta)}[1+O(\eta/t)].    \tag{10}
\]

还需要全动作尾界。由于 \(H_n\le\|g_n\|_\infty\)，(9) 中 H≤C 在全部实轴成立。全部 c 的倾斜矩统一有界（用 q_c≤2e^(−z)、M_c≥1）。因此 (9) 给
\(K_cv_{t-1}/v_t=M_c/M_a+O((t\eta)^{-1})\)，而 M_c/M_a 有统一正下界。对充分大 t，全部动作的 payoff Doob 概率核遂满足
\[
 d\Pi_c(z)=\frac{q_c(z)v_{t-1}(r-z)}{K_cv_{t-1}(r)}dz
 \le C\widehat q_c(z)dz\le Ce^{-(1-\delta_0)z}dz.              \tag{11}
\]
这是包括无限左移的归一化界，不是仅在有界 z 上的密度近似。

## 5. 正侧中间带的全动作单步比较

**定理。** 固定 epsilon>0、alpha=1+epsilon、源界 S<∞。存在 d,A,t0，使匹配 a=a_N 在 N≥t、t>t0 且
\[
 A\sqrt{\log t/t}\le\eta\le d\varepsilon,\quad r=t\Lambda_a'(\eta)
                                                                    \tag{12}
\]
时，对所有 c≥0 都有 \(K_cI_{t-1,N}(r)\le I_{t,N}(r)\)。这里
\[
 U_t=\alpha J_t+s_tf_t,\quad V_t=v_te^{w_t},\quad w_t=U_t/v_t,
 \quad I_t=V_t+\Lambda tb^2/N,
\]
0≤s_t≤S 非减，\(b=m(a)-a\ge0\)、\(\ell=m(a)-2b\)、\(|a-\mu|+b\le K_0\eta\)。本项目匹配序列满足这些条件；b的符号不能删去。

**相位导数与大跳。** 精确式
\[
 J_t/f_t=\frac{p_a}8\{t(y-\ell)^2+(C_t/t)(H_t-1)\}
\]
其中 H_t=∑_{j=1}^t 1/j、J_0=0，J_t=K_aJ_{t−1}+(p_a/8)f_t(r)(r−t ell)²/t²。辅助系数由精确恒等式 f_t C_t=t g_(t+2,t)+a²t p_a(1−p_a)g_(t,t−2)(r−a) 定义；g_(j,k) 是 Gamma(j,1)+a Bin(k,p_a) 的密度（见第36轮式(5)）。这给 C_t/t 及其前三阶 eta 导数有统一界。证明是在中心 Fourier 积分中乘以
\(G_1(\theta)=(1-\theta)^{-2}\)、\(G_2(\theta)=e^{a\theta}/(1-p_a+p_ae^{a\theta})^2\)；乘子及所需局部导数一致有界，外频由指数因子控制。

由 (2)–(4)，\(h_N^{-1}/(1+\Delta_N)\) 的 j 阶 eta 导数为 O(eta^(1−j))、j≤3。因 y−ell=O(eta)，在 eta/2到2eta的局部带有
\[
 |w_r|\le C\eta^2,\qquad |w_{rr}|\le C\eta/t,\qquad
 |w_{rrr}|\le C/t^2.                                        \tag{13}
\]
调和项的导数由 t eta²≥A²log t 吸收；固定有界源也由充分大 t 吸收。

对全部可达前态 r'≤r，还有 \(0\le w_{t-1}(r')\le Ct\eta^3\)。左于参考均值时，v≥c、左倾斜给 \(f(r')\le Ct^{-1/2}e^{-c(r'-tm(a))^2/t}\)；因此二次源≤C(t^(−1/2)+b²sqrt(t))，辅助密度混合的 sup≤C/sqrt(t) 控制调和项。右于均值时，中心局部 CLT 在 |u|≤c sqrt(t) 给
\(f/v\le C(\eta'+t^{-1/2})\)，对所有 0≤eta'≤2eta 成立（在终端负 u 方向积分长度 c/(eta'+t^(−1/2))）。代回 J/f 得上述全局 w 界。

取 L=c0 t eta，使 z≤L 时前态仍在 (13) 的稍宽带。用 (11)，局部振荡 O(eta²(1+z)) 给二阶指数矩 O(eta⁴)；z>L 部分有归一化上界
\(\operatorname{poly}(t,\eta)\exp[-c t\eta+C t\eta^3]=o(\eta^4)\)。所以对所有动作，虽 w 的绝对水平不小，仍有
\[
 \log E_{\Pi_c}e^{w_{t-1}}\le E_{\Pi_c}w_{t-1}+C\eta^4.        \tag{14}
\]

**精确参考源。** Poisson 递推及 s 非减给
\[
 E_{\Pi_a}w_{t-1}-w_t=-\gamma_t/v_t,\qquad
 \gamma_t=U_t-K_aU_{t-1}\ge\alpha\frac{p_a}8 f_t(y-\ell)^2.
\]
置 \(D_{ex}=p_a\eta(1-\eta)(y-\ell)^2/8\)。因 b≥0，
\(y-\ell=2b+\nu_a\eta+O(\eta^2)\ge c\eta\)，所以 D_ex≥c eta³。
又 \(|h_N^{-1}-\eta(1-\eta)|\le1/N\)，故 (4) 给
\[
 \gamma_t/v_t\ge\alpha D_{ex}-C\eta/t.                       \tag{15}
\]
有限 call 的误差是 O(eta²/N)，已被 eta/t 控制。

**附近动作及精确系数。** 由 (9)，payoff Doob 核在 |c−a|≤K eta 下的一阶中心矩变化仍为 O(eta²)：A1乘方差变化只增加 O(1/t)≤C eta²；二阶矩变化 O(eta)。对 (13) 作物理位移 Taylor，结合 (10)、(14)–(15)，得
\[
 \log(K_cV_{t-1}/V_t)
 \le\Lambda_c(\eta)-\Lambda_a(\eta)-\alpha D_{ex}
       +C(\eta/t+\eta^4).                                  \tag{16}
\]
令 a=mu+u eta、c=mu+v eta，u、v 在固定有界集。显式 M_c 的联合 Taylor 给
\[
 [\eta^3]\{\Lambda_c-\Lambda_a-D_{ex}\}
 =-\frac{\mu-1}{8\mu}(\mu-2v)^2\le0.                        \tag{17}
\]
这是本轮 Wolfram 精确核验的恒等式；解析函数在紧参数集的四阶余项为 C eta⁴。MGF 的全局最优动作因 m(c) 的严格最大值而趋于 mu，驻点方程再给 c−mu=O(eta)，故
\(\sup_c(\Lambda_c-\Lambda_a)\le D_{ex}+C\eta^4\)。式 (17) 的控制中心是 c−mu，不是 c−a。

**远动作排除。** 为使 (16) 真正覆盖全部 c，记 d_c=|c−mu|。在 K eta≤d_c≤d0 中，均值二次损失及方差 Lipschitz 界给 \(\Lambda_c-\Lambda_a\le-c\eta d_c^2\)（先取 K 大）。由 (9) 的相同中心矩计算，payoff 比值误差≤
\[
 C\{d_c^2/(t\eta)+d_c/t+\eta/t\}.
\]
Doob 一阶矩变化≤C(d_c²+eta²)，二阶矩变化≤C(d_c+eta)。式 (13) 的相位期望误差≤
\[
 C\eta^2(d_c^2+\eta^2)+C\eta(d_c+\eta)/t+C/t^2.
\]
这些项由 eta d_c² 吸收：使用 d_c≥K eta、t eta²≥A²log t，先取 K 大，再缩小 eta 上界并增大 t0。离开固定动作邻域后，m(c)≤mu−c_* 一致成立（c→0或∞时 m→1）；MGF 损失−c eta 支配桥误差 O(1/(t eta))及相位误差 O(eta²)。此处 (1) 使全部 c 的加权矩统一有界，没有把远动作换成参考控制。

由 (16)–(17)，附近动作满足
\[
 \log(K_cV_{t-1}/V_t)\le-c\varepsilon\eta^3+C(\eta/t+\eta^4)<0
\]
（先取 d 小，再取 A、t0 大）；远动作也已排除。最后
\(K_cI_{t-1}-I_t=K_cV_{t-1}-V_t-\Lambda b^2/N\le0\)。证毕。

## 6. 拼接范围与尚未证明的部分

结合 r036 自然管及 r040 正侧局部剪切，常数稍作缩窄即可覆盖正侧剩余区域；mu 参考 cutoff 与匹配 a 的倾斜坐标差由率连续性和中心带覆盖处理。该点态全动作比较允许后续反馈选择，不依赖参考桥路径占据假设。

早期左侧中间带仍需证明真实 clipped cost 的导数 buffer 与全动作最小化，不能以已在 r037 闭合的晚时段 floor 定理替代。故本轮没有得到全空间锐上界，亦未确定真实 C_n 的 log n 系数。

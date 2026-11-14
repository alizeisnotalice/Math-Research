# 零能量附近的维数一致相对费用

本稿使用已审 r721—r723 的指定实际角向源、穿点矩形与共同动作中位数最大点。它补足 r723 的能量趋零缺口；不处理无界中心化频率能量，不给一般强极大端点的系数。记 $c=c_n$、$r=r_n$、$q=(1+e^c)^{-1}$，$Z=W+cJ$，$M(c)=1+cq$，$v(c)=1+c^2q(1-q)$。$S_l$ 是 $l$ 个共同动作增量之和，密度为 $f_l$。沿用
\[
\varphi_c(z)=\frac{1-q+qe^{icz}}{1-iz},\qquad
\psi_c(z)=e^{-iM(c)z}\varphi_c(z),\qquad p_*=e^{-\mu}/\mu.
\]
已证 $c_n=\mu+O(n^{-1})$、$r_n=nM(c_n)+O(1)$，全动作负分位 Hessian 为 $H=(h-b)I+b\boldsymbol1\boldsymbol1^T=p_*I+O_{\rm op}(n^{-1})$，其中 $h=qf_{n-1}(r-c)/f_n(r)$。以下只需充分大 $n$。

## 1. 结论与精确微能量斜率

存在绝对常数 $\epsilon_0>0,C<\infty,n_0$，对 $n\ge n_0$ 以及所有实切向频率 $\sum_i\xi_i=0$、$E=\sum_i\xi_i^2\le\epsilon_0$，全轴最优二阶系数满足一致展开
\[
\mathcal K_n(\xi)=\frac14(\kappa_n-a_n)E+O(CE^2),\qquad
\kappa_n=\frac{h^2}{h-b}\left(\frac{nc-r}{n-1}\right)^2=O(n^{-2}),\qquad a_n\to\mu. \tag{1}
\]
因此缩小 $\epsilon_0$、增大 $n_0$ 后，
\[
\mathcal K_n(\xi)\le-\frac{\mu}{16}E \quad(0<E\le\epsilon_0,n\ge n_0). \tag{2}
\]
结合 r723 的能量环带定理：对每个固定 $K>0$，存在 $N_K,c_K>0$，使全部 $n\ge N_K$、$0<E\le K$ 的切向频率均有 $\mathcal K_n(\xi)\le-c_KE$。零向量的形式二阶系数为零。公共频移不变性把该结论推广到任意非恒定频率向量，其 $E$ 取中心化能量 $\sum_i(\xi_i-\bar\xi)^2$。

## 2. 复切向球上的统一 Fourier 界

取固定 $0<\delta<1/8$，在复切向球
\[
\mathcal B_\delta=\{\zeta\in\mathbb C^n:\sum_i\zeta_i=0,\ \sum_i|\zeta_i|^2\le\delta^2\}
\]
上工作。$\varphi_c$ 的分子零点全部位于 $\operatorname{Im}z=-1$，分母极点为 $z=-i$：这是由 $q/(1-q)=e^{-c}$ 直接得到的。因此 $\psi_c$ 在 $|\operatorname{Im}z|<1$ 无零且解析，可取 $\ell_c(z)=\log\psi_c(z)$、$\ell_c(0)=0$。当 $c$ 位于 $\mu$ 的固定紧邻域时，$\ell_c''$ 在整条闭带 $|\operatorname{Im}z|\le2\delta$ 上一致有界；写成
\[
\ell_c(z)=-iM(c)z+\log(1-q)+\log(1+e^{-c}e^{icz})-\log(1-iz)
\]
即可验证，两个对数的分母均一致远离零。对每个实 $u$ 沿复线段作 Taylor 展开，切向性消去一次项，给出**全实轴**的界
\[
\left|\prod_j\psi_c(u+\zeta_j)\right|
\le e^{C\sum_j|\zeta_j|^2}|\psi_c(u)|^n
\le e^{C\delta^2}(1+u^2)^{-n/2}. \tag{3}
\]
没有随 $n$ 增长的虚部损失。条件相位均值仍由中心化 Fourier 反演给出
\[
m(\zeta)=\frac1{2\pi f_n(r)}\int e^{-iu(r-nM(c))}\prod_j\psi_c(u+\zeta_j)\,du. \tag{4}
\]
因为 $f_n(r)\asymp n^{-1/2}$，(3) 给 $|m(\zeta)|\le C$，且在被积函数绝对值上插入 $|u|^a$ 后、除以 $f_n(r)$，积分不超过 $Cn^{-a/2}$，$a=0,1,2$。这些常数独立于 $n$ 与 $\zeta\in\mathcal B_\delta$。

记 $m_{i,t}(\zeta)$ 为设置 $Z_i=t$、其余坐标条件于和 $r-t$ 后的相位均值。$c$ 面的 Fourier 比值为
\[
R_i(u)=\frac{e^{ic(u+\zeta_i)}}{\varphi_c(u+\zeta_i)}
=\frac{e^{i(c-M(c))(u+\zeta_i)}}{\psi_c(u+\zeta_i)}.
\]
因为 $c-M(c)=O(1/n)$、$\psi_c'(0)=0$，在固定小 $|u|$ 区间有
\[
|R_i(0)-1|\le C(|\zeta_i|^2+|\zeta_i|/n),\quad
|R_i(u)-R_i(0)|\le C(|u||\zeta_i|+u^2+|u|/n).
\]
区间外 $|R_i(u)|\le C(1+|u|)$，可由 (3) 吸收。用 $f_{n-1}(r-c)/f_n(r)=1+O(1/n)$ 归一化，得到
\[
m_{i,c}(\zeta)=R_i(0)m(\zeta)+e_i(\zeta),\quad
|e_i(\zeta)|\le C(|\zeta_i|/\sqrt n+1/n). \tag{5}
\]
特别地，$m_{i,c}$ 对全部 $i$ 一致有界。$0$ 面的比值为 $1/\varphi_c(u+\zeta_i)$，同样有界并满足一致多项式尾控制，故 $m_{i,0}$ 也一致有界。

相位均值的总动作梯度在驻点处为 $g_i(\zeta)=h[m(\zeta)-m_{i,c}(\zeta)]$。由 (5)，
\[
\|g(\zeta)\|_2\le C\left(\delta^2+\frac{\delta}{n}+\frac{1+\delta}{\sqrt n}\right)\le C. \tag{6}
\]
这一步必须在 $\ell^2$ 中合并每轴误差；仅分别证明各面有界会错误地产生 $\sqrt n$。对固定 $n$，这些函数都是整个复频率空间上的整函数：条件于 $S=r$ 或任一面后各坐标有界，故条件指数矩在任意紧复集一致收敛。复球内的 Fourier 公式由原积分直接延拓，也可由解析唯一性得到。

## 3. 精确线性项与无维数的解析余项

对实切向单位向量 $e$，条件交换性给 $E Z_j\mid S_n=r=r/n$，因此 $\partial_z m(ze)|_{z=0}=0$。在 $c$ 面，
\[
E[e\cdot Z\mid Z_i=c,S_n=r]
=e_i\left(c-\frac{r-c}{n-1}\right).
\]
于是严格有
\[
g(0)=0,\qquad Dg(0)e=-ih \eta_n e,\qquad \eta_n=\frac{nc-r}{n-1}=O(1/n). \tag{7}
\]
由 (6) 对 $z\mapsto g(ze)$ 应用 Banach 值 Cauchy 估计，在 $|z|\le\delta/2$ 得
$\|g(ze)-zDg(0)e\|_2\le C|z|^2$，常数与 $e,n$ 无关。因此也有 $\|g(\xi)\|_2^2\le C(E/n^2+E^2)$。

为获得更精确的统一四阶余项，将面费用作全纯复化。写 $w_0=(1-q)f_{n-1}(r)/f_n(r)$，定义
\[
\mathcal D(\zeta)=w_0[1-m(-\zeta)m_{i,0}(\zeta)-m(\zeta)m_{i,0}(-\zeta)+m(\zeta)m(-\zeta)]
\]
\[
\hspace{18mm}+h[1-m(-\zeta)m_{i,c}(\zeta)-m(\zeta)m_{i,c}(-\zeta)+m(\zeta)m(-\zeta)]. \tag{8}
\]
对实频率它等于 $D_\xi/f_n(r)$；由解析唯一性，它与所选 $i$ 无关，且是偶函数、置换不变函数。也可直接用 $(\partial_s+1)[f_n(s)(1-m_s(\zeta)m_s(-\zeta))]/f_n(r)$ 定义同一函数。上述复球界使 $|\mathcal D|\le C$。其常数项为零，置换不变二次型在切向空间上必为 $a_n\sum_i\zeta_i^2$（置换不变矩阵只有一个对角元、一个非对角元）。偶性和一维 Cauchy 估计给
\[
\mathcal D(\xi)=a_n E+O(CE^2), \tag{9}
\]
统一于所有方向与维数。

进一步定义 $\mathcal S(\zeta)=g(-\zeta)^T H^{-1}g(\zeta)$。它全纯、偶且在复球上一致有界；实频率时就是 $g^*H^{-1}g$。由 (7)，其切向二次项为 $h^2\eta_n^2/(h-b)$，因为 $H^{-1}$ 在切向空间为标量 $1/(h-b)$。同样的偶函数 Cauchy 估计给 $\mathcal S(\xi)=\kappa_n E+O(CE^2)$。与 (9) 相减即为 (1) 的展开，且没有需逐轴求和的四阶常数。

## 4. 二次面费用的正极限及有限和公式

取 $e=(1,-1,0,\ldots,0)/\sqrt2$，在 (8) 选择第三轴作面。全层与两个面上 $T=(Z_1-Z_2)/\sqrt2$ 的均值均为零；展开 $m(ze)=1-z^2ET^2/2+O(z^3)$ 及面均值可见，每个方括号的二次系数精确为该面上的 $ET^2$，全层的二阶矩抵消。该纯代数抵消另经本机 Wolfram 核验，代码与原输出见 `face_quadratic_check.wl`、`face_quadratic_check_output.txt`。因此令 $N=n-1$，
\[
a_n=w_0\,\widetilde v_N(r)+h\,\widetilde v_N(r-c),\qquad
\widetilde v_N(x)=\frac{N}{N-1}\operatorname{Var}(Z_1\mid S_N=x). \tag{10}
\]
这里使用条件总和固定与交换性：$\operatorname{Cov}(Z_1,Z_2)=-\operatorname{Var}(Z_1)/(N-1)$。

在上述两个有界中心偏移处，任意固定两个保留坐标的条件化密度比 $f_{N-2}(x-Z_1-Z_2)/f_N(x)$ 趋一，且全局一致有界。前者由局部中心极限定理，后者由 Gamma 卷积密度峰值和 $f_N(x)\asymp N^{-1/2}$。因此以单步统一指数矩支配，条件二阶矩收敛至独立未条件二阶矩，$\widetilde v_N(x)\to v(\mu)=\mu$。又 $w_0\to1-q(\mu)$、$h\to q(\mu)$，故 $a_n\to\mu$。

(10) 还可用有限 Gamma–Binomial 和直接计算。给定 $S_N=x$，总 Bernoulli 数 $j$ 的归一化权为
\[
\omega_j(x)\ \propto\ \binom Nj q^j(1-q)^{N-j}\gamma_N(x-cj),\qquad
\gamma_N(t)=\frac{e^{-t}t^{N-1}}{(N-1)!}1_{t>0}.
\]
置 $t_j=x-cj$。条件于 $j,S_N=x$，指数向量在总和 $t_j$ 的单纯形上均匀，Bernoulli 子集在大小 $j$ 的子集中均匀且与该指数向量独立，故
\[
E[Z_1\mid S_N=x]=x/N,\qquad
E[Z_1^2\mid j,S_N=x]=\frac{2t_j^2}{N(N+1)}+\frac{2ct_jj}{N^2}+\frac{c^2j}{N}. \tag{11}
\]
将 (11) 对 $\omega_j$ 加权、减去 $x^2/N^2$ 即可代入 (10)。这是精确可复现公式；本稿不把待运行的数值当作证明。

最后，$a_n\to\mu$、$\kappa_n=O(n^{-2})$ 与统一 $O(E^2)$ 直接推出 (2)。r723 在任意固定环带 $\epsilon_0\le E\le K$ 给统一负常数，除以 $K$ 后与 (2) 拼合，得到第 1 节的整个有界能量球结论。剩余问题是无界中心化能量及从该实际局部族到一般端点的推论；两者均未由本稿解决。

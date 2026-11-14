# 全轴响应：精确组公式、曲率算子极限与弥散频率

本稿沿用 r721 的真实深层角向指标源与实际穿点盒。所有结论是该输入—见证族的局部二阶结论；一般端点的对数系数仍未闭合。记
\[
q(c)=(1+e^c)^{-1},\quad Z_c=W+cJ,\quad M(c)=1+cq(c),\quad
v(c)=1+c^2q(c)(1-q(c)),
\]
其中 $W\sim\mathrm{Exp}(1)$、$J\sim\mathrm{Bern}(q(c))$ 独立。$S_N$ 是 $N$ 个同动作增量之和，密度、分布函数、尾分别记 $f_N,F_N,\bar F_N$，共同动作中位数为 $r_N(c)$。令 $c_n$ 为使 $r_n(c)$ 最大的共同动作、$r=r_n(c_n)$、$q=q(c_n)$。以下讨论充分大 $n$，此时 $c_n$ 为内点，各轴分位梯度为零。置 $\mu=1+e^{-\mu}$、$p_*=e^{-\mu}/\mu$。

## 1. 精确全轴与组频率公式

取切向频率 $\sum_i\xi_i=0$，写
\[
m=E[e^{i\xi\cdot Z}\mid S_n=r],\qquad
m_{i,u}=E[e^{i\xi\cdot Z}\mid Z_i=u,S_n=r].
\]
设 $H=-\nabla_c^2r_n(c_1,\ldots,c_n)$ 为异构动作分位的负 Hessian。对称基准给
\[
H=(h-b)I+b\boldsymbol1\boldsymbol1^T,
\quad h=\frac{q f_{n-1}(r-c_n)}{f_n(r)},\quad
b=\frac{q^2[1/2-\bar F_{n-2}(r-2c_n)]}{f_n(r)}. \tag{1}
\]
注意 $b$ 的分子亦为 $q^2[F_{n-2}(r-2c_n)-1/2]$；改用 CDF 时不能反号。驻点满足 $F_n(r)=F_{n-1}(r-c_n)=1/2$。相位均值的动作梯度为
\[
g_i=h(m-m_{i,c_n}).
\]
若 $\lambda_\perp=h-b>0$、$\lambda_\parallel=h+(n-1)b>0$，则所有 $n$ 个动作均可随相位响应，最优二阶系数为
\[
4C_{\rm all}=\frac{\sum_i|g_i|^2}{h-b}
-\frac{b|\sum_i g_i|^2}{(h-b)[h+(n-1)b]}-\frac{D_\xi}{f_n(r)}, \tag{2}
\]
\[
\frac{D_\xi}{f_n(r)}
=\frac{(1-q)f_{n-1}(r)}{f_n(r)}(1-2\Re(\bar m m_{i,0})+|m|^2)
+h(1-2\Re(\bar m m_{i,c_n})+|m|^2). \tag{3}
\]
式 (3) 对任一 $i$ 相同，来自 r721 的两面方差恒等式。若频率按取值分成组、组大小 $N_a$，则只需组梯度 $g_a$，在 (2) 中用 $\sum_aN_a|g_a|^2$ 与 $\sum_aN_ag_a$。这保留了所有零频轴的响应。

这里给出用于独立数值核验的完全有限公式。对 $N$ 个频率含 $a$ 个 $+t$、$b$ 个 $-t$、$z$ 个零，$a+b+z=N$，令
\[
w_j=q^j(1-q)^{N-j}[u^j](1+e^{itc}u)^a(1+e^{-itc}u)^b(1+u)^z.
\]
写 $\gamma_N(x)=e^{-x}x^{N-1}1_{x>0}/(N-1)!$，$\chi_{a,b,z}(x)$ 为 $W$ 在 $\sum W=x$ 条件下的相位特征函数，则
\[
A_{a,b,z}(s)=\sum_{j=0}^N w_j\gamma_N(s-jc)\chi_{a,b,z}(s-jc),\quad
m=A_{a,b,z}(r)/f_N(r). \tag{4}
\]
对 $x>0$，有
\[
\begin{aligned}
\chi_{a,b,0}(x)&=e^{-itx}\,{}_1F_1(a;N;2itx),\\
\chi_{a,a,0}(x)&={}_0F_1(a+1/2;-t^2x^2/4),\\
\chi_{1,1,N-2}(x)&={}_1F_2(1;N/2,(N+1)/2;-t^2x^2/4),\\
\chi_{0,1,N-1}(x)&={}_1F_1(1;N;-itx).
\end{aligned} \tag{5}
\]
前两式来自 Beta$(a,b)$；第三式将 $D_1-D_2$ 的偶矩代入级数：$E(D_1-D_2)^{2k}=(2k)!\,\Gamma(N)/\Gamma(N+2k)$。最后一式来自 $D_2\sim\mathrm{Beta}(1,N-1)$。删去某轴以后，(4)—(5) 中的 $N$ 必须改为删后的总轴数。具体 $m_{i,c}=e^{i\xi_i c}A_{-i}(r-c)/f_{N-1}(r-c)$、$m_{i,0}=A_{-i}(r)/f_{N-1}(r)$。以上是校准公式，不是浮点结果。

## 2. 全 Hessian 的算子范数极限

**命题。** 在上述共同中位数极大点，
\[
c_n=\mu+O(n^{-1}),\quad h=p_*+O(n^{-1}),\quad
\lambda_\parallel=p_*+O(n^{-1}),\quad b=O(n^{-2}),\quad
\|H-p_*I\|_{\rm op}=O(n^{-1}). \tag{6}
\]
因此充分大 $n$ 时，(2) 的全 $n$ 轴逆矩阵存在且一致有界。

给出可微渐近所需的论证，而不从逐项极限推算子极限。先由统一 Berry–Esseen 和单步矩的统一有界性，$r_n(c)=nM(c)+O(1)$，且最大点 $c_n\to\mu$（端点 $c=0,\infty$ 的 $M=1<\mu$）。在 $\mu$ 的一个固定紧邻域，记
\[
\kappa_3(c)=2+c^3q(c)(1-q(c))(1-2q(c)),\quad a(c)=-\kappa_3(c)/(6v(c)).
\]
对 $c$ 作至多二次参数微分后，仍有统一展开
\[
r_n(c)=nM(c)+a(c)+O_{C^2}(n^{-1/2}). \tag{7}
\]
这里 $O_{C^2}$ 指余项及其前两阶 $c$ 导数。为核验 (7)，使用中心化特征函数
\[
\psi_c(u)=e^{-iM(c)u}\frac{1-q(c)+q(c)e^{icu}}{1-iu},\quad
\log\psi_c(u)=-v(c)u^2/2-i\kappa_3(c)u^3/6+O_{C^2}(u^4).
\]
对 $F_n(nM(c)+x,c)$ 的 Fourier 反演，$x$ 在固定有界区间，令 $u=y/\sqrt n$。在 $|y|\le n^{1/12}$ 逐项展开得到
\[
F_n(nM+x,c)=\tfrac12+\frac{x-a(c)}{\sqrt{2\pi nv(c)}}+O_{C^2}(n^{-1}). \tag{8}
\]
同一误差对 $x$ 的前两次及与 $c$ 的所需混合微分成立（微分增添 $u$ 的幂）。在中间频率，$|\psi_c(u)|\le e^{-\eta u^2}$；在远频率，公共指数因子给 $|\psi_c(u)|\le(1+u^2)^{-1/2}$。两次参数微分至多添 $C n^2(1+|u|)^4$ 的因子，故这两段积分及导数比任何 $n$ 的负幂小。CDF 的 $1/u$ 在零处由奇部相消处理；在中频边界以 $n^{1/12}$ 截断同样有效。将 (8) 乘 $\sqrt n$ 后用隐函数定理、密度的同一展开，得到 (7)。

由驻点 $0=r_n'(c_n)=nM'(c_n)+O(1)$ 与 $M''(\mu)=-p_*<0$，得 $c_n=\mu+O(1/n)$。对称性给严格恒等式
\[
\lambda_\parallel=-r_n''(c_n)/n=-M''(c_n)+O(1/n)=p_*+O(1/n).
\]
为单独核验密度精度，在一阶中心化特征函数展开后，余项为 $O_{C^2}((y^4+y^6)e^{-\eta y^2}/n)$；反演的尺度因子 $n^{-1/2}$ 给绝对 $O(n^{-3/2})$ 密度误差。三阶奇项在有界偏移处的相对贡献仅为 $O(1/n)$。这里是直接反演密度，不是对式(8)的CDF误差形式求导。因此，密度的相同局部展开在有界偏移处给 $f_{n-1}(r-c_n)/f_n(r)=1+O(1/n)$，故 $h=p_*+O(1/n)$。再从 $\lambda_\parallel=h+(n-1)b$ 得 (6)。这里不需要假设整个异构动作空间已全局最优。

本机 Wolfram 的代数微分检查验证了 $M''(\mu)=-p_*$，并另得 $a'(\mu)=0$；后者不用于 (6)，故无需为更高阶定位增加假设。可复现代码与原输出见同目录 `curvature_check.wl`、`curvature_check_output.txt`。

## 3. 全轴响应的两个严格渐近结论

**固定支撑。** 对固定 $d$、固定切向频率 $(t_1,\ldots,t_d,0,\ldots,0)$，现在允许所有 $n$ 轴响应，最优二阶系数仍趋于 r721 的式 (5)：
\[
4C_{\rm all}\longrightarrow \prod_{j\le d} A(t_j)-1+
\sum_{i\le d}B(t_i)\prod_{j\ne i}A(t_j)<0 \tag{9}
\]
（非零模式），其中 $A=|\phi_\mu|^2$、$B=p_*|\phi_\mu-e^{it\mu}|^2$。关键新增估计是每个支撑外轴 $g_i=O_d(1/n)$。为见此，对固定个数的保留坐标，密度条件化因子由 (8) 的密度版本等于 $1+O((1+z^4)/n)$，$z$ 是保留增量的和；可先在 $z\le n^{1/12}$ 展开，其余利用一致指数尾及全局密度比上界截断。因而无条件与条件特征函数之差为 $O(1/n)$，删去一个零频轴、固定其值为 $c_n$ 后亦然，故 $m-m_{i,c_n}=O(1/n)$。支撑外梯度总平方为 $O(1/n)$；(6) 使 $g^*H^{-1}g=p_*^{-1}\|g\|^2+o(1)$，r721 的有限支撑计算即给 (9)。收敛对固定紧频率集合一致。

**弥散频率。** 允许支撑随 $n$ 增长，假设
\[
\sum_i\xi_{n,i}=0,\quad \max_i|\xi_{n,i}|\to0,\quad
\sum_i\xi_{n,i}^2\to\tau^2<\infty.
\]
则所有 $n$ 轴同时响应的系数满足
\[
 C_{\rm all}(\xi_n)\longrightarrow -\tfrac14(1-e^{-\mu\tau^2}). \tag{10}
\]
特别是偶数 $n$ 的 balanced $n/2$ 个 $+t/\sqrt n$ 与 $n/2$ 个 $-t/\sqrt n$，极限为 $-(1-e^{-\mu t^2})/4$，对固定非零 $t$ 严格为负。

证明：记 $A_n(s)$ 为带相位的密度，$\varphi_c(u)=(1-q+qe^{icu})/(1-iu)$。Fourier 反演给
\[
A_n(r)=\frac1{2\pi}\int e^{-iur}\prod_j\varphi_{c_n}(u+\xi_{n,j})\,du.
\]
令 $r=nM(c_n)+x_n$，$x_n=O(1)$，$u=y/\sqrt n$。因频率和为零，线性项消失，且
\[
\sum_j(u+\xi_{n,j})^2=y^2+\sum_j\xi_{n,j}^2,\quad
\sum_j|u+\xi_{n,j}|^3\le C(|y|^3/\sqrt n+\max_j|\xi_{n,j}|\sum_j\xi_{n,j}^2).
\]
小频率中心化特征函数的对数展开与支配收敛得
$m=A_n(r)/f_n(r)\to e^{-\mu\tau^2/2}$。支配是统一的：当所有移位落在小邻域时，乘积模不超过 $e^{-\eta\sum_j(u+\xi_j)^2}\le e^{-\eta n u^2}$；当 $|u|$ 离零有固定距离，最终所有 $\xi_j$ 均小，公共指数因子给指数衰减。后者再分有界及无界区间即可控制全部反演尾。

对第 $i$ 个 $c_n$ 面，其未归一化相位密度的 Fourier 被积函数是在上述乘积上再乘
\[
R_i(u)=\frac{e^{ic_n(u+\xi_i)}}{\varphi_{c_n}(u+\xi_i)}
=\frac{e^{i(c_n-M(c_n))(u+\xi_i)}}{\psi_{c_n}(u+\xi_i)}.
\]
由于 $c_n-M(c_n)=O(1/n)$，小频率下
$|R_i(u)-1|\le C[(u+\xi_i)^2+n^{-1}|u+\xi_i|]$。
以刚才的高斯支配积分，远频率用 $1/|\varphi_c(u)|\le C\sqrt{1+u^2}$（$c$ 近 $\mu>0$），得到
\[
|m_{i,c_n}-m|\le C(n^{-1}+\xi_i^2+|\xi_i|/n),\qquad
\|g\|^2\le C(n^{-1}+\sum_i\xi_i^4)\longrightarrow0. \tag{11}
\]
归一化时用 $f_{n-1}(r-c_n)/f_n(r)=1+O(1/n)$。常数可依赖 $\sup_n\sum\xi_i^2$，不依赖 $i,n$。$0$ 面同理：其 Fourier 比值为 $1/\varphi(u+\xi_i)$，线性项未取消，但仍给 $m_{i,0}-m=O(n^{-1/2}+|\xi_i|)=o(1)$。因此 (3) 的两面系数之和趋一，两面相位均值趋 $m$，得到 $D_\xi/f_n(r)\to1-e^{-\mu\tau^2}$。结合 (6)、(11) 即证 (10)。

## 4. 尚未闭合的范围

全 Hessian 的算子范数问题已由 (6) 解决；(9)、(10) 均允许所有轴响应。仍未证明任意有限 $n$、任意切向频率的 $D_\xi/f_n(r)\ge g^*H^{-1}g$，也未覆盖随维数增加且频率能量无界的全部模式，或有有限大频率与弥散尾并存的一般统一估计。逐轴的 (3) 只能分别控制同一个预算，直接相加会重复支付；不能由这些式子独自宣称全频率 PSD，更不能宣称一般强极大端点已闭合。

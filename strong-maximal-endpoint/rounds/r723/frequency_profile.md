# 有界中心化能量的全频率、全轴二阶判别

沿用已审 r721—r722 的实际角向指标源、穿点矩形及共同动作中位数最大点 $c_n$。写 $r_n$ 为中位数、$q_n=(1+e^{c_n})^{-1}$，$M(c)=1+cq(c)$，$v(c)=1+c^2q(c)(1-q(c))$。记 $\mu=1+e^{-\mu}$、$p_*=e^{-\mu}/\mu$。本稿使用 r722 已证
\[
c_n=\mu+O(n^{-1}),\quad r_n=nM(c_n)+O(1),\quad
H_n=p_*I+O_{\rm op}(n^{-1}),\quad
h_n=q_nf_{n-1}(r_n-c_n)/f_n(r_n)=p_*+O(n^{-1}). \tag{1}
\]
其中 $H_n$ 是所有独立动作的负分位 Hessian。定义
\[
\varphi_c(t)=\frac{1-q(c)+q(c)e^{ict}}{1-it},\quad
\psi_c(t)=e^{-iM(c)t}\varphi_c(t),\quad \psi=\psi_\mu,
\]
\[
A(t)=|\psi(t)|^2=|\varphi_\mu(t)|^2,\qquad B(t)=p_*|\psi(t)-1|^2.
\]
以下所有二阶系数 $\mathcal K_n(\xi)$ 指 r721 的**指定实际见证族**在所有轴均允许侧比随相位微调时的最优对数评分系数。它们不是一般极大算子端点常数。

## 1. 公共频移的精确不变性

在固定分位层 $S=s$ 上，$m_{\xi+a\boldsymbol1}(s)=e^{ias}m_\xi(s)$，因此 $|m|^2$ 及 $D_\xi=(\partial_s+1)[f_n(1-|m|^2)]$ 不变。动作梯度取**总导数**：
\[
\nabla_c\{e^{iar(c)}m_\xi(r(c),c)\}
=e^{iar(c)}[g_\xi+ia m_\xi\nabla_c r(c)].
\]
基准的所有分位梯度为零，故 $g_{\xi+a\boldsymbol1}=e^{iar_n}g_\xi$，进而 $g^*H_n^{-1}g$ 不变。因此
\[
\mathcal K_n(\xi+a\boldsymbol1)=\mathcal K_n(\xi). \tag{2}
\]
实际深层相位均匀化只需频率不是常向量：固定深度的单纯形线性相位只要有两个系数不同就有绝对连续分布；加 $a\boldsymbol1$ 只使该层相位旋转 $as$。故对非恒定频率，上述不变性同样对应 r721 的真实输入与盒族的二阶公式。纯常频率的实际源没有这项角向均匀化，不在此物理结论中。可用唯一均值零代表 $\xi-\bar\xi\boldsymbol1$，其能量是 $\sum_i(\xi_i-\bar\xi)^2$。

## 2. 统一有界能量定理

**定理。** 对每个固定 $K<\infty$，在全部 $\sum_i\xi_i=0$、$E(\xi)=\sum_i\xi_i^2\le K$ 的实频率向量上，一致有
\[
4\mathcal K_n(\xi)=\prod_{i=1}^n A(\xi_i)-1+
\sum_{i=1}^nB(\xi_i)\prod_{j\ne i}A(\xi_j)+o_K(1). \tag{3}
\]
记 $r_*=e^{-\mu}$、$\kappa=(1-2r_*)/\mu>0$，则
\[
4\mathcal K_n(\xi)\le \exp\!\left[-\frac{\kappa E(\xi)}{1+K}\right]-1+o_K(1). \tag{4}
\]
特别对每个 $0<\delta\le K$，全部中心化能量位于 $[\delta,K]$ 的频率，在充分大维数中有统一严格负二阶系数；支撑无需固定，全部 $n$ 轴均可响应。

在零向量处，式(3)–(4)按代数连续延拓取 $\mathcal K_n(0)=0$；实际相位均匀化结论仅用于非零中心化频率。公共频移不变性也只针对驻点处的二阶系数，不宣称完整非线性评分不变。

### 条件 Fourier 估计

置 $x_n=r_n-nM(c_n)=O(1)$，$P_n=\prod_i\psi_{c_n}(\xi_i)$。在 $c$ 的固定邻域内 $q(c)<1/2$，因此 $\varphi_c$ 在实轴无零，$\psi_c$ 有从零开始连续选择的复对数 $\ell_c$。在任何固定紧实区间，其所有所需导数一致有界。因
\[
\ell_c(0)=\ell'_c(0)=0,\quad \ell''_c(0)=-v(c),
\]
且 $\sum_i\xi_i=0$，Taylor 展开给
\[
\sum_i\ell'_c(\xi_i)=O_K(1),\quad
\sum_i\ell''_c(\xi_i)=-nv(c)+O_K(1),\quad
\sum_i|\ell'''_c(\xi_i+u)|\le C_Kn\quad(|u|\le1). \tag{5}
\]
第二式要展开到一次项并用切向性抵消；只用 Lipschitz 界会丢失这一精度。$|P_n|\ge e^{-C_KK}>0$，因为 $\log|\psi_c(t)|\ge-C_Kt^2$。

带相位密度的反演在中心化后是
\[
J_n=\frac1{2\pi}\int_{\mathbb R}e^{-iux_n}\prod_j\psi_{c_n}(u+\xi_j)\,du,
\qquad m=J_n/f_n(r_n). \tag{6}
\]
中心化的额外常数相位精确为 $e^{iM(c_n)\sum\xi_j}=1$。令 $u=y/\sqrt n$，(5) 给
\[
\prod_j\psi_{c_n}(u+\xi_j)
=P_n\exp[-v(c_n)y^2/2+O_K((|y|+|y|^3)/\sqrt n+y^2/n)].
\]
由此及下述统一尾控制可得 $m=P_n+O_K(n^{-1/2})$，并有绝对积分矩界
\[
\frac1{f_n(r_n)}\int |u|^a\left|\prod_j\psi_{c_n}(u+\xi_j)\right|du
\le C_K n^{-a/2},\qquad a=0,1,2. \tag{7}
\]
说明尾控制：取一个固定足够小的 $d>0$，至多 $K/d^2$ 个频率超过 $d$。删去这些因子（其模至多一），其余因子在 $|u|\le d$ 时有 $|\psi_c(u+\xi_j)|\le e^{-\eta(u+\xi_j)^2}$。被删频率之和绝对值至多 $K/d$；配方后得到 $C_K e^{-\eta n u^2/4}$，对充分大 $n$ 成立。对 $|u|$ 离零固定距离，改删 $|\xi_j|>d/2$ 的有限多个因子，公共指数分母在有界区间给指数衰减；在 $|u|>2\sqrt K+2$，用 $|u+\xi_j|\ge|u|/2$ 控制无界尾。这亦允许额外有限次多项式权。中心区可先取 $|y|\le n^{1/12}$ 展开，再以此高斯界截断，故误差对全部能量球一致。这里 $f_n(r_n)\asymp n^{-1/2}$。

又因 $\partial_c\ell_c(t)=O_K(t^2)$，(1) 给
\[
P_n=P+O_K(n^{-1}),\qquad P:=\prod_i\psi(\xi_i),\qquad m=P+o_K(1). \tag{8}
\]
这一步保留了全部复相位；没有假定大频率原子之和为零。

### 所有动作梯度的合并误差

固定第 $i$ 轴于 $c_n$ 的未归一化相位密度，其 Fourier 被积函数等于 (6) 的被积函数乘
\[
R_i(u)=\frac{e^{ic_n(u+\xi_i)}}{\varphi_{c_n}(u+\xi_i)}
=\frac{e^{i(c_n-M(c_n))(u+\xi_i)}}{\psi_{c_n}(u+\xi_i)}.
\]
这里 $c_n-M(c_n)=O(1/n)$。面比值 $R(t)=e^{ict}/\varphi_c(t)$ 精确满足 $R'(0)=i(c-M(c))$；本机 Wolfram 核验及源码见 `face_ratio_check.wl`、`face_ratio_check_output.txt`，并核验了 $R''(0)|_{c=\mu}=\mu$ 与 $(\log\psi_c)''(0)=-v(c)$。在 $|u|$ 小时，$\psi'_c(0)=0$ 与紧区间导数界给
\[
|R_i(u)-R_i(0)|\le C_K(|u||\xi_i|+u^2+|u|/n). \tag{9}
\]
远区使用 $1/|\varphi_c(u+\xi_i)|\le C_K(1+|u|)$，可被上面的尾界吸收。结合 (7) 与 $f_{n-1}(r_n-c_n)/f_n(r_n)=1+O(1/n)$，得
\[
m_{i,c_n}=R_i(0)m+e_{n,i},\qquad
|e_{n,i}|\le C_K(|\xi_i|/\sqrt n+1/n),\qquad
\sum_i|e_{n,i}|^2=O_K(1/n). \tag{10}
\]
因此逐轴误差并未随轴数累积成常数量。$g_i=h_n(m-m_{i,c_n})$，而 $R_i(0)=1/\psi(\xi_i)+O_K(1/n)$。此外 $\sum_i|1-1/\psi(\xi_i)|^2\le C_K\sum_i\xi_i^4\le C_K$。用 (8)—(10) 得
\[
\left\|g-\left(p_*[\psi(\xi_i)-1]\prod_{j\ne i}\psi(\xi_j)\right)_{i=1}^n\right\|_2=o_K(1). \tag{11}
\]
由 (1) 的全算子逆极限，
\[
g^*H_n^{-1}g=\sum_i B(\xi_i)\prod_{j\ne i}A(\xi_j)+o_K(1). \tag{12}
\]

### 面方差与谱结算

选 $i=i_n$ 满足 $|\xi_i|\le\sqrt{K/n}$。(10) 给 $m_{i,c_n}=m+o_K(1)$。$0$ 面的 Fourier 比值为 $1/\varphi_{c_n}(u+\xi_i)$，其与一之差在中心区为 $O_K(|u|+|\xi_i|)$，故由 (7) 也有 $m_{i,0}=m+o_K(1)$。两个面密度系数分别趋 $q(\mu)$ 与 $1-q(\mu)$，于是精确两面方差公式给
\[
D_\xi/f_n(r_n)=1-|m|^2+o_K(1)=1-\prod_i A(\xi_i)+o_K(1). \tag{13}
\]
将 (12)—(13) 代入 $4C_n=g^*H_n^{-1}g-D_\xi/f_n(r_n)$ 得 (3)。r721 已由直接代数及本机 Wolfram 核验的单因子恒等式给
\[
A(t)+B(t)\le1-\kappa\frac{t^2}{1+t^2},\qquad A,B\ge0.
\]
非负乘积展开表明 $\prod A+\sum_i B_i\prod_{j\ne i}A_j\le\prod_i(A_i+B_i)$。再用 $1+x\le e^x$ 与 $\xi_i^2\le K$，得 (4)。

## 3. 原子与弥散尾：精确极限推论

对任意有界能量序列，按 $|\xi_i|$ 递减重排，必要时取子列，使每个固定坐标 $\xi_{n,j}\to t_j$ 且 $E(\xi_n)\to E$。则 $\sum_jt_j^2\le E$，置 $\tau^2=E-\sum_jt_j^2\ge0$。由 $\log\psi(t)=-\mu t^2/2+O(t^3)$、$B(t)=O(t^4)$，无限乘积与下述和收敛，且 (3) 给
\[
4C_n\longrightarrow e^{-\mu\tau^2}
\left[\prod_{j\ge1}A(t_j)+\sum_{i\ge1}B(t_i)\prod_{j\ne i}A(t_j)\right]-1\le0. \tag{14}
\]
证明尾交换时先留前 $J$ 个原子。其后最大频率的极限至多 $|t_{J+1}|\to0$，对数乘积余项有界于 $C\max_{i>J}|\xi_i|\sum_{i>J}\xi_i^2$，而 $B$ 项尾有界于 $C\max_{i>J}\xi_i^2\sum_{i>J}\xi_i^2$。依次令 $n\to\infty$、$J\to\infty$ 即得 (14)。即使 $\sum_jt_j$ 不收敛，中心化乘积 $\prod\psi(t_j)$ 仍绝对对数收敛；弥散尾补偿切向均值所需的相位已在 (6) 精确处理。

## 4. 留下的缺口

本稿统一覆盖全部有界**中心化**能量，不只是有限支撑或单一尺度，也没有固定响应轴数。它没有处理中心化能量趋无穷的任意频谱；公共频移不保证把该能量变成有界。能量趋零时 (4) 的负主项同时消失，本稿的加法 $o_K(1)$ 不提供相对误差，故不能断言每个极小非零频率在每个大维数均为负。没有证明全有限维频率 PSD，也没有由此确定一般强极大端点的对数系数。

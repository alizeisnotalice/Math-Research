# audit_orbit042：038 全轨道机制与粗解析 I 改进独审

## 审计结论与证据边界

**通过。** 在冻结的 alpha=1、16 条等权 body 方向、beta=1/2 模型及 029 已证的负矩界下，038 的完整群公式 (4.3) 和粗解析完整积分界 (7.4) 均成立。没有发现需要撤销这两条结论的常数、量词或换序错误。这里通过的是带精确 Gamma 表达式 `J_flag` 的解析结论，未运行新 I 数值、旧 flag/Gamma、038 的 SPD 实例、QMC 或空间搜索。

必须保留三项限定：完整 signed-permutation 平均先于正 face 范数；单 flip 的全局线性损失含 `s=kappa/(1+kappa)`；已有负矩只保证 Green 的 Lq / a.e. 结论，允许 canonical Green 在零测奇异集合上取无穷。这些限定与 038、029 的实际陈述兼容。

本轮唯一预登记符号/有理批实际 PID **44735**，0.1395 秒，36 项检查全部通过，退出码 0。该批支持代数系数；不等式与可积性依据下面的独立推导。无新增小数 I 预算，也无 C5 突破结论。

只读来源：

- `Math-Research-cube5-astra-cover/cube-dim-5-bounds/astra_orbit038/report.md`，SHA-256 `fd60fc7b14afd6d907532c9660b4deeea5b93e59ff1d88069ad61445f11cddb8`，已实际核对。
- `Math-Research-cube-dim-5/cube-dim-5-bounds/upper_alpha1_029/report.md`，SHA-256 `298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81`，已实际核对；只引入同一 Green 模型、群不变性及有限 flag 负矩定理，不重算其数值。

## 1. 固定对象、交叉能量及 Hessian

令 q=5/4，M 为 5 阶正定矩阵，d=det M，A=M^{-1}，并令

\[
 H_M(x)=d^{-1/2}(x^TAx)^{-2},\qquad
 a_M=\int_{S^4}H_M^q=\omega_4d^{-1/8},\qquad \omega_4=8\pi^2/3.
\]

这里使用 `a_M` 避免将球面积分与精度矩阵 A 混同。Gaussian 极坐标恒等式
\(\int_{S^4}(\theta^TC\theta)^{-5/2}d\sigma=\omega_4/\sqrt{\det C}\)
对每个 C>0 成立。Gamma 正积分与 Tonelli 给，对任意正交 U，

\[
 \rho_U=\frac{\langle H_M^{q/2},H_{UMU^T}^{q/2}\rangle}{a_M}
 =\mathbb E_{t\sim B(q,q)}\bigl[d\det(tA+(1-t)UAU^T)\bigr]^{-1/2}. \tag{1}
\]

原始交叉积分的系数为 \(d^{-q/2}\omega_4\)，除以 \(a_M\) 正好留下 \(d^{-1/2}\)；因此括号中的 d 必须存在。U=I 时被积函数为 1。

若 X(t)=tA+(1-t)UAU^T、E=A-UAU^T，则

\[
 \frac{d^2}{dt^2}\log\det X(t)
 =-\operatorname{tr}\bigl[(X(t)^{-1/2}EX(t)^{-1/2})^2\bigr]\le0.
\]

当 E 非零时严格小于零。两个端点 determinant 都是 d^{-1}，所以 (1) 给 0<rho_U<=1；rho_U=1 当且仅当 UMU^T=M。这里没有从 Hessian 擅自推出一个与 M 无关的定量严格常数。

## 2. Pair-gap 的 q 常数

对 1<q<2，令正数 u=me^t、v=me^{-t}，交换后可设 t>=0。对 t>0，

\[
 \partial_q^2\log\sinh(qt)=-t^2\operatorname{csch}^2(qt)<0.
\]

在 q=1、2 两端插值得
\(\sinh(qt)\ge2^{q-1}\cosh^{q-1}(t)\sinh t\)。因此
\(F_q(t)=\cosh(qt)-1-2^{q-1}(\cosh^qt-1)\)
满足 F_q(0)=0、F_q'(t)>=0，等价于

\[
 (u+v)^q\le u^q+v^q+(2^q-2)(uv)^{q/2}.
\]

零值用连续性补足。两 profile 的 q 范数相同，故

\[
 \int_{S^4}\left(\frac{H_M+H_{UMU^T}}2\right)^q
 \le a_M[1-\beta_q(1-\rho_U)],\qquad
 \beta_q=1-2^{1-q}=1-2^{-1/4}. \tag{2}
\]

所以 q 次积分损失是 beta_q，取范数时再乘指数 1/q=4/5；范数平方的指数为 2/q=8/5。beta_q 是 pair 常数，与稳定指数 beta=1/2 不同。

## 3. 完整群 Hölder、球面到 face 与形状期望的顺序

对任意有限正交群 Gamma，P_Gamma 是球面上的平均算子，h=H_M^{q/2}。Hölder 的两个正权重
\(a=2-2/q\)、\(b=2/q-1\)
满足 a+b=1、aq/2+bq=1。因此无需假设 L^{q/2} 是范数，就有

\[
 (\mathbb EX)^q\le(\mathbb EX^{q/2})^{2q-2}(\mathbb EX^q)^{2-q}.
\]

逐点用于轨道，再用球面 Hölder 的共轭指数 1/(q-1)、1/(2-q)，得

\[
 \int(P_\Gamma H_M)^q
 \le a_M^{2-q}\left(\int(P_\Gamma h)^2\right)^{q-1}
 =a_M r_\Gamma^{q-1}, \tag{3}
\]
\[
 r_\Gamma=\frac{1}{|\Gamma|}\sum_{U\in\Gamma}\rho_U
 =\frac{\|P_\Gamma h\|_2^2}{a_M}.
\]

最后一个等式使用 P_Gamma 自伴、幂等及球面积分旋转不变性；并非把平均平方误当平方平均。若 Gamma_0 是子群，P_Gamma P_Gamma_0=P_Gamma，故 r_Gamma<=r_Gamma_0，子群估计可安全代入上界。

现在且仅现在令 Gamma 为完整 signed-permutation 群，阶为 2^5*5!=3840。对 f=(P_Gamma H_M)^q，f 是群不变的 -5 次齐次函数。cube 每个面上令 theta=z/|z|，径向 Jacobian 为 d sigma=|z|^{-5}dS，因为 z dot n=1；而 f(theta)=|z|^5 f(z)，恰好抵消。10 个面、每面 16 个正交象限构成 160 个 congruent F_+ 面片，公共边界零测。因此

\[
 \int_{F_+}(P_\Gamma H_M)^q
 =\frac1{160}\int_{S^4}(P_\Gamma H_M)^q,
\]
\[
 \|P_\Gamma H_M\|_{L^q(F_+)}
 \le(\pi^2/60)^{4/5}d^{-1/10}r_\Gamma^{1/5}. \tag{4}
\]

单个 M 或单 flip 平均一般没有该 face 对称性，不能直接除以 160。

实际 body 模型中，正 stable 权重几乎处处严格正，16 个方向张成 R^5，所以 M>0 几乎处处。任意 signed permutation 置换 body 直线；vv^T 消除代表方向的整体符号。独立同分布权重及 T^{-1/2} 倾斜均保持此置换不变，故 R 的 M 分布在全部共轭下不变。有限轨道平均与非负期望交换合法，即使某点 Green 是无穷，也有

\[
 G/c_1=\mathbb E_R P_\Gamma H_M,\qquad c_1=(32\pi^3)^{-1}.
\]

029 已证 \(\mathbb E_Rd^{-1/10}\le J_{\rm flag}<\infty\)。因为 r_Gamma<=1，积分 Minkowski 对 (4) 的右端可积，给

\[
 I_1\le c_1^{5/4}\frac{\pi^2}{60}
 \left[\mathbb E_R(d^{-1/10}r_\Gamma^{1/5})\right]^{5/4}. \tag{5}
\]

可先对正形状积分作有限截断，再单调极限。此过程覆盖所有近退化 M，且未删除空间管、时间段或尺度；它是 Lq / a.e. 结论，并不声称 canonical Green 在 S1/S2 等零测奇异集逐点有限。不可将 determinant 与 r_Gamma 分离为独立随机变量。

## 4. 单 flip：协方差块与精度块的两个二次型

以下使用实际模型 M_ii=1。写

\[
 M=\begin{pmatrix}1&b^T\\ b&C\end{pmatrix},\quad
 s=b^TC^{-1}b\in[0,1),\quad
 A=M^{-1}=\begin{pmatrix}a&r^T\\r&D\end{pmatrix}.
\]

准确区分为

\[
 a=(1-s)^{-1},\quad \kappa=a-1=s/(1-s),\quad
 a-r^TD^{-1}r=1/M_{ii}=1,\quad r^TD^{-1}r=\kappa.
\]

因此协方差二次型是 s，精度二次型是 kappa；二者一般不相等。翻转第一坐标把 r 乘以 z=2t-1，Schur 公式给

\[
 d\det[tA+(1-t)SAS]
 =\frac{a-z^2r^TD^{-1}r}{a-r^TD^{-1}r}
 =1+(1-z^2)\kappa. \tag{6}
\]

若 t~Beta(q,q)，换元 Y=4t(1-t) 的两个分支相加，其密度正比于 Y^{q-1}(1-Y)^{-1/2}；所以 Y~Beta(q,1/2)、EY=5/7。函数
\(g(y)=1-(1+\kappa y)^{-1/2}\)
满足 g(0)=0 与 \(g''(y)=-3\kappa^2/[4(1+\kappa y)^{5/2}]\le0\)，故在 [0,1] 上 g(y)>=y g(1)。于是

\[
 1-\rho_i\ge\frac57\left(1-\frac1{\sqrt{1+\kappa}}\right)
 =\frac57(1-\sqrt{1-s})\ge\frac5{14}s
 =\frac5{14}\frac{\kappa}{1+\kappa}. \tag{7}
\]

不能把最后一个 s 换成 kappa；后一式在 kappa 无界时立即不可能。又因 tr C=4、C>0，C^{-1}>=I/4，因此 s>=sum_{j!=i}M_ij^2/4。

最后，Z=1-Y~Beta(1/2,q) 给
\(\rho_i=\sqrt{1-s}\,\mathbb E(1-sZ)^{-1/2}\)。二项正级数与 Beta 矩逐项积分准确给普通 Gauss 系数
\((1/2)_n^2/((7/4)_n n!)\)，没有额外 Gamma 归一化。038 式 (5.4) 的首项 `sqrt` 缺少反斜杠只是排版笔误。

## 5. 一般 sign split 与 rank 0/1/2 常数

精度矩阵按 I|J 分块为 E,F,L。标准化交叉块的平方
\(T=E^{-1/2}FL^{-1}F^TE^{-1/2}\)
具有谱 tau_l in [0,1)。非零谱数 k<=min(|I|,|J|)<=2，令 kappa_l=tau_l/(1-tau_l)。Schur 分解直接给

\[
 d\det[tA+(1-t)SAS]
 =\frac{\det(I-z^2T)}{\det(I-T)}
 =\prod_{\ell=1}^k(1+Y\kappa_\ell). \tag{8}
\]

块逆矩阵和 Jacobi 主子式恒等式分别给

\[
 K=\sum\kappa_\ell=\operatorname{tr}(M_{II}A_{II})-|I|,\qquad
 D_{I|J}=\prod(1+\kappa_\ell)=\frac{\det M_{II}\det M_{JJ}}d.
\]

乘积 >=1+YK，所以 (7) 同样适用 K：
\(\rho_S\le1-(5/14)K/(1+K)\)。每个因子还满足 1+Y kappa>=Y(1+kappa)，因此

\[
 \rho_S\le\min(1,C_kD_{I|J}^{-1/2}),\qquad
 C_k=\frac{B(q-k/2,1/2)}{B(q,1/2)}
 =4^{-k/2}\frac{B(q-k/2,q-k/2)}{B(q,q)}. \tag{9}
\]

积分条件是 k<2q=5/2，故本问题的 k=0,1,2 全部有限。具体 C_0=1，
\(C_1=B(3/4,3/4)/(2B(5/4,5/4))\)，
\(C_2=(q-1/2)/(q-1)=3\)。若实际 rank 低于 min block size，用较大 rank 常数仍安全，因为 0<Y<=1；但 rank 2 不能用 C_1。038 的两条 Beta 常数均正确。

## 6. 粗解析修正的完整闭合

令 S_i 是五个单坐标 flip，Q_i=(Id+S_i)/2 为对应的函数平均算子。由于 P_Gamma Q_i=P_Gamma，Jensen 与球面测度不变性给
\(\|P_\Gamma H\|_q\le\|Q_iH\|_q\)。这一步对每个 i 分别成立，再平均五个合法上界。写 \(w=\sum_{i<j}M_{ij}^2\)，由 (2)、(7) 得

\[
 \frac15\sum_i(1-\rho_i)
 \ge\frac15\frac5{14}\frac14\sum_i\sum_{j\ne i}M_{ij}^2
 =\frac{w}{28},
\]
\[
 \int_{S^4}(P_\Gamma H)^q\le a_M(1-\beta_qw/28).
\]

随后才能用完整群对称性除以 160。取 q 次根并用凹切线 (1-x)^{4/5}<=1-4x/5，得到

\[
 \|P_\Gamma H_M\|_{L^q(F_+)}
 \le(\pi^2/60)^{4/5}d^{-1/10}(1-\beta_qw/35). \tag{10}
\]

正定相关矩阵有 |M_ij|<1，故 0<=w<=10，括号为正。

独立核对 R 的二阶矩：一般稳定指数 0<beta<1、N 个等权变量，LT 为 exp(-t^beta)，倾斜 Z=ET^{-beta}=1/(N Gamma(1+beta))。不同变量的 Laplace 导数乘积给

\[
 E_R P_uP_v=\frac{1}{Z\Gamma(\beta+2)}
 \int_0^\infty\beta^2t^{3\beta-1}e^{-Nt^\beta}dt
 =\frac{2\beta}{N^2(1+\beta)}.
\]

对角再加二阶导数的 beta(1-beta)t^{beta-2} 项，得到
\(E_RP_u^2=(1-\beta)/(N(1+\beta))+2\beta/(N^2(1+\beta))\)。积分全为非负且有限，不使用未倾斜的正矩。beta=1/2、N=16 时交叉为 1/384、对角为 3/128；准确满足 \(N E P_u^2+N(N-1)E P_uP_v=1\)。

对每个坐标对，16 个 body 系数 epsilon_u=v_ui v_uj 平衡为 8 正 8 负，所以

\[
 E_RM_{ij}^2=N(E_RP_u^2-E_RP_uP_v)=1/3,\qquad E_Rw=10/3. \tag{11}
\]

Hadamard 给 0<d<=1。w 有界与 029 的有限负矩还保证 \(E_R(d^{-1/10}w)<\infty\)，所以以下减法不涉及 infinity-infinity。由 (10)、Minkowski、(11)，

\[
 \|G\|_q\le c_1(\pi^2/60)^{4/5}
 \left[E_Rd^{-1/10}-\frac{\beta_q}{35}E_R(d^{-1/10}w)\right]
 \le c_1(\pi^2/60)^{4/5}\left[J_{\rm flag}-\frac{2\beta_q}{21}\right].
\]

这里精确沿用已锁定的

\[
 J_{\rm flag}=4\Gamma(3/2)
 [\Gamma(6/5)/\Gamma(11/10)]^5(2661/64)^{-1/10}.
\]

所以审计通过的完整积分上界是

\[
 \boxed{I_1\le c_1^{5/4}\frac{\pi^2}{60}
 \left[J_{\rm flag}-\frac{2(1-2^{-1/4})}{21}\right]^{5/4}.} \tag{12}
\]

由于 J_flag>=E_Rd^{-1/10}>=1、0<beta_q<1/4，右端括号为正。相对旧同一公式的比值严格小于 1、严格大于 (41/42)^{5/4}>(41/42)^2>19/20。因此这确实改善完整 I 预算，但该固定粗公式不能提供 5% 的下降。它不是实际 I 的下界，也不能排除更强的完整轨道期望 (5)。

038 开头“更强”的表述不应读成：对任意有限群，Hölder 表达式必逐点优于所有 pair 表达式。两者应分别保持为合法上界，必要时取较小者；038 的唯一数值例已明确区分二元 Hölder 与 pair-gap，故这不是其 (12) 的证明缺口。

## 7. 复现、归档与锁定

唯一批只运行 `verify_coefficients.py`：Python 3.15.0、SymPy 1.14.0；未调用浮点特殊函数或 Arb，因为本审计没有新的数值端点。预登记 120 秒截止，实际所有阶段完成，最终 checkpoint 无活动句柄。

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/audit_orbit042/verify_coefficients.py \
--output-dir /absolute/path/to/new-audit-output
```

锁后复核请使用新输出目录，避免覆盖本轮实际收据。最小白名单为本报告、登记、验证器、结果、运行收据、checkpoint；其 hash 见 manifest.json。仅创建本独立 clone 的 audit_orbit042；未修改 038、主仓库、旧锁，未提交或推送。

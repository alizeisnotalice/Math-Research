# 一般强矩形端点的单中位数归约

设 $H_t(p)=\sup_E |\{M_{S,t}1_E>p\}|/|E|$，上确界取遍全部有限正测度的可测集。写 $\sigma^2=\mu=1+e^{-\mu}$、$g=(\mu-1)/8$、$a_3=(3-\mu)/6$、$c_Q=a_3+g$。令
\[
 \alpha=\limsup_{t\to\infty}\frac{\log H_t(1/2)-\mu t}{\log t}.
\]
则对项目定义的一般 Orlicz 端点，
\[
 \limsup_{N\to\infty}\frac{\log C_N-\mu N+d_N}{\log N}\le c_Q+\alpha . \tag{1}
\]
这里 $d_N=\sigma\sqrt N z_N$、$\bar\Phi(z_N)=\sigma z_N/\sqrt N$，$\Phi,\phi$ 为标准正态分布函数与密度。全坐标 $c=\mu$ 的深层实际指标见证配合第22轮中央 Edgeworth 估计，给固定足够大 $A$ 下 $H_t(1/2)\ge e^{\mu t-A}$；第91轮式 (6) 给粗一般上侧。因此 $0\le\alpha<\infty$。特别，若能证明 $\log H_t(1/2)\le\mu t+o(\log t)$，则 (1) 与一般下界一起确定系数为 $c_Q$。

## 每一步都保留真实中位数边界

定义 $h_t(x)=\inf\{p\in[0,1]:H_t(p)\le x\}$，约定 $H_t(1)=0,H_t(0)=+\infty$，并置 $D_t(b)=1-h_t(e^{\mu t+b})$。这是所有真实指标共用的非减亏损函数。取重排 $f^*(v)=\inf\{p:|\{f>p\}|\le v\}$。则
\[
 h_t(x)=\sup_E (M_{S,t}1_E)^*(x|E|).
\]
事实上当 $H_t(p)\le x$ 时右侧至多 $p$；当 $p<h_t(x)$ 时，某个 $E$ 的严格水平集超过 $x|E|$，故右侧至少 $p$。反向使用时，$h_t(x)<p$ 充分推出 $H_t(p)\le x$。

在 $t$ 维写 $E_z$ 为最后 $t-1$ 个坐标的截面，$G(z,y)=M_{S,t-1}1_{E_z}(y)$。Fubini 给 $|\{G>p\}|\le H_{t-1}(p)|E|$，故 $G^*(v)\le h_{t-1}(v/|E|)$。真实矩形给 $M_{S,t}1_E\le M_1G$，第11轮全局单轴重排不等式给 $(M_1G)^*\le T G^*$。其中
\[
 (Tk)(x)=\sup_c \mathbb E k(xe^{-Z_c}),\quad Z_c=W+cJ_c,\quad W\sim\mathrm{Exp}(1),\quad P(J_c=1)=1/(1+e^c),
\]
独立，$0\le c\le\infty$，无穷动作按极限取。可先取有界支撑源，再单调截断；对于此源 $G$ 有限正水平集、沿未用轴有界支撑，符合第11轮的纤维证明。即使原重排比较先以几乎处处陈述，由右连续的递减重排及 $Tk$ 在 $x>0$ 的连续性（其积分公式中 $r=e^{-c}\in[0,1]$ 上的紧最大值），也可逐秩使用。除以 $|E|$ 后取上确界得
\[
 D_t(b)\ \ge\ \inf_c\mathbb E D_{t-1}(b+Y_c),\qquad Y_c=\mu-Z_c\le\mu. \tag{2}
\]
每个 $D_t$ 同时满足它本维的几何中位数边界。给定 $\epsilon>0$，令 $L=\log N$, $n=\lfloor N/L^{40}\rfloor$, $s=(\alpha+\epsilon)L$。因 $\alpha\ge0$，充分大 $N$ 时对所有 $n\le t\le N$ 都有 $H_t(1/2)\le e^{\mu t+s}$，从而
\[
 D_t(b)\ge 1/2\quad(b\ge s). \tag{3}
\]
(2)、(3) 属于同一族物理包络，故可在任一步重新使用 (3)。

## 左半轴上的离散比较

写 $x=b-s$, $u=x/(\sigma\sqrt t)$, $f_t(x)=\phi(u)/(\sigma\sqrt t)$。取 $\lambda=1/32$，在 $u\le0$ 上置 $w(u)=2-u-e^u$；在 $[0,1/10]$ 保持此公式，随后作固定光滑截断延拓，使 $w$ 在正半轴非负、有界、各阶导数有界（例如与常数 1 光滑拼接，在 $1/5$ 前完成）。在负半轴 $w\ge1$ 且
\[
 \mathcal A w+\lambda w\le0,\qquad \mathcal A w=(w''-uw')/2.
\]
确实 $-2\mathcal A w=-u+(1-u)e^u$。在 $[-1,0]$ 上后一项至少 $e^{-1}$，且 $w\le3$；在 $u\le-1$ 时 $-2\mathcal A w\ge |u|$，$w\le3|u|$。取固定 $K>g$，定义
\[
 V_t(x)=\Phi(u)-[c_Q u^2+C_2+K\log(t/n)]f_t(x)
       -M L f_t(x)(n/t)^\lambda w(u). \tag{4}
\]
$C_2$ 固定大于 $2\mu+2$，$M$ 为充分大的固定常数。初始时第91轮 (6) 或第107轮式 (3) 给
\[
 D_n(s+x)\ge\Phi((s+x-C_H\log n)/(\sigma\sqrt{n+t_0}))-2n^{-4},\quad t_0=O(\log n).
\]
在 $-24\sqrt L\le u\le0$，右侧至少为 $\Phi(u)-C L f_n(x)-2n^{-4}$（均值定理，两个分位间距为 $O(L/\sqrt n+|u|L/n)$）。取 $M>C$。在更左处 $V_n\le\Phi(u)\le N^{-20}$。这给整个左半轴的初值比较。

以下验证 $-24\sqrt L\le u\le0$ 的一步比较。$P_1$ 表示中心方差 $\mu$ 的一步正态卷积。去掉 (4) 的末项时，依高斯乘积恒等式，其 $P_1V_{t-1}-V_t$ 除以 $f_t$ 恰为 $c_Q(u^2-1)/t+K\log(t/(t-1))$。对新增的 $F_t(x)=M L f_t(x)(n/t)^\lambda w(u)$，直接微分给
\[
 (\partial_t-\tfrac{\mu}2\partial_x^2)F_t
 =-ML f_t(n/t)^\lambda (\lambda w+\mathcal A w)/t\ge0 \quad(x\le0).
\]
在本管内，$P_1F_{t-1}-F_t=-(\partial_t-\tfrac{\mu}2\partial_x^2)F_t+O(f_t(1+L)^C/t^2)$。这是在当前点的时间一步 Taylor 和正态空间 Taylor 至三阶：两次时间导数、混合导数 $\partial_t\partial_x^2$、四次空间导数除以 $f_t$ 均为 $O((1+L)^C/t^2)$；正态余项位移 $v$ 可额外乘 $(1+|v|)^C e^{|v|/8}$ 后积分。由 $w$ 的延拓，该误差越过 $u=0$ 时仍有效。故 (4) 中负号的热校正至多损失该余项。

对非高斯动作，现列出第107轮式 (7) 的参数展开及本次使用条件。递减尾的写法换为此处的非减亏损，符号同该轮。逐次微分 (4) 给，以下导数取在 $(t-1,x)$，$P_L=(1+L)^C$ 可随式增大固定次数，
\[
 V'/f_t=1+O(P_L/\sqrt t),\quad V''/V'=-u/(\sigma\sqrt t)+O(P_L/t),\quad V'''/V'=(u^2-1)/(\mu t)+O(P_L/t^{3/2}).
\]
四阶 Taylor 的期望余项为 $O(f_t P_L/t^{3/2})$：$w$ 负半轴至多线性增长，导数有界，位移 $v$ 的高斯密度比至多 $C\exp(C\sqrt{L/t}|v|)$，$Y_c$ 一致有正负指数矩（例如阶数 $1/4$）。注意空间三阶展开的四阶系数来自整个 $V$，其主项为 $\Phi$。记 $\Delta=\mu-\mathbb E Z_c$, $v_c=\operatorname{Var}Z_c$, $\kappa_c=\kappa_3(Z_c)$。与正态比较，除以 $V'$ 的动作多项式为
\[
 \Delta+\tfrac12(v_c-\mu+\Delta^2)V''/V'
 +\tfrac16(-\kappa_c+3\Delta v_c+\Delta^3)V'''/V'.
\]
第22轮的矩估计 $|v_c-\mu|\le C\sqrt\Delta$ 和一致有界的三阶矩，与 $c=\mu$ 比较使近极小动作满足 $\Delta=O(P_L/t)$；此时 $\kappa_c=\mu(3-\mu)+O(\sqrt\Delta)$。为核对二次系数，$m(c)=1+c/(1+e^c)$ 在 $\mu$ 唯一非退化极大，$m''(\mu)=-(\mu-1)/\mu$, $v'(\mu)=\mu-1$。对小的两侧 $h$，$\inf_c[\Delta+(v_c-\mu)h/2]=-\mu g h^2+O(|h|^3)$，由非退化极小分支的 Taylor 公式得到。代入 $h=-u/(\sigma\sqrt t)$（该替换乘上 $O(\sqrt\Delta)$），得
\[
 \frac{\inf_c \mathbb E V_{t-1}(x+Y_c)-P_1V_{t-1}(x)}{f_t(x)}
 =-g u^2/t-a_3(u^2-1)/t+O(P_L/t^{3/2}).
\]
与高斯步相加，二次项相消，得到下界 $K\log(t/(t-1))-g/t-O(P_L/t^{3/2})>0$。全程 $t\ge n$ 大于 $L$ 的任何固定幂。

现在只沿 $x\le0$ 归纳实际 $D_t(s+x)$。在 $0\le x'\le\mu$，$\Phi(x'/(\sigma\sqrt{t-1}))-1/2\le (\mu+o(1)) f_{t-1}(x')$，故 $V_{t-1}(x')\le1/2$。由 $Y_c\le\mu$，从左半轴越过边界的所有前态都在这段，(3) 逐点给出支配，完全不用给越界概率收费。在 $u<-24\sqrt L$，直接用 $D_t\ge0$、$V_t\le\Phi(u)\le N^{-20}$。因此从初值出发，(2) 和管内非负残差给全时刻统一误差
\[
 D_t(s+x)\ge V_t(x)-2n^{-4}-N^{-20},\qquad n\le t\le N,\ x\le0. \tag{5}
\]
此处外带和中位数边界在每个时刻使用真实的 $D_t$，无需累计外带误差。

## 临界秩与原端点

取 $x_N=-d_N+(c_Q+\epsilon)L$，物理阈值 $p_N=1-d_N/N-BL/N$，$B$ 任意固定。此时 $u_N=-z_N+O(L/\sqrt N)$，$u_N^2=L+o(L)$，$f_N(x_N)\sim L/N$，且 $w(u_N)=O(\sqrt L)$。新增校正除以 $L f_N$ 至多 $O(L^{1/2-40/32})=o(1)$。由 $\log(N/n)=O(\log L)$、正态分位的局部展开和 (5)，
\[
 D_N(s+x_N)\ge d_N/N+(\epsilon-o(1))L f_N(x_N)>1-p_N.
\]
故 $\log H_N(p_N)\le\mu N-d_N+(c_Q+\alpha+2\epsilon)L$。第92轮的一般物理单指标归约 $C_N\le(1+o(1))[H_N(p_N)+Q_N N^{-P}]$ 允许对任一固定 $P>0$ 取相应固定 $B$。第22轮独立下界的展开从而使第二项不改变所需上界。令 $\epsilon\downarrow0$ 得 (1)。

本证明消去第110轮所需的一整段中央分位数据；任意真实源在中位数处的 $o(\log t)$ 上误差仍待证。现有反馈上界仍可与 (1) 取小。

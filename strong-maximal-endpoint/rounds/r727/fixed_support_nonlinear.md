# 固定有限角向支撑的非线性指标见证上界

本稿将 r725 的精确几何相位不等式连接到真实指标源的中位数，使用独立非活跃坐标的 Berry–Esseen 估计。允许全部 $n$ 个侧比依接收相位任意变化，允许相位轮廓具有 $O(\log n)$ 振幅；不要求轮廓可导或振幅足够小。结论限定于固定有限角向支撑的指定实际见证族，不是一般强极大 halo 的上界。

## 1. 陈述

固定整数 $d\ge1$，$n>d$、$N=n-d$。非零实频率 $\xi\in\mathbb R^n$ 只在前 $d$ 个坐标上可能非零，频率的数值及能量不受限制。令 $h$ 是实值连续周期函数，定义
\[
A=\inf_{a\in\mathbb R}\|h-a\|_\infty=\tfrac12(\max h-\min h).
\]
侧比场 $c(\theta)\in\mathbb R^n$ 可为任意有限实值 Borel 可测周期场。给定 $\theta$，所有坐标独立，
\[
Z_{c_i}=W_i+c_iJ_i,\quad W_i\sim\operatorname{Exp}(1),\quad
J_i\sim\operatorname{Bern}((1+e^{c_i})^{-1}),\quad M(c)=1+\frac{c}{1+e^c}.
\]
记 $r_\theta$ 为
\[
T_\theta=\sum_{i=1}^nZ_{c_i(\theta)}-h\left(\theta+\sum_{i\le d}\xi_iZ_{c_i(\theta)}\right)
\]
的唯一中位数，$K_ch(\theta)=E h(\theta+\xi\cdot Z_c)$。存在仅依赖 $d$ 的常数 $C_d$：对每个固定 $B<\infty$，若
\[
A\le BN^{1/3},\qquad N\ge C_d(1+B^3)^2, \tag{1}
\]
则统一于 $\theta$、全部动作值、全部上述频率及轮廓，有
\[
\left|r_\theta-\left(\sum_iM(c_i(\theta))-K_{c(\theta)}h(\theta)\right)\right|
\le C_d\left(1+\frac{A^3}{N}\right). \tag{2}
\]
从而精确相位泛函满足
\[
\boxed{\displaystyle
\frac{\langle e^{r_\theta}\rangle}{\langle e^{-h}\rangle}
\le\exp\left\{\mu n+C_d\left(1+\frac{A^3}{N}\right)\right\}.} \tag{3}
\]
这里 $\mu=1+e^{-\mu}$，尖括号为圆周 Haar 平均。特别地，$d$ 固定、$A=O(\log n)$ 时，全 $n$ 轴相位适应后的评分对数不超过 $\mu n+O_d(1)$。本稿不声称 $d$ 随维数增长时该常数仍然可用。

## 2. 所有实动作的统一中心矩

令 $q=(1+e^c)^{-1}$。对所有 $c\in\mathbb R$，
\[
1\le\operatorname{Var}Z_c=1+c^2q(1-q)\le1+4e^{-2}<2,
\]
因为 $q(1-q)\le e^{-|c|}$。此外
\[
E|c(J-q)|^3=|c|^3q(1-q)[q^2+(1-q)^2]\le27e^{-3},
\]
且 $E|W-1|^3=12/e-2$。由 $|a+b|^3\le4(|a|^3+|b|^3)$，得
\[
\sup_{c\in\mathbb R}E|Z_c-M(c)|^3\le4(12/e-2+27e^{-3})<20. \tag{4}
\]
因此即使某些动作很负、均值向负无穷移动，所需中心矩常数也不变化。

加常数到 $h$ 会使所有 $r_\theta$ 减去同一常数，并使 (3) 的分子、分母乘以同一因子。故可且以下假设 $\|h\|_\infty\le A$。对固定 $\theta$，写
\[
T_\theta=\sum_{i>d}Z_{c_i}+Y_\theta,\qquad
Y_\theta=\sum_{i\le d}Z_{c_i}-h\left(\theta+\sum_{i\le d}\xi_iZ_{c_i}\right).
\]
关键是 $Y_\theta$ 与全部 $N$ 个非活跃增量独立；相位对活跃坐标的任意非线性依赖完全保留在这个单独加项中。由 (4)，
\[
\operatorname{Var}Y_\theta\le4d+2A^2,\qquad
E|Y_\theta-EY_\theta|^3\le80d^3+32A^3. \tag{5}
\]
第一式使用 $\operatorname{Var}h\le A^2$；第二式可用 $|\sum_{i\le d}X_i|^3\le d^2\sum_i|X_i|^3$，无需对 $h$ 求导。

## 3. 将整个活跃块作为一个加项的中位数估计

令 $\sigma_\theta^2=\operatorname{Var}T_\theta$。独立性和 (4)—(5) 给
\[
N\le\sigma_\theta^2\le2N+4d+2A^2.
\]
对这 $N+1$ 个独立加项应用非同分布 Berry–Esseen 不等式，
\[
\sup_x\left|P\left(\frac{T_\theta-ET_\theta}{\sigma_\theta}\le x\right)-\Phi(x)\right|
\le \varepsilon_N:=\frac{C_d(N+A^3)}{N^{3/2}}. \tag{6}
\]
在 (1) 中增大同一个 $C_d$ 后，$\sigma_\theta\le C_d\sqrt N$，且 $\varepsilon_N$ 小于固定绝对常数。因为 $\Phi(x)-1/2\ge\phi(1)x$ 对 $0\le x\le1$ 成立，可在均值两侧各取 $C\sigma_\theta\varepsilon_N$，使分布函数分别小于、大于 $1/2$。故
\[
|r_\theta-ET_\theta|\le C\sigma_\theta\varepsilon_N
\le C_d(1+A^3/N),
\]
即 (2)。阈值 (1) 足够：$A^2/N\le B^2N^{-1/3}$ 及 $\varepsilon_N\le C_d(1+B^3)/\sqrt N$ 同时受控。

这里中位数确实唯一。任一个有限实动作的密度在 $(\min(0,c),\infty)$ 上严格为正；因此非活跃和的密度在其最小支撑点之上严格为正。与独立的 $Y_\theta$ 卷积后，$T_\theta$ 的分布连续，并在其最小支撑点之上严格递增。$Y_\theta$ 有有限下界，故中位数存在且唯一。对 $N=1$ 非活跃密度的跳跃不影响其分布函数连续；本定理实际使用充分大 $N$。分布函数对 $(\theta,x)$ 可测，因此唯一中位数 $r_\theta$ 也是 Borel 可测。

## 4. 几何相位不等式的直接结算

r725 的精确不等式对 $F=e^{-h}$ 给
\[
\left\langle\exp\left\{\sum_iM(c_i(\theta))-K_{c(\theta)}h(\theta)\right\}\right\rangle
\le e^{\mu n}\langle e^{-h}\rangle. \tag{7}
\]
对于无界但有限实值的 Borel 场，可先双边截断 $c_i^{(L)}=\max(-L,\min(c_i,L))$。每个相位处最终所有坐标均恢复原值；被积函数不超过 $\exp(\mu n+\|h\|_\infty)$，因为 $M(c)\le\mu$。支配收敛把有界场的 (7) 扩展到全部这种场。因此将 (2) 的上侧指数化、再使用 (7)，即得 (3)。同样的截断与 (2) 的统一上界说明评分泛函的无界场值可由有界场值逼近；这里不需要在每个维数预设统一动作紧集。

## 5. 连续一般轮廓确实来自真实指标源

固定上述 $n,h,\xi$ 和任意**有界 Borel** 动作场，取 $C_0$ 控制全部 $|c_i|$。定义实际指标源
\[
E_{T,h}=\{x\in[-1,1]^n:\ u_i=-\log|x_i|\ge0,\ \sum_i u_i>T+h(\xi\cdot u)\}.
\]
零坐标不影响体积。对接收点 $u_i\ge C_0+1$，用端点 $x_i,-e^{-c_i(\theta)}x_i$、$\theta=\xi\cdot u$ 的实际矩形；它位于立方体内。盒内源坐标精确为 $u+Z_c$，所以盒内指标平均超过 $1/2$ 当且仅当
\[
P(T_\theta>T-\sum_i u_i)>1/2
\quad\Longleftrightarrow\quad T-\sum_i u_i<r_\theta. \tag{8}
\]
这用的是连续严格递增的分布，而非对轮廓求导或小振幅展开。

总深度 $l=\sum_i u_i$ 的体积密度是 $2^ne^{-l}l^{n-1}/(n-1)!$，给定 $l$ 时 $u=lD$，其中 $D$ 在标准单纯形上均匀。由于 $\xi$ 非零且至少一个非活跃坐标为零，$\xi$ 不是常向量，所以 $\xi\cdot D$ 有 $L^1$ 密度。其 $l$ 倍模 $2\pi$ 的分布以**总变差**趋于 Haar：先将该 $L^1$ 密度逼近为有限个区间示性函数的线性组合，缩放后每个区间的绕圈密度误差为 $O(1/l)$；缩放、模映射及卷积均收缩总变差，再令逼近误差趋零。这可用于任意有界可测接收阈值，而不要求动作场连续。

令 $v=T-l$，将体积除以 $a_T=2^ne^{-T}T^{n-1}/(n-1)!$ 后，积分权为 $e^v(1-v/T)^{n-1}1_{v\le T}$。源阈值为 $v<-h(\theta)$，接收阈值为 $v<r_\theta$。固定场下这两个阈值均有界；上述总变差收敛可用于每个固定 $v$。低侧以 $e^v(1+|v|)^{n-1}$ 支配，高侧超过统一阈值即为零。接收点的额外限制 $u_i\ge C_0+1$ 只移除固定深度单纯形中 $O_n(1/l)$ 的相对体积。于是源与该所选接收集 $\Omega_{T,c}$ 满足
\[
\frac{|E_{T,h}|}{a_T}\longrightarrow\langle e^{-h}\rangle,\qquad
\frac{|\Omega_{T,c}|}{a_T}\longrightarrow\langle e^{r_\theta}\rangle. \tag{9}
\]
因此 (3) 是真实指标与实际矩形见证族的非局部深层评分上界。对每个维数先取 $T\to\infty$，再讨论 $n\to\infty$；无需 $h$ 的导数、非退化水平集或小参数假设。它限制该见证族能达到的评分，没有把整个指标 halo 替换成该接收集。

## 6. 固定支撑结论的范围

固定 $d$、$A=O(\log n)$ 时，全部动作的相位适应以及任意固定支撑频率大小，都不能通过这类深层指标见证产生正的额外 $\log n$ 项。这里保留全部活跃块内部相关性，独立性只用于非活跃坐标与整个活跃块之间。第7节将量化其中一部分增长支撑范围；相位依赖全部坐标时，本证明不给所需控制。一般强极大端点仍未闭合。

## 7. 追加推论：量化的增长支撑范围

前面的独立块论证还给出以下更大范围；这一段使用通用常数，替代固定 $d$ 的粗三阶矩估计。有绝对常数 $\eta>0,C,n_0$：若
\[
n\ge n_0,\qquad 1\le d\le n/2,\qquad A\le\sqrt n,\qquad
d^{3/2}+A^3\le\eta n^{3/2}, \tag{10}
\]
则所有前述动作、频率及轮廓统一满足
\[
|r_\theta-ET_\theta|\le C\left(1+\frac{d^{3/2}+A^3}{n}\right),\qquad
\log\frac{\langle e^{r_\theta}\rangle}{\langle e^{-h}\rangle}
\le\mu n+C\left(1+\frac{d^{3/2}+A^3}{n}\right). \tag{11}
\]
证明只需改善 (5)。指数变量有 $E(W-1)^4=9$，而
$E|c(J-q)|^4\le |c|^4e^{-|c|}\le256e^{-4}$，故
\[
\sup_c E|Z_c-M(c)|^4\le8(9+256e^{-4})=:C_4.
\]
指数三、四阶矩另经本机 Wolfram 实际核验，见 `centered_moments_check.wl` 和原输出。对独立中心化的活跃增量 $X_i$，精确展开给
\[
E\left(\sum_{i\le d}X_i\right)^4
=\sum_{i\le d}EX_i^4+6\sum_{i<j}EX_i^2EX_j^2
\le C_4d+12d(d-1)\le C d^2.
\]
由 Lyapunov 不等式，$E|\sum X_i|^3\le C d^{3/2}$。因此 (5) 的三阶矩改为 $C(d^{3/2}+A^3)$，总方差仍满足 $N\le\sigma_\theta^2\le2N+4d+2A^2$。在 (10) 下，$N\ge n/2$、$\sigma_\theta\le C\sqrt n$，Berry–Esseen 误差至多
\[
C\frac{n+d^{3/2}+A^3}{n^{3/2}}.
\]
先将 $\eta$ 取得足够小、再将 $n_0$ 取得足够大，该误差即可落入第 3 节固定的中位数反演区间。相同反演与几何相位结算给 (11)。因此条件中明确保留了小的绝对比例 $\eta$，未仅由 $d\le n/2$ 就推断 BE 误差足够小。

特别地，对序列
\[
d=o((n\log n)^{2/3}),\qquad A=o((n\log n)^{1/3}), \tag{12}
\]
(10) 自动最终成立，且 (11) 的附加费用为 $o(\log n)$。这包含 $A=O(\log n)$、$d=o((n\log n)^{2/3})$ 的情形。每个维数的真实源与接收体积仍按第 5 节先取深度极限；没有使用对维数一致的相位混合速率。该推论扩大了可控支撑范围，仍不覆盖依赖全部坐标的一般相位轮廓或一般强极大端点。

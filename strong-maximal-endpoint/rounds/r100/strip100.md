# 临界指标估计在增加自由坐标后的传递

令 $n=N-k$, $L=\log N$, $k\log L=o(L^2)$，$k\ge1$。$M_j$ 为非中心轴平行强极大算子。采用 r096 的 $p_N=1-\tau_N-BL/N$, $q_{0,N}=\mu N-d_N$, $d_N=N\tau_N$。令 $\mathcal L_n$ 为 r094 的真实无条件对数凸下集族，加上任意平移、非零对角伸缩与空集。若 $E\subset\mathbb R^k\times\mathbb R^n$ 可测且具有正有限体积，几乎每个截面 $E_z=\{w:(z,w)\in E\}$ 都属于 $\mathcal L_n$（允许零测修改），则
\[
 \log\sup_E\frac{|\{M_N1_E>p_N\}|}{|E|}=q_{0,N}+c_Q L+o_B(L).
\]
截面的形状、体积、中心及对角伸缩都可以依赖 $z$。充分性也适用于在下述低门槛有同一上界的其他截面族。

设 $\delta=L/N$, $t=p_N-\delta$, $v=t+A(k+L)L/N$，其中 $A$ 为稍后选择的常数。最终 $0<t<p_N<v<1$。写 $g(z,w)=M_n1_{E_z}(w)$。可取有理矩形使其联合可测。每个真实矩形给 $M_N1_E\le M_k(g(\cdot,w))(z)$。由
\[
 g\le t+(v-t)1_{g>t}+(1-v)1_{g>v},
\]
严格成功必使同一个前块平均的两项之一超过 $\delta/2$。先按 $z$ 使用指标弱常数，再对 $w$ 积分并用 Tonelli，得到
\[
 \frac{|\{M_N1_E>p_N\}|}{|E|}
 \le H_k\!\left(\frac\delta{2(v-t)}\right)H_n^{\mathcal L}(t)
   + H_k\!\left(\frac\delta{2(1-v)}\right)H_n(v).                 \tag{1}
\]
每层使用原来截面的实际体积，积分恰为 $|E|$。只有第一项需要截面假设。

由 $\bar\Phi(z_N)=\sigma z_N/\sqrt N$ 及 Mills 界，$z_N\asymp\sqrt L$、$\phi(z_N)\asymp L/\sqrt N$。隐式求导给 $z'_N=O(1/(Nz_N))$，从而 $\tau_n-\tau_N=O(k\sqrt L/N^{3/2})=o(L/N)$，$q_{0,n}=q_{0,N}-\mu k+o(1)$。因此 $t$ 可由维数 $n$ 的两个固定临界系数门槛夹住。r094 给 $\log H_n^{\mathcal L}(t)\le q_{0,n}+c_Q L+o_B(L)$。

对任意 $P\ge2$, GMS 的一维强型界和逐坐标支配给 $H_k(s)\le s^{-P}\exp(k b_P)$，其中 $b_P=P\log c_P=\mu+O(1/P)$ (r001)。令 $P=\max(2,\sqrt{k/\log(k+L)})$，在第一项 $s=1/[2A(k+L)]$，于是 $\log H_k(s)\le\mu k+O_A(\sqrt{k\log(k+L)}+\log(k+L))=\mu k+o(L)$。这抵消了后块基准中的 $-\mu k$。

在第二项，$1-v\asymp\sqrt{L/N}$；取 $P=2$，有 $\log H_k(\delta/[2(1-v)])\le2k\log(1+\sqrt2)+L+O(1)$。为支付此项使用一般指标界 r091 (6)。在测试秩 $q=q_{0,n}-4(k+L)$，该式对应的正态自变量为 $-z_n+O((k+L)/\sqrt N)$，其中粗常数 $C_H$ 和 $t_0=O(L)$ 已包含在误差内。因为 $z_n(k+L)/\sqrt N=o(1)$，整个短区间上密度 $O(L/\sqrt N)$。故此秩的正态 CDF 至少 $\tau_n-C(k+L)L/N$。选固定 $A$ 足够大，使 $1-v+2n^{-4}$ 小于该下界。由正态分位的单调性，$\log H_n(v)\le q_{0,n}-4(k+L)$。代入 (1)，第二项至多 $\exp(q_{0,N}-2L)$；第一项给所需上系数。所有界对满足同一指定增长序列的 $E$ 一致。

完整 $N$ 维等权双曲十字属于这一族：固定任意前 $k$ 个坐标后，截面为有坐标截断的对数半空间、整箱或空集。r093 的深度极限和 r092 的确定向量展开给匹配下系数。这个下界保留所有 $N$ 个坐标的基准。

对于任意指标源，(1) 仍成立，只须将第一项的截面价换成 $\int|\{M_n1_{E_z}>t\}|dz/|E|$。一般上界的未闭合项因此是低门槛的真实截面平均；反复剥去小块时误差如何累计仍待证明。

## 对任意源的截面必要条件

固定 $\gamma>0$。若某列一般指标源的完整价达到 $\exp(q_{0,N}+(c_Q+\gamma)L+o(L))$ 或更高，则任取上述大小的坐标前块，(1) 的高门槛项仍可支付。用概率测度 $d\pi(z)=|E_z|dz/|E|$，记 $h(z)=|\{M_n1_{E_z}>t\}|/|E_z|$（空截面略去）。有
\[
 \int h\,d\pi\ \ge\ \exp(q_{0,n}+(c_Q+\gamma)L-o(L)).
\]
由 r091 (6) 的正态逆秩在 $t$ 处估计，$h(z)\le\exp(q_{0,n}+C_H L+O_B(1))$。于是对每个固定 $0<\epsilon<\gamma$，集合 $h(z)>\exp(q_{0,n}+(c_Q+\gamma-\epsilon)L)$ 占据的**源体积分数**至少 $\exp((c_Q+\gamma-C_H)L-o(L))$。证明是从积分下界中减去门槛，再除以上述逐截面上界。这些高价截面对大 $N$ 必须超出 r094 几何类。坐标集合预先固定，各集合对应的异常截面可以不同。这个必要条件也适用于任意源的微小高价子序列，只需在同一子序列使用相同对数余量。

一维强型根方程为 $(P-1)c_P^P-Pc_P^{P-1}-1=0$。它给出上文使用的 $b_P$ 展开。[Grafakos 与 Montgomery-Smith 的原始论文](https://doi.org/10.1112/S0024609396002081)

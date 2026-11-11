# 第 103 轮鞅传递的独立核验

核验 `moment-audit.md` 及所引 r011 §1–2、r092 (7)–(8)、r094 (5)、(11)–(15) 和 r100 截面类。对每个给定整数列 $k=o(\sqrt{N\log N})$，本稿证明了上、下匹配的 $q_{0,N}+c_Q\log N+o_B(\log N)$。上界对该序列的全部可测物理源一致，要求几乎每个后缀截面属于 r100 的无条件对数凸下集及其给定的几何变换；截面可随前缀任意变化。

## 扩大的局部逆窗

令 $n=N-k,L=\log N,\rho=\sqrt{L/N}$。两个固定临界系数夹住 $p_N$，r092 给 $s_0=q_{0,n}+c_Q\log n+o_B(L)$。即使把 $s_0+x$ 写作 $q_{0,n}+h$，在 $|x|\le W,\ W\ge L,\ \rho W=o(1)$ 内也有 $|h|=O_B(L+W)$。将基准 Legendre 值在 $q_{0,n}$ 展开，鞍点约为 $z_n/(\sigma\sqrt n)$，指数误差为 $O_B(\rho(L+W)+(L+W)^2/n)=o(1)$；其三阶余项 $n\rho^3=O(\rho L)$ 已计入。对于竞争向量 $F_S(q)\le2\tau_n$，r094 (11) 及 Jensen 缺损给

\[
 c n(s-s_*)^2+A+P/2\le C+\tfrac12\log(1+P/L).
\]

此处基准率与 $\log(2\tau_n)$ 的 $-z_n^2/2$ 相消，粗前因子的 $\tfrac12\log L$ 与 Mills 项 $-\log z_n$ 相消。故三种缺损全为 $O(1)$，$s\asymp s_*\asymp\rho$。r094 (12)–(15) 的平滑坐标比例、三阶误差与精确倾斜 Mills 比较适用，给 $F_S(q)\ge e^{-C/\sqrt n}F_{c_q}(q)$。对非竞争者，基准 iid 的 Edgeworth 或相对 Mills 估计给 $F_{c_q}(q)=(1+o(1))\tau_n$，自动满足同一不等式。上述 $o(1)$ 对指定窗中的 $q$ 一致。

比较 $s_0$ 与 $q=s_0+x$ 时，必须对两个端点均使用同一个动作 $c_q$。因 $s_0$ 为等权确定向量的分位上确界，$F_{c_q}(s_0)\ge1-p$。两端绝对 CDF 展开误差为 $O(1/n)$，不用乘窗口长度。匹配高斯密度相对于 $1/\beta$ 的误差为 $O(\rho L+\rho|x|)$，三阶修正相对导数为 $O((1+z_n^3)/\sqrt n)=O(\rho L)$。积分允许 $x$ 取任意符号；乘 $e^{-C/\sqrt n}$ 另损 $O(\tau_n/\sqrt n)$。在秩尺度上这两种绝对误差分别为 $O(1/L),O(1/\sqrt L)$，可并入常数。这给稿件 (1) 右端。由于 $|x|/(\beta\tau_n)=O(\rho W)=o(1)$，该阈值始终严格位于 $(0,1)$。r094 (5) 对实际方向的权重取 $\lambda_i=nw_i$ 后，$q$ 是所需全体确定向量命中分位的共同上界。截面逐一积分再用右连续的广义逆，得到整个 $N$ 维空间的重排，不需对截面几何参数作可测选择。

## 前缀、严格成功与鞅

真实全维框的平均由同一前缀框对 $g(z,y)=M_n1_{E_z}(y)$ 的平均支配。r011 全环境重排可单步施行在任一前缀轴上；保序迭代 $k$ 次。可先取有界紧盒截断 $g_m\uparrow g$，坐标极大算子沿该极限递增，$g_m^*\le g^*$。若实际 halo 体积大于 $a$，在 $(a,|\mathcal H|)$ 内选正常秩 $a'$，迭代重排在该处严格超过 $p$；由 Bellman 有限时域的近最优策略及 $a\le a'$，存在可预测策略给 $E_\pi g^*(ae^{-S_k})>p$。仅需近似取优，无需等号策略存在。零测截面修改不改变盒平均。

以 $M_k=S_k-\sum E[Z_i\mid\mathcal F_{i-1}]$ 表示鞅差和，则 $S_k\le\mu k+M_k$ 路径成立，且 $EM_k=0$, $EM_k^2\le Ck$。Bernoulli 位移 $c$ 的条件发生率为 $1/(1+e^c)$，故中心化双侧矩生成函数在一个固定零邻域统一满足 $\log E(e^{t\Delta M_i}\mid\mathcal F_{i-1})\le Ct^2$。迭代、Chernoff 得文中 Bernstein 型尾界。重排单调性保证把 $S_k$ 换成上界可给期望上界。

以 $x=\epsilon L-M_k$ 取 $W=A(\sqrt{kL}+L)$，$A$ 可以依赖固定 $\epsilon$。在所证范围 $\rho W=o(1)$，窗外概率可压至 $N^{-4}$；用 $g^*\le1$ 和 Cauchy–Schwarz 支付包括截断线性项在内的 $o(1)$ 秩误差。在窗内取期望再估计二次项，$Ex^2\le Ck+\epsilon^2L^2$。从而期望相对于 $p$ 的秩增量上界是

\[
 -\epsilon L+O_{B,\epsilon}(1+\rho L(\sqrt k+L)+\rho(k+L^2))+o(1).
\]

关键恒等式 $\rho k/L=k/\sqrt{NL}\to0$，并且 $\rho\sqrt k\to0$，$\rho L\to0$，故上式对充分大的 $N$ 为负。$q_{0,n}+\mu k-q_{0,N}=d_N-d_n$。由 $z'_t=O(1/(t z_t))$ 有 $d'_t=O(\rho)$、$|\tau'_t|=O(\rho/N)$（$n\le t\le N$），故两项维数迁移分别为 $O(k\rho)=o(L)$ 和 $O(k\rho/N)=o(L/N)$。固定 $\epsilon$ 得上极限后再令其趋零；完整 $N$ 维等权十字属于该类并给匹配下界。

## 幂矩前因子及有效范围

截断幂矩的整体 Gamma 权确为 $e^P\Gamma(P+1)/P^P\sim\sqrt{2\pi P}$。$W=A(\sqrt{PL}+L)$ 在 $P\rho=o(1)$ 时使 $\rho W^2=o(L)$，窗外的 Gamma 浓度可支付粗截距 $C_BL$。故该路线所得上界保留 $\tfrac12\log P$；当 $k=N^{o(1)}$ 可取 $P=\max(L,k)$ 并吸收这一项。鞅证明绕过该前因子，使锐系数覆盖全部 $k=o(\sqrt{NL})$。若 $k\asymp\sqrt{NL}$，单是 $d_N-d_n$ 就是正的 $\Theta(L)$：由 $d'_t\sim(\sigma/2)\sqrt{\log t/t}$ 可见，现有无符号二次误差也为 $O(L)$。本文的估计可以容许这一量级的额外对数成本，不能凭此给临界常数倍列的锐系数。

本机 Python 单进程标量复核：$P=16,64,256$ 时 Gamma 前因子除以 $\sqrt{2\pi P}$ 分别为 $1.0052212393,1.0013029208,1.0003255737$。取 $k=\sqrt{NL}/\log L$，$L=20,80,160$ 时，$\rho k/L$ 与 $1/\log L$ 双精度吻合；这些数值仅核对尺度，以上范围来自解析不等式。

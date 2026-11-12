# 第216轮　平方根熵界近等号迫使实际见证共享同一族密度层

2026-11-12。固定中央门槛 \(0<p<1\)，\(u=\bar\Phi^{-1}(p)\)。沿用第215轮的任意有限成功族及实际归属。设沿维数序列，联合量在平方根精度达到已证上包络：
\[
 V+D_X=\mu n+\sigma\sqrt n\,\varphi(u)/p+o(\sqrt n).
\]
这在第213轮构造中确实成立。不要求额外指定几何深度函数。用全局实际源密度定义同一批嵌套子集
\[
 F_z=\{x\in E:\log(|U|f_X(x))>\mu n+\sigma\sqrt n\,z\}.
\]
对任何固定 \(z>u\)，在均匀接收点的负责盒 \(R\) 中，
\[
 \frac{|R\cap F_z|}{|R|}\longrightarrow\bar\Phi(z),\qquad
 \frac{|R\cap F_z|}{|E\cap R|}\longrightarrow\bar\Phi(z)/p
 \quad\text{于 }L^1(P_Y).
\]
对 \(z\) 在 \((u,\infty)\) 内的任意紧区间，盒内整条曲线的上确界误差亦于接收概率下趋零。结论关联不同见证中同一个全局源层，强于仅有源边际抽样的尾分布；量词是按实际归属的接收概率，不涵盖每个接收点的任意见证。

**普适单侧覆盖费用。** 置 \(a_R(z)=|F_z\cap R|/|R|,\ r_R=|E\cap R|/|R|>p\)。因 \(|F_z|\le|U|e^{-\mu n-\sigma\sqrt n z}\)，对任意固定 \(h>\bar\Phi(z)\)、\(h<1\) 有
\[
 \Pr_Y(a_R(z)>h)\le H_n(h)|F_z|/|U|
 \le \exp\{-\sigma\sqrt n(z-\bar\Phi^{-1}(h))+O_h(\log n)\}.
\]
这里只使用第140轮中央粗界，没有使用近等号。因此 \(\mathbb E(a_R(z)-\bar\Phi(z))_+\to0\)。

**熵近等号使平均尾取满。** 写 \(Z_n=(\log(|U|f_X(X))-\mu n)/(\sigma\sqrt n)\)。由实际抽样，
\(\Pr(Z_n>z)=\mathbb E[a_R(z)/r_R]\)。上段和 \(r_R>p\) 给出逐点上极限 \(\min(1,\bar\Phi(z)/p)\)。正尾一致可积：在第211轮中尾证明直接使用 \(|F_z|/|U|\le e^{-\mu n-\sigma\sqrt n z}\)，取 \(q=2n^{-4}+\bar\Phi((\sigma\sqrt n z-A\log n)/(\sigma\sqrt{n+t_0}))\)，对 \(z\ge1\) 取 \(A\) 充分大，得到 \(q/p+n^{-3}\) 的尾上界，积分到 \(z=K\sqrt n\)。再用 \(H_n(q)\le24^n/q^2\)，取 \(q=e^{-\sigma\sqrt n(z-K\sqrt n)/3}/2\)，\(K\) 足够大，余尾积分为 \(O_p(n^{-1/2})\)。这些估计均无需给 \(|U|/|E|\) 下界。于是任取 \(z_0<u\)，\(\mathbb E\max(Z_n,z_0)\) 的上极限至多为 \(\mathbb E[G\mid G>u]=\varphi(u)/p\)。假设给 \(\mathbb E Z_n\to\varphi(u)/p\)，因此从 \(z_0\) 起的尾积分饱和。逐点上极限、单调性及极限上包络的连续性迫使每个 \(z>u\) 处的尾趋于 \(\bar\Phi(z)/p\)：若某点沿子列少一个固定量，右边一小段上的积分也会少一个固定量，其他区间由上极限及尾可积性不能弥补。

**从源边际回到同一实际盒。** 令 \(q_z=\bar\Phi(z)\)。由 \(p\Pr(Z_n>z)\le\mathbb E a_R(z)\) 得均值下极限至少为 \(q_z\)，上段单侧覆盖界于是给 \(a_R(z)\to q_z\) 于 \(L^1\)。同理 \(a_R/r_R\le a_R/p\) 且两边均值之差趋零，得第二个结论。因 \(q_z>0\)，这还推出 \(r_R\to p\) 于接收概率下。对于紧区间，先在有限网格上并估各点误差，再用曲线单调及正态尾一致连续插值即可。由任意固定小 \(\epsilon>0\) 下 \(P_X(F_z\setminus F_{z+\epsilon})\to(q_z-q_{z+\epsilon})/p>0\)，还得 \(\log(|U|/|F_z|)=\mu n+\sigma\sqrt n z+o(\sqrt n)\)：在此差集上密度至多 \(e^{\mu n+\sigma\sqrt n(z+\epsilon)}/|U|\)，最后令 \(\epsilon\downarrow0\)。

一般端点系数差仍为 \((\mu-1)/8\)。平方根熵近等号是本定理的假设，近中央体积本身尚未推出该假设。下一步可研究这些全局源层在移动门槛 \(q_z\) 附近是否继承近极值体积，或在联合量未饱和时能否把其缺口用于逐门槛估计。

# 第218轮　按联合熵缺口定量结算源层移动

2026-11-12。固定 \(0<p<1\)、目标紧区间 \(J\Subset(0,p)\) 及允许损失的接收比例 \(0<\theta<1\)。沿用第216轮的全局密度层 \(F_t\)、\(V,D_X\)。令
\[
 K_n=\log n+[\mu n+\sigma\sqrt n\,\varphi(u)/p-(V+D_X)]_+,\qquad \delta_n=K_n/\sqrt n.
\]
若 \(\delta_n\to0\)，则对每个精确目标 \(q\in J\)，可选择源层 \(F_{z-h(q)}\)，其中 \(z=\bar\Phi^{-1}(q)\)。该层上严格超过 \(q\) 的原盒并 \(U_q\) 覆盖至少 \((1-\theta)|U|\)，且存在同一套可测选择满足
\[
 \sup_J h(q)\le C_{J,p}\sqrt{\delta_n/\theta}+n^{-1},
 \qquad \int_J h(q)dq\le C_{J,p}\delta_n/\theta+|J|/n .
\]
对数体积比至少为 \(B_n(q)-\sigma\sqrt n\,h(q)+\log(1-\theta)\)。特别在 \(K_n=O(\log n)\) 时，给定固定 \(\theta\)，平均对数损失为 \(O(\log n)\)，而预先指定门槛的最坏保证为 \(O(n^{1/4}\sqrt{\log n})\)。

**缺口控制曲线之间的面积。** 写 \(a_R(t)=|R\cap F_t|/|R|\)、\(Z=(\log(|U|f_X(X))-\mu n)/(\sigma\sqrt n)\)、\(S(t)=\Pr(Z>t)\)、\(S_*(t)=\min(1,\bar\Phi(t)/p)\)。取固定区间 \([A,B]\subset(u,\infty)\) 在两端严格包住目标分位区间。我们有
\[
 \int_A^B\mathbb E(\bar\Phi(t)-a_R(t))_+dt\le C_{J,p}\delta_n
\]
（可依 \(J,p\) 预先选定端点）。具体地，取 \(t_0<u\) 固定；第216轮的中尾及极远尾比较给 \(\int_{t_0}^\infty(S-S_*)_+dt=O_p(\log n/\sqrt n)\)。为明确精度，在 \(t\le u\) 无正差，从 \(u\) 到 \(K\sqrt n\) 用第140轮选门槛 \(2n^{-4}+\bar\Phi((\sigma\sqrt n\,t-C\log n)/(\sigma\sqrt{n+t_0'}))\)，其中 \(t_0'=O(\log n)\) 为时间膨胀量。附加的 \(O(n^{-3})\) 概率误差可积，正态尾相对原尾的积分绝对差是 \(O_p(\log n/\sqrt n)\)（同一正态变量作仿射耦合，比较尾面积）。余段按第216轮极远尾估计为 \(O_p(n^{-1/2})\)。因 \(\mathbb E\max(Z,t_0)\ge\mathbb EZ\)，上包络的均值为 \(\varphi(u)/p\)，故 \(\int_{t_0}^\infty(S_*-S)_+dt\le C_p\delta_n\)。

对 \(t\in[A,B]\)，由 \(|F_t|/|U|\le e^{-\mu n-\sigma\sqrt n t}\)，取门槛 \(\bar\Phi(t-C'\log n/\sqrt n)\)，一般中央粗界使超过它的盒负责比例至多 \(n^{-2}\)。于是 \(\mathbb E(a_R(t)-\bar\Phi(t))_+=O_{J,p}(\log n/\sqrt n)\) 一致。结合 \(\mathbb E a_R(t)\ge pS(t)\)，即得所称面积界。

**把面积转成可严格成功的层。** 对每个 \(q\in J\)，接收坏比例 \(b_q(t)=\Pr_Y(a_R(t)\le q)\) 随 \(t\) 非减。对 \(A\le t<z\)，若 \(b_q(t)>\theta\)，则 \(\mathbb E(\bar\Phi(t)-a_R(t))_+\ge\theta(\bar\Phi(t)-q)\)。坏比例超过容许值的 \(t\) 在 \([A,z]\) 形成末段，长度记为 \(l(q)\)。由于 \(\varphi\) 在 \([A,B]\) 上有正下界，积分得 \(\theta l(q)^2\min_{[A,B]}\varphi/2\le C\delta_n\)。特别是坏段最终不会延伸至 \(A\)。另一方面，交换 \(t,q\) 积分，由 Markov 界
\[
 \int_J l(q)dq\le\theta^{-1}\int_A^B\int_{q<\bar\Phi(t),\,q\in J}\Pr(a_R(t)\le q)dq\,dt
 \le C\delta_n/\theta .
\]
取 \(h(q)=l(q)+1/n\) 避开端点平局，保证至少 \(1-\theta\) 的实际接收归属给严格成功盒。所有选择均由单调可测的坏比例决定。最后由源层体积上界直接得到对数比下界；无需把重新抽样的熵误差误当作已有定量速率。

一般系数差仍为 \((\mu-1)/8\)。本轮量化了第217轮的体积传递，没有量化新联合熵的饱和速率。门槛平均与逐点之间仍有平方根损失。下一步需检验层曲线的真实矩形结构是否排斥将全部缺口集中于窄门槛段，否则面积控制本身不足以给逐点对数精度。

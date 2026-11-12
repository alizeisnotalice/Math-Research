# 第318轮　超临界双曲壳的平方根项与对数级恢复门槛

2026-11-12。设 $T=n^4$，$S=E_T\setminus E_{T+h}$，固定 $0<p<1$。置 $\lambda=h/\sqrt{\mu n}$，假设 $\lambda\ge\lambda_0>2\Phi^{-1}((1+p)/2)$。定义 $u(\lambda)$ 为
\[
 \Phi(u+\lambda)-\Phi(u)=p
\]
的较大根。对所有这样的 $h$（无上界）一致有
\[
 \log\frac{|\{M_S1_S>p\}|}{|S|}
 =\mu n+\sqrt{\mu n}\,u(\lambda)+O_{p,\lambda_0}(1).
\]
因此宽度是有限超临界常数倍 $\sqrt n$ 时，虽然恢复线性指数率，仍比完整中央基线低严格正的平方根项。以下证明不覆盖切点 $\lambda=2\Phi^{-1}((1+p)/2)$。

**均值损失控制全部动作的方差偏离。** 沿用[第317轮](../r317/report.md)截后单轴真实深度的矩 $(v_i,m_i)$。在两条动作曲线上一致有 $|v_i-\mu|\le C\sqrt{\mu-m_i}$：穿零均值在 $r_*=e^{-\mu}$ 有唯一非退化最大值（其二阶导在该点为 $-1/[r_*(1+r_*)]$），方差光滑且 $v(r_*)=\mu$；邻域外由紧性，单侧动作由 $m_i\le1$。写 $D=\sum(\mu-m_i)$、$v=n^{-1}\sum v_i$，得 $|v-\mu|\le C\sqrt{D/n}$。

**移动壳端点的统一分位。** 令 $a=h/\sqrt n$，记 $w(v,a)$ 为 $\Phi((w+a)/\sqrt v)-\Phi(w/\sqrt v)=p$ 的较大根。$w(\mu,a)=\sqrt\mu\,u(\lambda)$。在 $v$ 的一个固定 $\mu$ 邻域内，对 $a\ge\sqrt\mu\lambda_0$，根一致有界且左端点导数（对 $w$ 的导数）严格为负，绝对值有统一正下界。这对有界参数由远离切点和紧性得到；$a\to\infty$ 时较大根一致趋于 $\sqrt v\,U(p)$，$U=\Phi^{-1}(1-p)$，由正态尾在有界左端点上的一致收敛及单调根夹逼即得。方程对 $v$ 的偏导也一致有界，因为 $|t|\varphi(t)$ 有界。因此 $|w(v,a)-w(\mu,a)|\le C|v-\mu|$；命中概率误差 $O(n^{-1/2})$ 给较大根增加至多 $O(n^{-1/2})$（同样论证对稍微移动的 $p$ 一致）。

若 $D\ge\epsilon n$，Chebyshev 已把最大可达深度位移限制为 $\mu n-D+O_p(\sqrt n)$。其余情形取 $\epsilon$ 足够小，使 $v$ 位于上述邻域。由第317轮的统一三阶矩，Berry–Esseen 将任意真实壳区间的命中概率从上方控制在相应正态区间概率加 $O(n^{-1/2})$。故截后基准深度到 $T$ 的位移至多
\[
 \mu n-D+\sqrt n\,w(v,a)+O(1)
 \le \mu n+\sqrt n\,w(\mu,a)-D+C\sqrt D+O(1).
\]
最后两损失项之和有绝对常数上界（常数允许依赖 $p,\lambda_0$）。用第315轮外部面求和，此位移给全空间对数体积相同的上界，源分母由 $|E_T|$ 换为 $|S|$ 仅付 $O(1)$，因为 $h\gtrsim\sqrt n$ 且 Gamma 尾比随 $h$ 递减。

反向在立方体内部各轴采用 $r_*$。均值和方差分别恰为 $\mu n,\mu n$。取接收深度使所需左端点为 $\mu n+\sqrt{\mu n}\,u(\lambda)-C$，并在该深度附近取长度一的带。选固定 $C$ 足够大，正态命中率的增加压过 Berry–Esseen 误差，整带真实平均严格超过 $p$。其 Gamma 体积相对 $|E_T|$ 的对数等于所取位移加 $O(1)$，统一性来自位移为 $O(n)$ 和 $T=n^4$。

**恢复到 $o(\log n)$ 的尺度。** 若 $\lambda\to\infty$，则
\[
 U(p)-u(\lambda)\sim \bar\Phi(U(p)+\lambda)/\varphi(U(p)).
\]
由方程知 $\bar\Phi(u)-\bar\Phi(U)=\bar\Phi(u+\lambda)$。$u\to U$，均值定理先给 $0<U-u\le C\bar\Phi(U-1+\lambda)$，故 $\lambda(U-u)\to0$，再用正态 Mills 比即得等价式。因此距 $B_n(p)$ 为 $o(\log n)$ 当且仅当 $\sqrt n\bar\Phi(U+\lambda)=o(\log n)$。例如 $U+\lambda=\sqrt{\log n-3\log\log n+c}$、$c$ 固定时，对数损失除以 $\log n$ 趋于 $\sqrt\mu\,e^{-c/2}/(\sqrt{2\pi}\varphi(U))$。本例把线性指数率恢复与对数精度恢复区分开。

一般上下界系数差仍为 $(\mu-1)/8$。下一步检验临界切点处根合并的尺度，并研究一般源是否存在可类比的上下两端命中约束。

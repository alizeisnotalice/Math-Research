# 中央门槛的一个统一中位数界

取独立 $Z_i=Y_i+c_i B_i$，其中 $Y_i$ 为单位指数，$B_i$ 的成功概率 $p_i=(1+e^{c_i})^{-1}$，$c_i\ge0$，端点 $c_i=\infty$ 用纯指数。设 $S=\sum w_i Z_i$，$V=\operatorname{Var} S>0$，$a_i=w_i/\sqrt V$，$w_i\ge0$，$\delta=\max a_i$。存在绝对常数 $\delta_0,c>0$，只要 $\delta\le\delta_0$，就有
\[
 P(S\le ES)\ge 1/2+c\sum a_i^3.
\]
于是对 $\sum w_i=1$ 的这类法向，所有确定动作的中位数都不超过 $\mu$。这个假设允许单个权重远大于 $1/n$。

各阶中心矩有统一界，因为 $c_i^j p_i$ 有界。方差在 $[1,C]$ 内，所以 $1/C\le\sum a_i^2\le1$。第三累积量是 $2+c_i^3p_i(1-p_i)(1-2p_i)\ge2$。记 $\gamma=\sum a_i^3\kappa_3(Z_i)$，$b_4=\sum a_i^4$。中心标准化和的特征函数 $f(t)$ 在 $|t|\le c_0/\delta$ 上满足
\[
 \log f(t)=-t^2/2-i\gamma t^3/6+O(b_4t^4).
\]
每个中心化因子在原点附近统一远离零，统一四阶矩给该对数 Taylor 余项。选 $c_0$ 小，使余项实部小于 $t^2/4$，再对纯虚三次相位使用二阶余项，就能以 $C e^{-t^2/4}(b_4t^4+\gamma^2t^6)$ 控制 $f(t)-e^{-t^2/2}(1-i\gamma t^3/6)$ 的绝对值，且 $\gamma^2\le C b_4$。

全轴上指数分量给 $|f(t)|\le\prod(1+a_i^2t^2)^{-1/2}\le(1+\delta^2t^2)^{-1/(2C\delta^2)}$，最后一步为对数凹性在 $[0,\delta^2]$ 的弦界。将此上界在 $|t|>c_0/\delta$ 对 $dt/|t|$ 积分得到 $O(e^{-c_1/\delta^2})$：换元 $s=\delta |t|$，分出一半幂次提取端点指数，另一半在 $[c_0,\infty)$ 可积且有统一界。高斯多项式尾也满足同样类型估计。这些尾部是 $O(b_4)$，因为 $b_4\ge\delta^4$。

最后在零点用 Gil–Pelaez 反演（分布有密度），得 $P(S\le ES)=1/2+\phi(0)\gamma/6+O(b_4)$。因 $b_4\le\delta\sum a_i^3$，选择充分小的绝对 $\delta_0$ 即得严格余量。由 $ES\le\mu\sum w_i$ 得中位数断言。若有限体积无条件对数凸下集几乎所有法向均满足 $\max_i w_i/(\sum w_i^2)^{1/2}\le\delta_0$，r094 的真实角向积分便给完整中央门槛 $H(E,1/2)\le e^{\mu n}$。其余法向与一般指标源的全空间体积仍待结算。

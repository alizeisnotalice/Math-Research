# 强连通分量大小随维数增长时的中位数误差

沿用[线性锚点证明](linear_anchors.md)的 $Y_b=-\log|P_b(V)|$、$Y=\sum_bY_b$ 和独立 $Z\sim\mathrm{Exp}(1)$。下面的 $C,C_0,C_1,\ldots$ 均为绝对常数，不依赖矩阵、$n$ 或最大块大小 $K$（$1\le K\le n$）。对每个维数，式 (1) 的右侧可以取为
\[
 \mu n+\exp(C K^2\log(2K)). \tag{1}
\]
因此 $K^2\log(2K)=o(\log\log n)$ 时，中位数上界达到 $\mu n+o(\log n)$。该结论仍需整个见证规则采用一个指定的线性锚点。

**矩的绝对常数。** 对 $k$ 变量多仿射 $F$，若 $\|F\|_\infty=1$，原证明逐轴归纳中的一维常数可取 4，每次递推至多再加 4，故 $P(|F|\le e)\le4k e^{1/k}$，$0<e<1$。取 $H=\|P_b\|_\infty\ge1$，$U=-\log|P_b/H|\ge0$，于是 $U$ 被 $k\log(4k)+kZ$ 随机支配。特别 $\mathbb EU\le k(\log(4k)+1)$。对 $j=2,3$，令 $a_k=k(1+\log(4k))$，由指数变量的固定阶矩及 $Y_b=U-\log H$ 有
\[
 \mathbb E|Y_b-\mathbb EY_b|^j\le C_j a_k^j. \tag{2}
\]

**好块的方差。** 若 $\mathbb EY_b\ge\mu k/2$，则 $\log H\le a_k$。各轴偏导至多 $2H/(1+r)$（固定其他坐标后，仿射函数两端差不超过 $2H$）。取 $\delta=(1+r)/(8kH)$；在 $[1-\delta,1]^k$ 上 $|P_b(v)-1|\le1/4$，故 $Y_b\le\log(4/3)<1/3$，与 $\mathbb EY_b\ge\mu k/2\ge k/2$ 相差至少 $k/6$。该角盒的概率为 $(8kH)^{-k}$，所以
\[
 \operatorname{Var}Y_b\ge(k/6)^2(8kH)^{-k}\ge\exp[-C_0 K^2\log(2K)]. \tag{3}
\]
这里用 $k\le K$ 和 $\log H\le a_k$，吸收了 $K=1$ 的固定系数。

**从概率误差转回原始秩。** 写 $m=\mathbb E(Y+Z)$、$V=\operatorname{Var}(Y+Z)\ge1$、$B=\sum_b\mathbb E|Y_b-\mathbb E Y_b|^3+\mathbb E|Z-1|^3$。因 $\sum_b k_b=n$，(2) 给 $V\le C n K(1+\log(4K))^2$、$B\le C n K^2(1+\log(4K))^3$。若 $\mathbb EY\le3\mu n/4$，Chebyshev 给上尾严格小于 $1/2$ 的一个秩为 $m+2\sqrt V$。以 $-\mu n/4+C\sqrt{nK}(1+\log(4K))\le C'K(1+\log(4K))^2$ 完成平方，即得 (1)，包括 $n$ 较小的情形。

若 $\mathbb EY>3\mu n/4$，逐块均值上界 $\mu k_b$ 给至少 $n/(2K)$ 个好块。因此 $V\ge n e^{-C_0K^2\log(2K)}/(2K)$，故 $B/V\le\exp(C_1K^2\log(2K))$。为使秩界在 Berry–Esseen 概率误差较大时仍成立，置 $q=C_{\mathrm{BE}}B/V^{3/2}$，$C_{\mathrm{BE}}>0$ 为固定绝对常数，$\phi$ 为标准正态密度。若 $q\le\phi(1)/2$，取标准化位移 $z=2q/\phi(1)\le1$，则 $\Phi(z)-1/2\ge z\phi(1)=2q$，Berry–Esseen 保证在秩 $m+z\sqrt V$ 上上尾严格小于 $1/2$。秩费用是 $O(B/V)$。若 $q>\phi(1)/2$，则 $\sqrt V\le (2C_{\mathrm{BE}}/\phi(1)) B/V$，直接取 $m+2\sqrt V$ 使用 Chebyshev，费用仍为 $O(B/V)$。两种情形均无需另设起始维数。由 $m\le\mu n+1$ 和线性换元的真实中位数比较得 (1)。

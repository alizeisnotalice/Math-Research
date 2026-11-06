# 软指数轮廓下界审计

日期：2026-11-06。结论：取 \(r=e^{-\mu}\)、\(\sigma^2=\mu\)，并对 \(n\ge1024\) 设
\[
 z=\sqrt{\log n-2\log\log n},\qquad q=n\mu-\sigma\sqrt n\,z,
\]
则所审计的软轮廓论证可给出
\[
 \log C_n\ge n\mu-\sigma\sqrt{n(\log n-2\log\log n)}-\log(2n).
\]
没有发现倾斜符号、密度界、严格性或截断顺序的错误。

令 \(W=\mu-Z_r\)，其累积量
\[
 K(\lambda)=\lambda\mu-\log(1+\lambda)
 +\log(1+r^{1+\lambda})-\log(1+r)
\]
满足 \(K(0)=K'(0)=0\)、\(K''(0)=\mu\)，且
\[
 K'''(\lambda)=-{2\over(1+\lambda)^3}
 -\mu^3p_\lambda(1-p_\lambda)(1-2p_\lambda)<0
\]
对 \(\lambda\ge0\)。故 \(K(\lambda)\le\mu\lambda^2/2\)。在指数倾斜下，\(\sum W_i\) 是平移负 Gamma 的混合，密度上确界至多
\((1+\lambda)/\sqrt{2\pi(n-1)}\)。以
\(\lambda=z/(\sigma\sqrt n)\) 积分而非仅用 Chernoff，得
\[
 \delta=P(\sum W_i\ge\sigma\sqrt n z)
 \le {1+\lambda\over\lambda\sqrt{2\pi(n-1)}}e^{-z^2/2}.
\]

对 \(g(t)=e^{(t-q)/n}1_{\{t\ge q\}}\)，Jensen 和坏事件切除给
\(Eg(T)\ge e^{\sigma z/\sqrt n}-\delta>1\)。其充分条件的比值为
\[
 R={ (1+\lambda)\sqrt{n/(n-1)}\log n\over
 \sqrt{2\pi}(\log n-2\log\log n)}<1
\]
在 \(n\ge1024\) 成立：\(\lambda<.083\)、其余两个因子分别小于 \(1.001\)、\(2.28\)，而 \(\sqrt{2\pi}>2.5\)。

精确模数积分为
\[
 J(g)=e^{-q}I_n,\qquad
 I_n=\sum_{k=0}^{n-1}{(n-1)_k n^{-k}\over(1-1/n)^{k+1}}
 \le {n^2\over n-1}\le2n.
\]
最后，对每个固定 \(n\) 以 \(g_M=\min(g,M)\) 截断。单调收敛选出有限 \(M\) 仍有严格 \(Eg_M(T)>1\)，且 \(J(g_M)\le J(g)\)，所以可合法调用固定维数的 \(L\to\infty\) 下界归约；这一步不声称维数一致收敛率。

Wolfram 15 精确求和在 \(n=1024\) 给 \(I_n\approx40.7962<2048\)；MATLAB R2026a 给 \(R\approx0.948035<1\)。这些运行仅作算术交叉核验，严格结论来自不等式链。

# r016：固定 \(c\) 下尾分位数的精确导数路线

日期：2026-11-06。本笔记审计一个避免“两个极小 CDF 相减”的固定比例比较方法。结论范围仍只限于等比例族 \(Z_c=Y+cJ_c\)。

## 1. 精确分位数导数恒等式

记 \(F_{n,c}\)、\(f_{n,c}\) 分别为 \(S_{n,c}=\sum_{i=1}^nZ_c^{(i)}\) 的 CDF 和密度，
\[
 M_c(\lambda)=\mathbb Ee^{-\lambda Z_c}
 =\frac1{1+\lambda}\frac{1+e^{-(1+\lambda)c}}{1+e^{-c}},
 \qquad p_c=\frac1{1+e^c}.
\]
直接微分给出
\[
 \partial_cM_c(\lambda)=p_cM_c(\lambda)-p_ce^{-\lambda c}. \tag{1}
\]
CDF 的 Laplace 变换是 \(M_c(\lambda)^n/\lambda\)。逐项反演 (1) 因而给出
\[
 \partial_cF_{n,c}(q)=np_c\{F_{n,c}(q)-F_{n-1,c}(q-c)\}. \tag{2}
\]
若 \(q_c(\beta)\) 由 \(F_{n,c}(q_c)=\beta\) 定义且密度为正，则
\[
 q_c'(\beta)=np_c\frac{F_{n-1,c}(q_c-c)-F_{n,c}(q_c)}{f_{n,c}(q_c)}. \tag{3}
\]
这一步是精确恒等式，且符号正确：分位数导数来自 \(-\partial_cF/f\)。

## 2. 用一个密度的局部 Taylor 展开替代微小 CDF 比较

置 \(N=n-1\)、\(a=q_c-c\)、\(Z=Z_c\)。由卷积
\(F_{n,c}(q_c)=\mathbb EF_{N,c}(a-(Z-c))\)，有
\[
 F_{n,c}(q_c)-F_{N,c}(a)
 =(c-m(c))f_{N,c}(a)+\frac{\mathbb E(Z-c)^2}{2}f_{N,c}'(a)+\mathcal R. \tag{4}
\]
令 \(\lambda>0\) 唯一满足倾斜后 \(N\) 项和的均值为 \(a\)。在本轮窗口
\(c=\mu+O(\sqrt{\log n/n})\)、\(\beta=\sigma\sqrt{\log n/n}\) 中，已有相对中偏差的粗分位数控制足以推出
\[
 \lambda=\frac1\sigma\sqrt{\frac{\log n}{n}}(1+o(1)). \tag{5}
\]

## 3. 倾斜 Fourier 引理：直接证明所需的可微密度控制

在 \(e^{-\lambda S_{N,c}}\) 倾斜下，连续指数部分的速率为 \(1+\lambda\)，离散 Bernoulli 参数仍在 \((0,1)\) 的固定紧集内。居中后一步特征函数的模满足
\[
 |\varphi_{c,\lambda}(t)|
 \le\left(1+\frac{t^2}{(1+\lambda)^2}\right)^{-1/2}. \tag{6}
\]
故 \(N\) 次幂在 Fourier 反演中可积，并且对 \(c\) 的紧邻域、\(0\le\lambda\le2\sqrt{\log n/n}\) 一致有
\[
 \|g_N\|_\infty=O(N^{-1/2}),\quad
 \|g_N'\|_\infty=O(N^{-1}),\quad
 \|g_N''\|_\infty=O(N^{-3/2}). \tag{7}
\]
这里 \(g_N\) 是倾斜后居中和的密度。式 (7) 直接由 Fourier 积分
\(\int |t|^j|\varphi|^Ndt\) 得到，不需要 Edgeworth 定理。

在中心，统一 Taylor 展开为
\[
 \varphi_{c,\lambda}(t)=1-v_{c,\lambda}t^2/2+O(|t|^3).
\]
因此在固定小邻域内
\[
 |\varphi_{c,\lambda}(t)^N-e^{-Nv_{c,\lambda}t^2/2}|
 \le CN|t|^3e^{-cNt^2}; \tag{8}
\]
邻域外由 (6) 指数小。积分 (8) 得
\[
 g_N(0)=(2\pi Nv_{c,\lambda})^{-1/2}+O(N^{-1}). \tag{9}
\]
对导数，Gaussian 主项是奇函数而积分为零；在 (8) 的差上乘 \(|t|\) 后积分为
\(O(N^{-3/2})\)。故 \(g_N'(0)=O(N^{-3/2})\)。因此，恢复到原密度后
\[
 \frac{f_{N,c}'(a)}{f_{N,c}(a)}=\lambda+O(N^{-1}). \tag{10}
\]

对 (4) 的 Lagrange 余项，(7)、(9)--(10) 给出，在线段 \(s\) 连接 \(a\) 与 \(a-(Z-c)\) 上，
\[
 \frac{|f_{N,c}''(s)|}{f_{N,c}(a)}
 \le C e^{\lambda(s-a)}(\lambda^2+N^{-1}). \tag{11}
\]
由于 \(c\) 在紧集且 \(Z_c\) 有一致指数矩，(11) 的右端可积，故
\[
 \mathcal R=f_{N,c}(a)\,O(\lambda^2+N^{-1}). \tag{12}
\]
同一密度 Taylor 展开还给出
\[
 \frac{f_{n,c}(q_c)}{f_{N,c}(q_c-c)}=1+O(\log n/n). \tag{13}
\]

## 4. 导数与 \(O(\log n)\) 改善

将 (4)、(10)--(13) 代入 (3)，并用
\(m(\mu)=\mu\)、\(v(\mu)=\mu\)、\(p_\mu=(\mu-1)/\mu\)，得到对
\(c=\mu+t\rho\)、\(\rho=\sqrt{\log n/n}\)、固定 \(t\) 一致的主项
\[
 q_c'(\beta_n)=-np_\mu\left\{(c-\mu)+\frac\mu2\lambda\right\}+O(\log n)+o(\sqrt{n\log n}). \tag{14}
\]
由 (5)，\(\mu\lambda/2=(\sigma/2)\rho+o(\rho)\)。从 \(0\) 到 \(t\rho\) 积分 (13) 得
\[
 \frac{q_{\mu+t\rho}(\beta_n)-q_\mu(\beta_n)}{\log n}
 \longrightarrow-
 \frac{p_\mu t^2}{2}-\frac{p_\mu\sigma t}{2}. \tag{15}
\]
在 \(t=-\sigma/2\) 处，(14) 的值为
\[
 \frac{p_\mu\sigma^2}{8}=\frac{\mu-1}{8}. \tag{16}
\]
所以只要第 3 节的统一 Fourier 局部引理完整写出，就能严格证明一个**具体的有界位移等比例候选**比 \(c=\mu\) 高
\((\mu-1)\log n/8+o(\log n)\)。它不需要全局优化所有 \(c\)，也不需要比较两个 CDF 的同阶 FGL 余项。

## 5. 当前证明状态

式 (1)--(4) 为完全代数的精确步骤。式 (6)--(13) 已给出充分的直接 Fourier 证明骨架：只需在正式稿中逐项写出统一 Taylor 常数和 Fourier 分区积分。它不依赖极小 CDF 差分或普通 Edgeworth 余项。剩余工作是把该骨架写成可审计的引理并核对所有端点参数，而非另寻更强的概率定理。

# r018：固定 \(c=\mu\) 的指数斜坡软根

日期：2026-11-06。本笔记只处理固定比例 \(c=\mu\) 的指数斜坡
\(h(u)=e^{u/n}\mathbf1_{u\ge0}\)，并核对软根的 \(O(\log n)\) 修正。记
\(L=\log n\)，\(d_n=\sigma\sqrt n z_n\)，其中
\[
 \overline\Phi(z_n)=\frac{\sigma z_n}{\sqrt n}. \tag{1}
\]

## 1. 基础量的实际 Wolfram 核验

对 \(Z_\mu=Y+\mu J\)，Wolfram 15 的符号微分和 25 位数值计算给出
\[
 \kappa_1=\mu,\qquad\kappa_2=\mu,\qquad
 \kappa_3=\mu(3-\mu)=2.20092204118594\ldots. \tag{2}
\]
指数斜坡的模数积分有精确式
\[
 I_n=n e^{n-1}\frac{\Gamma(n,n-1)}{(n-1)^n}
 \sim\sqrt{\frac{\pi n}{2}}. \tag{3}
\]
例如 \(I_n/\sqrt{\pi n/2}\) 在 \(n=100,1000,10000\) 分别为
\(1.05959,1.01742,1.00538\)。

## 2. 软根方程

令 \(S_n=\sum_{i=1}^nZ_\mu^{(i)}\)，并设
\[
 A_n(q)=\mathbb E\left[e^{(S_n-q)/n}\mathbf1_{\{S_n\ge q\}}\right]. \tag{4}
\]
它严格递减，因为
\[
 \partial_qA_n(q)=-A_n(q)/n-f_n(q)<0. \tag{5}
\]
故存在唯一根 \(q_{\rm soft}\) 使 \(A_n(q_{\rm soft})=1\)。写 \(q=n\mu-d\)。未截断部分满足
\[
 e^{-q/n}\mathbb Ee^{S_n/n}-1=\frac dn+O\left(\frac Ln\right). \tag{6}
\]
另一方面
\[
 0\le F_n(q)-\mathbb E\left[e^{(S_n-q)/n}\mathbf1_{\{S_n<q\}}\right]
 \le\frac1n\mathbb E(q-S_n)_+. \tag{7}
\]
对固定 \(c=\mu\)，负倾斜累积量满足 \(K''(\lambda)\le\mu\)。Chernoff 积分在
\(d\asymp\sqrt{nL}\) 给出右端为 \(O(\sqrt L/n)\)，因而至少为 \(O(L/n)\)。所以根方程可化为
\[
 F_n(n\mu-d)=\frac dn+O\left(\frac Ln\right). \tag{8}
\]

## 3. 足够精确的一项 CDF Fourier--Edgeworth 展开

在 \(z=O(\sqrt L)\) 的窗口，直接对标准化特征函数作四阶 Taylor 展开，并在 Gil--Pelaez/Fourier CDF 反演中保留三阶项，可得到
\[
 F_n(n\mu-\sigma\sqrt n z)
 =\Phi(-z)+\frac{\kappa_3}{6\mu^{3/2}\sqrt n}(1-z^2)\phi(z)
 +O(n^{-1}). \tag{9}
\]
这里余项是**绝对** \(O(n^{-1})\)，并且在 \(z\le2\sqrt L\) 一致。一个直接充分证明是把
\[
 \varphi(t/\sqrt n)^n
 -e^{-t^2/2}\left(1+\frac{\kappa_3(it)^3}{6\mu^{3/2}\sqrt n}\right)
\]
除以 \(|t|\) 后积分：小 \(t\) 区由四阶 Taylor 余项给出 \(O(n^{-1})\)，外区利用指数变量的显式特征函数衰减。这个误差远小于判定 \(O(\log n)\) 分位数位移所需的 \(L^2/n\)。普通 FGL 不能替代 (9)，因为它的相对误差恰与三阶项同阶。

## 4. 根的位置与得分

令
\[
 d=d_n-CL,
 \qquad z=\frac d{\sigma\sqrt n}.
\]
由 (1) 与 Mills 比，
\[
 \phi(z_n)=\frac{\sigma z_n^2}{\sqrt n}(1+o(1))
 =\frac{\sigma L}{\sqrt n}(1+o(1)). \tag{10}
\]
故正态尾因 \(-CL\) 的 \(d\) 位移增加
\[
 \overline\Phi(z)-\overline\Phi(z_n)
 =\frac{CL^2}{n}(1+o(1)). \tag{11}
\]
式 (9) 的三阶项在 \(z_n\) 处为
\[
 -\frac{\kappa_3}{6\mu}\frac{L^2}{n}(1+o(1)). \tag{12}
\]
而 (8) 右端相对 \(d_n/n\) 仅有 \(O(L/n)=o(L^2/n)\) 的变化。因此 (11)--(12) 强制
\[
 C=\frac{\kappa_3}{6\mu}=\frac{3-\mu}{6},
\]
即
\[
 q_{\rm soft}=n\mu-d_n+\frac{3-\mu}{6}\log n+o(\log n). \tag{13}
\]

由 (3)，固定 \(c=\mu\) 的斜坡认证得分为
\[
 \log\frac{e^{q_{\rm soft}}}{I_n}
 =n\mu-d_n-\frac\mu6\log n+o(\log n). \tag{14}
\]
若 r017 的实际斜坡参数移动定理成立，\(c\) 的最优小位移再贡献
\((\mu-1)\log n/8+o(\log n)\)。合并后的候选 \(\log n\) 系数是
\[
 -\frac\mu6+\frac{\mu-1}{8}=-\frac{\mu+3}{24}. \tag{15}
\]

根定义确保 \(A_n(q_{\rm soft})=1\)。为满足严格超水平见证，可在最终构造中把阈值替为
\(q_{\rm soft}-n^{-1}\)；由 (5) 这使平均严格大于 1，而得分的对数损失为 \(o(1)\)。

## 5. 证明状态

(2)--(8)、(10)--(12) 的代数和尺度比较已核对。要把 (13) 升格为正式定理，需把 (9) 的 Fourier CDF 反演余项写成完整引理，包括分区、四阶矩常数和外区积分。该引理比 FGL 强，但对固定 \(c=\mu\) 的显式 Gamma--Bernoulli 特征函数是直接可证的；本笔记不把它当作已引用的外部定理。

## 6. 常数项：乘法型固定 \(\mu\) 斜坡渐近

这一节把 (13) 细化到 \(o(1)\)。设
\[
 a=\frac{\kappa_3}{6\mu}=\frac{3-\mu}{6},\qquad
 q_{\rm test}=n\mu-d_n+a z_n^2+b,
 \qquad d=d_n-a z_n^2-b. \tag{16}
\]
由 (1) 的 Mills 关系，\(\phi(z_n)/(\sigma\sqrt n)=(z_n^2+1+o(1))/n\)。将 (9) 在
\(z=d/(\sigma\sqrt n)\) 处 Taylor 展开，三阶 Edgeworth 项精确抵消 \(a z_n^2\) 的正态位移，得到
\[
 F_n(q_{\rm test})
 =\beta_n+(b+a)\frac{\phi(z_n)}{\sigma\sqrt n}+O(n^{-1})+o(n^{-1}). \tag{17}
\]
这里 normal shift 的二次余项和 Edgeworth 项的位移均为 \(o(n^{-1})\)。因此需要的是 (9) 的**绝对**
\(O(n^{-1})\) Fourier 余项，而不是微小尾的相对近似。

未截断指数在同一尺度满足
\[
 e^{-q_{\rm test}/n}\mathbb Ee^{S_n/n}-1
 =\frac dn+\frac{\mu}{2n}+\frac{d^2}{2n^2}+o(n^{-1}). \tag{18}
\]
加权下尾与普通 CDF 的差 \(D_n\) 满足
\[
 0\le D_n\le\frac1n\mathbb E(q_{\rm test}-S_n)_+
 =O\left(\frac{\sqrt L}{n}\right)=o\left(\frac{z_n^2}{n}\right). \tag{19}
\]
固定 \(c=\mu\) 的负倾斜 Chernoff 积分已足以证明 (19)；若改用倾斜中心密度积分还能得到更小量，但不影响常数项。

根条件 \(A_n(q)=1\) 等价于普通 CDF等于 (18) 的右端再加 \(D_n\)。乘以 \(n\)，用
\(d_n^2/n=\mu z_n^2\) 及 \(d_n(a z_n^2+b)/n=o(1)\)，两侧 \(z_n^2\) 系数比较为
\[
 b+a=\frac\mu2-a.
\]
所以
\[
 b=\frac\mu2-2a=\frac{5\mu-6}{6}. \tag{20}
\]
若把 \(b\) 替换为 \(b\pm\varepsilon\)，(17) 改变
\(\varepsilon\phi(z_n)/(\sigma\sqrt n)\asymp\varepsilon L/n\)，而 (18)--(19) 及 Fourier 余项都是
\(o(L/n)\)。由 \(A_n\) 的严格单调性完成符号夹逼，故
\[
 q_{\rm soft}=n\mu-d_n+a z_n^2+\frac{5\mu-6}{6}+o(1). \tag{21}
\]

结合 (3)，固定 \(c=\mu\) 家族的严格斜坡下界具有乘法型渐近
\[
 \frac{e^{q_{\rm soft}}}{I_n}
 \sim \sqrt{\frac2\pi}\,
 \exp\left(\frac{5\mu-6}{6}\right)n^{-1/2}
 \exp\{n\mu-d_n+a z_n^2\}. \tag{22}
\]
最终将 \(q_{\rm soft}\) 下移任意 \(o(1)\) 的正量（例如 \(n^{-1}\)）即可把平均从等号改为严格大于 1，而不改变 (22)。

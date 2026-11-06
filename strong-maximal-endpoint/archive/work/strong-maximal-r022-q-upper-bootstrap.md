# r022：\(Q_n\) sharp-\(\log n\) 上界的 bootstrap 推导

日期：2026-11-06。本文只研究独立比例向量、单调轮廓常数 \(Q_n\)，不涉及真实一般输入常数 \(C_n\)。目标是补齐 r021 的三个分析缺口。下文把需要在正式稿中单列证明的统一 Fourier 命题明确列出；在这些命题成立时，其余推导是闭合的。

设 \(m_i=m(r_i)\)、\(v_i=v(r_i)\)、\(D=\sum(\mu-m_i)\)、\(M=n\mu-D\)、\(V=\sum v_i\)。令
\(s\) 是轮廓预算 \(J\le e^{-s}\) 的对数阈值，\(d=n\mu-s\)，\(y=M-s=d-D\)。

## A. 统一倾斜密度引理（待正式写出）

若粗上界已迫使 \(D=o(\sqrt{n\log n})\)、\(0<y\asymp\sqrt{n\log n}\)，令 \(\lambda>0\) 使负倾斜后的和均值为 \(s\)。则
\[
 \lambda\sim\frac1\sigma\sqrt{\frac{\log n}{n}},\qquad
 g_\lambda(s)\asymp n^{-1/2},\qquad\sup_tg_\lambda(t)=O(n^{-1/2}). \tag{A1}
\]
这里 \(g_\lambda\) 是倾斜和的密度。因每一坐标均含参数为 1 的指数连续部分，其特征函数乘积有统一的指数型 Fourier 外区衰减；内区 Taylor 与 \(V\asymp n\) 给出 (A1)。于是
\[
 f_S(s-u)\le C f_S(s)e^{-\lambda u},\qquad u\ge0, \tag{A2}
\]
并且
\[
 \int_{-\infty}^s e^{t-s}f_S(t)dt=O(f_S(s)),\qquad
 \mathbb E(s-S)_+=O(F_S(s)/\lambda),\qquad
 F_S(s)\asymp f_S(s)/\lambda. \tag{A3}
\]

## B. 轮廓预算 bootstrap

低侧点态预算 \(h(t)\le e^{t-s}\) 与 (A3) 给出
\[
 \mathbb E(1-h(S))_+\ge F_S(s)-O(f_S(s)). \tag{B1}
\]
先用粗高侧反演和 \(\mathbb E|S-M|=O(\sqrt n)\)，可得
\(\mathbb E(h(S)-1)_+=O(\sqrt{\log n/n})\)。又全局密度界给出
\(f_S(s)=O(n^{-1/2})\)，所以 (B1) 先推出
\[
 F_S(s)=O(\rho),\qquad \rho=\sqrt{\frac{\log n}{n}}. \tag{B2}
\]
由 (A3)，此时 \(\mathbb E(s-S)_+=O(1)\)。在中段精确反演 Orlicz 函数，
\[
 h(t)-1\le\frac{t-s}{n}+O\left(\frac{(t-s)^2}{n^2}\right),
\]
加上正矩母函数控制的远尾，得到
\[
 \mathbb E(h(S)-1)_+\le\frac{y}{n}+O\left(\frac{\log n}{n}\right). \tag{B3}
\]
结合 (B1)--(B3) 及 \(f_S(s)=O(\lambda F_S(s))=O(\log n/n)\)（由 (A1)--(A3) 和 \(F=O(\rho)\) 升级），有
\[
 F_S(s)\le\frac yn+O\left(\frac{\log n}{n}\right). \tag{B4}
这一步保留了 \(y=d-D\)，是 sharp 系数的必要条件。

## C. 三角阵列 Edgeworth 与均值缺口压缩

在 (B2) 的参数范围，直接在独立坐标特征函数乘积上作三阶 Taylor/Fourier CDF 反演，可望得到一致绝对展开
\[
 F_S(s)=\Phi(-y/\sqrt V)+\frac{K_3}{6V^{3/2}}(1-y^2/V)\phi(y/\sqrt V)+O(n^{-1}). \tag{C1}
\]
指数连续部分给出外区衰减，所有四阶矩一致有界，故 (C1) 与 r018 的固定参数证明属于同一 Fourier 计算；它不是 FGL 的替代叙述。

将 (C1) 代入 (B4)，并在正态根附近反演，得到
\[
 s\le M-d_{\rm root}(V)+a\log n+o(\log n),
 \qquad a=\frac{3-\mu}{6}. \tag{C2}
\]
由 \(|V-\mu n|\le C\sqrt{nD}\)，
\[
 M-d_{\rm root}(V)
 \le n\mu-d_n-D+C\sqrt{D\log n}+o(\log n). \tag{C3}
\]
与 r019 下界比较即先迫使 \(D=O(\log n)\)。此后 (C1) 中累积量的向量误差
\(K_3-n\mu(3-\mu)=O(\sqrt{nD})\) 为 \(o(\log n)\)。

## D. 用全局负 Laplace 最小化结算方差补偿

令 \(K_c(\lambda)=\log\mathbb Ee^{-\lambda Z_c}\)。r015 的单步全局最小化（最小元满足
\(c_\lambda=\mu-\mu\lambda/2+O(\lambda^2)\)）给出
\[
 \min_c\{K_c(\lambda)-K_\mu(\lambda)\}
 =-\frac{\mu(\mu-1)}8\lambda^3+O(\lambda^4). \tag{D1}
\]
把每个坐标的 Taylor 累积量式相加，并用刚得到的 \(D=O(\log n)\)，得到
\[
 D+\frac\lambda2(V-\mu n)
 \ge-\frac{n\mu(\mu-1)}8\lambda^2+o(\log n)
 =-\frac{\mu-1}{8}\log n+o(\log n). \tag{D2}
\]
将 (D2) 代入 (C2)，即得到
\[
 \log Q_n\le n\mu-d_n+
 \left\{\frac{3-\mu}{6}+\frac{\mu-1}{8}\right\}\log n+o(\log n)
 =n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n). \tag{D3}

## E. 状态与剩余风险

Bootstrap 的代数、(D1) 的符号和 (D2) 的补偿系数均已核对。要把 (D3) 标为定理，必须实际写出并验证 (A1)--(A3)、(C1) 的统一常数，以及“粗上界迫使 \(D=o(\sqrt{n\log n})\)”的起始引理。特别地，不能只以 \(e^tf_S(t)\) 单调性控制 (B1) 的积分，也不能以 FGL 相对余项替代 (C1)。在这些引理未完成前，(D3) 是完整的证明蓝图，而非已证结论。

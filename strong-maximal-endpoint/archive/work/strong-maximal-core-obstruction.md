# 有界近阈值核心：两个严格受限结论

日期：2026-11-06。此记录不声称一般上界；它隔离了任何近阈值稳定性定理必须处理的几何内容。

## 1. 一般有界输入到示性 halo 的归约

设 \(0\le f\le1+\varepsilon\)，其中 \(0<\delta<1\)，并令
\(E_\delta=\{f>1-\delta\}\)。逐点有
\[
 f\le(1-\delta)+(\varepsilon+\delta)1_{E_\delta}.
\]
所以任意平均大于 1 的矩形满足
\[
 { |R\cap E_\delta|\over|R|}>{\delta\over\varepsilon+\delta}.
\]
记 \(T_n(\alpha)=\sup_{0<|E|<\infty}|\{M_S1_E>\alpha\}|/|E|\)，则
\[
 |\{M_Sf>1\}|
 \le {T_n(\delta/(\varepsilon+\delta))\over1-\delta}
 \int\Phi_n(f).\tag{A}
\]
最后一步只用 \(\Phi_n(f)\ge f>1-\delta\) 于 \(E_\delta\)。

这是严格的、但通常偏松的归约。比如选 \(\delta=\sqrt\varepsilon\) 后，示性 halo 阈值为
\(1/(1+\sqrt\varepsilon)\)。因此，若要由有界核心单独得到维数一致的尖锐上界，至少需要控制随维数变化的近 1 Tauberian 常数；单纯的 \(L^p\) 常数没有包含这部分矩形几何。

## 2. 常幅单盒不能产生主指数底数

令 \(Q=[0,1]^n\)，\(f=(1+\varepsilon)1_Q\)，并设
\(A=\log(1+\varepsilon)\)。一维非中心极大函数为
\[
 m(x)=M_1 1_{[0,1]}(x)=
 \begin{cases}
 1,&x\in[0,1],\\
 x^{-1},&x>1,\\
 (1-x)^{-1},&x<0.
 \end{cases}
\]
矩形独立性给出精确恒等式
\(M_S1_Q(x)=\prod_{i=1}^nm(x_i)\)。令 \(y_i=-\log m(x_i)\)。每一坐标在 \(y=0\) 有质量 1，并在 \(y>0\) 有密度 \(2e^y\)。故 halo 体积精确为
\[
 V_n(A)=1+\sum_{k=1}^n{n\choose k}{2^k\over(k-1)!}
 \int_0^A e^s s^{k-1}\,ds.\tag{B}
\]
这直接保留了同一矩形的所有坐标，不是坐标复合上界。

由 \(\int_0^Ae^ss^{k-1}ds\le e^AA^k/k\) 及
\({n\choose k}\le n^k/k!\)，
\[
 V_n(A)\le\sum_{k\ge0}{(2e^AnA)^k\over(k!)^2}
 \le\exp\{2\sqrt{2e^AnA}\}.\tag{C}
\]
因此常幅单盒的归一化比满足
\[
 \log { |\{M_Sf>1\}|\over\int\Phi_n(f)}
 \le2\sqrt{2e^AnA}-A-(n-1)\log(1+A).
\]
特别地，若 \(A\asymp\sqrt{\log n/n}\)，该式为
\(-\Theta(\sqrt{n\log n})+O(n^{1/4}(\log n)^{1/4})\)，远小于 \(n\log B\)。
所以“一个近阈值单盒及其直接 halo”不能是主指数底数 \(B^n\) 的机制；任何一般核心上界的近饱和者必须利用多尺度或相关矩形结构。

## 本机数值核对

MATLAB R2026a 对 \(n=32\)、\(A=\sqrt{\log32/32}\) 直接计算 (B)，得到
\(\log V_n=7.17048638\)，而 (C) 的右端对数为 \(10.82040872\)，并且
\(\log\Phi_n(e^A)=9.14856795\)。这仅核对了闭式公式的量级；(A)--(C) 本身是解析不等式。

# r014：隐式 Gaussian 根的上下界、严格性与记号边界

日期：2026-11-06。本笔记只整理隐式根路线。它区分三种常数：

\[
 Q_n\le C_n\le R_n. \tag{1}
\]

其中 \(Q_n\) 是已构造的单层双曲指数斜坡下界，\(C_n\) 是题目中的真实强极大算子常数，\(R_n\) 是经全局 rising-sun/重排归约后得到的自适应一维复合轮廓上界。本文的上界结论是关于 \(R_n\)；只有在该独立归约已严格证明时，才可由 (1) 转移给 \(C_n\)。本文绝不把 \(R_n\) 的上界误写成 \(Q_n\) 的等式。

设 \(L=\log n\)，\(\mu=\log B\)，\(\sigma^2=\mu\)，并以唯一正根定义
\[
 \overline\Phi(z_n)=\frac{\sigma z_n}{\sqrt n},\qquad
 d_n=\sigma\sqrt n\,z_n. \tag{2}
\]

## 1. 指数斜坡的严格下界

固定 \(r=B^{-1}\)。令 \(T=\sum_{i=1}^nZ_r^{(i)}\)，则 \(\mathbb ET=n\mu\)、
\(\operatorname{Var}T=n\mu\)。取
\[
 h_n(u)=e^{u/n}\mathbf1_{\{u\ge0\}},\qquad q=n\mu-d.
\]
在双曲层边界的锚定矩形平均是
\[
 A_n(d)=\mathbb E\left[e^{(T-q)/n}\mathbf1_{\{T\ge q\}}\right].
\]
未截断指数的 Jensen 不等式和 \(e^{(T-q)/n}\le1\) 于 \(T<q\) 给出严格的可用下界
\[
 A_n(d)\ge e^{d/n}-\Pr(T<q). \tag{3}
\]
因此只要右端大于 1，接收集合的内部有严格平均余量。

令 \(d=d_n+C_0L\)，其中 \(C_0>0\) 稍后固定。Fan--Grama--Liu 相对中偏差在此独立同分布情形给出
\[
 \Pr(T<q)=\overline\Phi\left(z_n+\frac{C_0L}{\sigma\sqrt n}\right)
 \left[1+O\left(\frac{L^{3/2}}{\sqrt n}\right)\right]. \tag{4}
\]
由 (2)、Mills 比 \(\phi(z_n)=\sigma z_n^2/\sqrt n\,(1+O(L^{-1}))\)，以及
\(z_n^2\asymp L\)，把 (4) 在 \(z_n\) 一阶展开可知，尾概率因 \(C_0L\) 的位移减少
\[
 \frac{C_0L}{\sigma\sqrt n}\phi(z_n)(1+o(1))
 =\frac{C_0L^2}{n}(1+o(1)). \tag{5}
\]
而 (4) 的绝对相对误差至多 \(O(L^2/n)\)。故取 \(C_0\) 大于该误差常数后，充分大的 \(n\) 有
\[
 \Pr(T<q)<\frac{d_n}{n}<e^{d/n}-1. \tag{6}
\]
结合 (3)，这给出所需的**严格** \(A_n(d)>1\)，并非仅边界等号。

对每个固定 \(n\)，令 \(L_0\to\infty\) 并取
\(f_{L_0}(x)=h_n(S(x)-L_0)\mathbf1_{[-1,1]^n}(x)\)。指数斜坡的模数极限是
\[
 I_n=\int_0^\infty e^{-(1-1/n)v}(1+v/n)^{n-1}\,dv\le2n. \tag{7}
\]
严格余量使 \(E_{L_0-q}\) 的内部成为 \(\{M_Sf_{L_0}>1\}\) 的见证；固定 \(n\) 后的截断极限于是给出
\[
 \log Q_n\ge q-\log I_n
 \ge n\mu-d_n-O(L). \tag{8}
\]

## 2. 自适应轮廓的隐式根上界

本节假设全局归约已将问题化为 \(R_n\)，并使用已核验的停止策略 Fan--Grama--Liu 相对中偏差结论。对全部自适应策略、全部
\(0\le d\le2\sigma\sqrt{nL}\)，它给出
\[
 \Pr\{S_n\le n\mu-d\}
 \ge(1-O(\kappa_n))\overline\Phi\left(\frac d{\sigma\sqrt n}\right)-n^{-3},
 \quad
 \kappa_n=O(n^{-1/8}L^{5/8}). \tag{9}
\]

设某可行轮廓的模数为 \(J=e^{-s}\)，并令 \(d=n\mu-s\)。轮廓单调性给出
\(\Phi_n(h(t))\le e^{t-s}\)。在 \(0\le t-s\le n^{3/4}\)，反演此不等式得到
\[
 h(t)-1\le(1+O(n^{-1/4}))(t-s)/n. \tag{10}
\]
为控制中段期望，令 \(M=S_n-\sum_i m(r_i)\)。对每个 \(r\) 和 \(\lambda\ge0\)，
\[
 \log\mathbb E e^{-\lambda(Z_r-m(r))}\le\frac{V_{\max}\lambda^2}{2},
 \qquad V_{\max}=1+4/e^2. \tag{11}
\]
这是因为负倾斜累积量的二阶导数递减，且
\(v(r)\le1+r(\log r)^2\le V_{\max}\)。迭代 (11) 并积分 Chernoff 尾得到
\[
 \mathbb E(-d-M)_+
 \le\frac{nV_{\max}}d\exp\left(-\frac{d^2}{2nV_{\max}}\right). \tag{12}
\]
由于 \(S_n-s\le d+M\)，(12) 将中段超额写成
\[
 \mathbb E(S_n-s)_+
 \le d+\mathbb E(-d-M)_+=d(1+o(\kappa_n)) \tag{13}
\]
于临界窗口。这里 \(\mu/(2V_{\max})>1/8\)，故 (12) 相对 \(d\) 的误差小于 \(\kappa_n\)。

为明确 (14) 中未写出的远尾项：在 \(t>s+n^{3/4}\) 处分割
\(h(t)<e^{4t/n}\) 与其补集。对前一部分，取 \(\theta=cn^{-1/4}\) 并用一致展开
\(\log\mathbb Ee^{\eta Z_r}\le\eta\mu+C\eta^2\)，可得贡献 \(e^{-c_1\sqrt n}\)。对后一部分，任意自适应策略的和密度都满足
\[
 f_{S_n}(t)\le\frac{2^ne^{-t}t^{n-1}}{(n-1)!},
\]
因为每一步条件密度不超过 \(2e^{-z}\)。又
\(\Phi_n(h)\ge h(4t/n)^{n-1}\)，故这一部分至多
\[
 J\frac{2^n(n/4)^{n-1}}{(n-1)!}=e^{-c_2n}
\]
于临界窗口。两者均为 \(o(d\kappa_n/n)\)。

把第 10 节的低值带取为 \(t\le s-2L\)，并用 (10)、(13) 及 \(n^{3/4}\) 外的低/高振幅远尾估计，得到可行性强制
\[
 \Pr\{S_n\le s-2L\}
 \le\frac dn(1+O(\kappa_n)). \tag{14}
\]
将 (9) 施于 \(d+2L\)（该平移对相对误差更小），可得
\[
 \overline\Phi\left(\frac{d+2L}{\sigma\sqrt n}\right)
 \le\frac dn(1+O(\kappa_n)). \tag{15}
\]

令 \(H(d)=\overline\Phi(d/(\sigma\sqrt n))/(d/n)\)。由 (2)，\(H(d_n)=1\)，且
\[
 (\log H)'(d_n)=-\frac{z_n}{\sigma\sqrt n}(1+o(1))
 =-\asymp\sqrt{\frac Ln}. \tag{16}
\]
因此 (15) 强制
\[
 d\ge d_n-O\left(\kappa_n\sqrt{\frac nL}\right)-O(L)
 =d_n-O(n^{3/8}L^{1/8}), \tag{17}
\]
从而
\[
 \log R_n\le n\mu-d_n+O(n^{3/8}L^{1/8}). \tag{18}
\]

结合 (8)、(18) 及 (1)，并在全局归约已证的前提下，最终可安全表述为
\[
 n\mu-d_n-O(L)\le\log Q_n\le\log C_n\le\log R_n
 \le n\mu-d_n+O(n^{3/8}L^{1/8}). \tag{19}
\]

## 3. Gaussian 根的显式展开

令 \(a=\log(\sigma\sqrt{2\pi})\)。Mills 展开代入 (2) 后，令 \(w=z_n^2\)，得到
\[
 w+2\log w=L-2a+2\log\left(1-w^{-1}+3w^{-2}+O(w^{-3})\right).
\]
迭代解得
令 \(H=\log L+a\)。Wolfram 的符号级数计算（以下系数也可由逐项代回上式核验）给出
\[
 \begin{aligned}
 z_n={}&\sqrt L-\frac H{\sqrt L}
 +\frac{-H^2/2+2H-1}{L^{3/2}}\\
 &+\frac{-H^3/2+4H^2-7H+9/2}{L^{5/2}}
 +O\left(\frac{H^4}{L^{7/2}}\right). \tag{20}
 \end{aligned}
\]
等价地，隐式根中心项满足
\[
 \begin{aligned}
 n\mu-d_n={}&n\mu-\sigma\sqrt{nL}
 +\sigma(\log L+a)\sqrt{\frac nL}\\
 &+\sigma\left\{\frac12H^2-2H+1\right\}
 \frac{\sqrt n}{L^{3/2}}\\
 &+\sigma\left\{\frac12H^3-4H^2+7H-\frac92\right\}
 \frac{\sqrt n}{L^{5/2}}
 +O\left(\frac{\sqrt nH^4}{L^{7/2}}\right).
 \end{aligned} \tag{21}
\]

实际运行的 Wolfram 15 符号输入为

```wl
w = 1/t + b0 + b1 t + b2 t^2;
eq = Series[w + 2 Log[w] -
  (1/t - 2 a + 2 Log[1 - 1/w + 3/w^2 - 15/w^3]),
  {t, 0, 3}] // Normal // Expand;
sol = Solve[Thread[Coefficient[eq, t, #] == 0 & /@ Range[0, 2]],
  {b0, b1, b2}];
Series[Sqrt[w /. sol[[1]]], {t, 0, 3}]
```

它返回 \(b_0=-2H\)、\(b_1=4H-2\)、\(b_2=(2H-3)^2\)，并给出 (20) 的
\(L^{-5/2}\) 系数。先前只写到 \(L^{-3/2}\) 时的余项必须是
\(O(H^3/L^{5/2})\)，不能写成平方；上式通过显式加入该项避免了这一歧义。

Wolfram 15.0.0 以 80 位精度直接解 \(\log(\operatorname{Erfc}(z/\sqrt2)/2)=\log\sigma+\log z-L/2\)，并与 (20) 比较：

| \(L\) | 精确 \(z_n\) | (20) 近似 |
|---:|---:|---:|
| 100 | 9.42969333311276 | 9.42965605837283 |
| 1,000 | 31.3708609726983 | 31.3708625392463 |
| 10,000 | 99.8974458448172 | 99.8974458633889 |

数值只核验代数展开；(20)--(21) 的严格余项来自 Mills 级数的标准可控余项。

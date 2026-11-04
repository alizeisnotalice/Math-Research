# 第十五轮独立审查：幂权外层端点

## 审查结论

一维候选成立。对

\[
d\mu_\tau(x)=|x|^{\tau-1}\,dx,\qquad0<\tau<1,
\]

外层 Lebesgue 平均和 \(\mu_\tau\) 光晕满足

\[
\lim_{t\downarrow0}\frac{C^{\rm out}_{c,w_\tau}(1-t)-1}{t^\tau}=1.
\tag{15.1}
\]

并且在

\[
0<\tau<\frac{5-\sqrt{17}}2
\]

有精确二阶式

\[
C^{\rm out}_{c,w_\tau}(1-t)
=1+t^\tau+t^{2\tau}+o(t^{2\tau}).
\tag{15.2}
\]

追加的 \(n\) 维 \(\ell^\infty\)-径向幂权
\(d\mu(x)=\|x\|_\infty^{\tau-n}dx\)、\(0<\tau<n\)，也可安全保留：

\[
\lim_{t\downarrow0}
\frac{C^{\rm out}_{c,w}(1-t)-1}{t^{\tau/n}}=1.
\tag{15.3}
\]

这不自动推广到 Euclidean 径向权；下述重排和孔构造都使用 \(\ell^\infty\) 球即方体。

## 一维小见证

令 \(H=\{M_c\mathbf1_E>1-t\}\)，其中
\(0<|E|<\infty\)、\(0<\mu_\tau(E)<\infty\)。为消除可能的重叠，取

\[
\begin{aligned}
H_s&=\{x\in H\setminus E:\ \exists\text{ 严格见证 }Q(x,r),\
r\le C|x|\},\\
H_l&=(H\setminus E)\setminus H_s .
\end{aligned}
\tag{15.4}
\]

设 \(A_k=\{2^k\le|x|<2^{k+1}\}\)、\(L=2(1+C)\)、
\(D_k=[-L2^k,L2^k]\)。若 \(x\in A_k\) 的见证满足 \(r\le C|x|\)，
则 \(Q(x,r)\subset D_k\)。故一维无权中心弱型/扰动界给出

\[
|(H_s\setminus E)\cap A_k|
\le\frac{2t}{1-t}|E\cap D_k|.
\tag{15.5}
\]

在 \(A_k\) 上 \(w_\tau\le2^{k(\tau-1)}\)。对固定 \(y\ne0\)，

\[
\sum_{k:y\in D_k}2^{k(\tau-1)}
\le\frac{L^{1-\tau}|y|^{\tau-1}}{1-2^{\tau-1}}.
\]

Tonelli 因而给

\[
\frac{\mu_\tau(H_s)}{\mu_\tau(E)}
\le
\frac{2L^{1-\tau}}{1-2^{\tau-1}}\frac{t}{1-t}.
\tag{15.6}
\]

取 \(C=t^{-a}\)、\(0<a<1\)，这为
\(O(t^{1-a(1-\tau)})=o(t^\tau)\)。零点没有短见证，且不影响此估计。

严格集的可测性也可闭合：固定半径的 Lebesgue 平均连续；严格见证可用更小的
有理半径近似。因此 \(H_s\) 可写成可数个开条件的并。

## 一维长见证

对 \(x\in H_l\)，每个严格见证均有 \(r>C|x|\)。令 \(R\) 为所有这些
**实际见证半径**的上确界（空集时 \(R=0\)）。记

\[
A_C=(1+C^{-1})^\tau+(1-C^{-1})^\tau.
\]

因为 \(r>|x|\)，

\[
\mu_\tau(Q(x,r))\ge\frac{r^\tau A_C}{\tau}.
\]

而 \(|Q(x,r)\setminus E|<2tr\)，递减重排给

\[
\mu_\tau(Q(x,r)\setminus E)\le\frac{2(tr)^\tau}{\tau}.
\]

因此

\[
\mu_\tau(E)\ge\frac{r^\tau(A_C-2t^\tau)}{\tau}.
\tag{15.7}
\]

当 \(t\) 小到 \(A_C>2t^\tau\) 时，\(R<\infty\)。取一列实际长见证
\((x_j,r_j)\) 使 \(r_j\uparrow R\)。由于 \(C>2\)，充分大的 \(j\) 满足

\[
r_j-|x_j|>r_j(1-C^{-1})>R/C,
\]

故 \(Q(x_j,r_j)\) 包住 \([-R/C,R/C]\)，而 \(H_l\) 包含于该核心且与
\(E\) 不交。于是 \(|H_l|\le2tR\)，并在 (15.7) 中令 \(r_j\uparrow R\)，得到

\[
\frac{\mu_\tau(H_l)}{\mu_\tau(E)}
\le\frac{2t^\tau}{A_C-2t^\tau}.
\tag{15.8}
\]

这里若把 \(R\) 错换成 (15.7) 给出的粗上界，核心未必能被一个实际见证方体包含；
必须采用上述实际半径上确界。

由 (15.6)、(15.8) 得一阶上界。取 \(C=t^{-a}\) 后
\(A_C\to2\)，故其主项为 \(t^\tau\)。

## 一维下界与二阶项

取

\[
E_t=[-1-t,-t]\cup[t,1+t].
\]

对 \(0\le x<t\)，令 \(r=1+t-x\)。该方体中的 \(E_t\) 长度为
\(2-2x\)，故平均

\[
\frac{1-x}{1+t-x}>1-t
\quad\Longleftrightarrow\quad t(t-x)>0.
\]

由对称性，\((-t,t)\) 模端点包含于严格光晕，因而

\[
C^{\rm out}_{c,w_\tau}(1-t)-1
\ge\frac{t^\tau}{(1+t)^\tau-t^\tau}.
\tag{15.9}
\]

实际 Wolfram Language 15 符号计算返回

\[
\frac{A_C}{2}
=1-\frac{\tau(1-\tau)}{2C^2}+O(C^{-4}).
\tag{15.10}
\]

为使 (15.8) 的 \(C^{-2}\) 误差在 \(t^{2\tau}\) 后消失，需
\(2a>\tau\)；为使 (15.6) 为 \(o(t^{2\tau})\)，需

\[
a<\frac{1-2\tau}{1-\tau}.
\]

这与 \(a>\tau/2\) 同时可取，当且仅当
\(\tau<(5-\sqrt{17})/2\)。此时 (15.8) 和 (15.9) 都展开为
\(t^\tau+t^{2\tau}+o(t^{2\tau})\)，证明 (15.2)。本轮实际 Wolfram 输出还给出
\((5-\sqrt{17})/2=0.4384471871911697\ldots\)，仅作代数交叉核验。

## \(\ell^\infty\) 径向的 \(n\) 维推广

令 \(0<\tau<n\)、\(\rho(x)=\|x\|_\infty\)、\(d\mu=\rho^{\tau-n}dx\)。写

\[
\mu(Q(0,r))=K_{\tau,n}r^\tau,\qquad
K_{\tau,n}=\frac{2^n n}{\tau}.
\tag{15.11}
\]

在环带 \(2^k\le\rho(x)<2^{k+1}\) 上重复短见证论证，使用 \(n\) 维无权中心
弱型界，得到

\[
\frac{\mu(H_s)}{\mu(E)}
=O_{\tau,n}\bigl(tC^{\,n-\tau}\bigr).
\tag{15.12}
\]

对长见证 \(r>C\rho(x)\)，尺度不变性和
\(\rho(x)/r\le C^{-1}\) 给

\[
\mu(Q(x,r))\ge K_{\tau,n}r^\tau A_C,\qquad A_C\to1.
\]

若 \(|Q(x,r)\setminus E|<t(2r)^n\)，\(\ell^\infty\) 的递减重排给

\[
\mu(Q(x,r)\setminus E)\le K_{\tau,n}r^\tau t^{\tau/n}.
\]

以实际长见证半径上确界和同样的核心包含论证，遂有

\[
\frac{\mu(H_l)}{\mu(E)}
\le\frac{t^{\tau/n}}{A_C-t^{\tau/n}}.
\tag{15.13}
\]

取 \(C=t^{-a}\)、\(0<a<1/n\)，(15.12) 为 \(o(t^{\tau/n})\)，而
(15.13) 的首项为 \(t^{\tau/n}\)。

下界取 \(E=Q(0,R)\setminus Q(0,h)\)，

\[
\frac hR=\frac{t^{1/n}}{1+t^{1/n}}.
\]

当 \(t\) 足够小时，对孔内 \(x\) 取 \(r=R-\rho(x)\)。则 \(Q(x,r)\subset Q(0,R)\)
并包含整个内孔（充分条件为 \(h<R/3\)），且其缺失比例严格小于 \(t\)。
故孔内模零属于严格光晕，给

\[
C^{\rm out}_{c,w}(1-t)-1
\ge
\frac{(h/R)^\tau}{1-(h/R)^\tau}
\sim t^{\tau/n}.
\tag{15.14}
\]

(15.12)--(15.14) 证明 (15.3)。

## 审查意见

一维一阶、指定 \(\tau\) 范围内的二阶、以及 \(\ell^\infty\) 径向高维一阶结论都可
保留。主文须明确长类取 (15.4) 的不重叠定义，并让 \(R\) 是实际长见证半径上确界；
高维结论必须保留 \(\ell^\infty\) 径向权与方体重排，不能无证明替换为 Euclidean
径向权。

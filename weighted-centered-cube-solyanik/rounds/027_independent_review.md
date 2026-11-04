# 第 27 轮独立审查：三尺度上界与任意范数角向权极限

审查时scratch式号R27对应当前paper/small_homogeneity_weak.tex中的SW；本文保留审查时的一维分解以供独立核对，当前正文只保留统一一般证明。

## 裁决

一维三尺度论证成立，并给出
\[
\lim_{\tau\downarrow0}\tau W_\tau=1
\]
（这里 \(W_\tau\) 是 Lebesgue 平均、\(|x|^{\tau-1}dx\) 测量的外层弱型范数）。
推广到任意有限维范数球、固定本质有界且正的角向剖面
\(h_\tau=\rho^{\tau-n}a(\theta)\) 也成立：
\[
\lim_{\tau\downarrow0}\tau W_{\tau,K,a}
=\frac{n\,\mathbb E_{\sigma_K}a}{\operatorname*{ess\,inf}_{\sigma_K}a}.
\tag{R27.1}
\]
这些是完整弱型范数结论，和第 26 轮有限原子下界不同；其依据是对任意输入的
统一上界与独立的平滑原子下界。

## 一维上界

令 \(A=\int f(x)|x|^{\tau-1}dx<\infty\)，可取 \(f\ge0\)，并固定
\(0<\eta<1\)。对每个严格 Lebesgue 见证 \(Q=x+[-r,r]\) 分类为
\[
r<|x|/2,\qquad |x|/2\le r\le|x|/\eta,\qquad r>|x|/\eta.
\]
三类的存在性并集覆盖严格 halo；可用有理半径与固定半径平均的连续性取得可测版本。
\(L^1(\mu_\tau)\) 输入自动局部 Lebesgue 可积，因为 \(w_\tau\) 在每个有界区间
有正的本质下界。

局部类不穿过原点，区间上的权本质比最多为 \(3^{1-\tau}\le3\)。
Lebesgue 平均大于一遂给内禀平均大于 \(1/3\)。中心区间族对任意 Radon 测度的
弱型常数为 \(2\)，所以
\[
\mu_\tau(H_{\rm loc})\le6A.
\]
中间类令 \(C=1+\eta^{-1}\) 及
\[
T(x)=|x|^{-1}\int_{|y|\le C|x|}f(y)\,dy.
\]
见证给 \(T(x)>1\)，而 Tonelli 精确给
\[
\int_{\mathbb R}T(x)|x|^{\tau-1}dx
=\frac{2C^{1-\tau}}{1-\tau}A.
\]
本机 Wolfram 已直接核验该内层积分。中心类中
\[
A>2r(r+|x|)^{\tau-1}
\ge2r^\tau(1+\eta)^{\tau-1},
\]
从而严格 halo 被包含在一个以零为中心、其 \(\mu_\tau\) 质量至多
\(\eta^\tau(1+\eta)^{1-\tau}A/\tau\) 的区间中。故对每个输入和每个阈值
（由齐次性先取阈值一）统一有
\[
W_\tau\le
6+\frac{2C^{1-\tau}}{1-\tau}
+\frac{\eta^\tau(1+\eta)^{1-\tau}}{\tau}.\tag{R27.2}
\]
先固定 \(\eta\) 令 \(\tau\downarrow0\)，再令 \(\eta\downarrow0\)，得
\(\limsup\tau W_\tau\le1\)。第 25 轮的单原子下界给相反的 liminf，故极限成立。
这里没有把任何有限原子搜索当作上界。

## 任意范数与角向剖面

设 \(V=|K|\)，cone probability 为 \(\sigma_K\)，
\(0<a_*\le a\le a^*<\infty\) a.e.，且 \(E_a=\mathbb E_{\sigma_K}a\)。
对 local/middle/central 三类，完全相同的几何与 Tonelli 计算给
\[
\begin{aligned}
\mu(H_{\rm loc})&\le B_K(a^*/a_*)3^{n-\tau}A,\\
\mu(H_{\rm mid})&\le
\frac{n2^nE_a}{a_*(n-\tau)}C^{n-\tau}A,\\
\mu(H_{\rm cen})&\le
\frac{nE_a}{a_*\tau}\eta^\tau(1+\eta)^{n-\tau}A.
\end{aligned}\tag{R27.3}
\]
其中 \(B_K\) 是范数球族的测度一致中心弱型常数（已有 BIP/Aldaz 路线保证其有限）。
中间项的关键恒等式是
\[
\int_{\rho(x)\ge\rho(y)/C}\rho(x)^{\tau-2n}\,dx
=\frac{nV}{n-\tau}C^{n-\tau}\rho(y)^{\tau-n}.
\]
这里若已经把极坐标的 Jacobian 合并，才是径向积分
\(nVE_a\int_{\rho(y)/C}^\infty s^{\tau-n-1}ds\)。本机 Wolfram 15.0
分别在 \(n=2,3\) 给出 \(2\pi C^{2-\tau}y^{\tau-2}/(2-\tau)\) 与
\(4\pi C^{3-\tau}y^{\tau-3}/(3-\tau)\)，与此式相符。将 (R27.3) 乘以 \(\tau\)、
先固定 \(\eta\) 再令 \(\eta\downarrow0\)，给出上界
\(nE_a/a_*\)。

下界不需要角向连续性。对任意 \(\varepsilon>0\)，低角集合
\(\{a\le a_*+\varepsilon\}\) 在 \(1<\rho<2\) 的锥区域有正 Lebesgue 测度；
在其中取 \(h_\tau\) 的非零 Lebesgue 点 \(z\)，设 \(s=\rho(z)\)。以 Lebesgue
质量 \(V(rs)^n\) 的窄平滑原子置于 \(z\)，其 strict halo 包含
\((r-1)s\,\operatorname{int}K\)，而成本趋于
\(Vr^ns^\tau a(z)\)。于是
\[
W_{\tau,K,a}\ge
\frac{nE_a}{\tau(a_*+\varepsilon)}
\frac{(r-1)^\tau}{r^n}.
\]
右端在 \(r=n/(n-\tau)\) 最大，故
\[
\tau W_{\tau,K,a}\ge
\frac{nE_a}{a_*+\varepsilon}
\frac{\tau^\tau(n-\tau)^{n-\tau}}{n^n}.
\]
先 \(\tau\downarrow0\)、再 \(\varepsilon\downarrow0\)，得到匹配下界。
这说明 (R27.1) 对任意固定 \(K,a\) 成立。

需要在正式稿中明确：\(a_*\) 是 cone measure 下的本质下确界，选取的是低角
锥区域的 Lebesgue 点而非假设面上最小值或角向连续性。


## 可选一致收敛速率

设 \(L=nE_a/a_*\)，且 \(a^*/a_*\le R\)。选
\(\eta=\tau^{1/(n+1)}\)、\(0<\tau<\min\{1,n/2\}\)。由 (R27.3) 除以 \(L\) 后，
中间项为
\[
\frac{2^n\tau}{n-\tau}(1+\eta^{-1})^{n-\tau}=O_n(\eta),
\]
中央项为
\[
\eta^\tau(1+\eta)^{n-\tau}=1+O_n(\eta).
\]
local 项也为 \(O_{n,R}(\tau)\)，故可吸收。注意若 \(B_K\) 指原始测度一致
弱型常数，完整 local 项是
\[
B_K\,3^{n-\tau}R\,\tau/n,
\]
而非省略 \(3^{n-\tau}\) 的式子；用 \(B_K\le3^n\) 后仍是统一
\(O_{n,R}(\tau)\)。

优化的下界除以 \(L\) 为
\[
L_0(\tau)=\frac{\tau^\tau(n-\tau)^{n-\tau}}{n^n}.
\]
其对数展开为
\[
\log L_0=\tau\log\tau-\tau(\log n+1)+O_n(\tau^2),
\]
所以 \(1-L_0=O_n(\tau|\log\tau|)=O_n(\eta)\)。从而
\[
\left|\frac{\tau W_{\tau,K,a}}{L}-1\right|
\le C_{n,R}\tau^{1/(n+1)}
\]
可对全部范数球与所有允许随 \(\tau\) 变化、但 contrast 至多 \(R\) 的角向剖面
一致成立。低角 Lebesgue 点构造逐个 \(\tau\) 取 \(\varepsilon\downarrow0\)，不需要
角向族的统一代表选择。

## 最终源文件复核（2026-11-05）

逐行复核当前稿 `paper/small_homogeneity_weak.tex` 与
`paper/single_atom_lower_bound.tex` 后，未发现需要修改正文的数学问题。

- 正文 (SW.6) 已采用正确的 \(dx\) 幂次 \(\rho(x)^{\tau-2n}\)，并在下一行明确了
  极坐标后的 \(s^{\tau-n-1}\)；本记录旧版的 \(dx\) 幂次笔误现已更正，完全不影响
  (SW.1)--(SW.9)。
- (SW.8) 的显式常数正确：以 \(\eta=\tau^{1/(n+1)}\)，局部、中尺度、中央三项分别
  被 \(9^nR\eta/n\)、\(2^{2n+1}\eta/n\)、
  \(1+(2^n-1)\eta\) 控制。故
  \(C_{n,R}=(9^nR+2^{2n+1})/n+2^n-1\) 正是所列上界。低角集合与有理半径的写法也
  消除了可测分类和“面上取得最小角”的隐含假设。
- (SW.9) 的极限顺序与 (OW.2) 一致：先固定 \(\tau\) 取 \(\lambda\downarrow\lambda_c\)，
  再令 \(\tau\downarrow0\)。在一维、等指数双线性情形
  \(\rho_c\lambda_c=1/2\)，故留数 \(\sim1/(2\tau)\)，符号正确。
- 单原子严格 halo 的端点被排除，平滑的 Fatou 回迁保留严格余量；(TA.3) 的三条有理
  不等式经本机 Python `Fraction` 复验。该段只说单原子族极大，并未误称为 \(W_\tau\)
  的全局精确值；随后双原子定理的范围也仍止于“严格改进单原子族”。

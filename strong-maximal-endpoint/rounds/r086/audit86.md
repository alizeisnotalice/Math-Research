# r086 任意替代矩形结论独立审计

结论：[free_rectangles86.md](free_rectangles86.md) 在其明确范围内通过审计。它把 r085 的 canonical-witness 熵障碍扩展到每个相对方向各有一个任意可测真 source piece、receiver 可重选该方向任意正边长的情形；不涉及一般 endpoint、一般残差分解或 r083 固定物理角点 column 架构。

## 端点统一化与异常边资格

任意真矩形可写为 `Y_i=x_iU_i`，其中 `U_i` 独立、均匀于端点 `1,d_i` 的区间。令 `W_i=(-log|U_i|)_+`。即使边越过零、短于到零的边或远超 canonical 尺度，cube cutoff 与绝对坐标单调性仍给
\[
F(Y)\le h_N\!\left(t(x)+\sum_iW_i\right).
\]
因 `h_N` 单调且 `Phi_N>=id`，
\[
h_N(s)\le e^s j_N.
\]

端点分类、均值和尾部均正确：
\[
m(-r)=1-\frac{r\log r}{1+r}\ (0\le r\le1),\qquad
m(-r)=\frac2{1+r}\ (r>1),
\]
其余 non-crossing 情形的均值至多一；前式在 `r_*=e^{-mu}=mu-1` 唯一取最大 `mu`。所有端点均满足
\[
\Pr(W_i>s)\le2e^{-s},\qquad
\log\mathbb Ee^{bW_i}\le bm_i+\frac{2b^2}{(1-b)^3}.
\]
这里不需要错误地要求短区间有统一密度界；尾界已足够。

因此若 `sum m_i<=(mu-v)N`，Chernoff、`t=o(N)` 与 `j_N=exp(-mu N+o(N))` 使该矩形的 `F` 平均严格小于一，统一排除真成功。固定 `a` 后，每个非 `d_i=-r_i, |r_i-r_*|<=a` 的边有正的均值缺口，故异常数为 `o(N)`。

新增的定量资格估计也正确。令
\[
q=-\log j_N-T_N-\log4.
\]
在 `sum W_i<=q` 上有 `F(Y)<=1/4`，故严格成功强迫
`P(sum W_i>q)>3/8`。代入矩估计并取 `b=N^{-1/2}` 得
\[
0\le\mu N-\sum_i m_i<B_N,
\]
其中原稿所列的 `B_N` 正是 Chernoff 代数结果。以
\[
g_a=\min\{\mu-1,\mu-m(r_*-a),\mu-m(r_*+a)\}>0
\]
得到 `|I|<=(B_N)_+/g_a`。由 r085 的标量精度，
`B_N=O(sqrt(N log N))`，所以异常数是统一的
`O_a(sqrt(N log N))`；没有假称固定维数的数值起效点。

## 异常坐标收费与标签

异常坐标用 noncentered one-dimensional maximal composition `P_I` 支配，不必枚举尺度。单侧 Hardy rearrangement 给单侧 `L^2` 范数至多 2，noncentered 范数至多 4，故
\[
\sum_\sigma\int(P_IG_\sigma)^2
\le16^{|I|}\sum_\sigma\int G_\sigma^2
\le2\,16^{|I|}\sum_\sigma J_N(G_\sigma).
\]
最后一步只用 cap two 和实际 Orlicz fee。正则边扩大到 `[-(r_*+a)x_i,x_i]` 的体积损失至多 `e^{ND_a}`，故原成功推出 `A_{I,d}>e^{-ND_a}`。

`d` 固定异常坐标的绝对方向，正则标签为 `-sign x_j`。尽管 `P_I` 在异常坐标中采用完整 noncentered operator，它只提供上界；余下 `G_{sigma(x;I,d)}` 标签在后续符号变换中仍保持。使用有理长度、局部 `L^1` 连续性、Fubini 和非零线性变换处理 supremum/零集是充分的。对固定 `a`，pair 数为 `exp(o(N))`。

## 模式熵、Jacobian 与真成功集

在正则坐标，`p_a=(r_*+a)/(1+r_*+a)<1/2`，且 `H_2(p_a)>=H_2(p)`。固定 `U` 后
`(X_I,X_J)->(Y_I,Y_J)=(X_I,X_JU_J)` 的密度至多 `e^{S_J}/|E|`；给定 `Y,d`，正则 sign pattern 到 `sigma(X;I,d)` 是双射。模式与 `S_J` 无需独立，联合事件直接以模式边缘概率上界即可。因此
\[
|E|^{-1}\int_EA_{I,d}^2
\le O(e^{-cN})+
\frac{e^{(\mu-H_2(p)+3\delta)N}}{|E|}
\sum_\sigma\int(P_IG_\sigma)^2.
\]
代入 `|E|/J(F)=e^{mu N+o(N)}`、实际 budget、及上式平方收费，指数为 `lambda-H_2(p)+3delta+o(1)`。先定 `delta` 和独立于最终 `a` 的 `c`，再取 `a` 使 `2D_a<c/2`，最后求和 `exp(o(N))` 个 `(I,d)`，避免了常数循环。Markov 的平方形式遂给所有 `M^sigma G_sigma` 真成功集在 shell 中指数小。

## 熵率、匹配与限制

必要成本的正确熵是
\[
p=\frac{r_*}{1+r_*}=\frac{\mu-1}{\mu},\qquad
H_2(p)=\mu-1+\log\mu,
\]
即 `sign U` 负号模式熵，不能误用 r085 的未定 receiver-sign 计数。r085 Hamming-ball mask 保留每个原 canonical rectangle，费用正是该熵率；允许换矩形只会降低 infimum。因此对固定 `theta in (0,1]`，
\[
\lim N^{-1}\log b_N^{\rm free}(\theta)=H_2(p).
\]

本机 MATLAB 公式 sanity check 得
`mu=1.27846454276107`、`r_*=0.278464542761074`、
`p=0.2178117057198`、`H_2(p)=0.524124324657529`、`m(r_*)=mu`；PID `11609` 已退出并由 `ps` 验证。数值不作为证明依据。

必须保留的范围：结论只约束“每个 relative direction 单独支付一个真实 source piece”的 `M^sigma G_sigma`，不排除截源后用另一种架构（尤其 r083 固定物理角点跨方向复用）成功；也不否定直接以物理体积支付该 shell，且不导出一般 `B,D` 分解或新的 `C_N` 上界。

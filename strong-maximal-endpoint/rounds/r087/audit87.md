# r087 固定物理角点 shell 效率独立审计

结论：`corner87.md` 的 restricted fixed-physical-corner 效率界通过审计：
\[
\log K_E\le-\frac12N\log N+\frac12N\log\log N+O(N).
\]
它是 r084 tiny slab 之外、覆盖 r085 近 `Q_N` 物理 shell 的更强障碍；但仍只否定真固定物理角点 columns 在该大 shell 上以 `exp(O(N))J(F)` 的 source fee 支付正比例，不否定把此 shell 按物理 residual 支付，也不改变一般 `C_N` 结论。

## 1. 从实际矩形到 r022 的 fixed-vector 资格

对一个严格成功的真实矩形，写 `U_i` 在 `1,d_i=a_i/x_i` 之间均匀，`W_i=(-log|U_i|)_+`。由 cube cutoff 和径向单调性，
\[
F(xU)\le h_N(t(x)+\sum_iW_i).
\]
若 `d_i=-r_i, 0<=r_i<=1`，这正是 r022 的 increment `Z_{r_i}`。其余端点可由独立 `Z_0~Exp(1)` 单调耦合支配：short inward `0<d<1` 有
\((e^{-s}-d)_+/(1-d)\le e^{-s}\)，outward `d>1` 的正截断为零，long crossing `d=-r,r>1` 的尾为 `2e^{-s}/(1+r)<=e^{-s}`。因此严格实际成功强迫一个确定的 `tilde r in [0,1]^N` 满足
\[
\mathbb E h_N(t+\sum_iZ_{\widetilde r_i})>1.
\]

移位 profile 的 modular 至多 `e^t j_N`。由于 `t<=sqrt(log N)`、`log j_N=-log Q_N+o(log N)`，它满足 r022 已验收的 near-`Q_N` fixed-constant test；所需是 r022 内部的 action-deficit 结论，而不是只引用最终 `log Q_N` 展开。于是
\[
D(\widetilde r)=\sum_i[\mu-m(\widetilde r_i)]\le C\log N.
\]
这一步的量词正确：vector 只为这一 receiver/矩形服务，r022 正是允许任意确定 vector 和非减 profile 的标量定理，并未把 receiver-dependent actions 当作反馈或一个共同源替换。

固定小 `e` 后，不规则或 `|r_i-r_*|>e` 的坐标有固定缺口；`m''(r_*)<0` 给其余坐标的二次缺口。因此
\[
|I|=O_e(\log N),\qquad
\sum_{j\notin I}(r_j-r_*)^2=O_e(\log N).
\]
这里保留所有原始 irregular sides（包括正截断后可有 atom 的情况）是必要且原稿已做到。

## 2. 反射、outer box 与 signed Hardy (L^2)

固定角点后，对 regular coordinates 的反射只是一套随 `a,I` 固定的坐标变换；它同时作用于 source/receiver，保持 `F` 的坐标偶性、`J(G)` 和 Lebesgue measure。正文“same true source”的措辞宜理解为这个固定坐标系中的同一 source，而不是不经反射逐点替换 source。

令 `A_j=|a_j|`。若 regular `a_j=0`，该类为空；否则反射后
\[
x_j=A_j/r_j\in[A_j/r_+,A_j/r_-],\qquad O_j=[-A_j,A_j/r_-].
\]
直接计算给
\[
\frac{|O_j|}{A_j+x_j}\le
R_e=\frac{r_+(1+r_-)}{r_-(1+r_+)},\qquad
\frac{u_j-\ell_j}{|O_j|}=
P_e=\frac{2e}{r_+(1+r_-)}.
\]
实际 regular box 包含于 `O`，非负性使真实成功推出 `T(x_I)>R_e^{-n}`，方向无误。

异常变量上的 `H_{a_I}^{(I)}` 是真正固定角点 tensor Hardy。按 relative-to-`a_I` 的 source orthants 分割后，接收 orthants 是直和；每一维 half-line Hardy 的 `L^2` norm 是 2，故完整 signed tensor norm 为 `2^k`，没有再次叠加 orthant 个数。Jensen/Fubini 给
\[
\int T^2\le\frac{4^k}{|O|}\int_{O\times\mathbb R^k}G^2
\le\frac{2\,4^k}{|O|}J_N(G).
\]
这正确使用了 cap two 的 `G^2<=2G<=2Phi_N(G)`，而非未证明的 maximal `L^1` 界。

对于每一个固定 `I`，G 保持同一个函数；仅 success subset 被分组。将式 (8) 对 `I` 求和只是 Lebesgue 上界的相加，不产生 source reuse 的免费假设。角点无限、`a_j=0` 和 coordinate-plane null sets 均已被 uniform bound、空类约定和有理端点/局部 `L^1` 连续性处理。

## 3. 球体积与组合计数

在 regular 参数球中使用 `x_j=A_j/r_j`，其 Jacobian 至多 `r_-^{-2n}prod A_j`，而 `|O|=(1+1/r_-)^n prod A_j`。因而所有 `A_j` 消去，且
\[
\frac{|X_{J_0}|}{|O|}
\le C_e^n\left(\frac{\log N}{n}\right)^{n/2}.
\]
参数 cube 截断只能减小该球体积。代入 (8) 的主要对数项为
\[
-\frac n2\log\frac n{\log N}+O(N)+O(k).
\]
因 `k=O_e(log N)`，
\[
\log\sum_{k\le C_e\log N}\binom Nk4^k=O((\log N)^2),
\]
且以 `n=N-O(log N)` 替换只损失同阶项。故没有漏掉的 `N log N` 组合项，最终为
\[
-\frac12N\log N+\frac12N\log\log N+O(N).
\]

## 4. 覆盖/残差结论和范围

对 countable fractional cover，
\[
|E|-D_E=\int_E\min(1,\Gamma)
\le\sum_iw_i|E\cap\{H_{a_i}G_i>1\}|
\le K_E\sum_iw_iJ(G_i)
\]
是 Tonelli 的直接应用，不要求 successes 不重叠。结合 `|E|/J(F)=Q_N exp(o(log N))` 给出文中的 (2)--(3)。全覆盖或支付固定正比例确实需要 superexponential source budget。

范围必须保留：该障碍正适用于 r083 的真实固定 physical-corner、每列 `G<=F`、constant countable weights 的跨 relative-direction 复用列架构；在这个指定 `E` 上，即使作任意跨方向复用，支付固定正比例仍需超指数 budget。它不否定 r083 的一般 fixed-corner 常数 `K_H`，不排除 `D_E/|E|->1` 的直接物理 residual 付款或其他非列架构，也不改变一般 `C_N` 上下界。

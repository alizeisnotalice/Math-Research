# r085 方向残差结论独立审计

审计对象为当前 `core85.md` 的第 1--2 节，并对任务说明中尚未写入该文件的“薄壳--背向坐标--低方向重合”拟议结论核对其所需前提。结论分两部分：

* 当前第 1 节的**固定单侧方向**定理可成立，且它不与双侧非中心极大算子的反馈模型矛盾；第 2 节的有限方向计数也成立。
* 薄壳方向重合需要把 `K_{H,N}`、Gamma 一侧 Hardy 常数和 r028 的一般固定比例标量 `Q_N` 明确区分，也需要从固定维数的 `L -> infinity` 极限作合法的随 `N` 对角选择；core85 后续第 3--4 节已经补入这些论证，最终结论见第 4 节。

## 1. 单一方向的精确常数：当前证明路线正确

令 `sigma in {+,-}^N`，并允许接收点处各边长自由选择。该算子是单一正交方向的真实 corner-rectangle maximal operator；它不是固定物理锚点本身。关键区别是：

* 对**单侧**一维算子 `P^sigma`，rising-sun 给
  \[
  s|\{P^\sigma f>s\}|=\int_{\{P^\sigma f>s\}}f.
  \]
  全空间层蛋糕/Hardy--Littlewood 因而给
  \((P^\sigma f)^*\le Hf^*\)。
* 对同一坐标中**两侧的并**，会出现交集/正尾项，才得到 r011/`strong-maximal-adaptive-upper.md` 的反馈算子 `T`，而非 `H`。因此不得把这里的一侧结论误推广到双侧 `M_i`。

在有界紧支撑输入上，逐坐标应用一侧结论可得
\[
 M^\sigma F\le P_N^{\sigma_N}\cdots P_1^{\sigma_1}F,
 \qquad (M^\sigma F)^*\le H^N F^*.
\]
这里不需声称前一个不等式取等；反复把非负的边长 supremum 拉入积分已足够。当前稿关于中间函数的尾部/未使用坐标紧支撑的说明，足以避免把非 `L^1` 的中间函数错误地直接送进全局 Fubini。对紧支撑纤维和固定正水平，rising-sun 集是有界开集，故其连通分支至多可数且端点有限；逐分支的端点等式可相加。这是“sum the component equalities”的充分严格化。

严格阈值也处理正确：先在 `s'>s` 得到
`Hf*(D_{s'})>=s'>s`，再令 `s' downarrow s`，不能直接以等号阈值代替 strict halo。

令
\[
h_u(s)=F^*(ue^{-s}),\qquad S\sim\Gamma(N,1),
\]
则
\[
H^NF^*(u)=\mathbb E h_u(S),\qquad
\int_0^\infty e^{-s}\Phi_N(h_u(s))ds
=u^{-1}\int_0^u\Phi_N(F^*(v))dv.
\]
这正是 r028/r083 的 `r=0` Gamma 标量问题，因此得到上界 `K_{H,N}^{[A]}`。反向构造也对：把 r028 的坐标偶对称低输入、其 receiver 集都限制到同一正交象限，会把 source modular 和 receiver 体积同乘 `2^{-N}`；`[0,x]` 在该象限内，故全后向单侧方向保留严格成功与同一比值。反射给任意 `sigma`。这证明的是每个**固定**方向的 cap-`A` 最优常数恰为 `K_{H,N}^{[A]}`。

需在正式稿中保留的定义性细节：`M^sigma` 的边长取有理正数即可保证可测，随后用局部 `L^1` 下矩形积分关于端点的连续性恢复全体边长；一般有限模输入的盒截断则用单调收敛处理。这些细节目前已基本写明。

## 2. 方向多重度与 source-list 计费：正确但须固定记号

方向集合是有限集 `{+,-}^N`，故
\[
\int m_F=\sum_\sigma |\{M^\sigma F>1\}|
\le2^NK_{H,N}^{[A]}J_N(F),
\]
再以 Markov 得到高多重度集的式 (1)。任意固定子集 `U` 的 union bound 同理为 `|U|KJ(F)`；这里不涉及不可数 directions 或可测选择。

若 `Lambda(y)` 是 source point 的可测方向列表，则
\[
\sum_\sigma J_N(F1_{\{\sigma\in\Lambda\}})
=\int\#\Lambda(y)\Phi_N(F(y))dy,
\]
故式 (2) 正确。文字上应将“receiver-dependent fixed rectangle lengths and anchors”改为：**长度/由长度确定的角点可随 receiver 变，但 `Lambda` 是 source-dependent 的分配**；否则容易被读成未证明的 receiver-dependent source sharing。

`Q_N` 在第 2 节必须明确定义并注明是否固定同一 cap `A`。只要它是所需比较中满足 `Q_N>=K_{H,N}^{[A]}` 的标量基准，
`q=ceil(2^NK/Q_N)` 的残差代数正确；若欲使用整数阈值，也应说 `q>=1`。

## 3. 薄壳/背向坐标命题的初始审计清单

任务说明中的下列链条不能仅凭当前两节宣称：薄壳 `t in [0,sqrt(log N)]` 上成功方向至多有 `o(log N)` 个背离原点坐标，继而
`m<=exp(o((log N)^2))`，以及大方向清单只能捕获指数小 receiver 比例。它至少需要以下独立引理和量词。

1. 对单调径向 profile `h`、正交象限内 receiver，背向原点的坐标平均可被其在 receiver 坐标的值上界，故可删除；须先说明向外区间不穿越零（或截到零后外部 source 为零）。余下 `n=N-k` 个向内坐标经跨零规约后给一个固定 ratio vector 核；仅其 `r_i=0` 特例才是 Gamma-`n` Hardy 平均。

2. 在向内区间已规约为 `[-r_i x_i,x_i]` 后，剩余坐标的核是 r028 的一般**固定 ratio vector**核（可跨原点），不是 Gamma/one-sided 核。由某点成功直接把移位 profile 纳入该 fixed-vector 标量定义，随后对 ratio vector 取上确界，得到
   \[
   1\le Q_n^{[A]}J_n(h(t+\cdot))
   \le Q_n^{[A]}e^tJ_N(h).
   \]
   strict feasibility gives the lower cost `J_n>=1/Q_n`; the second inequality uses `Phi_n<=Phi_N`. The Gamma constant occurs only in the special all-one-sided (`r=0`) case. It cannot be inserted before the vector optimization here, and `Q_n` is not `K_{H,n}`.

3. 若目标是从该式推出 `k=o(log N)`，必须给出一列 source profiles 的定量近最优性，例如
   \[
   \log(Q_N^{[A]}J_N(h_N))=o(\log N),
   \]
   并在所需的 `n=N-k` 区间证明
   \[
   \log Q_{N-k}^{[A]}-\log Q_N^{[A]}=-\mu k+o(\log N)
   \]
   的**统一**版本（至少足以排除 `k>=epsilon log N`）。对每个固定 n 的渐近、或没有统一性的 `o(1)`，均不能推出该结论。随后才可用
   \(\sum_{k\le o(\log N)}\binom Nk=\exp(o((\log N)^2))\)。

4. r028 的低输入使用先固定维数、再令空间参数 `L -> infinity` 的极限。它本身不给随 `N` 的薄壳质量。需要显式选择一个对角序列 `L_N`，并同时控制 modular、receiver volume 及 `0<=t<=sqrt(log N)` 的占比误差；不能交换这两个极限。历史 r069 还明确指出：全局近 `Q` 效率不自动给每个 receiver 的 local band/capture qualification。因此这里若只讨论特定单调 radial near-`Q` 源，应把这个限制写成假设，不能升级为一般近极值输入结论。

上述是对初版尚未写出后续段落时的必需条件清单；下一节记录补充后的实际核对结果。

## 4. core85 后续补充后的最终复核

随后加入的第 3--4 节补足了上述待证环节，审计结论更新为接受其**限定性**结论。

第 3 节现在只针对显式构造的一列有界、单调径向 `h_N`，并写明
`log(Q_N j_N)=o(log N)`；没有从一般 global near-optimality 推出局部性质。它还先固定 `N`，再选有限的 `L(N)` 同时控制 shell volume 和 modular 到相对误差 `1/N`。这是合法对角选择，得到 `|E|/J(F)=Q_N exp(o(log N))`，并未交换 `L,N` 极限。

向外坐标的单调删除和跨零的向内区间规约为 `[-r_i x_i,x_i]` 是正确的。每一个成功 receiver 随后只产生一个固定 ratio vector，故 r028 的 fixed-vector scalar bound 可以直接使用；不存在把 receiver-dependent ratios 当作 feedback 的步骤。式 (4) 正确地由该 vector scalar definition 给出 `J_n>=1/Q_n`，再以 `Phi_n<=Phi_N` 获得右端。这里不应先插入 Gamma/one-sided 常数：跨原点 ratio kernel 一般不是 Gamma，且 `Q_n` 不等同 `K_{H,n}`。对每个固定 `epsilon>0`，`k=floor(epsilon log N)` 时 `-mu k` 压过 shell 的 `t<=sqrt(log N)` 及近最优误差；既有两个维数的 `o(log N)` 展开之差仍为 `o(log N)`，而 `d'` 估计控制移动维数的二阶项。因此可对角地取得 `k_N=o(log N)`，二项式方向计数及 sign-orthant 体积计数成立。

第 4 节的精确熵率匹配定理也通过审查。这里须用
\[
p=\frac{r_*}{1+r_*}=\frac{\mu-1}{\mu},\qquad
H_2(p)=\mu-1+\log\mu,
\]
即 `sign U` 为负的模式熵；这不同于仅数“receiver sign 尚未确定的坐标”时的 `2r_*/(1+r_*)`。典型模式各自的概率至多
`exp(-(H_2(p)-delta)N)`；固定 `U` 下 `X -> Y` 的 Jacobian 密度至多 `e^S/|E|`，且对每个固定 `y`，
`b -> -sign(y)b` 是方向标签的双射。于是实际 modular 预算直接给出
\[
\mathbb E\,G_{-\operatorname{sign}X}(Y)
\le e^{-cN}+\frac{e^{(\mu-H_2(p)+2\delta)N}}{|E|}
\sum_\sigma J_N(G_\sigma).
\]
结合 `|E|/J(F)=e^{\mu N+o(N)}`，预算指数率小于 `H_2(p)` 时，保留的 canonical-witness 平均指数小；由非负性和 Markov，同一见证仍成功的 receiver 比例亦指数小。

反向 Hamming-ball 构造同样正确：每个 source point 保留距 `-sign y` 至多
`Np_N+2 sqrt(N log N)` 的方向，实际 modular 倍数正是该球的基数，Stirling 给出 `H_2(p)`；对任意原 receiver，丢失只在 binomial 尾，Hoeffding 使平均损失至多 `2N^-8`，小于原来的严格余量 `N^-2`。故它保留全部原始 witness rectangles，并为每个固定 `theta in (0,1]` 给出同一指数上界。

需保留三项范围说明：

1. 这是“用一个 reusable source piece per relative direction 保留所展示的 `R_x`”的准确指数率；它不排除截源后以**另一**矩形成功。
2. `lim N^-1 log b_N(theta)=H_2(p)` 不意味着临界预算 `exp(H_2(p)N+o(N))` 自动保留某个固定正比例；这种临界 liminf 需要额外的次指数精度。
3. 它没有否定将整个 shell 按物理体积作为 residual 支付，也不导出一般 `B,D` 分解或新的 `C_N` 上界。

综上，当前 `core85.md` 的固定方向精确常数、方向残差、多重度、低重合 shell，以及 canonical-witness retention entropy-rate 结论均可接受，以上范围限制不可删。

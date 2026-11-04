# r005：阈值误差的定量形式及一阶矩障碍

本轮只研究真实上界证明中的 S3/S14 交界成本，不展开周期下包络。已读 `r004-global/global_sharp_candidate.md` 与 `r004-global-review/sharp-charge.md`。沿用已审核通过的 r004 有限覆盖、深分配、私有点见证和矩形交叠结论。未改论文或代码仓库，未执行 Git。

**本轮得到的明确命题：** 将分位削去量从 \(\sqrt v\) 改为 \(2^{n+2}\delta\)，可把 r004 的上界误差定量写成

\[
C_{n,C_0}a\{\delta^2+(v+\delta)F(L)+t\},\qquad \delta=Lv.
\tag{P}
\]

因此，若某个有界阈值的例外体积是 \(O_n(vS)\)，则目标
\(O_n(t^{1+2/n}+t^2)\) 余项确实成立。另一方面，仅有 r004 已证明的一阶矩约束，不能把 (P) 普遍改善到 \(O(av^2+at)\)：下面给出有限尾分布反例，其所有阈值及这些阈值证书的凸平均都保留 \(\Theta(av/\log(1/v))\) 的误差成本。这是现有全局阈值账本的确定障碍，**不是**真实中心方体常数的反例。

## 1. 定量阈值命题

固定 \(n\ge2\)、\(C_0<\infty\)。令 \(H=2^n\)、\(a=Ht\)、\(v=a^{1/n}\)。取 r004 的有限见证覆盖数据，并假设

\[
M\ge a|E|(1-C_0v).
\]

记

\[
d_i=e_i+1-\frac{m_i}{ar_i^n}\ge0,\quad
z_i=d_i/v,\quad
F(L)=\frac1S\sum_{z_i>L}r_i^n.
\]

已有账本给

\[
\sum_i d_i r_i^n\le C_{n,C_0}vS,
\qquad S/|E|=1+O_{n,C_0}(v).
\tag{1}
\]

**命题 1。** 存在仅依赖 \(n\) 的 \(\delta_n>0\)，使对任意 \(L\ge1\)、\(\delta=Lv<\delta_n\)，

\[
\boxed{
\frac M{|E|}
\le a-nav+
C_{n,C_0}a\bigl[\delta^2+(v+\delta)F(L)+t\bigr].}
\tag{2}
\]

这里 \(L\) 可以依赖输入与有限族；常数不依赖 \(L\)。

### 证明：用线性分位量代替平方根分位量

规则指标仍定义为 \(d_i\le\delta\)。选

\[
B_n=2^{n+2},\qquad b=B_n\delta.
\]

取 \(\delta_n\) 足够小，使原局部填充、私有簇局部化、深分配、见证包含 \(G_i\) 的条件都成立，并保证 \(2nb<1/2\)。S10 的相同削去步骤得到

\[
\frac{m_i^*}{ar_i^n}
\ge1-C_n(\delta+b+\delta^n)
\ge1-C_n'\delta.
\tag{3}
\]

小半径点仍用其本身局部输入的外部 halo 界删除，所占比例是 \(O_n(\delta^n)\)。没有将小半径点按 \(o(1)\) 估计后误作更高阶量。

S13 中唯一使用 \(b=\sqrt v\) 的预算步骤，是需要

\[
bm_j>\delta ar_i^n.
\]

相交邻盒的半径比至少为 \(1/2\)，且 \(m_j\ge(1-\delta)ar_j^n\)，所以新的选法给

\[
bm_j\ge B_n\delta(1-\delta)2^{-n}ar_i^n
=4\delta(1-\delta)ar_i^n
>\delta ar_i^n
\]

只需 \(\delta<3/4\)。其严格余量为 \(\delta(3-4\delta)ar_i^n\)。因此一个临界方向的交叠仍然消失，不要求 \(\delta/b\to0\)；足够大的固定比值 \(b/\delta\) 已经够用。

由 (3)、\(R_i\ge(1-\delta)r_i\)，规则盒的体积增益定量为

\[
|X_i|\ge g_i/a+(n-C_n\delta)v r_i^n.
\tag{4}
\]

两个临界方向的交量及规则–例外交量与分位量无关，保留原始界

\[
\sum_{i<j}|X_i\cap X_j|
\le C_n(\delta^2S+\delta S_b),
\quad S_b=F(L)S.
\tag{5}
\]

把 (4)、(5) 代入体积账本，例外盒仍贡献其全部 \(g_i/a\) 首项，得到

\[
V\ge M/a+nvS
-C_n\bigl[v\delta+\delta^2+(v+\delta)F(L)\bigr]S.
\tag{6}
\]

利用 \(V/|E|\le1+O_n(t)\) 和 (1)，乘以 \(a/|E|\) 消去 \(V,S\)。把由 \(S/|E|-1\) 产生的 \(O(av^2)\) 合并，且 \(L\ge1\) 意味着 \(v^2\le v\delta\le\delta^2\)，就得到 (2)。所有小量条件只要求 \(\delta<\delta_n\)，所以常数对 \(L\) 一致。证毕。

### 直接可证明的中间结论

若固定 \(K,A<\infty\)，并存在一个 \(1\le L\le K\) 满足

\[
F(L)\le Av,
\tag{7}
\]

则 (2) 立即给

\[
\boxed{
M/|E|\le a-nav+O_{n,C_0,K,A}(av^2+at).}
\tag{8}
\]

这包括所有选中指标都满足 \(d_i\le Kv\) 的子类，也允许 \(O(vS)\) 体积的任意大局部亏损例外。换回 \(t\)，余项为

\[
O_{n,C_0,K,A}(t^{1+2/n}+t^2).
\]

要将 (8) 升格为真实 \(C_n\) 的统一上界，尚需证明：对所有与该上界有关的输入，总能选取覆盖使 (7) 成立，并使 \(K,A\) 仅依赖 \(n\)。r004 的一阶矩结论本身没有提供这一步。

## 2. 一阶矩不能普遍消去阈值对数：严格有限反例

**命题 2。** 存在有限概率分布，满足 \(1\le z\le1/v\)、\(\mathbb E z<2\)，但其尾分布 \(F(L)=\mathbb P(z>L)\) 使

\[
\inf_{L\ge0}
\left\{v^2L^2+v(1+L)F(L)\right\}
\asymp\frac v{\log(1/v)}.
\tag{9}
\]

因此不存在仅从 \(\mathbb E z\le2\) 推出、对所有输入通用的 \(O(v^2)\) 阈值成本界。

### 有限分布及完整证明

取任意整数 \(J\ge2\)、\(v=2^{-J}\)。在 \(z_k=2^k\)、\(1\le k\le J\) 处赋予质量

\[
w_k=\frac1{J2^k}.
\]

其余质量赋予 \(z=1\)。因为

\[
\sum_{k=1}^Jw_k=\frac{1-v}J<1,
\quad
\mathbb E z=2-\frac{1-v}J<2,
\]

这确为满足所述矩条件的有限分布，且没有用 \(z=0\) 的背景质量制造退化。

若 \(L<1\)，则 \(F(L)=1\)，所以成本至少为 \(v\)。若 \(1\le L<2^J\)，令唯一整数 \(m\in\{0,\ldots,J-1\}\) 满足

\[
2^m\le L<2^{m+1}.
\]

精确尾分布为

\[
F(L)=\frac{2^{-m}-2^{-J}}J,
\]

包括 \(L=2^m\) 时的严格尾约定。因此

\[
LF(L)\ge2^mF(L)
=\frac{1-2^{m-J}}J\ge\frac1{2J}.
\]

该区间的成本至少为 \(v/(2J)\)。若 \(L\ge2^J=1/v\)，则 \(v^2L^2\ge1\)。所以对全部 \(L\ge0\)，

\[
v^2L^2+v(1+L)F(L)\ge\frac v{2J}.
\tag{10}
\]

反之选 \(L=1\)，得到

\[
v^2+\frac{2v(1-v)}J\le\frac{3v}J,
\]

因为 \(J2^{-J}\le1\)。这证明 (9)，且 \(J=\log_2(1/v)\)。

这个分布也使任意固定有界阈值的尾质量约为 \(1/J\)，远大于 \(v\)，所以不能由一阶矩自动推出 (7)。

若把不同阈值的证书按任意概率测度 \(\mu\) 作凸平均，则

\[
\int\left[v^2L^2+v(1+L)F(L)\right]d\mu(L)
\ge\frac v{2J}.
\]

因此，**对既有全局阈值证书作线性混合**不会改善此障碍。这里没有声称所有逐盒连续权重或所有空间多层构造都属于这种凸平均；那些做法有机会使用本反例舍弃的空间关联信息。

## 3. 结论的准确边界

命题 2 只反驳“已知一阶矩足以把现有统一阈值成本降到 \(O(v^2)\)”的推论。它不是实际立方体族的构造，不认证给定尾分布能同时实现 r004 的全部几何约束，更不反驳真实 \(C_n\) 存在所需下一阶余项。

本轮确定的障碍位于信息层：若后续仍只保留 \(\sum_i d_i r_i^n=O(vS)\)，然后只按一个共同阈值削去例外，或线性混合这些证书，\(v/\log(1/v)\) 已是该账本的最优一般量级。继续改变 \(T\) 的幂、机械 Taylor 展开或把硬阈值改写成这些证书的软平均，都不会解决目标。

命题 1 比 r004 的固定 \(b=\sqrt v\) 定量更强：规则部分在 \(\delta=O(v)\) 时只损失 \(O(v^2S)\)，而非 \(O(v^{3/2}S)\)。因此真正剩下的是高亏损邻盒交界的空间计费，或比一阶矩更强的局部结构。

## 4. 下一条明确攻关命题：邻边自适应分位

以下为**尚未证明的下一命题**，不列为本轮成果。它改变空间接口而非再次优化全局阈值。

在半径可比、扩大半盒可能相交的邻接图上，每个顶点的邻居数仅依赖 \(n\)。可尝试令

\[
b_i\asymp_n v+\max_{j\sim i}d_j,
\]

使每个规则–规则唯一临界面预算逐邻边闭合。若还能得到适当的局部平方成本

\[
\sum_i d_i^2r_i^n=O_n(v^2S),
\tag{11}
\]

则有界邻居次数会把 \(\sum_i b_i r_i^n\) 控制为 \(O(vS)\)，双薄方向的邻边成本预计由 \(\sum_i d_i^2r_i^n\) 控制，最终给 \(O(v^2S)\)。不过 (11) 不是现有一阶矩的推论，也未证明可从二阶近极值自动导出；规则–例外的局部平方收费还需完整几何论证。后续若做此路，应直接证明这两个接口，不能以“连续权重”命名替代它们。

## 5. 实际 Mathematica 核验

实际运行本机 `/usr/local/bin/wolframscript`，版本 `15.0.0 for Mac OS X ARM (64-bit) (May 15, 2026)`。未设置数学运行时限，没有启动空 CPU 等待任务。可复现核心输入：

```wl
Clear[J,k,m,v,dd,n];
wsum=FullSimplify[Sum[1/(J 2^k),{k,1,J}],
 Assumptions->Element[J,Integers] && J>=2];
FullSimplify[1-wsum+Sum[2^k/(J 2^k),{k,1,J}],
 Assumptions->Element[J,Integers] && J>=2]
tail=FullSimplify[Sum[1/(J 2^k),{k,m+1,J}],
 Assumptions->Element[{J,m},Integers] && J>=2 && 0<=m<=J-1];
FullSimplify[2^m tail,
 Assumptions->Element[{J,m},Integers] && J>=2 && 0<=m<=J-1]
FullSimplify[(2^(n+2) dd)(1-dd)/2^n-dd,
 Assumptions->n>=2 && dd>0]
FullSimplify[(v^2+2v(1-v)/J)/(v/J)]
```

实际输出包括：

```wl
wsum = (1-2^(-J))/J
mean = (-1+2^(-J)+2J)/J
tail = (-2^(-J)+2^(-m))/J
2^m tail = (1-2^(-J+m))/J
budgetMargin = (3-4dd)dd
costRatioAtL1 = 2+(J-2)v
```

这是符号转录核验。命题 1 的严格论证是重新核对分位预算、矩形体积和全局收费；命题 2 的严格反例是上面的有限分布及全部阈值分支，不以软件枚举替代证明。

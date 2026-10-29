# 增长网格上的凸混合比较：概率预算与好/坏尺度合并

先实际运行 `check_budget.py`：使用 Fraction 做参数预算，并完成16个有限原子模型（均匀和固定非均匀权）检查，随后给出本证明。它是受限 source family 与有限增长网格的比较，不是一般 $C_n$ 上界。

## 1. 前提与记号

沿用接受基线：固定 $0<p<1,w\ge4$，$q=K+1=2^d$，

\[
v=\frac{9p^2n}{29K^2}\ge16\log n,\quad
v\le B<2v,\quad\epsilon_{box}=K\delta=1/(wv),\quad
K/v^s\longrightarrow\infty\quad(\text{每个固定 }s>0).
\]

设 $R\le v^r$，r固定。各 bank 用同一域、同一查询和同一 Fourier 补偿 B_L；各自真实周期来源质量1且严格正。令

\[
Y_{b,L}=1+a_L(H_{b,L}+B_L),\quad a_L=p/(KL),\quad
A_{b,L}(X)=\prod_{i=1}^nY_{b,L}(X_i),
\]

所以在同一 Haar 空间 $P$ 上 $A_{b,L}>0$、$EA_{b,L}=1$。这些是实际平均形成的 anchor 密度，并非原 input density 倾斜。

令有限网格 $\Lambda\subset I=[10/9,11/9]$，$J=|\Lambda|\le C_JKv$。按另一部分的 joint-safe 判据分成 $\Lambda_s$ 与 $\Lambda_u$。本证明使用如下精确接口：

1. safe点的全部所需单端点映射秩为2b，因而各表同一 $D_L=\int H_b^2$、各异表同一 $C_L=\int H_bH_c$，且
\[
\Delta_L=D_L-C_L=\frac{2q}{K}\operatorname{var},\qquad
\operatorname{var}=\frac{B^2-1}{6B^2}-\frac{\epsilon_{box}(B-1)}{6B}.
\]
2. 坏点数 $m=|\Lambda_u|\le C_uR^2B^2(v+1)$，常数独立于n。

显式计数可使用 $N=B-1$、$H\le1+R^2N^2$，网格间距至少 $a/(Kv)$ 时
$m\le\min\{J,H(\lfloor4v/a\rfloor+1)\}$。source family存在性、bad lag计数与秩判据见同轮 [joint/grid证明](../round050_joint_grid/PROOF.md) 的(4)、(7)、(9)；本证明不把这些条件扩为任意既定banks或连续所有尺度。

## 2. safe点的统一亲和度指数

$V=\sum c_m\le40/9$，故 $|Y_{b,L}-1|\le\kappa_L=a_L(2+40/9)$，最终统一有 $\kappa_L\le1$。此外

\[
\operatorname{var}=\frac16[1-B^{-2}-\epsilon_{box}(1-B^{-1})]\ge1/7
\]

最终成立；例如 B≥v、w≥4、v≥4 时，损失项≤$v^{-2}+1/(4v)\le1/8<1/7$。

共同 B_L 在 $Y_b-Y_c$ 中准确抵消。49轮的 Hellinger 推导于是逐 safe L 给

\[
\int\sqrt{A_{b,L}A_{c,L}}\,dP\le e^{-E_L},\quad
E_L=\frac{n a_L^2\Delta_L}{4(1+\kappa_L)}.
\]

由 q/K≥1、var≥1/7、κ≤1、L≤11/9，

\[
E_L\ge\frac{np^2}{28K^2L^2}
=\frac{29v}{252L^2}
\ge\frac{261}{3388}v>\frac v{16}\ge\log n.
\tag{1}
\]

其中 $261\cdot16=4176>3388$，系数无舍入。这里没有删除HB，没有Gaussian迁移，也不需要bank独立性。

## 3. 均匀或固定非均匀权的safe overlap费用

取确定性权重 $\pi_b>0$、$\sum\pi_b=1$，在所有尺度使用同一组权重；允许它们随n变化。零权bank可先删除。令 $G_L=\sum_b\pi_bA_{b,L}$，$\varepsilon=1/v<1$。在每个safe L设

\[
E_{b,L}=\left\{\sum_{c\ne b}\frac{\pi_cA_{c,L}}{\pi_bA_{b,L}}\le\varepsilon\right\}.
\]

这些事件两两不交，且在其上
$\pi_bA_{b,L}\le G_L\le(1+\varepsilon)\pi_bA_{b,L}$。
记 $Q_{b,L}=A_{b,L}P$、$Q_{mix,L}=G_LP$、Bad_L为所有好事件并集的补集。平方根次加性及Markov给

\[
Q_{b,L}(E_{b,L}^c)
\le\varepsilon^{-1/2}\sum_{c\ne b}\sqrt{\pi_c/\pi_b}\ e^{-E_L}.
\]

因此

\[
Q_{mix,L}(Bad_L)
\le\varepsilon^{-1/2}\sum_{b\ne c}\sqrt{\pi_b\pi_c}\ e^{-E_L}
\le(R-1)\sqrt v\,e^{-E_L}=:\delta_L.
\tag{2}
\]

最后用 $(\sum\sqrt\pi)^2\le R$。无需权重有统一正下界；小权bank的单独付款大，但混合加权后有上述统一界。权重若依赖位置、样本或随机赢家，上述固定凸混合结构与质量1不再自动成立，不在本命题内。

全部safe节点费用满足

\[
\delta_s:=\sum_{L\in\Lambda_s}\delta_L
\le C_JKv(R-1)\sqrt v\,n^{-1}
=O_p\!\left(\frac{v^{r+1/2}}K\right)=o(1).
\tag{3}
\]

等式使用 $n=29vK^2/(9p^2)$，最后用接受基线 $K/v^s\to\infty$。各 L 有自己的 anchor 测度；下文只逐 L 用它付款后求和，绝不把这些不同测度当同一个概率律。

## 4. 好/坏尺度必须在事件层面合并

记 $M_s=\max_{L\in\Lambda_s}G_L$，$M_u=\max_{L\in\Lambda_u}G_L$；空集最大值取0。**实际完整尺度最大值是 $\max(M_s,M_u)$，不是 $M_s+M_u$。**

对任意严格阈值λ>0，safe见证事件的49轮论证给

\[
\lambda P(M_s>\lambda)
\le(1+\varepsilon)\sum_b\pi_b
\left\|\max_{L\in\Lambda_s}A_{b,L}\right\|_{1,\infty}
+\delta_s.
\tag{4}
\]

具体地，除去所有 $\{G_L>\lambda,Bad_L\}$ 后，某个bank必有
$\max_s A_{b,L}>\lambda/[(1+\varepsilon)\pi_b]$；有限bank并集及阈值缩放给第一项。坏事件逐 L 满足
$\lambda P(G_L>\lambda,Bad_L)\le Q_{mix,L}(Bad_L)$，从而给第二项。

对不受控节点，只用其质量1及Markov：

\[
\lambda P(M_u>\lambda)
\le\sum_{L\in\Lambda_u}\lambda P(G_L>\lambda)
\le\sum_{L\in\Lambda_u}EG_L=m.
\tag{5}
\]

现在使用同一个λ的精确事件并集

\[
\{\max_{L\in\Lambda}G_L>\lambda\}
=\{M_s>\lambda\}\cup\{M_u>\lambda\}.
\]

将(4)、(5)相加，再对λ取sup，得到

\[
\boxed{
\left\|\max_{L\in\Lambda}G_L\right\|_{1,\infty}
\le(1+\varepsilon)\sum_b\pi_b
\left\|\max_{L\in\Lambda}A_{b,L}\right\|_{1,\infty}
+m+\delta_s
\le(1+\varepsilon)T_*+m+\delta_s,}
\tag{6}
\]

其中 $T_*=\max_b\|\max_{L\in\Lambda}A_{b,L}\|_{1,\infty}$。

这一步**没有**使用错误的weak-L¹函数和三角不等式。对非负函数和 $F+H$，一般必须作 $\theta\lambda,(1-\theta)\lambda$ 阈值split并支付相应因子；这里因为原对象恰为两个尺度max的max，同阈值union直接成立，不需要split，也没有给T_*多乘常数。若改用函数和去支配原对象，反而会引入不必要损失。

## 5. 渐近量级及可选更强写法

坏节点数

\[
m=O(v^{2r+3}),\qquad
\frac m{\sqrt n}=O_p\!\left(\frac{v^{2r+5/2}}K\right)=o(1).
\]

结合(3)、ε=1/v，(6)给所需受限机制结论

\[
\boxed{\|\max_\Lambda G_L\|_{1,\infty}
\le(1+o(1))T_*+o(\sqrt n).}
\]

其实此有限网格还可用纯质量预算 $T_*\le\max_b\int\max_L A_{b,L}\le J$，从而
$\varepsilon T_*\le J/v=O(K)=o(\sqrt n)$；故也可写
$\|\max_\Lambda G_L\|_{1,\infty}\le T_*+o(\sqrt n)$。
这只用了节点数上界J，不是给individual product source证明了 $O(\sqrt n)$ 弱型上界；T_*本身仍可能比√n更大。

## 6. 已运行的算术/有限模型检查与限制

`check_budget.py` 在 n10³²、p.99、K2⁴⁷−1、B2048、w64 的实际参数下，精确计算 var≥1/7、κ<1、$E_{min}\ge118.3020891$。取保守 $J\le\lceil Kv\rceil+1$、R4，费用以有理 $J(R-1)\lceil\sqrt v\rceil2^{-\lfloor E_{min}\rfloor}$ 包络，约 $7.80\times10^{-17}$；只用(1)中的 $n^{-1}$ 则约 $2.59\times10^{-13}$。a=1时显式粗坏点帽411847268335，约为√n的 $4.12\times10^{-5}$。

同时保存 R≈v、R≈v² 的纯算术预算；在当前n，R≈v的粗坏点帽已约6.07√n，并不小。没有把渐近 $o(\sqrt n)$ 冒充当前n对全部多项式R的有限认证，也没有宣称已构造这些大R来源。

16个原子模型使用36原子的共同概率空间、6尺度、3bank正质量1函数，分别用均匀权和固定权[1/9,4/9,4/9]，以Fraction检查好事件不交、混合坏集质量、坏节点mass1预算、同λ事件union及完整比较。这里用抽象有限模型核对概率接口，不替代真实GF joint-safe条件；脚本没有重算或修改旧轮来源。

最小交付为 `check_budget.py`、`RESULTS.json`、本证明。依赖Python标准库与已有NumPy（仅生成固定seed排列）：

```sh
python3 check_budget.py
```

源SHA与精确分数保存在结果中。实际预算小数仅作参考；有限Fraction核对是精确的，渐近证明依赖明列的source/网格前提。本命题不包含连续所有尺度，不给单表T_*的新上界，不排除单个product source本身的进一步增长，也没有证明或反驳一般 $C_n$ 的正对数改进。

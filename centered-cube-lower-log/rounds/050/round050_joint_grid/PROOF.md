# 第50轮：真实 joint GF family 的存在与增长网格的安全剔除

本稿在本机运行 `probe.py`、保存 `RESULTS.json` 后提出。已证明：任意固定多项式数量的bank最终可在同一有限域中同时选择；对适当分离的有限增长网格，删除至多多项式数量的坏节点后，所有二阶端点对都具有满秩。这里只服务于指定真实凸混合source family，不是一般 C_n 上界，也不覆盖连续全部尺度。

## 1. 单bank预算与joint family存在

固定 q=2^d、K=q−1、一个域表示和本原元g。取 B=2^b≥4，并写

\[
N=B-1,\quad T=N(N+1)/2,\quad D=N(N-1),
\]
\[
H(B)=N+T(T-1)+2D^2+4D<3B^4.
\]

使用已接受单bank有限域引理的同一个随机选择：独立均匀取b个域元素作为候选基。坏概率≤H(B)/q；四项分别支付线性不独立、含平方的乘法Sidon冲突、g邻接比率、四个短位移。每个多项式坏事件用次数/q支付：N个线性零式、`binom(T,2)`个非恒零二次式、D²个非恒零二次式及4D个线性式。故该概率界不要求先给候选基加条件。

现在独立均匀选择R个候选基块 `X_{a,0},…,X_{a,b−1}`。对每个非零二进制向量u，记
`L_{a,u}=Σ_i u_i X_{a,i}`。对不同bank a,c 和非零u,w，

\[
P(L_{a,u}=L_{c,w})=1/q.
\]

理由是等式为全部未加条件变量中的非零线性式，至少有一个变量系数为1；固定其他变量时，该变量只有一个解。每对bank共N²个等式。因此联合坏概率至多

\[
\boxed{\mathfrak B(q,B,R)
=\frac{R H(B)+\binom R2(B-1)^2}{q}.} \tag{1}
\]

不使用条件事件之间的独立性。若(1)<1，至少有一组同时各自合格且

\[
W_a\cap W_c=\{0\}\quad(a\ne c), \tag{2}
\]

其中W_a是第a个基的F₂线性张成。因为各bank内部基向量线性独立已保证每个W_a中的非零域元素对应非零u，而上述所有跨bank等式均被排除，(2)恰成立。

这是存在性命题，不是R个标签块联合独立的断言；尤其不要求Rb≤d，也不承诺枚举或随机搜索大型R family的实际成本。各bank之后均用相同p、δ、M和真实Fourier项构造已有正密度source。

### 多项式数量的渐近结论

沿用基线参数 `v→∞`、v≤B<2v、`K/v^s→∞` 对每个固定s>0。对任意固定r≥0及常数C，若整数R≤C v^r，则

\[
\mathfrak B\le\frac{48C v^{r+4}+2C^2v^{2r+2}}q\longrightarrow0.
\]

故充分大n时joint family存在，特别覆盖题设R≤v^r。有限n上预算≥1只是此union bound不够，不是不存在。

## 2. 精确定义opposite坏lag集合

写W_a*=W_a\{0}。在Z/KZ中定义

\[
\mathcal D_a=\{\log_g(u/w):u,w\in W_a^*,\ u\ne w\},
\]
\[
\mathcal E_{ac}=\{\log_g(u/w):u\in W_a^*,\ w\in W_c^*\}\quad(a\ne c),
\]
\[
\boxed{\mathcal D_* =\{0\}\cup\bigcup_a\mathcal D_a
\cup\bigcup_{a\ne c}\mathcal E_{ac}.} \tag{3}
\]

单bank Sidon性给|D_a|=N(N−1)，而跨bank只需粗界|E_ac|≤N²，不要求跨比率互异。所有有序bank对均已包括，故 `D_*=-D_*`；0额外处理同bank两端点恰在同一粗格的退化。

\[
\boxed{H_*:=|\mathcal D_*|
\le\min\{K,\ 1+R N(N-1)+R(R-1)N^2\}
\le1+R^2N^2.} \tag{4}
\]

“joint-safe”在本稿只指下述下述三类不同标签块配对的**2b秩**，不是全部4b、Rb行联合满秩。48轮已有44联合秩亏1而二阶仍可精确计算的真实例子。

## 3. mod K距离、floor/ceil和所有相位

对实数x和集合D_*定义

\[
\operatorname{dist}_K(x,\mathcal D_*)
=\min_{d\in\mathcal D_*,\ m\in\mathbb Z}|x-d-mK|.
\]

称尺度L安全，如果 `dist_K(KL,D_*)>2`。固定任意θ∈[0,1)，两端点粗格偏移为

\[
k_+(\theta)=\lfloor\theta+KL/2\rfloor,\qquad
k_-(\theta)=\lfloor\theta-KL/2\rfloor.
\]

对x=KL，写 `theta−x/2=m+u`、0≤u<1，则

\[
k_+-k_-=\lfloor u+x\rfloor
\in\{\lfloor x\rfloor,\lceil x\rceil\}. \tag{5}
\]

x整数时两者均为x；x非整数时两个整数与x的距离均<1。相位断点取指定floor定义仍满足(5)，没有例外。由距离的三角不等式和安全缓冲，(5)中的任一k均不在D_*模K中。2的缓冲是保守选择，不必取最小缓冲。

迹配对非退化，所以端点trace行秩等于相应域子空间之和的维数。于是：

- **不同bank，同一符号端点：** 子空间是g^{k_±}W_a及g^{k_±}W_c，由(2)交为0，秩2b；不靠尺度安全性。
- **不同bank，对端点：** 若 `g^{k_+}u=g^{k_-}w` 非零，则k_+−k_-属于E_ca⊂D_*；反向符号用D_*取负对称性。安全性排除该情形，秩2b。
- **同bank，两端点：** 非零交对应k=0或k∈D_a，均被D_*排除，秩2b。

这些结论对所有θ的选定floor表示成立。在相位边界，端点平滑函数本身连续，故标签写法的边界选择也不改变任何Haar积分。

因此在每个安全尺度，48轮非零态纤维计数和θ分段积分公式可对任意bank对使用。具体地，若

\[
\epsilon=K\delta=1/(wv),\qquad
\operatorname{var}=\frac{B^2-1}{6B^2}
-\frac{\epsilon(B-1)}{6B},
\]

则同尺度 `HH_diag−HH_cross=2q var/K`。这是二阶接口，完全不把共享Fourier项当成独立项。

## 4. 有限增长网格的坏节点数

设K≥7，I=[10/9,11/9]，Λ⊂I为任意有限网格，任意两点最小间距至少

\[
\Delta_L\ge a/(Kv),\qquad a>0\text{ 为固定常数}.
\]

首先

\[
J:=|\Lambda|\le\left\lfloor\frac{Kv}{9a}\right\rfloor+1. \tag{6}
\]

将满足 `dist_K(KL,D_*)≤2` 的节点标为坏，安全集取严格>2。这一闭坏集包括距离恰为2的边界，不会漏掉ceil/floor临界节点。

对固定lag d，坏邻域由lifted闭区间

\[
[(d+mK-2)/K,(d+mK+2)/K],\quad m\in\mathbb Z,
\]

与I相交组成，每段长≤4/K。**在I内至多一段：** 若两个相邻lift均与I相交，则在KL坐标中必有 `K−4≤K/9`，推出K≤9/2，与K≥7矛盾。此论证也处理circle在0处wrap的情况，不能把模距离直接误写成普通绝对值。

一个长≤4/K的闭区间内，间距至少a/(Kv)的点至多

\[
\left\lfloor\frac{4v}{a}\right\rfloor+1;
\]

加1不能省略：当长度正好是间距的整数倍，两端点可同时属于网格。因此

\[
\boxed{m_{\rm bad}\le\min\left\{J,
H_*\left(\left\lfloor\frac{4v}{a}\right\rfloor+1\right)\right\}
\le H_*(4v/a+1).} \tag{7}
\]

结合(4)，`m_bad=O_a(R²B²(v+1))`。若R=O(v^r)、B=Θ(v)，则m_bad=O(v^{2r+3})。基线 `n/K²=Θ_p(v)`给√n=Θ_p(K√v)，所以

\[
\boxed{m_{\rm bad}=o(\sqrt n).} \tag{8}
\]

(8)是渐近结果；有限参数有时粗计数已覆盖整个网格，必须保留这种失败情形。

本计数只删除有限网格节点。它既不控制坏区间上的连续上确界，也不证明节点间的模连续性或高阈值事件可离散化。不能据此宣布连续全尺度结论。

## 5. 坏节点的mass-one付款：限定于指定凸混合family

令每个bank的真实n维输入为乘积密度，A_{b,L}(x)是其中心立方体尺度L平均，故对Haar查询 `∫A_{b,L}=1`。取**固定、非数据依赖**权重π_b≥0、Σπ_b=1，考虑先作乘积后作全局凸混合的真实family

\[
G_L=\sum_b\pi_b A_{b,L},\qquad \int G_L=1.
\]

这不是每坐标重新抽bank的乘积混合。对有限坏节点，Markov和同阈值事件并集给

\[
\lambda P(\max_{L\in\Lambda_{\rm bad}}G_L>\lambda)
\le\sum_{L\in\Lambda_{\rm bad}}\lambda P(G_L>\lambda)
\le m_{\rm bad}.
\tag{9}
\]

因此该坏节点最大值的weak profile上确界≤m_bad，(8)使这一费用为o√n。这里没有使用错误的函数和weak-L¹三角式；只有max事件的同阈值并集。

其余安全节点的affinity/混合比较需另用对应命题。本稿不限制individual source本身的增长，也不对一般输入、一般C_n给上界。

## 6. 实际Fraction/Decimal预算和失败记录

参数取p=.99、n=10^m，m=32/42/60；选择最大d使 `v=9p²n/(29(2^d−1)²)≥16 log n`，B为不小于v的最小二幂。脚本用

`log10=2 Σ_{j≥0}(9/11)^{2j+1}/(2j+1)`

的16项有理部分和及几何余项给严格log10夹逼，同时认证当前d满足及下个d不满足，因此d选择不依赖Decimal对数误差。得到d=47/63/93，v约1535.66/3575.49/3101.24，B=2048/4096/4096。

| n | R | 联合坏预算(近似显示) | (1)<1？ |
|---|---:|---:|---|
| 10³² | 4 | 1.12195197 | 否 |
| 10³² | ceil(v) | 430.864585 | 否 |
| 10³² | ceil(v²) | 744251.980 | 否 |
| 10⁴² | 4 | .000274286 | 是 |
| 10⁴² | ceil(v) | .245223208 | 是 |
| 10⁴² | ceil(v²) | 1025.19564 | 否 |
| 10⁶⁰ | 4 | 2.55449×10⁻¹³ | 是 |
| 10⁶⁰ | ceil(v) | 1.98109×10⁻¹⁰ | 是 |
| 10⁶⁰ | ceil(v²) | 6.92520×10⁻⁷ | 是 |

“否”只表示预算不够。48轮实际四bank在10³²参数下已经构造成功；本轮另以F₂消元精确检查其六对同端点union rank均22，即满足新增条件(2)。这是现成具体family，不与预算1.12195矛盾。

对网格a=1，以(7)取保守坏点帽：10⁴²、R=ceil(v)的坏点/√n≤.00306690；10⁶⁰、R=ceil(v²)≤.00001924188。10³²中R=ceil(v)或ceil(v²)的粗坏点帽已达到全网格规模，约2.40139√n；这些失败/不够强的结果均保留于JSON。未假设每个有限n都足够大。

脚本还对K=255的逆对称toy lag集合，精确检查模距离、闭坏区间计数、含相位断点的floor/ceil关系；该toy仅检验组合步骤，不冒充新GF证书。

复现：

```sh
/opt/homebrew/opt/python@3.15/bin/python3.15 probe.py
```

仅使用标准库Fraction/Decimal和已有48轮manifest的四组基，不生成大型R的labels或不可承受的全集。严格性来自有理预算/计数和上述证明；Decimal只是显示。最终范围：joint GF family存在、joint-safe有限增长网格剔除及该特定凸混合family的坏点费用；连续尺度和一般C_n仍未闭合。

# 第48轮：四个真实 GF bank 的同查询二阶接口

结论：已构造并逐个精确认证四个 d47、b11 bank，保留同一 Haar 查询、同一 δ、p=.99、n=10³² 和 M40 Fourier 补偿。目标三尺度的真实有限 source 线性扰动相关系数分别被认证包含在

| L | 任意两个不同 bank 的 Pearson 相关系数（向外9位十进制） | 共享 B 方差占总基线方差（数值参考） |
|---|---|---:|
| 10/9 | [0.129939147, 0.130043816] | .129996662 |
| 7/6 | [0.210814227, 0.210942039] | .210884050 |
| 11/9 | [0.280925645, 0.281048845] | .280992029 |

前一列来自精确有理 HH、已接受谱恒等式及尾界、Gauss–Walsh 的有理 HB 包络。端点用 Fraction 精确向外取整，不把 NumPy 三角函数小数当区间端点。后一列只是实际有限谱数值评价。两者均为 H+B 或线性 d 的二阶信息，不是 log-score 协方差、Gaussian 迁移结论或弱型下界。

## 1. 构造、冻结和实际成本

`MANIFEST.json` 含模多项式、本原元、四组有序基、种子及全部个体认证统计；字段 `banks[i].basis` 可插拔。域为 `GF(2^47)`，模 `X^47+X^5+1 = 0x800000000021`，g=2。完整域认证沿用46轮已认证结果，复制保存在 manifest；不生成 K=140737488355327 个状态。

先实际执行 d8,b2 小例与 d47,b8 成本例，保存 `COST.json`；之后才写 `FROZEN_PLAN.json`，固定四种子

```text
35347851986945, 35347851991041, 35347851995137, 35347851999233.
```

每 bank 最多16候选，取首个个体合格者，不按跨表秩择优。四组均第0候选成功，随后全部按显式 basis 重验，总约8.02秒。每次 C 预分配约67MB；顺序运行四bank，不同时保存四套大比率数组。

每 bank 的精确检查均为：rank11；含平方多重集的2,096,128个无序乘积全异；4,188,162个非平凡有序比率全异；S∩gS为空；g^±1、g^±2均不在S。算法采用复制自46轮、哈希相同的 uint64 C checker，只排序 O(B²) 项，不做 O(B⁴) 事件枚举。

共同参数按已接受定义取

\[
v=\frac{9p^2n}{29K^2}\approx1535.6608350596682,
\quad\epsilon=K\delta=\frac1{64v},\quad\delta=\frac1{64Kv}.
\]

Fraction 精确核验 ε<1/B；四个正宽峰来源保持质量1且有共同正背景。代码还核验
`1−p−(20p/(3K))H_40>0`，故原密度的严格正性成立。

## 2. 实测秩：44联合秩不是二阶计算的必要条件

对每个 L 及 θ∈[0,1)，两端点的粗格偏移为

\[
k_\pm(\theta)=\lfloor\theta\pm KL/2\rfloor\pmod K.
\]

各 trace 块的域行是 `g^{k_±} beta_{b,i}`。按使偏移改变的有理 θ 边界分段，用域乘法和 F₂ 消元计算秩；迹配对非退化使该域行秩就是标签映射秩。边界集是 Haar 零测集，不影响积分。

冻结尺度为 10/9、7/6、11/9、整数1及1±1/K、1±2/K。共17个 θ 相位区间、102个 bank 对的四端点联合秩：44计84例，43计6例，整数退化22计12例。亏1的六例已逐项保存在 `RANKS.json`，未被忽略或换基消除。

然而全部408个**跨 bank 单端点对**的秩都是22；所有非整数尺度、每个 bank 自身两端点秩也都是22。二阶 H_A H_B 展开只有四个单端点乘积，所以这些2b秩足够给精确二阶式；无需假设4b=44。

## 3. 精确 θ 分段积分及非零纤维修正

记单端点标签函数

\[
f_y(r)=\operatorname{clip}\!\left(\frac{r-(y+1/2)/B+\epsilon/2}{\epsilon},0,1\right)-r,
\quad 0\le r<1.
\]

令均匀标签均值 `m(r)=B^{-1}Σ_y f_y(r)`，零标签函数为 f₀。每个跨bank单端点标签映射秩22，故全域状态上两个标签独立均匀；实际 source 只遍历非零域，给任意两端点函数 F,G 的精确计数

\[
E_{j<K}FG=\frac qK E_{y,z}F(y)G(z)-\frac1K F(0)G(0). \tag{1}
\]

零态项不能删除。对一个 bank 内两个不同端点，测得秩22时式 (1) 同样适用。Haar 下每个 H 的均值0，所以后述二阶矩就是协方差。

### 3.1 单端点的标签方差积分

直接在 ramp 的两端点 `U_y±ε/2` 分段，平方积分给

\[
\int_0^1 f_y(r)^2dr
=\frac13-U_y+U_y^2-\frac\epsilon6+\frac{\epsilon^2}{12},
\quad U_y=(y+1/2)/B.
\]

这里未平滑阶跃的平方积分为1/3−U_y+U_y²；将阶跃换成线性 ramp 的局部改变量是 `−ε/6+ε²/12`。中点平均满足 `E U_y=1/2`、`E U_y²=1/3−1/(12B²)`，故

\[
\int E_y f_y^2=\frac16-\frac1{12B^2}-\frac\epsilon6+\frac{\epsilon^2}{12}.
\]

定义居中单 ramp 周期函数

\[
f_w(t)=\operatorname{clip}((t-1/2+w/2)/w,0,1)-t,
\quad t\in[0,1),\quad 0<w<1.
\]

将 r 分成 B 个标签格可准确写出

\[
m(r)=B^{-1}f_{B\epsilon}(\{Br\}).
\]

因而

\[
\int m(r)^2dr
=\frac1{12B^2}-\frac\epsilon{6B}+\frac{\epsilon^2}{12}.
\]

两式相减，得到本轮使用的精确量

\[
\boxed{\operatorname{var}:=\int(E_y f_y^2-m^2)
=\frac{B^2-1}{6B^2}-\frac{\epsilon(B-1)}{6B}.} \tag{2}
\]

### 3.2 两端点均值差和零态差

令

\[
J(s,w)=\int_0^1[f_w(\{t+s\})-f_w(t)]^2dt.
\]

程序 `J` 只在 ramp 端点、其平移及周期连接处分段。每段差函数为线性；若端值为 a,b、段长ℓ，则平方积分准确为 `ℓ(a²+ab+b²)/3`。所以全程 Fraction 有理运算，不作数值求积。

两个物理端点的粗局部差为 KL。利用上一节 m 的周期表示，以及 f₀ 与 f_ε 只差一个平移及常数，得

\[
\boxed{j_{\rm mean}=\int[m(r_+)-m(r_-)]^2d\theta
=B^{-2}J(BKL,B\epsilon),}
\]
\[
\boxed{j_{\rm zero}=\int[f_0(r_+)-f_0(r_-)]^2d\theta
=J(KL,\epsilon).} \tag{3}
\]

J 第一参数按模1解释。这个压缩只用小于十个 ramp 段，不枚举2^44个标签输出。

### 3.3 真实 HH 的最终公式

所有非整数的测得秩条件给，对任意不同 bank A,B，

\[
\boxed{HH_{AB}=\frac qKj_{\rm mean}-\frac1Kj_{\rm zero},}
\]
\[
\boxed{HH_{AA}=\frac qK(2\operatorname{var}+j_{\rm mean})
-\frac1Kj_{\rm zero}.} \tag{4}
\]

前式将四项全域乘积组合成均值差的平方；后式把同bank两个单端点二矩减去交叉项，得到两倍(2)加(3)。最后两式的零态均为全零 trace 标签，故同一个 j_zero 同时适用于所有bank。

实际数值参考如下，精确 Fraction 原值在 `COVARIANCE.json`：

| L | HH_AA | HH_AB（六个 bank 对相同） |
|---|---:|---:|
| 10/9 | .333329903478 | 3.95520975×10⁻⁸ |
| 7/6 | .333329915251 | 5.13258552×10⁻⁸ |
| 11/9 | .333329921138 | 5.72127340×10⁻⁸ |

另以 GF(2⁸)、两个b2标签块、255个真实非零态和每态的 Fraction θ 分段积分，独立直接求出 HH_AB、HH_AA，两者与(4)逐分数相等。详见 `SMALL_COVARIANCE_CHECK.json`。小例只检验秩到二阶积分接口，不宣称该小例有目标Sidon参数。

## 4. 共享 Fourier 项、严格有理包络与相关系数

保持同一真实共享项

\[
B_L(x)=\sum_{m=1}^{40}c_m\sin(\pi mL)\cos(2\pi mx),
\quad c_m=\sqrt{2b_m/3}.
\]

Haar 正交性准确给

\[
BB(L)=\int B_L^2=\frac13\sum_{m=1}^{40}b_m\sin^2(\pi mL).
\]

直接 NumPy 有限和给三目标值 .0498064465218、.0890793837581、.1302670548184。严格包络不依赖这些小数：已接受的定义为
`F_M=1+Σ_{m≤M}b_m sin²(πmL)`，确实**包含常数1**。在 I=[10/9,11/9]，`F_∞=γL²`、γ=27/29，且 `0≤F_∞−F_M≤48/(29M³)`。故

\[
\frac{\gamma L^2-1-48/(29M^3)}3\le BB(L)
\le\frac{\gamma L^2-1}3, \tag{5}
\]

下端再与0取大，均为有理数。

绝不假设 HB 独立。固定 θ 时 B_L 的循环 DFT 均值0，因为40<K。由46/47轮的 Gauss–Walsh 协方差估计、`||H||_W≤24`，得每 bank

\[
|HB_A|=\left|\int H_A B_L\right|
\le\frac{\sqrt q}K\,24\,V_L,
\quad V_L=\sum_{m\le40}c_m|\sin(\pi mL)|.
\]

为得到纯有理包络，取

\[
\bar\lambda=\frac{\lceil\sqrt q\rceil}K,
\quad\bar V_L=\min\left\{\frac{40}9,
\frac{440}{63}\operatorname{dist}(L,\mathbb Z)H_{40}\right\},
\quad e_L=24\bar\lambda\bar V_L. \tag{6}
\]

第二项来自 `c_m≤20/(9m²)`、`|sin(πmL)|≤πm dist(L,Z)` 和 `π<22/7`。`ceil√q` 用整数 isqrt 加1获得；q=2^47不是平方数。所有量都是确切有理数。

由(4)至(6)，设 BB∈[b_−,b_+]，则真实 G_AA 和 G_AB 满足

\[
G_{AA}=HH_{AA}+BB+2HB_A\in[D_-,D_+],
\quad D_\pm=HH_{AA}+b_\pm\pm2e_L,
\]
\[
G_{AB}=HH_{AB}+BB+HB_A+HB_B\in[C_-,C_+],
\quad C_\pm=HH_{AB}+b_\pm\pm2e_L.
\]

代码精确检查 D_->0，所以 `sqrt(G_AA G_BB)` 也在[D_−,D_+]。正分子下端除D_+、负分子下端除D_-；正分子上端除D_-、负分子上端除D_+，得到符号正确的有理相关区间，再与[−1,1]相交。表首十进制端点由这些 Fraction 向外9位取整。

这些同样是线性扰动
`d_{A,L}=p(H_A+B)/(KL)` 的 Pearson 相关界，因为共同正比例系数消去。`n_scaled_linear_d_variance_baseline_numeric` 只是 n·a_L²(HH_AA+BB) 的数值参考，尚未包含未知但有界 HB；不是 n坐标 log-score 的精确方差。

## 5. 整数节点与邻域

在 L=1，物理区间长度为一整周期，所有bank Y=1，H=B=0，方差0。因此相关系数**未定义**，输出 `null`，没有赋值0。

在 L=1±1/K、1±2/K，测得非整数同bank秩22、跨bank单端点对秩22，且 KL 是整数，故(3)两项恰为0，HH_AB=0。HH_AA约.333329863926；共享BB数值约2.91×10⁻²⁸或1.17×10⁻²⁷。此处不使用 I 内的 plateau 恒等式，而用正交性给有理
`0≤BB≤bar V_L²/2`，再结合(6)得到合法相关包络；高精度有理端点在JSON，显示到9位为[−.000000001,.000000001]。

这说明“共同 Fourier 背景必然使不同bank高度相关”不是这些真实源与节点邻域的普遍事实；同时节点外确实存在约.13至.28的共享相关。不能由此选择任意Gaussian相关参数、宣称Gaussian尾迁移、或推出连续尺度极大事件的收益。

## 6. 复现和证明范围

在本目录按顺序执行：

```sh
/usr/bin/clang -O3 -Wall -Wextra -std=c11 gf_check.c -o gf_check
/opt/homebrew/opt/python@3.15/bin/python3.15 build_banks.py
/opt/homebrew/opt/python@3.15/bin/python3.15 ranks.py
/opt/homebrew/opt/python@3.15/bin/python3.15 covariance.py
```

`build_banks.py` 重放保留既有冻结计划及成本记录，核对固定参数；运行用时和输出哈希可能随重放更新。所有源码、最小manifest、原始个体检查、rank、精确小例与二阶结果均在本目录，不改旧轮。

已证：这四个具体source在八个固定尺度上的个体GF条件、精确rank、HH、HB/BB与线性相关系数包络。条件式(4)只在实测秩前提下使用；整数情况另算。未证：所有连续尺度的rank性质、log-score协方差、joint Gaussian/高尾迁移、多表最大事件下界及 C_n 新log因子。

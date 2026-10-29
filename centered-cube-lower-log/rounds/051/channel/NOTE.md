# 第051轮：真实来源的固定半弧 gate channel

本轮先实际计算45轮所定义的正密度，而非自由 Gaussian 通道。小源的全部峰窗端点作分段 Gauss 积分；随后只对d47实际来源作32个有理查询 pilot。得到可核对的真实 posterior 数值，并独立证明任意固定半周期弧的严格几何约束。

## 1. 实际测量

来源始终是

\[
\rho=1-p+p\mu_{K,\delta}+(p/K)D_M,
\quad D_M=\sum_{m\le M}\pi m c_m\cos(2\pi mx).
\]

d10使用原labels、p1/4、M32；d12使用原labels、p.99、M40。两者B8、v1024、w64，均严格正；d12的三角下界为.00310418。它们是可分辨的小来源诊断，B8不符合大参数增长族B≥v，不冒充其有限维证书。

两个固定macroclass是半弧w₀=1_[0,1/2)，以及w_c=1_([0,1/4)∪[3/4,1))。所有尺度共用同一个ρ及固定gate。完整小源结果见RESULTS；d12摘要为：

|L|r²，w₀|Corr(p_L,logY_L)，w₀|r²，w_c|Corr(p_L,logY_L)，w_c|
|---|---:|---:|---:|---:|
|10/9|.00851706|.00312248|.00851882|.246686|
|7/6|.01586978|.00156606|.01587400|.358860|
|11/9|.02325802|.00152668|.02326526|.461343|

d12的π分别为.499879120879及.500276204254。p_L的query均值与π吻合约1e−16。4点与8点分段Gauss的r²差约1e−18、相关差不超过约5e−14；这只是嵌套计算的一致性，**不是严格求积误差证书**。断点包含全部原峰两端、平移端点及gate边界，未用不分辨正宽峰的普通均匀网格。本机NumPy的longdouble实际只有53位尾数（已记入JSON）；小源最窄δ约3.73e−9仍远大于其机器epsilon，高d超小δ从不走这一路径。

小源Var_Q(p_L)围**真实class质量π**积分。高d pilot则把有限32查询的Y加权均值当作经验中心，计算经验variance；两者不能混称同一个认证量。32点高d的Corr波动明显，例如L7/6 centered half约.661，不能据此给全Haar相关界。高d pilot仅验证可执行的真实posterior接口；下面的r²定理比这些pilot估计更可靠。

## 2. 精确几何通道与协方差恒等式

在圆周上取周期区间核k_L（L>1时保留绕圈的重数），定义真正联合概率律

\[
Q_L(dx,dy)=k_L(x-y)\rho(y)\,dx\,dy.
\]

其source边缘为ρ(y)dy，query边缘为Y_L(x)dx，其中Y_L=A_Lρ。令C=w(y)∈{0,1}、π=∫ρw，N_L=A_L(ρw)，则

\[
p_L(x)=Q_L(C=1\mid x)=N_L(x)/Y_L(x),\quad
E_Qp_L=\pi.
\]

因此门控后的真实source密度为ρw/π，其查询平均为N_L/π=Y_Lp_L/π，质量精确为1；没有跳过π归一化费用。

对任意平方可积query函数g(x)，tower property与对称区间核给

\[
\operatorname{Cov}_Q(C,g)=\operatorname{Cov}_Q(p_L,g)
=\int\rho(y)(w(y)-\pi)A_Lg(y)\,dy.
\]

另有Cov(C,p_L)=Var(p_L)，故

\[
r^2=\frac{\operatorname{Var}_Qp_L}{\pi(1-\pi)}\in[0,1],\quad
\operatorname{Corr}_Q(C,p_L)=r,
\quad\operatorname{Corr}_Q(C,\log Y_L)
=r\operatorname{Corr}_Q(p_L,\log Y_L).
\]

这里r与后一个相关系数均由同一真实来源和核决定，不能像自由Gaussian探路那样独立指定。

若固定这一来源后取iid坐标，对(C,p_L,logY_L)的中心化和作通常有限方差多元CLT，极限Gaussian满足上述协方差关系；非退化时它构成C→p_L→logY_L的Gaussian Markov链。这只是**固定来源的CLT**；不提供随n变化源族的统一中偏差、罕见尾或多尺度迁移。

## 3. 任意固定半弧的严格r²约束

以下对本来源所有标签、任意固定平移的半周期弧，以及所有L∈I=[10/9,11/9]同时成立。置V=Σc_m≤40/9，原函数

\[
F_\rho(t)=t+\frac pK\left[\mathcal P_\delta(t)
+\frac12\sum_m c_m\sin(2\pi mt)\right],\qquad F_\rho'=\rho.
\]

因为|Pδ|≤1，令

\[
D=\frac pK(1+V/2),\quad
\kappa=\frac p{KL_{min}}(2+V),\quad
|F_\rho(t)-t|\le D,\quad |Y_L-1|\le\kappa.
\]

半弧π满足|π−1/2|≤2D。长度L<3/2的物理窗口最多交到两个lifted半弧段；每段的source质量与长度差≤2D，因此若p₀=A_Lw为uniform来源posterior，

\[
|N_L-p_0|\le4D/L_{min},\qquad
\|p_L-p_0\|_\infty\le e_p:=\frac{4D/L_{min}+\kappa}{1-\kappa}.
\tag{1}
\]

此处κ<1；所有d47实际参数都满足。因为p、p₀∈[0,1]，以及∫Y=1，

\[
\left|\operatorname{Var}_Qp_L-\operatorname{Var}_{dx}p_0\right|
=\left|\int Yp_L^2-\pi^2-(\int p_0^2-1/4)\right|
\le\kappa+2e_p+3D.
\]

再用π(1−π)=1/4−(π−1/2)²，得到完全显式的

\[
\boxed{|r^2-r_0^2|\le
\mathcal E:=\frac{\kappa+2e_p+3D+4D^2}{1/4-4D^2},}
\tag{2}
\]

只需D<1/4。该证明没有做高d数值积分，也不依赖GF随机性、Gaussian近似或δ相对于浮点分辨率的大小。

uniform来源下令ℓ=L−1∈(0,1/2)。从整周期后的剩余ℓ窗口与半弧的交长度O出发，E O=ℓ/2、E O²=ℓ²/2−ℓ³/3，故

\[
r_0^2(L)=\frac{\ell^2}{L^2}(1-4\ell/3).
\]

它在I上递增，因为导数为
$2\ell(1-2\ell-2\ell^2/3)/(1+\ell)^3>0$。
三个目标值是23/2700、1/63、76/3267。脚本用**独立的Fraction分段线性平方积分**核对这三个闭式，而不只代入同一公式。

d47、p.99、V上界40/9的(2)由Fraction精确求值为

`1232096962299927464329893605746397909 / 871122859317548356107121824529730061012906945456`

约1.41438e−12，并严格核验76/3267+E<1/36。因此

\[
\boxed{r<1/6\quad\text{对d47实际源、任意固定半弧及全部 }L\in I.}
\]

沿增长源族，E=O(1/K)→0。故此半弧channel的r=1/3自由参数不符合真实几何。结论**不覆盖其他macroclass、位置依赖门或全部尺度**，也没有证明Corr(p,logY)统一远离1。低d数值的Corr随gate平移大变，仍应保留同一核的联合协方差费用。

## 4. 实现、数据与复现

Fρ同时给ρ区间积分；将gate的有限区间分段并周期延拓即可得N_L。高d计算先用Fraction作整数粗格/局部phase及所有端点运算，再用046 oracle 的70位mpmath原函数，不将global float用于解析约7.23e−20的δ。高d只有32个分层有理查询，不是n10³²维MC或Haar误差认证。

归档输入仅含labels：`inputs/source_d10.npz`745字节、`inputs/source_d12.npz`2133字节；原NPZ包含的atoms/defects不需要复制。d10 labels已与045 accepted快照逐项一致；原文件、labels字节和本轮快照哈希见inputs/MANIFEST。d12不声称能由缺失basis重建，因此直接保留这2133字节必要数据。

支持`CUBE_INPUT_DIR`覆盖输入目录、`CUBE_ORACLE_DIR`指向046 oracle目录；mpmath依046接口使用正常安装或`CUBE_MPMATH_VENDOR`。默认不查找旧numerics目录。高d manifest可用`--manifest`显式传入46轮参数/基。复现：

```sh
python3 channel.py
python3 channel.py --high-pilot --manifest ../round046_gf/MANIFEST.json
```

最小归档：`channel.py`、`RESULTS.json`、`HIGH_PILOT.json`、`NOTE.md`、`inputs/MANIFEST.json`及两份labels NPZ。最初成本PILOT已被正式RESULTS覆盖，不列入归档。外部依赖是046 source_oracle.py及高d manifest，哈希在结果中；普通运行还需NumPy和mpmath。70位是任意精度浮点，不是区间算术；严格(2)及r<1/6认证只用已证明的解析界和Fraction。

没有得到新门控弱型下界、实际极大事件概率或C_n改进。root的固定软分层计数比较处理不同的标签模型；这里测得的是明确固定空间半弧的真实通道，不能混为一类。

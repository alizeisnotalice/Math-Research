# Astra014：Cauchy 模型的完整 Gaussian shape 面积分预算

本批证明同一 16 条未单位化 body 方向的任意 `0<alpha<2` 形状积分接口；将已锁 I022 的五类概率重新代入 `beta=1/2`，得到

\[
 I_1\le U_1=7.2224324476276829882\ldots\,10^{-5}
 <\frac{36112162239}{500000000000000}.
\]

这覆盖全部时间、完整正单位面和退化形状尾，不依赖截时网格。用此预算使条件成本低于 32，需要另证
\[
 g>5(U_1/32)^{4/5}=0.00015198362980537152698\ldots.
\]
这个等式使用精确 Gamma 表达式 `U1`；采用归档的有理向上端时，充分门槛略升，应重新代入该端点。没有新全域 `g`，没有 `C5<32`。下一步低阶 Cauchy Young 证书有研究价值，但本批不生产矩表或空间证书。

## 1. 模型、常数及完整时间积分

固定 Lebesgue 测度、Fourier 逆变换因子 `(2pi)^-5`、中心 cube `[-1,1]^5`，以及
\[
 V=\{(1,\pm1,\pm1,\pm1,\pm1)\},\qquad
 \widehat p_t(\xi)=\exp\{-t\sum_{v\in V}|v\cdot\xi|^\alpha\}.
\]
以下 `0<alpha<2`，`beta=alpha/2`，`r=(5-alpha)/2`，`q=5/(5-alpha)`。取独立正 stable 变量 `S_v`，Laplace 变换为 `exp(-lambda^beta)`。写
\[
 A=\sum_v S_vvv^T,\quad T=\sum_vS_v,\quad P_v=S_v/T,\quad M=A/T.
\]
`A,M` 几乎处处正定：所有权正且方向张成五维。`M` 的对角元均为 1。条件 Gaussian 协方差为 `2t^(1/beta)A`，故
\[
 p_t(x)=(4\pi)^{-5/2}t^{-5/(2\beta)}
 \mathbb E\left[(\det A)^{-1/2}
 e^{-x^TA^{-1}x/(4t^{1/\beta})}\right].
\]
非负 Tonelli 及 `u=t^(-1/beta)` 给，对每个 `x!=0`，允许两边为无穷，
\[
 G_\alpha(x)=\kappa_\alpha\mathbb E\left[
  (\det A)^{-1/2}(x^TA^{-1}x)^{-r}\right],\qquad
 \kappa_\alpha=\frac{\beta\Gamma(r)}{2^\alpha\pi^{5/2}}.
 \tag{1}
\]
尺度分离恰产生 `T^-beta`，不是其它幂。负 stable 矩的 Laplace 积分给
\[
 Z_T:=\mathbb ET^{-\beta}=\frac1{16\Gamma(1+\beta)}.
\]
令形状分布 `R_beta` 是密度 `T^-beta/Z_T` 下 `M` 的推前分布；这里不声称 `T` 与 `M` 独立。记
\[
 d=\det M,\quad Q=x^TM^{-1}x,\quad B=x^T\operatorname{adj}(M)x=dQ,
\]
\[
 H_\alpha(M,x)=d^{-1/2}Q^{-r}=d^{(4-\alpha)/2}B^{-(5-\alpha)/2}.
\]
于是
\[
 \boxed{G_\alpha=c_\alpha\mathbb E_{R_\beta}H_\alpha},\qquad
 c_\alpha=\frac{\Gamma((5-\alpha)/2)}
 {2^\alpha\pi^{5/2}\,16\Gamma(\alpha/2)}.
 \tag{2}
\]
特别 `c1=1/(32pi^3)`，`c_(3/2)=3/(2^(15/2)pi^(5/2))`；程序分别检查了这两个归一化。下文用 `b=beta/5` 表示负矩指数，以大写 `B` 表示伴随二次型，避免同名混淆。

## 2. 轨道、球面与 cube 面的准确系数

对固定正定 `M`，`H_alpha^q=d^(-q/2)Q^(-5/2)`，对 `x` 为负五次齐次。球面 Gaussian 径向积分给
\[
 \int_{S^4}Q^{-5/2}\,d\sigma=|S^4|d^{1/2},\qquad |S^4|=8\pi^2/3.
\]
对于任意负五次齐次非负函数 `F`，径向投影把球面积积分化成 cube 边界通量积分：
\[
 \int_{S^4}F=\sum_{\text{10 faces}}\int_F F(x)\,(x\cdot n)\,dS
 =\sum_{\text{10 faces}}\int_F F(x)\,dS.
\]
最后一步因单位 cube 每面 `x.n=1`。这个等式也可用两种极坐标参数对任意正厚度齐次 annulus 积分证明，因此没有关于 `F` 在面上连续延拓的假设。

令 `O` 为全部 signed permutation 群。每个原方向在群作用下成为另一个方向或其负向，因此 `R_beta` 在 `M->UMU^T` 下不变。对正面片 `F+={(1,y):y in[0,1]^4}`，有限轨道平均把全部 `160=10*16` 个面片等权覆盖。凸性和上述积分给
\[
 \left\|\frac1{|O|}\sum_{U\in O}H_\alpha(M,Ux)\right\|_{L^q(F_+)}
 \le \left(\frac{\pi^2}{60}\right)^{1/q}d^{(1-q)/(2q)}
 =\left(\frac{\pi^2}{60}\right)^{1/q}d^{-\beta/5}.
\]
再对 `R_beta` 用 Minkowski；下一节先证明右侧有限，因而该步对非负截断取极限也合法。记 `J_beta=E_R d^-b`，则
\[
 \boxed{I_\alpha:=\int_{[0,1]^4}G_\alpha(1,y)^q\,dy
 \le c_\alpha^q\frac{\pi^2}{60}J_\beta^q.} \tag{3}
\]
这里保留了形状轨道；没有对单个 `M` 的一个面片擅自使用完整球面积。

## 3. 完整尺度抵消及 flag 熵优化

因 `5b=beta`，有准确恒等式
\[
 J_\beta=\frac{\mathbb E[T^{-\beta}(\det A/T^5)^{-b}]}{Z_T}
 =16\Gamma(1+\beta)\mathbb E(\det A)^{-b}. \tag{4}
\]
采用 I022 已锁的纯整数几何分解：
\[
 \det A=16^5\sum_f\theta_f^0F_f,\qquad \sum_f\theta_f^0=1.
\]
每个 `F_f` 为五个正线性 stable 因子的乘积；其中两个为 singleton，其余三个概率向量保存在 `flag_types.json`。五块的支持互不相交且覆盖全部 16 个变量。原整数身份为
\[
 \det A=\frac1{60}\sum_{J,X,Y}(S_iS_j)X_SY_SL_S;
\]
每个五基含十个初始二子集，两次商方向选择分别重复三次、两次，产生系数 `60`。原批的整数 Schur 补及平行类分组已经核验该身份，本批消费其锁定结果，不把旧 `beta=.75` 的概率幂或尺度常数搬来。

对任意正 flag 权 `w_f`、和为 1，加权 AM–GM 与随后普通 Hölder 给
\[
 \mathbb E(\det A)^{-b}
 \le16^{-\beta}e^{bD(w\Vert\theta^0)}
 \prod_f(\mathbb EF_f^{-b})^{w_f}. \tag{5}
\]
Hölder 的权是 `w_f`，不是在形状概率下错误地认为各块独立。只有原 `S` 产品概率内同一 flag 的五个块互相独立。

若块的概率向量为 `p`，其线性 stable 和的 Laplace 指数为 `lambda^beta sum p_i^beta`，故
\[
 \mathbb E\left(\sum_i p_iS_i\right)^{-b}
 =\frac{\Gamma(1+b/\beta)}{\Gamma(1+b)}
   \left(\sum_i p_i^\beta\right)^{-b/\beta}.
\]
设 `R_f(beta)=product_blocks(sum_i p_i^beta)`，singleton 因子为 1。负矩均有限，因而式 (5) 无退化可积性缺口。写
\[
 Z_\beta=\sum_f\theta_f^0R_f(\beta)^{1/\beta},\qquad
 w_f^*=\theta_f^0R_f(\beta)^{1/\beta}/Z_\beta.
\]
准确 KL 恒等式
\[
 D(w\Vert\theta^0)-\frac1\beta\sum_fw_f\log R_f
 =D(w\Vert w^*)-\log Z_\beta
\]
证明这是全部正权的最优值，而非有限网格。合并 (4)–(5) 得
\[
 \boxed{J_\beta\le
 16^{1-\beta}\Gamma(1+\beta)
 \left[\frac{\Gamma(6/5)}{\Gamma(1+\beta/5)}\right]^5
 Z_\beta^{-\beta/5}.} \tag{6}
\]
五类型合并时保留类型内原 `theta^0` 比例，故 type 权优化准确等价于全部 4560 个 flag 权优化。

## 4. beta=1/2 的新准确类型计算与统一成本

五类型质量依次为 `(1/8,3/16,3/16,1/8,3/8)`。直接对原概率取平方根，得到
\[
 R_f(1/2)^2=(81/8,27,36,36,64),
\]
\[
 Z_{1/2}=2661/64,\qquad
 w^*=(27,108,144,96,512)/887.
\]
因此式 (6) 成为完全显式的
\[
 J_{1/2}\le4\Gamma(3/2)
 \left[\frac{\Gamma(6/5)}{\Gamma(11/10)}\right]^5
 (2661/64)^{-1/10}
 =2.04469055665376492848\ldots.
\]
一次 384-bit Arb 批准确包围该 Gamma 表达式并产生有理外端，见 `budget.json`。同一五概率输入的 `beta=3/4` 控制得到熵优化预算 `I_(3/2)<=1.2967687959491083460...e-5`。

统一传递成本是
\[
 \mathcal B_\alpha(g)=(5/\alpha)^q I_\alpha/g^q.
\]
它由齐次 Green 的层蛋糕积分和截断水平 `c=alpha*g/5` 得到；`g` 必须是同一模型、同一单位 cube 的认证 essential 下界。令 `h=g/c_alpha`，预算成本也可写成
\[
 (5/\alpha)^q(\pi^2/60)(J_{\beta,\rm bound}/h)^q.
\]
所以比较两个 alpha 不应只比 `g`：

| 固定模型 | 完整 `I` 预算 | 达到成本 32 所需 `g` | 所需 `g/c_alpha` |
|---|---:|---:|---:|
| alpha=1 | 7.22243244763e-5 | 1.51983629805e-4 | .150798287252 |
| alpha=3/2，同一 flag 熵预算 | 1.29676879595e-5 | 1.11757723479e-4 | .117965752134 |

这些 `I` 是上界表达式，不能从它们的大小推出真实 `I` 的大小；必要 `g` 仅是使用该预算的充分门槛。预算不构成模型排除或新 C5 上界。

## 5. S1/S2 发散与完整 Lq 有限性同时成立

Cauchy 模型不能沿用四维的“仅在线上发散”。若 `x=av+bw`，`a,b!=0`，补三方向组成五基，将其余十一条 Cauchy 变量限制在 `O(t)` 区间；此事件概率为固定正数。前两个基密度下界 `ct`，三个补基密度下界 `c/t`，卷积和 Jacobian 因而给
\[
 p_t(x)\ge c_x/t,\qquad G_1(x)=\infty.
\]
同理一条方向线上 `p_t>=c t^-3`，原点由缩放发散。`11110=(11111+1111(-1))/2` 就在二方向平面内。它未被数值运行。

这些 `S1/S2` 集与四维面相交是低维集合，不妨碍 (3)–(6) 的完整 `L^(5/4)` 面积分有限性；我们的证明直接覆盖它们附近的全部可积尾，没有删除管状邻域。新点
\[
 x_*=(1,7/8,3/4,5/8,1/2)
\]
在所有 pair plane 之外，因为 `av+bw` 每坐标绝对值只可能是 `|a+b|`、`|a-b|`，而 `x_*` 有五种不同绝对值。其 Green 值有限性也可由至少三条非小 Cauchy 参数和五基卷积的局部 `p_t<=C t` 支配证明，具体同已锁 `upper/round001_singular_strata.md`。这里不作整个线外连续的错误断言。

## 6. 唯一非奇异点的实际诊断及边界

预登记经过主修订，仅对 `x_*` 用 seed `202610260141`、32768 个 scrambled Sobol 样本。Kanter 生成 16 个 `beta=1/2` 正 stable 变量；直接对式 (1) 已准确积分掉的时间表达式求样本均值。未设时间截断，未运行其它点或种子。实际结果为

- Green 样本均值 `0.0002996073565858675`；代入本模型预算的**单点样本成本代理**为 `13.699521644680962`。
- `E exp(-S_1)` 的样本检查为 `.36787803734136515`，准确目标 `e^-1=.36787944117144233`。
- `E T^(-1/2)` 的样本检查为 `.0705167278445729`，准确目标 `1/(16Gamma(3/2))=.07052369794346955`。
- 全部 32768 样本成功；无非正 Cholesky、无特征值失败、无丢弃、无裁剪、无矩阵修复。详细条件数、分解残差及八个块的描述数据保存在 `diagnostic.json`。
- 最小样本特征值 `1.2834917478e-11`，最大条件数 `3.8956230206e11`，最大 Cholesky 绝对重构残差 `3.3306690739e-16`；小残差不等于这些病态样本的前向误差已受控。八块均值范围 `.00028500756` 至 `.00031289747`，最大单样本贡献 `.0327184000`。

这是有限浮点/QMC 诊断，没有随机误差或舍入误差包围；块均值不是置信区间。一个点的有限均值不是 `g`，也不是该点的严格上下界。该数值支持继续低成本解析证书筛选，既不能证实达到所需全域门槛，也没有排除 alpha=1。

## 7. 复现与下一步的确切缺口

在本独立 clone 中运行：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper014/verify_budget.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper014/direct_shape.py
```

Python/SymPy 的五类型平方根简化为有理数，Arb 384-bit 的向外常数比较属于有限严格核验；NumPy/SciPy 样本不属于证明。源码、登记、最小输入、两个结果和本报告为白名单，无大缓存、无新图矩。完整几何输入锁定来源与 SHA256 在 `flag_types.json`，本批文件在 `manifest.json`。

后续如获新批登记，可用 alpha=1 的
\[
 H=d^{3/2}B^{-2}\ge3\sqrt d\,Z^2-2BZ^3,
 \qquad Z=a(1+\eta W).
\]
必须重算 beta=1/2 的图矩，不能复用 beta=3/4 的数值。一个需避免的符号陷阱：直接对 `sqrt d` 作两节点三次 Hermite 是**上界**，因为其四阶导数为负；可对 `d^(-1/2)` 作三次 Hermite 下界再乘 `d`，或先用 `d` 乘该函数的单切线下界。`B` 与形状的联合矩必须保留。新 f 为线性时空间每个 `X_i=x_i^2` 的次数至多 2，正面最多 `3^4=81` 个 Bernstein 系数；这只是结构规模说明，本批未运行该证书。

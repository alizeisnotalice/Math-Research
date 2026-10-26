# 153：时间单调性增强的锚点曲率上界

**所提公式成立。** 固定五维Lebesgue、16条未单位化body Cauchy方向、权1共同时间、逆Fourier因子(2π)^−5，BBᵀ=16I。对每个有限ε>0、实锚点c，置
\[
H_\epsilon(c)=\int_\epsilon^\infty p_t(c)dt,\quad
K_\epsilon(c)=\int_\epsilon^\infty p_t(c)t^{-2}dt,\quad
Q_\epsilon(c)=K_\epsilon(c)/(16H_\epsilon(c)),\quad
A=1/(16\epsilon^2),\quad r=H_{2\epsilon}(c)/H_\epsilon(c).
\]
则
\[
\boxed{Q_\epsilon(c)\le\min\{\tfrac23A,\ A(\tfrac7{10}-\tfrac8{15}r)\}.}
\tag{153.1}
\]
若有同点同模型可靠端点0≤L≤H2ε(c)、Hε(c)≤U，U>0，则
\[
\boxed{Q_\epsilon(c)\le q_{\rm new}
=\min\{q_{\rm old},\ \tfrac23A,\ A(\tfrac7{10}-\tfrac8{15}L/U)\},}
\tag{153.2}
\]
其中q_old可省略；如保留，它必须是同一anchor、同ε、同物理Euclidean二次式的已资格上界。线性项仍是原Hε的梯度；不混合H2ε梯度。

本轮只证明lemma、准备精确有理接口并核新符号/人工比值。0真实H、梯度、log、点曲率或角点计算；不消费公共g115新Delta，不更新g/Cn。

## 1. 时间概率、归一化与反向Chebyshev

029已证正Gaussian混合，A_G为不依赖t的正定随机矩阵：
\[
p_t(c)=(4\pi)^{-5/2}t^{-5}E[(\det A_G)^{-1/2}
\exp\{-c^TA_G^{-1}c/(4t^2)\}].
\]
所以β_c(t)=t⁵p_t(c)非减。此处指数5来自物理维数；不是16个Cauchy参数的维数。任何常数Fourier/密度归一化均在后面的时间比值中抵消，ε没有被首坐标scale16t重新缩放。

对0<a<b≤∞，以t^−5dt归一化成概率。β非减、f(t)=t^−2非增，双重积分的协方差恒等式给
\[
\frac{\int_a^bp_t(c)t^{-2}dt}{\int_a^bp_t(c)dt}
\le\frac{\int_a^bt^{-7}dt}{\int_a^bt^{-5}dt}.
\tag{153.3}
\]
正时间p_t(c)>0，且p_t≤Ct^−5，故各分母正有限、β有界；b=∞仍可用有限区间极限/DCT。无需时间导数或任何Cauchy绝对矩。

Kε正有限，因为Kε≤ε^−2Hε<∞。Hε与H2ε均处处正有限，且Hε−H2ε是直接定义的正有限band积分，不出现∞−∞。因此0<r<1，所有分段条件概率均合法。

## 2. 三个准确常数与质量混合

直接积分t^−5、t^−7得
\[
\frac{\int_\epsilon^\infty t^{-7}}{\int_\epsilon^\infty t^{-5}}=\frac2{3\epsilon^2},\qquad
\frac{\int_\epsilon^{2\epsilon}t^{-7}}{\int_\epsilon^{2\epsilon}t^{-5}}=\frac7{10\epsilon^2},\qquad
\frac{\int_{2\epsilon}^\infty t^{-7}}{\int_{2\epsilon}^\infty t^{-5}}=\frac1{6\epsilon^2}.
\]
最后一项也等于(2/3)(2ε)^−2。对两个真实时间band分别用(153.3)，其后验权分别为1−r与r，故
\[
\frac{K_\epsilon}{H_\epsilon}
\le\epsilon^{-2}\left[\frac7{10}(1-r)+\frac16r\right]
=\epsilon^{-2}(\frac7{10}-\frac8{15}r).
\]
再除16即得(153.1)。两band交点2ε是时间Lebesgue零集，没有重复收费或遗漏尾域。

另一个有用的一致性是**真实r≥1/16**：对全尾t^−5概率，β与1_{t≥2ε}同向单调，所以其β加权尾概率不小于原概率(2ε)^−4/ε^−4=1/16。原点p_t(0)=p_1(0)t^−5给等号，并使Qε=(2/3)A；故2/3是全空间统一常数的准确边界，不需要重算p_1(0)。

在r=1/16处，7/10−(8/15)r恰为2/3。因此对任何有效λ=L/U，
\[
\min\{\tfrac23A,A(\tfrac7{10}-\tfrac8{15}\lambda)\}
=A[\tfrac7{10}-\tfrac8{15}\max\{1/16,\lambda\}].
\]
如果保存质量下端很弱（λ<1/16），统一时间界仍提供改进。

## 3. 与128/129的全空间锚点Jensen接口

075固定t的物理坐标log-Hessian界给
\[
p_t(c+u)\ge p_t(c)\exp\{\nabla\log p_t(c)\cdot u-|u|^2/(16t^2)\}.
\]
按固定锚点概率p_t(c)dt/Hε(c)做指数Jensen，128/129已严证
\[
\log H_\epsilon(c+u)\ge\log H_\epsilon(c)
+\nabla\log H_\epsilon(c)\cdot u-Q_\epsilon(c)|u|^2.
\tag{153.4}
\]
于是可直接把Q换成(153.2)的任意可靠上端。它对所有实u成立，不是沿空间线段冻结一个局部Hessian后验。首坐标位移为0时，|u|²就是四个物理尾坐标的平方和。

新因子1/16来自BBᵀ=16I与Taylor的1/2；t^−5时间权来自维数5。方向若单位化或时间权改变，不能原样复制常数。此处没有改变模型或物理度量。

Hε处处正、光滑，空间求导与完整正时间积分的交换仍由075有界score和t^−5尾支配。S1/S2上完整G可为∞，但本证明只用正时间p、H和K，全部奇异集合、零/重复投影、闭边界均保留。u=0给平凡等号；ε=0、非有限ε或非正U不属本lemma。

## 4. 有向端点与可用范围

因为r=H2ε/Hε≥L/U且(153.1)的混合项对r递减，代入L/U是安全上界。**分母必须是Hε上端U，不能换成其下端。** 由严格正band质量，可靠实际输入自动满足0≤L<U。L=0是允许的“没有非平凡质量证据”；L≥U则与正时间模型不相容，应查源，不clip成1。

若有多个真实质量下界，可取max。例如EPS128锚点有原H64整叶祖先证明，以及同点已资格H64缓存下端，可取两者max。不同anchor也可提供已严证的H64(c)下界，但必须实际覆盖同一个c；数值同点缓存必须核完整精确坐标/ε/模型/hash。只用这些H64质量，不混入其gradient或Delta。

仅由非减凸代数得到的混合项范围是[A/6,7A/10]；与统一项取min后在[A/6,2A/3]。但再与**任意已资格q_old**取min，不能断言最终q_new≥A/6，更不能断言q_new≥A/4。

对EPS64、L32=0，A=256，得到统一q≤512/3。对EPS128，A=1024；只要L64/U128>27/32，混合项就小于256。人工λ=9/10给q=5632/25<256。

因此未来153 consumer必须删除旧129专用的`256<=q`或`A/4<=q`断言，**不能把新值clip回旧下限**。可核的新统一范围是0<q_new≤2A/3，并按(153.2)和源身份逐字段验证。原129范围仍正确描述其旧公式；不修改冻结旧报告或旧证书。q非负即可保持角点安全函数凹性，16闭角点规则不需要256下限。本轮不运行该角点规则。

若采用129同一λ的简单cutoff界A(1−3λ/4)，其与新混合项之差准确为A(18−13λ)/60>0（0≤λ≤1）；这只是公式族的代数比较，不是重跑某个旧点的曲率。

## 5. H下界的域限制与当前gstar

已读017主接纳及主证明：全部111201闭叶仅使用H64或H128本身；正时间cutoff单调性和保存正余量给**整个闭F上H128≥gstar**。`primary_global_acceptance.json`明确区分：H128 floor只在complete_F，完整G floor才通过G的−4齐次性扩展到Q。

- EPS128用L64=g105，必须有该c的真实H64祖先/同点等证据。旧完整g105本来含H128叶，且旧H64存在严格障碍；不能把global G或新的global H128 floor反向变成整个F的H64 floor。
- EPS64无H32证据时取L=0即可；若将来有真实同点H32下端，另按同一资格使用。
- future EPS256，2ε=1/128，若anchor c∈F，可直接用017已接纳的H128≥gstar作为L。若anchor在F外，不能仅凭G在Q的下界使用这个L；Hε不作固定ε齐次扩展。
- future H256的Delta/梯度必须仍来自H256同一anchor的合格记录；H128只提供质量下界。此处不授权或启动H256。

本轮不计算gstar的新比值或成本，当前正式C5≤30.121409937411276只引用主已接纳结果；153没有新g或Cn结论。

## 6. 有限新符号与接口检查

本轮新登记后运行本机Python/SymPy，PID **68840**、exit0、约 **0.0988秒**。核全尾/两band的二阶逆时间矩、8/15混合系数、1/16尾比例、旧新公式差及27/32阈值；这些不同于150已完成的一阶逆时间矩，没有重跑150 symbols。

`penalty.py`只提供精确Fraction算术接口，不读取任何H缓存、target、Delta、gradient或角点。三个纯人工接口例（EPS64无L、EPS128λ=9/10、保留一个更小的假设旧上界）和四个非法参数拒绝例已通过。人工例只核接口，未声称这些端点或旧q在真实模型某点成立。

实际点q、H、gradient、log与角点求值数均为0。全部准确符号结果、注册和源条件另由manifest锁定；旧150/139及其他冻结材料不变。只新增`upper_curvature153/`，无下界、Astra、MATLAB、提交或推送。

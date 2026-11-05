# 负剪切宏观孔与微孔的共同参数桥

2026-11-05。冻结送审候选。只研究明确的负共同剪切宏观构造及其内部微孔，给真实有限输入之间的比较；不识别C_n绝对九阶系数，不把三角收益和微孔常数未经结算直接相加。固定n≥3。

## 1. 同一设计、同一完整负算子

记u=t^(1/n)、v=2u、a=(2u)^n、g=(1+2u)⁻¹，
B_n=binom(n,2)、K_n=B_n n²/27。
κ(u)在固定正紧区间，Λ(u)有界。所有隐含常数可依n及这些预选区间，不依u、micro间距、孔数或最后安全C；涉及C的最终质量余项会明确注明。令
\[
 e=\kappa u^3,\quad T=(1+e)I-eJ,\quad d=\det T,\quad
 w=ug+\Lambda u^6,\quad\rho=1-2w,\quad C_w=[-w,w]^n,\quad M=a\rho^n.
\]
每点使用原规则s=w−||x||∞、r(x)=ρ+s。置b=(2w)^(n−2)，
\[
 F_w(s)=\frac a b[(\rho+s)^n-\rho^n],\quad f=F_w'(0),\quad
 \delta=e/f=O_n(u),\quad
 R_{w,V}(x)=b^{-1}\sum_{k\ne0}|(Tk+(C_w\setminus V))\cap Q(x,\rho+s)|,
\]
\[
 \mathcal B_w(V)=\partial C_w\cup\{R_{w,V}\ge F_w(s)\},\quad
 U_0(w)=\mathcal B_w(\varnothing),\quad
 U_{j+1}(w)=\mathcal B_w(U_j(w)),\quad K(w)=U_4(w).       \tag{1}
\]
全部测试集使用同一设计M及原半径，不将|C_w\V|重设成自身质量。设旧宏观参数安全：
\[
 D_\Lambda-|K(w)|>0,\qquad D_\Lambda=(2w)^n-a\rho^n.     \tag{2}
\]
下文给非空供应，尤其R072的负宏观族满足该条件。

定义
\[
 L_0=\frac{B_n\kappa^3}{8n^2g},\qquad
 \Lambda^\mu=\Lambda-aL_0+Cau.                         \tag{3}
\]
C是最后选择的大固定安全常数。所有新w^μ、ρ^μ、M^μ、B^μ、K^μ均按(3)重新构造，然后才放置微孔。没有将旧host搬到新孔，也没有在旧算子中代入新质量。

本稿目标为对每个小实u构造负剪切纯宏输入与负剪切宏微输入，使其完整halo比值满足
\[
 \mathscr R_{\rm neg,micro}-\mathscr R_{\rm neg,macro}
 =\left[\frac{B_n\kappa_0^3}{256n}+o(1)\right]av^{n+6},
 \qquad \kappa(u)\to\kappa_0.                          \tag{4}
\]
κ₀=8n/3时系数为2K_n。三维再与同一旧参数的正纯宏族直接比较，得到(6+o(1))av⁹。

## 2. 负算子的完整几何与任意目标小常数

为避免把正剪切合同改名，先列负号实际接口及理由。若任意原见证碰Tk+C_w，则
||Tk||∞≤||x||∞+r+w=1。又
||Tk||∞≥[1−(n−1)e]||k||∞，故小e时全部||k||∞≥2整数被排除。
对{-1,0,1}^n中的k，T k=(1+e)k−e(Σk_l)1；混合正负同时严格进入要求Σk_l>1且<−1，不可能。只有
k=±1_A、|A|=j≥2可有正质量：
支撑中心±[1−(j−1)e]，非支撑中心∓je。支撑1仅边界接触。

每个非支撑方向，对任意源x∈C_w和目标原宽z_l∈[−w,w]，
|x_l−(Tk)_l−z_l|≤2w+ne<ρ。故非临界位移被**精确整宽包含**消去，不是小误差。
单目标支撑cap宽≤(j−1)e，归一化总R≤C_ne²。因此所有坏点s≤C_ne²/f=O(u⁴)≪w−ne。
在最近面符号及全局最浅轴q的图中，d_l=s+y_l、y_q=0，
\[
 0\le s\le H_w(y):=w-\max_l y_l.
\]
活动目标支撑方向必须与源最近面相同，穿透前缀为
[0,((j−1)e−y_l)_+]，不依s；所有非支撑目标宽完整保留。
真实输出是全部目标质量先相加后的高度
min{H_w,F_w⁻¹(R_V(y))}。odd支持、三面图及clip全部保留。

固定一个目标点打开Tonelli，支撑j的可见源底面积，当q在支撑内、外时分别≤C_ne^(j−1)w^(n−j)、C_ne^jw^(n−j−1)。
除bf，核常数≤C_nδ(e/w)^(j−2)及多一个e/w。有限清单求和，得到对任意可测V,W
\[
 |\mathcal B_w(V)\triangle\mathcal B_w(W)|
       \le L_u|V\triangle W|,\qquad L_u\le C_n\delta=O_n(u)<1/2. \tag{5}
\]
这里不要求V,W是规则轮廓，也不要求微孔帽内均匀；新增odd支撑进入同一常数。F′≥f及真实clip的1-Lipschitz性不会放大核。

令V₀=(2w)^(n−2)e³/f。名义j2帽的全最浅轴积分给
|U₀|=2B_nV₀+O_n(uV₀)：固定一朝向，pair内两个最浅轴的积分各为b e³/(2f)，合V₀；两相反朝向合2V₀，再对B_n个轴对求和。pair外最浅轴多e/w，源clip相对O(eδ/w)，高支持j≥3多至少(e/w)²，F曲率相对O(eδ/ρ)，都更小。
由(5)，|U₄−U₀|≤C_nu|U₀|，故
\[
 \frac{|K(w)|}{a}
 =nL_0u^5g^{n-2}+O_n(u^6),\qquad
 f=4nu^2g[1+O_n(u^5)].                                \tag{6}
\]
这是负号全图自身的领先体积，不假设它逐点等于正号坏集。

有限目标cap交量随前缀端点连续，F严格增，故每张图高度连续；K闭，边界为有限图分界和高度图，n维零测，且|int K|=|K|。K为偶层，有B_w(K)⊂K。

## 3. 跨w稳定：仅pullback及体积，不是物理对称差

用最近面坐标
\[
 x_l=\sigma_l(w-d_l),\quad \sigma_l\in\{\pm1\},\quad0\le d_l\le w
\]
将集合指标函数嵌入固定有限符号直和⊕_σ ℝ_+^n，在各自[0,w]^n域外补零，记
\[
 (I_w A)(\sigma,d)=\mathbf1_{[0,w]^n}(d)\,
           \mathbf1_A(\sigma_1(w-d_1),\ldots,\sigma_n(w-d_n)).
\]
这是拉回后的指标函数，不是对原物理坐标集合做L¹比较。
除坐标中分面零测集外每点表示唯一，各chart雅可比绝对值1。因此
\[
 \big||A|-|A'|\big|
 \le\|I_wA-I_{w'}A'\|_1.                               \tag{7}
\]
不声称原x坐标下|A△A′|同阶：角面整体位移Δw可能产生大至a u⁷|ΔΛ|的物理对称差。预算只需(7)的体积差，所有新host将在新K重做。

固定u、κ，仅改变有界Λ。对相近w,w′，有
\[
 \|I_w\mathcal B_w(V)-I_{w'}\mathcal B_{w'}(V')\|_1
 \le C_n\delta_{\max}\|I_wV-I_{w'}V'\|_1
       +C_na u^5\,\frac{|w-w'|}{w_{\min}}.             \tag{8}
\]
完整理由如下。

* 在共同intrinsic target坐标内，负支撑j的临界前缀[0,((j−1)e−y_l)_+]及其面对符号独立w；非临界±je仍精确吞入整宽。输入差先由(5)的同一可见底面核收费，系数C_nδ_max。
* 其他target坐标的完整域[0,w]改变，差条至多占完整名义cap的C_n|Δw|/w_min；b⁻¹=(2w)^(2−n)的改变同样只有此相对比例。临界前缀始终远离中分面，不引入新的奇目标。
* 在0≤s≤C_nu⁴，F_w(s)/s为正，|∂_w log(F_w(s)/s)|≤C_n/w；这由ρ=1−2w及有限乘积展开直接得出。F′≥f，所以固定名义R的高度改变≤名义高度×C_n|Δw|/w。s=0由连续延拓，不作0/0操作。
* 源底域与H_w(y)也改变。若临界y_l≤(n−1)e，则它远离w，不触该端域；只有自由底坐标的端条变化。各名义帽对该自由变量独立，所以Fubini按长度收费C_n|Δw|/w。clip变化只在自由坐标近w的带；涉及该坐标的目标帽在此恒零，其余帽与它无关，沿该坐标积分仍是同一名义体积比例，不是将|Δw|乘整个底域。
* 这些名义全图体积和≤C_nV₀=O_n(a u⁵)，有限支持j及全部最浅轴求和只改C_n。因而各参数成本合为(8)第二项，不要求测试删除集边界规则。

从空集开始对固定j=0,…,4迭代(8)，缩小u使C_nδ_max<1/2，得到
\[
 \|I_wU_j(w)-I_{w'}U_j(w')\|_1
 \le C_na u^5|w-w'|/w_{\min}.
\]
由w′−w=u⁶(Λ′−Λ)，
\[
 \boxed{\big||K(w')|-|K(w)|\big|
       \le C_na u^{10}|\Lambda'-\Lambda|.}             \tag{9}
\]
这不是跨w运输micro集合的许可，只比较旧／新宏观体积。所有常数在预选紧参数盒上一致。

## 4. 固定新host的一步保护引理

先固定任意当前负算子B_w和偶层K=U₄(w)，设H⊂K是可测微孔，所有微孔的短见证位于K中。定义
\[
 V=K\cup\mathcal B_w(K\setminus H).
\]
则
\[
 H\subset K\subset V,\qquad
 \mathcal B_w(V\setminus H)
 \subset\mathcal B_w(K\setminus H)\subset V,\qquad
 |V|\le|K|+L_u|H|.                                    \tag{10}
\]
前两个关系由反单调性；最后用B_w(K)⊂K和(5)，
V\K⊂B_w(K\H)\B_w(K)，体积≤L_u|H|。
目标上包恰是
(C_w\V)∪H=C_w\(V\H)，因此(10)保护全部宏观见证，包括负odd支撑，不只是一个指定角。既有host始终在K⊂V，宏孔以后避开整个V，所以所有短见证host不变；没有先强迫H却允许宏孔重新进入host的问题。

## 5. 在新K内放置全部短见证host

以下全部用(3)的新参数，省略μ仅表示当前新对象。先构造K=U₄(w)，再放H；不使用旧K的任何点态包含。
K闭、边界零测且|intK|=|K|。K的canonical代表可包含整∂C_w；本文不对它的字面周长作任何小量断言，只用n维零测边界与每个固定u的内部体积极限。固定当前u，取
\[
 t'=t(1-u)<t,\quad v'=2(t')^{1/n},\quad a'=(v')^n=a(1-u).
\]
对dyadic间距ℓ，只保留q∈ℓℤ^n且完整host Q(q,ℓ)⊂intK的点，置
\[
 \omega=\frac{\ell v'}{2(1+v')},\qquad
 H=\bigcup_q(q+[-\omega,\omega]^n).
\]
host可以重叠；计量使用不交基本胞q+[-ℓ/2,ℓ/2)^n。令其有效总体积为L_ℓ。每个有效基本胞在K中，所以L_ℓ≤|K|；K内距补集大于3ℓ/2的点由有效胞覆盖，所以L_ℓ→|K|。可选第一个满足ℓ≤u²⁰且L_ℓ≥(1−u)|K|的间距。有限q清单给
\[
 |H|=\frac{a'}{(1+v')^n}L_\ell=a|K|[1+O_n(u)],\qquad
 \frac{|H|}{a}=n aL_0u^5g^{n-2}+O_n(a u^6).            \tag{11}
\]
这里H不用在任意宏观帽里均匀分布，只使用实际总质量。

每个微孔点z∈q+[-ω,ω]^n取明确短半径
\[
 r_\mu(z)=\ell-\omega-\|z-q\|_\infty.
\]
小u时r_μ≥ℓ−2ω≥2ω>0，含自身完整微孔。任意其他格孔的某轴中心距至少ℓ，目标点到z该轴距≥ℓ−ω−||z−q||=r_μ，所以至多零体积接触。又Q(z,r_μ)⊂Q(q,ℓ−ω)⊂Q(q,ℓ)，整个见证位于host。缺陷率≤
\[
 \left(\frac{\omega}{\ell-2\omega}\right)^n
       =(v'/2)^n=t'<t.                                \tag{12}
\]
这是每个点、包括微孔边界点的严格原阈值见证。

不同宏胞的外包C_w平移有固定正间隔：||T k||∞≥1−(n−1)e>2w。微见证在自身C_w内的host，不能碰其他宏胞的孔。应用(10)取V，未来宏孔全部在C_w\V，故host内没有宏孔。host不随保护操作移动或缩小，(12)因此在最终输入仍成立。

## 6. 新设计预算：先定范围，再选安全常数

设旧参数Λ在[-B,B]，L₀在固定正区间[L_*,L^*]。先把允许Λ区间扩大到
[-B−L^*−1,B+L^*+1]，在该盒中固定§2–3和领先量(6)的全部常数。它不依赖尚未选择的安全C。之后选择C，最后要求Cu≤min(1,L_*/2)。于是
|Λ^μ−Λ|≤a(L^*+1)，新参数留在预选盒内，且Λ^μ<Λ。这个顺序避免把误差常数反过来依赖未知C。

由binomial展开，D_Λ/a=(g+Λu⁵)^n−(g−2Λu⁶)^n，故在固定盒上一致
\[
 \partial_\Lambda(D_\Lambda/a)
       =n u^5g^{n-2}+O_n(u^{10}).
\]
由(3)，
\[
 \frac{D_{\Lambda^\mu}-D_\Lambda}{a}
 =-n aL_0u^5g^{n-2}+n C a u^6g^{n-2}+O_n(a u^{10}).  \tag{13}
\]
(9)只作体积比较，给
\[
 \big||K^\mu|-|K|\big|\le C_na^2u^{10}.                \tag{14}
\]
新H按新K重新构造，(11)给|H|=n a²L₀u⁵g^(n−2)+O_n(a²u⁶)。
一步保护的代价≤L_u|H|≤C_nu|H|=O_n(a²u⁶)。因此
\[
 \begin{aligned}
 D_{\Lambda^\mu}+|H|-|V|
 &\ge D_\Lambda-|K|
       +a^2u^6[n Cg^{n-2}-C_0].
 \end{aligned}                                       \tag{15}
\]
先前预选盒使C₀与C无关；O(a²u¹⁰)已包含其中。g→1，先选大固定C，再小u，右边严格为正。旧安全余量没有被拿来购买新增微孔；即使旧余量极小，(15)仍新增正空间。

特别M^μ−|H|>0，因为M^μ≈a、|H|=O(a²u⁵)。
(15)等价于
|C_w^μ\V|>M^μ−|H|。在此开集内取有限内部紧盒总体积超过该目标，再连续分位切得紧宏孔A，质量恰为M^μ−|H|。总孔
\[
 P^\mu=A\cup H,\qquad |P^\mu|=M^\mu
\]
是实际新设计质量，不是在旧质量上额外加H。

对每个宏孔点x∈A，其原新半径r=ρ^μ+s含整个本胞P^μ，且外国P^μ由
(C_w^μ\V)∪H=C_w^μ\(V\H)上包。由(10)，x∉B^μ(V\H)，所以
\[
 M^\mu+\sum_{k\ne0}|(Tk+P^\mu)\cap Q(x,r)|
 \le M^\mu+b^\mu R^\mu_{V\setminus H}(x)
 <M^\mu+b^\mu F^\mu(s)=ar^n.                          \tag{16}
\]
所有odd支持和非临界位移都在同一个完整算子中，不是多个局部帽分别安全。
微孔点沿用(12)，其host在K^μ⊂V，不受A污染。
因此宏、微的全部实际孔点都有原阈值严格见证。

## 7. 质量增量与一般维数的真实micro收益

T、e、u、κ在此次Λ变化中固定，故determinant d也固定。
由ρ^μ−ρ=2aL₀u⁶−2Ca u⁷、ρ=g+O(u⁶)，
\[
 M^\mu-M
 =2n a^2g^{n-1}L_0u^6+O_{n,C}(a^2u^7).               \tag{17}
\]
d/(d−m)对m导数为d/(d−m)²=1+O_n(a+u⁶)，所以周期比值的增量同为(17)。
若κ→κ₀，
\[
 \frac{M^\mu-M}{av^{n+6}}
 \longrightarrow\frac{nL_0(0)}{32}
 =\frac{B_n\kappa_0^3}{256n}.                         \tag{18}
\]
κ₀=8n/3时恰2K_n。
虽然新增保护的体积费用O(u)|H|在设计预算中需支付，但它由+Ca u那项支付后，对实际M只损O(a²u⁷)，比(17)晚一阶。未假定微孔均匀或保护免费。

## 8. 共同有限化，涵盖全部半径及所有host

旧负宏观族取任意质量M的紧孔P⊂C_w\K；新混合孔是上节P^μ，质量M^μ。
分别按同一负T周期重复。所有宏观半径<1，所有微半径≤ℓ≤u²⁰<1，且||T⁻¹||∞<2。
统一取
\[
 N=\lceil u^{-(n+12)}\rceil,\qquad \Omega_N=T([-N,N]^n)
\]
并保留全部边界切孔。周期积分精确给两输入质量
(2N)^n(d−M)、(2N)^n(d−M^μ)。

除O_n(N^(n−1))个边界宏胞外，所有长、短见证都保留。边界孔总量≤C_naN^(n−1)，常数不依micro间距、micro孔数量或host数量：每个宏胞的全部总孔质量≤C_na。不能逐微孔数额外乘一个边界损失。

输入密度点几乎处处进入halo。对凸Ω_N外任一中心、任意正半径，分离半空间及中心方体对称性给输入占比≤1/2，小t时不在严格halo。因此两组有限输入的**完整极大算子比值**分别为
\[
 \mathscr R_{\rm neg,macro}
       =\frac d{d-M}+O_n(a/N),\qquad
 \mathscr R_{\rm neg,micro}
       =\frac d{d-M^\mu}+O_n(a/N).                    \tag{19}
\]
由于a/N=O_n(a²u¹²)=o_n(av^(n+6))，得到(4)。
统一N按一般n的micro尺度选，不把R072的u⁻¹²有限化在大n错误地称为足够深。
每个充分小实u均有这样的有限输入，不限子列；宏孔、H、所有host及P^μ都在取N之前已固定。

## 9. 三维triangle与micro在一个实际族中的兼容

从R072的旧共同参数出发：正纯宏偶层有余量≥c a u⁸，κ→8n/3。负U₄按同参数安全，且仍有≥c′a u⁸余量。实际参数供应可取R072源(22)；若正U₂安全，正U₄⊂U₂，R072固定层比较使负U₄也安全。因此(2)的旧负输入不是空假设。

将上面全过程应用于这个旧负族；新负macro+micro用Λ^μ、总M^μ及新K^μ全部重建，旧正纯宏族仍用原Λ、M。
用同一深N有限化，记d_±=detT_±。两个周期比值差有精确恒等式
\[
 \begin{aligned}
 \frac{d_-}{d_--M^\mu}-\frac{d_+}{d_+-M}
 ={}&\frac{M(d_+-d_-)}{(d_--M)(d_+-M)}\\
 &+\frac{d_-(M^\mu-M)}{(d_--M^\mu)(d_--M)} .
 \end{aligned}                                       \tag{20}
\]
第一项才是旧相同M下的triangle，第二项才是同负T下经实际保护后的新质量增量。这是在一个新负混合输入已构造、预算与全部见证已核之后的代数结算，不是把两个独立模型的常数先相加。

n=3时d_+−d_-=4e³，κ→8，第一项/(av⁹)→4；(18)第二项/(av⁹)→2K₃=2。共同有限误差o(av⁹)，所以对实际有限输入
\[
 \boxed{\mathscr R_{\rm neg,micro}
       -\mathscr R_{\rm pos,pure}=(6+o(1))av^9\quad(n=3).} \tag{21}
\]
旧正族本来的九阶安全常数保留在共同参考部分；(21)不是C₃绝对九阶系数。

一般n≥4，(4)仍给相对旧负纯宏的2K_n av^(n+6)收益。若与正纯宏同时写，须保留(20)的**精确triangle项**再加micro收益；不能把triangle的未知高阶余项误写成o(av^(n+6))。

## 10. 依赖、可复核节点与证据范围

本稿使用R072的完整负邻孔／核计数，但§2逐式重建其负号几何；§3另证跨w的intrinsic稳定性，未把正号合同改名。R060第8节作为实际micro模式参照，真正必要的host、不变K、odd支持、归一化、总质量与有限化在本稿全部展开。R073二维的更细保护系数和R074全输入上界不在此证明的依赖中。

| 节点 | 实际接口 | 关键量词／误差 |
|---|---|---|
| 负完整算子 | 每个原r、全部unipolar j≥2、同新设计M | 非临界2w+ne<ρ为精确包含 |
| 小常数 | 任意可测V,W，含foreign H | L_u≤C_nδ=O(u)，不需帽内密度 |
| 跨w | I_w的Jacobian1及共同域外补零 | 只pullback／体积差a u¹⁰ΔΛ，不是物理对称差 |
| 新K与host | 先新Λ→新K→H与全部hosts→V | 旧host不搬；V⊃K保持所有短见证 |
| 设计预算 | D_new+H−V而不是旧M+H | leading抵消，保护费a²u⁶由最后C au支付 |
| 新实际质量 | 宏质量M_new−H、总孔质量M_new | 全长见证(16)，全短见证(12)严格 |
| 有限化 | 所有比较族共同N=u^−(n+12) | a/N=o(av^(n+6))，不依micro数 |
| n3联合 | 先真实构造，再用同一分式(20) | (6+o)av⁹仅是实际族比较 |

已实际读取本机Mathematica 15.0.0 ARM的以下脚本及成功日志：

| 文件 | SHA256 |
|---|---|
| scratch/r077-computation/joint.wl | 5e7623c4453677266050b3df078f80bec94d556b2e50b831b0966fb963300341 |
| joint.log | 34aa4d41bd43ad99658210fe0508312ec41cfb2d0585149390f5e7a071135350 |
| joint_oldmass.wl | c90f86dadb8128bf1a9debf6856e01e1aefc8ec38b12b1639de02f5730e8dfdd |
| joint_oldmass.log | b8a49a112715c189e5e7704449f4de10c1c3d9ebb892b2c5bec2e4248227fe39 |
| joint_fraction.py | f4f536f2012c7f6262f730ab12776c5f0497653561535dfe72be53eb5287f59d |
| joint_fraction.json | 6556a825dbde06962d0fcf7a01c8e39ef271868e75f434ba41e5dab4b6d3843f |

general-n预算导数与首项抵消残差0；micro归一化与2K_n之差0。
第一份joint脚本用det@M_new+mass@d_+的等价分拆；第二份专门对正文(20)的det@M_old+mass@d_-核验，精确残差0，n3两个归一极限4、2，合6。
独立Fraction三个有理u样本核同一代数分式和分母／质量正性；其Λ=C=0仅为scalar样本，**不认证宏观安全预算，不是本文实际输入的有限参数示例**。小数只为阅读精确有理数结果。
首次joint运行误用了受保护的D并有替换优先级错误，joint-first-invalid.log原样保留，不作通过证据；修复后两份成功WL日志如上。这些计算不能代替§2–6的实际集合证明。

源材料按实际字节绑定：
R060实际micro原证明SHA c7bcd28338661fa782029c87c883467d7ea6dd812d3b16301a02c74598cc9ccc；
R072 signed比较SHA f8f655fda9be633b93aa2d5bfdd20960176ae93bb1f30f74aefb587915f8b91d。
读取时仓库HEAD e4b3773e839c1bdcdaa70c05c6c2767394eb009c，工作树干净，未修改仓库／旧源。
本稿在完成独立审查前不作确认定理发布；即便本桥成立，其结论仍只是实际构造族间的收益，不是C_n绝对九阶定理或一般micro最优性。



# R053 · R052准确保留tuple的显式门槛（条件于R051全链接受）

## 0. 结论身份与固定输入

**Conditional draft。** 本报告已把下列固定tuple的全部eventual条件量化到L=log n≥2^32；但只有R051泛参数机制全链正式接受后，才可激活其operator结论。没有从标量基点提前宣布新的C_n界。

固定准确参数为

a=73/8，s=293/16，d=1023/2^20，ε0=2^(−20)，
εJ=εS=εT=εD=1/1024，εA=1/64。

读入的是R052 `results.json`的第一条完整分数certificate，不从显示小数重构，也不重跑有限列表。其独立保留点审计已经接受精确算术身份。

若R051全链接受，则对全部L≥2^32有

C_n≤Φ:=F(A_n(t),B_n(t)−γn)<T−2^(−16)n，                       (1)

T=nL−nlogL，x=(1−σ)L，σ=sd²，t=n/x。本文另外证明新Φ对旧完整cost的比较，保持各旧定理自己的有效域。

冻结R051源为geometry `upper/round051/flexible_budget/PROOF.md`，SHA-256
`aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6`。

本轮不改该源及任何旧证明，不做旧网格、生产扫描、七维、MATLAB或Git/网络操作。

## 1. 完整准确budget，不省scalar buffer

以R051原公式重建，P0=3579/4096，qmin=1−3d/4，并有

u_lo=3579/8200，u_hi=8796093022208/7674592741893，M*=4，
m*=64449/134315936，b*=11148919/38310400，
Γ=718536680631/20582788938137600，
γ=(1−εD)Γ=735063024285513/21076775872652902400，
σ=306632997/17592186044416。

所有domain条件a>3、s>2a、0<d≤1/16、d<1/(32a)、σ<1/2、0<εJ,εA,εD<1、εS,T>0、b*>0都用Fraction重新核对。M*=max(4,2u_hi/qmin)的分支确实为4。

完整净预算是

c_raw=γ−σ/(1−σ)
 =6468305130094850733787980147/370780099533215320770594052505600
 >c=2^(−16)。                                                  (2)

必须保留σ/(1−σ)，不能把不同x下的正cost遗漏，也不能在此用粗糙2σ代替该准确价格。固定εD给实际缺额Γn与所用γn之间明确正余量εDΓn。

## 2. L0=2^32的全部几何域

记L0=2^32，n=e^L，h=ε0/(Ld)，ζ=n^(−3)，κ=1+1/t−ζ，θ=e^(−aLd²)，k=nhθ。目标r=1+u、sourceρ=tq，u,v=1−q∈I=[d/4,3d/4]，matching为|r−q−d|<h。

Taylor给n≥L²/2及n≥L^8/8!。因σ<1/2，x∈[L/2,L]，t≥L/2≥16且t<n；d≤1/16给r<2。

在L0核验h/d<1/8，此后递减；且

ζt/h≤2d/(ε0n²)≤1/2。

所以h+ζt<d/2、κ>1、β<2。固定target的baseline q0=r−d+h与整个bridge [q0−2h,q0]都在(1/2,1)内；R(q)=q+d−h≥r−2h>1，R≤r<2，ρ≥t/2>4r。每个actual source q本来也有q≥qmin>1/2。

于是g=d−h+ζtq满足0<g≤d、g≤R/2；全部exact J/A/B公式、局部Lipschitz derivative及pinned source-face公式的正域明确成立。coordinate bound只用这些域及x∈[L/2,L]，不依赖原d=2^(−24)或某个旧soft exponent。

## 3. Good probability至少1−εA，而非沿用旧5/8

baseline R=r时，free target坐标准确EJ=1，且

|E logJ|≤6g³/ρ²，Var(logJ)≤12g³/ρ²，|logJ_pin|≤2g/ρ。

因此对每张target面、全部r统一有

|EΣlogJ|≤[24d³L²+4dL]/n，
Var(ΣlogJ)≤48d³L²/n。

用n≥L^8/8!，上包为递减powers

M(L)=8![24d³/L^6+4d/L^7]，V(L)=48·8!d³/L^6。

准确基点证书给M(L0)≤1/32以及

V(L0)≤(3/32)²·εA/3。

Chebyshev遂给baseline log条件|ΣlogJ|≤1/8的坏率≤εA/3。

又K_w=1+Bin(n−1,2d/r)，mean≤1+2nd、variance≤2nd。nd≥2时

P(K_w>4nd)≤8/(9nd)。

基点用n≥L²/2核验该上包≤εA/3，尾部随L递减。

原corner-captured target部分每面比例≤(3/4)^(n−1)：d≤1/16和L≥64保证1/r−1/t>1/2。n≥33时(3/4)^32<εA/3，基点n≥L0²/2远大于33。

三项并集坏率≤εA，故每张面保留≥1−εA。于是准确得到

λ(A)≥(1−εA)n log[(1+3d/4)/(1+d/4)]≥m*n。                      (3)

这里m*使用真实d/[2(1+3d/4)]，没有换成较粗nd/4或5/8比例。

## 4. Matching product及individual J的显式精度

R051的已给coordinate导数界在第2节正域可积分跨kink。K_w≤4nd及n/ρ≤2L给matching log变化≤1280Ldh=5/4096。结合baseline 1/8，

|ΣlogJ|≤517/4096。

用e^(−z)≥1−z与e^z≤1/(1−z)得到原预算准确需要的

P0≤∏J≤1/P0。

此外，piecewise闭式在全部T≤2d给A(T)≤4T/3、B≤g²/(3R)。因此uniform有

|J−1|≤2d/ρ≤4dL/n≤4d·8!/L^7。

基点该量≤εJ，尾部严格递减，故对所有target坐标、全部matched/baseline q，

1−εJ≤J≤1+εJ。                                                 (4)

这是逐坐标uniform bound，不是对n个坐标作一次失败率并集估计。源面ell_*≥1−d/2、1≤N_i≤2及(4)给准确u_lo n≤U≤u_hi n。

inverse-row概率下E(−logG)≤d+2Ld²≤3Ld²（Ld≥1已核）。θ=e^(−aLd²)故retention≥1−3/a。完整soft row下界为(1−3/a)u_lo k，上界2u_hi k/qmin；加column≤4k即M*k。

## 5. Source cap与τ两个不同的指数尾

原source corner质量≤nd e^(−Ld/32)。除k给

Cap/k≤2^20 Ld² e^(−ηL)，η=d/32−ad²>0。

Taylor4给全尾上包

Cap/k≤24·2^20 d²/(η^4L³)。                                    (5)

L0的准确Fraction比较小于εS；尾部按L^(−3)递减。因此原generic最终row下界R_λ≥[(1−3/a)u_lo−εS]k在该显式域成立。使用b≤1支付cap删除，未用inverse-weighted mass冒充原质量。

原预算给B0>8n、B0≤nL、A0≤ne^(−σL)，且γ<1，所以τ=√[A0(B0−γn)]为实数。其相对row尺度满足

τ/k≤2^20 dL^(3/2)e^(−νL)，ν=(s/2−a)d²=d²/32>0。

在L0，νL0=1046529/8192>64，prefactor=1023·2^48<2^58。Taylor64于是给

τ/k<2^58·64!/2^384<εT。                                       (6)

对应的全尾power为L^(−125/2)，不是只核验一个指数样本。

(3),(5),(6)及M*给每个原历史的ambient缺额≥m*b*n/M*=Γn>γn。fixed εD的缺额余量完全保留。

## 6. Finite-a、同history及一般L1全部尺度

用mathfrak a表示原small-cube dilation，避免与soft exponent a混淆。对每个mathfrak a≥8n^4，β<2及κ>1给

κ+2/(mathfrak a+1)<αβ，α=mathfrak a/(mathfrak a+1)。

actual target半径<3、actual source半径>4，故所有source先于同一个target；center guard使其实际被bQ_i捕获。source b(z)≤1，所以一次原eligibility给weighted μ-row≤τ，不出现θ^(−1)因子。

在A/P⊂D∞上λ_a=α^(−n)λ准确成立，energy、degree、target mass按正倍率缩放；缺额Γn不会减小。数学证明使用同一实际increment measureμa，不假设它服从完整壳概率分布。

对任意一般f∈L1与全部真实共同半径见证，原紧水平集有限覆盖、selection及witness支付保持不变。固定n和τ，只在所有mathfrak a≥8n^4已成立的数值上界中令mathfrak a→∞，再用内正则性。若R051泛链接受，这给(1)的Cn≤Φ，不要求selector或历史收敛，也没有有界支撑/尺度截断。

## 7. 完整scalar余项小于准确slack

当前buffer为σ=sd²，必须使用

B0≤T+[σ/(1−σ)]n+1。

因此

Φ≤T−c_raw n+nR(L)，
R(L)=1/n+e^(−σL)+2√L e^(−σL/2)。

L0处σL0>2^16。n≥L²/2及Taylor16给

R(L0)≤2^(−63)+16!/2^256+2^17·16!/2^240
 <c_raw−2^(−16)。                                              (7)

三个全尾上包分别按L^(−2)、L^(−16)、L^(−31/2)递减。故所有L≥L0均有Φ<T−2^(−16)n。该结论使用完整正scalar cost，不是仅把gamma当作最终负系数。

## 8. 新Φ与已接受旧完整cost的比较

以下先是公式cost比较；新operator身份仍条件于R051接受。新scalar上包自L≥2^32起成立，比较可取L≥Lcmp=2^40。

R049端点的完整cost使用budget 2^(−40)n；R042使用budget 2^(−31)nL^(−1/3)。后者在L≥2^40时也≤2^(−40)n。其相应x均在[L/2,L]，故两个旧完整cost都有公共下包

Φ_old≥T−n[4/L+2^(−40)]。

Lcmp处4/L+2^(−40)=5·2^(−40)<2^(−16)，以后递减。因此新Φ严格小于这些旧完整cost。

更旧列明的R033/R034/R035/R036/R040固定t_old族，在各自有效域有D_old≤2n√(logL/L)，故完整cost≥T−n[4/L+2√(logL/L)]。该loss在logL>1时递减。Lcmp处logLcmp=40log2<40<7²，故

4/Lcmp+2√(logLcmp/Lcmp)<2^(−38)+14/2^20<2^(−16)。

这比较的是旧完整cost下包，未从display负系数推断完整根式排序。R034准确q取既有证书、理论q<1/8；不再重建小数或跑旧grid。

有效域分别保留：R049算子门槛是L≥2^52，所以作为两个算子界的比较取L≥2^52；R042仍为L≥2^108；R036的logL≥128也不能改成Lcmp即生效。R040未量化的参数N以及R034各Y²门槛各自保留。对R030/R031单位核心全t族，旧V≥U_old−2且U_old=T+(3l+1)n/L+O(nl²/L²)只给eventual V>T；本报告不为其O项虚构显式起点。

## 9. 验证与conditional状态

`verify_threshold.py`直接从R051原budget公式构造全部Fraction，并与R052保留certificate每个对应字段逐一相等；它不重跑R052枚举。脚本只核L0/Lcmp的最小exact gates与source hashes，所有无限尾由上文powers或解析导数覆盖。

在R051泛族full ACCEPT之前，本报告只是conditional threshold/cost theorem，不能将R052 scalar candidate直接当Cn新界。若主链接受且本阈值独审通过，则式(1)在全部L≥2^32激活；旧界仍按各自原有效域取minimum，不被本报告缩窄。没有修改作者源或旧receipt。

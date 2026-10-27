# R049 · R045 endpoint真实operator门槛：log n≥2^52

## 结论与来源

**门槛闭合。** 对所有满足L=log n≥2^52的整数n，已接受R045 bounded soft-inverse-column机制可在端点d=2^(−24)实施，并有

C_n≤F(A_n(t),B_n(t)−2^(−40)n)
 <nlog n−nloglog n−2^(−41)n，                                  (1)

t=n/[(1−2^(−43))L]。

这里是量化全部几何、典型集、inverse row、source cap、τ、有限a及一般L1桥接之后的operator门槛，不是只给scalar基点。没有沿用comparison的2^240门槛，也不声称2^52最小。

冻结机制源：`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md`，SHA-256
`add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9`。

本报告使用该已接受源中的准确coordinate公式及其非渐近常数，逐项证明其适用条件。本轮只写本目录，不改R045源，不做旧网格、生产计算、七维、MATLAB或Git/网络操作。

## 1. 固定参数及全部几何正域

固定

d=2^(−24)，ε0=2^(−20)，σ=32d²=2^(−43)，L0=2^52，
x=(1−σ)L，t=n/x，β=1+1/t，ζ=n^(−3)，κ=β−ζ，
h=ε0/(Ld)，θ=exp(−12Ld²)，k=nhθ，I=[d/4,3d/4]。

目标半径r=1+u、u∈I，source半径ρ=tq、q=1−v、v∈I；matching为|r−q−d|<h。

L≥L0时x∈[L/2,L]、x>1、n=e^L≥L²/2及n≥L^8/8!。所以

t≥L/2≥16，r<2，ρ≥t/2≥8>4r，n≥17，nd≥(L²/2)d≥8，Ld≥1。

并且

h/d=2^28/L≤1/8，
ζt/h=[Ld/(ε0 x n²)]≤2d/(ε0 n²)=2^(−3)/n²≤1/2。                 (2)

这里实际基点h/d=2^(−24)<1/8，故所有需要严格的R>1条件均有余量；并且h+ζt≤3h/2≤3d/16<d/2。κ=1+(1−ζt)/t>1且κ<β<2。

对固定source q，参考目标半径R=q+d−h、source gap g=d−h+ζtq。于是R>1，0<g≤d、g≤R/2，ρ≥4R，ell_*=1−g/(2R)≥3/4。

对固定目标r，baseline q0=r−d+h满足1/2<q0<1，R(q0)=r。整个matching bridge q∈[q0−2h,q0]内，均有q>1/2、1≤R≤r<2、ρ>4R、g≤d、g'=ζt≤d；T=(|w|−(κ−1)ρ)_+≤2d。上述范围即冻结源coordinate式(7)–(9)的全部域条件。

## 2. Baseline log J均值、方差与固定面

在baseline R=r下，冻结源的准确inverse coordinate为

J(w;q)=1+[2B(R,g)−A(R,g,T)]/(2ρ)，T=(|w|−(R−g))_+。

对free W uniform[−R,R]，Fubini准确给EJ=1。坐标界B≤g²/(3R)、A(R,g,T)≤(4/3)T（0≤T≤g）和ET²=g³/(3R)，给

E(J−1)²
 ≤[2g^4/(9R²)+8g³/(27R)]/ρ²≤3g³/ρ²。

同时|J−1|≤g/ρ≤1/8。对|u|≤1/8，|log(1+u)|≤2|u|且log(1+u)≥u−2u²。因此每个free坐标有

−6g³/ρ²≤E logJ≤0，Var(logJ)≤12g³/ρ²。

upper来自log凹性及EJ=1，不假设log均值为零。

固定target面有一个pinned坐标W=±r；其log绝对值≤2g/ρ。其余n−1坐标独立，故S_log=Σ_i log J(w_i;q0)满足统一界

|E S_log|≤[24d³L²+4dL]/n，
Var(S_log)≤48d³L²/n。                                         (3)

这里q0≥1/2使1/ρ≤2L/n。pinned坐标只贡献mean而不贡献variance；没有把全部n坐标当free uniforms。

用n≥L^8/8!把(3)上包为

M(L)=8![24d³/L^6+4d/L^7]，V(L)=48·8!d³/L^6。

两者随L>0递减。纯Fraction基点给

M(L0)≤1/32，V(L0)≤9/8192。

因此Chebyshev统一给

P{|S_log|>1/8}≤V(L)/(1/8−1/32)²≤1/8。                         (4)

这对每张面、全部target r均成立，不是只有某个采样半径成立。

## 3. Target计数、cap与总质量

在每张target面，K_w=1+Bin(n−1,2d/r)，其均值≤1+2nd、方差≤2nd。nd≥8给

P{K_w>4nd}≤[2nd]/(2nd−1)²≤8/(9nd)≤1/8。                     (5)

阈值2d<r已由第1节保证。target corner-captured部分每面比例≤(3/4)^(n−1)：

1/r−1/t≥1/(1+3d/4)−2/L>1/2

可由更粗的d≤1/16、L≥64检验，故ell_t(r)<3/4。n≥17时(3/4)^16<1/8。

从(4),(5)及cap并集界，每张面至少5/8保留；只需更弱的质量结论

a*=λ(A)≥nd/16。                                               (6)

这里完整target annulus质量≥nd/4。upper a*≤nd也成立。

## 4. 全matching bridge及完整soft row

冻结源coordinate导数界在第1节正域内成立：

|∂_q logJ(w;q)|≤64/ρ[d²+1_{|w|>r−2d}] a.e.

J局部Lipschitz，故可跨onset kink积分。K_w≤4nd、n/ρ≤2L及matching长度≤2h给

|Σ logJ(w_i;q)−Σ logJ(w_i;q0)|≤1280Ldh
 =1280ε0=5/4096<1/8。

结合baseline |S_log|≤1/8，得到每个good target与每个matched q上

1/2≤∏_iJ(w_i;q)≤2。                                          (7)

例如|ΣlogJ|<1/4且e^(1/4)≤4/3<2即可。

matched坐标单独也满足1/2≤J≤2：T≤2d；piecewise公式在T≥g分支也给A(R,g,T)≤(4/3)T，而B≤g²/(3R)，准确公式给J≥1−4d/(3ρ)>1/2与J≤1+d²/(3Rρ)<2。这补明不能从product界反推individual bounds。

完整inverse source壳系数为

U(w,ρ)=[1/(2ell_*)](∏J_i)Σ_i N_i/J_i，1≤N_i≤2。

由ell_*∈[3/4,1]及(7)，明确有n/8≤U≤8n。

对相应capture/G的inverse-row概率，free坐标正规化因子J≥1/2，且ell≥3/4。由−logell≤(4/3)(1−ell)、1/ell≤4/3及E_full(1−ell)=g²/(4Rρ)，

E_inv[−logell]≤8g²/(9Rρ)≤d²/ρ。

pinned source坐标的−logell_*≤d。因此V=−logG满足

E_inv V≤d+nd²/ρ≤d+2Ld²≤3Ld²，

最后只需Ld≥1。Markov在该**inverse-row概率**下给P_inv{G<θ}≤1/4。于是bounded soft source weightb=min(1,θ/G)的完整壳row系数介于3nθ/32与8nθ。

matching source log长度至少h、至多4h，所以去cap前的完整row介于3k/32与32k；该结论现在没有eventual前提。

## 5. Source cap删除的显式全尾

source原corner补集的壳指数满足

(n−1)(1/ρ−1/t)/2=(1−1/n)xv/(2q)≥xd/16，

因为n≥2、v≥d/4、q≤1。其captured壳系数≤(n/2)exp(−xd/16)，source log宽≤d，故总质量≤nd exp(−xd/16)≤nd exp(−Ld/32)。b≤1使每个row的删除量不超过这整个unweighted质量。

除k，利用12d≤1/64，得到

Cap/k≤2^20Ld² exp[−Ld/32+12Ld²]
 ≤2^20Ld² exp(−Ld/64)
 ≤24·2^44/(L³d²)，                                            (8)

末步是Taylor4的全域界e^u≥u^4/24。常数2^44来自2^20·64^4，不少任何d因子。

(8)右侧随L递减，L=L0的Fraction基点严格小于1/36<3/64。因此最终rows及columns有

R_λ(w)≥3k/64，R_λ≤32k，C_λ≤4k。                              (9)

column bound由b C_raw≤4nhθ=4k直接成立，不需要source concentration或source质量的另一次Markov。

## 6. τ/k的Taylor64门槛

令Δ=d/65536=2^(−40)，δ=Δn，σ=32d²=2^(−43)。原cap界在x∈[L/2,L]给

B0≥n[L−logL−4/L]>8n，B0≤nL，A0≤ne^(−σL)。

L0≥64足以使用前两个界；x≥L/2>e给B0≤1+n(L−logx)≤nL。δ<n，所以B0−δ严格正。选择

τ=√[A0(B0−δ)]≤n√L exp(−16d²L)。

因此

τ/k≤2^20 d L^(3/2) exp(−4d²L)。

Taylor64在u=4d²L给出全尾上包常数乘L^(−125/2)，严格递减。在L0，4d²L0=64且2^20d L0^(3/2)=2^74，故

τ/k≤2^74·64!/2^384<3/128。                                    (10)

最后一项为纯整数/Fraction检查。这里指数中的source θ已经准确扣除：−16d²L+12d²L=−4d²L，未漏真实source支付成本。

## 7. 能量缺额、有限a与一般L1全部尺度

由(6),(9),(10)，最终weighted kernel的max degree≤32k、actual μ-row≤τ，energy缺额引理给

Δ_missing≥[(3k/64−τ)a*]/(32k)
 ≥3nd/65536≥δ。

对每个a≥8n^4，有a>2t、c=(aβ−1)/(a+1)>1，且

αβ−κ−2/(a+1)=ζ−(β+2)/(a+1)>ζ−4/(8n^4)>0。

任何actual target increment遇到A时半径<α+r<3；source半径≥ρ>4，所以source在同一target之前被选。上述center guard把capture条件传回原bQ_i中心资格，b(z)≤1使其weighted实际row仍≤τ，不能乘任何θ^(−1)。

在A/P上λ_a=α^(−n)λ准确成立，degrees、energy和target质量按原正倍率缩放；缺额下界不会变差。由此δ对每个a≥8n^4的所有真实有限历史成立。没有引入实际μ的分布假设。

原小立方体的a倍仍是同一一般f∈L1的真实共同半径witness。对任意紧水平集取这些小立方体的有限子覆盖，使用量化后的旧selection与新post-selection energy；固定n,t,τ后仅对已证明的数值不等式令a→∞，再用内正则性，得到

C_n≤F(A0,B0−δ)。

所有原正尺度都保留，既不要求f紧支撑，也没有有限尺度operator或选择历史收敛前提。

## 8. 完整scalar cost与端点负线性项

原scalar预算给

B0≤T+64d²n+1，A0≤ne^(−σL)，
F(A0,B0−δ)≤T−Δn+2σn+nR(L)，
R(L)=1/n+e^(−σL)+2√L e^(−σL/2)。

这里2σ=64d²=2^(−42)，Δ/2=2^(−41)，故需要R≤2^(−42)，而不是只检查−δ。

用n≥L²/2以及Taylor16，L≥L0时R有三个递减power上包。在基点σL0=512、σL0/2=256，

R(L0)≤2^(−103)+16!/2^144+2^27·16!/2^128<2^(−42)。             (11)

全尾分别按L^(−2)、L^(−16)、L^(−31/2)递减，故(11)对全部L≥L0成立。于是完整root cost严格小于

T−Δn+2σn+2^(−42)n=T−2^(−41)n，

闭合(1)。新operator门槛是所有第1–8节条件的共同闭合结果，不是由(11)单独命名。

## 9. 最小exact验证与状态

`verify_threshold.py`只使用Fraction/integer在L0核验baseline mean/variance、target概率域、source cap、Taylor64、Taylor16、固定常数及源hash；同时记录解析尾的power。没有重新跑任何已通过网格、生产LP、数值积分或MATLAB。

本报告的所有无限尾由全域Taylor上包和递减powers覆盖；q/r范围内的界在解析中统一给出。验证receipt保留每个gate与冻结源hash。此显式门槛仍须独立threshold审计后归档，不改变R045原有“各fixed d eventual”的已接受结论。

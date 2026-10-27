# R055 · 简单同族tuple的standalone证书与L≥2^27门槛

## 身份及结论

本报告从R051原有公式重新构造一个指定tuple的完整准确证书，不运行新搜索、不重跑R052列表。它属于R052登记有限列表，但未保留在top64；本报告的真实性不依赖该tuple曾在生产枚举中如何运行。

固定

a=9，s=19，d=2^(−10)，ε0=2^(−20)，
εJ=εS=εT=εD=1/1024，εA=1/64。

条件于R051泛族全链接受，下列全部finite guards已量化到L=log n≥2^27：

C_n≤Φ_simple:=F(A_n(t),B_n(t)−γn)<T−2^(−16)n，                 (1)

T=nL−nlogL，x=(1−σ)L，σ=19/2^20，t=n/x。R051泛族接受身份由主报告统一激活；本有限门槛仍需独立threshold审计，不把主估或scalar基点提前当算子界。

该display系数与R053相同，但显式起点降低。两个完整cost可按各自有效域取minimum；不能仅因本门槛更小就宣称它的完整cost始终低于R053。

## 1. 独立准确budget

从P0=3579/4096、qmin=1−3d/4及R051公式直接得到

u_lo=3579/8200，u_hi=[P0(1−d/2)(1−εJ)]^(−1)，M*=4，
m*=63/131168，b*=151679/524800，
Γ=9555777/275347865600，
γ=(1−εD)Γ=9775559871/281956214374400，
σ=19/1048576，σ/(1−σ)=19/1048557，
c_raw=γ−σ/(1−σ)=4893063658542547/295647162275777740800>2^(−16)。

全部domain a>3、s>2a、d≤1/16、d<1/(32a)、σ<1/2、error domains及b*>0均为严格或非严格准确Fraction检查。M=4是重新计算max(4,2u_hi/qmin)的结果，没有借用某张旧表的branch。

在R052登记记号中此点为A=72、ρ=16/16、k=1024，即a=A/8、s=2a+ρ、d=k/2^20。该membership不是新的最优性声明，也没有从小数重构系数。

## 2. 统一几何及target好概率

令L0=2^27，h=ε0/(Ld)、ζ=n^(−3)、κ=1+1/t−ζ。Taylor给n≥L²/2及n≥L^8/8!。x∈[L/2,L]、t≥L/2≥16；target r<2，source及baseline/bridge q>1/2，ρ>4r。

L0处h/d=2^(−27)<1/8，且ζt/h≤2d/(ε0n²)≤1/2；所有尾部保留这些界。故R=q+d−h>1、0<g=d−h+ζtq≤d、g≤R/2及κ>1；整个q0±matching bridge满足accepted coordinate公式正域。

每个固定target face的baseline log和有

|EΣlogJ|≤M(L)=8![24d³/L^6+4d/L^7]，
Var(ΣlogJ)≤V(L)=48·8!d³/L^6。

在L0准确有M≤1/32、V≤(3/32)²εA/3。故log坏率≤εA/3。又K_w=1+Bin(n−1,2d/r)，n≥L²/2使8/(9nd)≤εA/3。原corner的face比例≤(3/4)^(n−1)≤(3/4)^32<εA/3，因n≥33且ell_t(r)<3/4。

所以每张面保留比例≥1−εA，给准确λ(A)≥m*n，而不是沿用旧5/8。所有M,V、Kw概率上包随L递减；corner条件由d≤1/16、t≥L/2及n≥33保持。

## 3. Product与individual J的准确精度

matching log变化≤1280Ldh=5/4096，故baseline 1/8给

P0≤∏J≤1/P0。

pointwise闭式统一给|J−1|≤2d/ρ≤4dL/n≤4d·8!/L^7。在L0该上包≤εJ，尾部继续递减。因此generic的u_lo,u_hi都可准确使用，不需对n个坐标的J误差再作并集概率。

ell_*≥1−d/2、1≤Ni≤2及inverse-row Markov（Ld≥1）给soft retention≥1−3/a。matching源log宽为[h,2h/qmin]，所以raw row在(1−3/a)u_lo k与2u_hi k/qmin之间，column≤4k。

## 4. Source cap与τ的有限门槛

本tuple有两个准确衰减率

η=d/32−ad²=23/2^20，ν=(s/2−a)d²=2^(−21)。

原corner source删除量满足

Cap/k≤2^20Ld² exp(−ηL)=L exp(−ηL)。

Taylor4给递减L^(−3)上包。在L0，ηL0=2944，故

Cap/k≤24·2^27/2944^4<εS。                                    (2)

对于实际τ=√[A0(B0−γn)]，原cap界给B0>8n、B0≤nL及A0≤ne^(−σL)，从而

τ/k≤2^20dL^(3/2)exp(−νL)。

L0处νL0=64，prefactor=2^(101/2)<2^51。Taylor64给

τ/k<2^51·64!/2^384<εT。                                       (3)

其全尾power为L^(−125/2)。没有数值指数外推。

现在row−τ≥b*k、maxdegree≤M*k及λ(A)≥m*n都成立，实际能量缺额≥Γn>γn；固定εD余量保持。

## 5. Finite-a和一般L1映射

原dilation mathfrak a≥8n^4保证κ+2/(mathfrak a+1)<αβ。actual target半径<3，source半径>4，故同一target的一次eligibility支付全部source；bounded source weight≤1使weighted μ-row不超过τ。

λ_a=α^(−n)λ在两侧准确成立，故finite-a缺额不减。原共同半径witness、actual disjoint increments、紧水平集有限覆盖与内正则性不变。固定n,t,τ后仅对已经有效的数值界令mathfrak a→∞，即可在R051全链接受的条件下得到Cn≤Φ_simple。没有新输入类、尺度截断或selector极限。

## 6. 完整scalar cost

准确使用σ/(1−σ)，得到

Φ_simple≤T−c_raw n+nR(L)，
R(L)=1/n+e^(−σL)+2√L e^(−σL/2)。

L0处σL0=2432>2^11，σL0/2=1216>2^10；Taylor16与n≥L²/2给

R(L0)≤2^(−53)+16!/2^176+2^15·16!/2^160
 <c_raw−2^(−16)。                                              (4)

对应全尾按L^(−2)、L^(−16)、L^(−31/2)递减，故(1)适用于全部L≥L0。这个结论已经支付完整正buffer及交叉根式，不是只显示γ。

## 7. 与R053准确raw cost保留minimum

R053的准确rawtuple具有更高c_raw，但本门槛更早。无需仅凭两个c_raw就断言完整cost排序；可以另外证明如下公式事实。

记γ_s,σ_s为本报告simple tuple，c_raw,r为R053准确rawtuple。原cap lower给

Φ_simple≥T−n[γ_s−σ_s+4/L]，

因为B_simple≥T−nlog(1−σ_s)−4n/L≥T+σ_s n−4n/L。

R053准确scalar上包为Φ_raw≤T−c_raw,r n+nR_raw(L)。在Lr=2^32，R053已证明

R_raw≤2^(−63)+16!/2^256+2^17·16!/2^240。

精确Fraction另核

c_raw,r−R_raw(Lr)>γ_s−σ_s+4/Lr。

左侧的误差upper随L递减，右侧4/L递减；可直接把两项正loss移到同一侧，得到全L≥Lr的Φ_raw<Φ_simple。该排序是full cost而不是display大小。

所以在2^27≤L<2^32时，simple tuple提供更早的同display证明；进入两个门槛交集后保留两完整cost的minimum，且上式给raw完整cost更小。所有operator身份仍依各自threshold独审及R051接受，不追溯改写R053原conditional生成状态。

## 8. 核验身份

`verify.py`仅从R051公式重建一个指定tuple与有限门槛，不重新枚举R052列表。receipt给完整有理certificate、全部domain、概率/cap/τ/scalar/finite-a gate及原source hashes。没有参数搜索、生产数值积分、MATLAB、Git或旧源修改。

R051 frozen source SHA-256：`aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6`。R053 frozen proof及原retained-candidate身份保持不变。本报告不声明连续参数最优、全列表新排名或文献纪录。

# R069 · R066首profile机制的全operator显式门槛

## 0. 状态与准确对象

本报告量化已接受R066的所有eventual条件，证明对L=log n≥L0=2^30的整数n有

C_n≤Phi66:=F(A_n(t),B_n(t)−gamma n)<T−(23/100)n，
T=nL−nlogL，F(A,B)=(sqrt A+sqrt B)^2。

新门槛仍待soft独立threshold审计，审毕前不写入已接受Cn账本。它不改变R066原来的eventual运行身份，不声称2^30最小，也不改变同一Phi66的核、缺额或scalar函数。

冻结R066作者源为geometry clone `upper/round066/target_profile/PROOF.md`，SHA256
`4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1`；完整审报告为uniform clone `audit/round066/target_profile/review.md`，SHA256
`5a13e5a4364a2fbcc4b9d89fb2536ebeded2569033041017a7355237fd110ea8`。

相对coordinate输入为roots clone `audit/round064/relative_coordinate/LEMMA.md`，SHA256
`c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`。
准确log价格读取本clone `audit/round066/exact_scalar/SCALAR_PROOF.md`，SHA256
`d058c87d90db9748ac05a6d42fd711f68546daa217eee5ef8106b9c6e1d689f0`；原八项receipt不重跑、不改写。

## 1. 固定tuple与已接受预算

准确固定

s=3/2，qmin=2/3，qmax=917/1000，delta=2^(-30)，
rlo=1+delta，rhi=(3/2)qmax−delta=2751/2000−delta，
sigma=83/1000，chi=9/200，x=(917/1000)L，t=n/x，
zeta=n^(-3)，epsilon=zeta t，beta=1+1/t，kappa=beta−zeta，
h=h_c/L，h_c=2^(-80)，theta=exp(−chi x)，k=nh theta/rhi。

e=2^(-20)用于epsJ/N/S/T/D/R和baseline log tolerance；epsRet=1/4096，epsA=1/65536。保留

P0=1−e−2^(-60)，Ulo=P0(1−e)/(1+e)，Uhi=(1+e)/[P0(1−e)]，
M=2(1+e)Uhi，b=2(1−epsRet)(1−e)Ulo−2e，
mL=637/2000，Gamma=mL*b/M，gamma=317/1000。

冻结full proof的准确证书已给(1−epsA)log(rhi/rlo)>mL及(1−epsD)Gamma>gamma，故Gamma>gamma且M>2。新target mean常数mu0=1/2−(5/2)log(6/5)<221/5000，profile转source mean和variance的准确常数也已经独审。这些固定预算不是本轮重跑对象。

## 2. 所有几何正域、source延拓与相对gap

令N(L)=L^8/8!，Taylor全域给n=e^L≥N(L)及n≥L²/2。N、N/L严格递增。L≥L0时x∈[L/2,L]、x>1，所以t≥N/L且t<n。准确基点gates给

h<min(delta/2,1/1000)，epsilon≤1/(xN²)<1/2，
qmax/(xN²)<2h/3，N/L>16，qmin N/L>4rhi>4，N≥400，N≥65。

h按L^(-1)递减；[1/(xN²)]/h按L^(-16)递减；其余lower递增。故全尾保持所有条件。

整个full source band q∈[qmin,qmax]上R=sq−h≥1−h>0，a0=(1−epsilon)q>0，g=(s−1+epsilon)q−h>0。并且

g−R/3=epsilon q−2h/3≤0。

因此ell_pin≥5/6，所有inactive source的G/b仍有正定义。局部R≥1估计只用于active sources：匹配到r≥rlo时R≥rlo−2h>1。对于固定target r，q_±=(r±h)/s全在原source band内，q_b=q_+使R=r，整个matching桥R∈[r−2h,r]，长度2h/s=4h/3。

于是局部域为q≥2/3、1<R≤r<3/2、rho=tq≥qmin N/L>4r、0<g≤R/3<1/2、epsilon q≤h/s、epsilon≤h、h≤1/1000。以上正是相对coordinate lemma的全部域，而非旧absolute g≤1/3域；g可以超过1/3。epsilon≤h来自epsilon qmax≤2h/3和qmax>2/3。beta<2、kappa>1亦成立。

两annuli均在norm<t内；target r<3/2，source rho>4，故A/P不交。Inactive source无匹配target，最终column为零，不能对其误用active R≥1。

## 3. 五个target坏事件，各付epsA/8

所有概率在每个固定full target面和每个r上计算；pinned坐标确定，其余n−1坐标独立uniform[-r,r]。实际未知history measure不被假定有此law。本轮不需要K_w事件：新的relative derivative无条件成立。

### 3.1 Baseline log J

冻结relative lemma给free坐标mean绝对值≤6g³/rho_b²、variance≤12g³/rho_b²，pinned log绝对值≤g/rho_b。使用g<1/2、rho_b≥n qmin/x及x≤L，有

|EΣlogJ|≤[6(1/2)^3 L²/qmin²+(1/2)L/qmin]/n
 ≤Mlog(L):=8![6(1/2)^3/(qmin²L^6)+(1/2)/(qmin L^7)]，
VarΣlogJ≤Vlog(L):=12(1/2)^3*8!/(qmin² L^6)。

基点准确核Mlog<e/2、Vlog<(e/2)² epsA/8，所有power递减。因此Chebyshev给baseline失败率<epsA/8。

### 3.2 两个source-face count

在q_-和q_+，Nsum有准确mean n+(n−1)a0(q_±)/r、variance≤n/4。以e*n/4为误差阈值，各失败率≤4/(e²n)≤4/(e²N)<epsA/8。

此处ell_pin'(q)=−(1−epsilon)h/(2R²)<0，而Nsum非减，因此比率Nsum/(2n ell_pin)非减。原full proof保留n−1的准确endpoint估计，在2h+1/n≤3e/4时给全matching的[1−e,1+e]。基点2h+1/N<3e/4，随后两幂递减，故无任何旧additive导数符号移植。

### 3.3 原target cap

准确基点1/rhi−L/N>1/2，随后N/L递增。因此每个正pinned target面的补集free fraction

1+1/(2t)−1/(2r)<3/4。

负pinned面补集为空。N≥65和(3/4)^64<epsA/8给逐面坏率<epsA/8，没有只在signed-face平均才成立的1/2。

### 3.4 新的Jphi profile事件

记H_i=Jphi(w_i;q_b)，此处只在baseline TARGET law取矩，不能换成inverse-source law。冻结Fubini身份给nE_free H≤mu0*x，其中mu0<221/5000；pinned项≤1/rho_b，且0≤H_i≤1/rho_b。因此

EΣH_i≤mu0*x+1/rho_b，
VarΣH_i≤n/rho_b²≤x²/(n qmin²)。

准确N qmin≥200000给1/rho_b≤x/(n qmin)≤x/200000。定义条件ΣH_i≤(221/5000+1/100000)x时，threshold距mean至少x/200000。逐target face Chebyshev遂给

P(profile失败)≤200000²/(n qmin²)≤90000000000/N<epsA/8。

该上包按L^(-8)递减。没有把pin变成free，也没有假定真实mu集中。

五项并集≤5epsA/8<epsA，所以完整good target mass至少(1−epsA)nlog(rhi/rlo)>mL*n。无需这五个事件相互独立。

## 4. 全matching传输、individual J与径向宽度

冻结a.e./AC derivative给|∂q logJ|≤256/rho；用n/rho≤(3/2)L和宽度4h/3，得到

|ΣlogJ(q)−ΣlogJ(q_b)|≤512Lh=2^(-71)<2^(-60)。

于是baseline tolerance e和P0给P0≤∏J_i≤1/P0。individual界单独满足

|J−1|<1/(3rho)≤L/(3qmin N)<e，

基点通过后按L^(-7)递减；不把n个individual errors相乘替代product桥。

Profile的a.e./AC derivative |∂q Jphi|≤8/rho，在同一matching区间给

|ΣJphi(q)−ΣJphi(q_b)|≤16xh。

基点h<1/1600000，全尾递减，故16xh≤x/100000，和profile条件合并得全q统一

ΣJphi(w_i;q)≤(2211/50000)x。

其作用是使每个captured source face的free mean精确为Σ_(l!=i)Jphi/J_l；并非另取无关worst target。

令u=h/r。原匹配source log宽与target radial weight的准确组合为

(r/rhi)log((r+h)/(r−h))=(2h/rhi)*log((1+u)/(1−u))/(2u)。

积分或非负级数给最后ratio在[1,1/(1−u²)]。基点h/rlo<1且1/[1−(h/rlo)²]<1+e，随后h/rlo递减。因此冻结的relative radial width e预算合法，完整rows的同一normalization k=nh theta/rhi被保持。

## 5. 每个source face的retention与完整度预算

冻结全链已证明每个captured signed source face上

E_inv V<1/5+(223/5000)x，Var_inv V<x/150，V=−logG。

其中mean使用上一节profile；variance才通过fullsource的非负二阶矩去capture。两种law没有交换。x≥1000给Qpin/x<1/5000，且

chi−223/5000−1/5000=1/5000。

因此先在每face应用Chebyshev、再混合概率，得到

P_inv(V>chi*x | face)≤5000²/(150x)。

准确L0处该值严格小于epsRet=1/4096；随后按L^(-1)递减。不能把within-face variance直接说成mixture variance。Fullsource soft shell rows因此保留至少1−epsRet的theta-scaled inverse row。

计数、product和individual J校准给Ulo*n≤U≤Uhi*n。结合径向宽度，raw rows下界2(1−epsRet)(1−e)Ulo*k、上界M*k。完整column仍是原精确target面耦合：T_i/(2ell_i)≤1、Σ pinned faces系数≤n∏ell_i(r)≤nG；(dr/r)*a_rad=dr/rhi，target ordinary匹配宽≤2h，故soft后column≤2k<Mk。Inactive列为零；没有漏掉一个面或另加放大系数。

## 6. Source cap、tau及全部scalar正域

宽source band质量≤Csrc*n，Csrc=(qmax−qmin)/qmin；准确Csrc*rhi<1。源原cap删减除实际k的上包为

(Csrc*rhi/h) exp[−((1−1/n)eta−chi)x]，eta=83/1834。

准确eta−chi=47/183400。基点N≥400，且纯Fraction核eta−chi−eta/400>1/10000。因此全尾指数coefficient>1/10000，令cS=(917/1000)/10000，得到全域Taylor10上包

S10(L)=10! Csrc*rhi/[h_c cS^10 L^9]。

L0处prefactor<2^110、cS L0>2^16，所以S10(L0)<10!/2^50<e。之后按L^(-9)递减。cap删除在fullsource inverse-face law估计之后支付，不能把cap条件律仍当独立坐标。

原accepted scalar积分界在x∈[L/2,L]、L≥64给

T−4n/L≤B_n(t)≤nL，A_n(t)≤ne^(−sigma L)。

下包可由I_n(t)≤2t/(n−1)和非负core得到B≥nlogt−2n/x≥T−4n/L；upper由B≤1+nlogt及x>4得到。在L≥16有logL<L/4，而L/4−4/L−1>0，所以B_n−gamma n>nL/2>0。该正域既付sqrt也付所有row阈值，不保留未知eventual项。

取同一个tau=sqrt[A_n(t)(B_n(t)−gamma n)]，则tau≤n sqrtL exp(−sigma L/2)。准确扣除theta的rate

nu=sigma/2−chi(1−sigma)=47/200000>0。

故tau/k≤(rhi sqrtL/h)e^(−nu L)，Taylor10全域上包

T10(L)=10!rhi/[h_c nu^10 L^(17/2)]。

L0平方根准确为2^15，prefactor<2^126、nu L0>2^17，所以T10(L0)<10!/2^44<e。随后按L^(-17/2)递减。支付epsS*k和同一epsT*k后，row excess≥b*k，所有度≤M*k，target质量>mL*n。因此energy缺额>Gamma*n>gamma*n；epsD留下已认证严格余量。

## 7. 所有真实finite histories与一般L1全尺度

这里逐项消去冻结源§8–9中的large-n guard。对每个algorithmic dilation a≥8n^4，由t<n和n>1有a>2t，故ccap=(a beta−1)/(a+1)>1。又beta<2，

alpha beta−kappa−2/(a+1)=zeta−(beta+2)/(a+1)>n^(-3)−4/(8n^4)>0。

在D上有限envelope为lambda_a=alpha^(-n)lambda；D⊂H_a的差[((beta+1)m−1)/(a+1)]>0对所有m≥1成立，不要求另一个维数阈值。真实uncaptured increments限制到D后仍受lambda_a支配，总finite缺额至少这个限制上的缺额。

实际target item遇A时ccap*r_j≤||u−w||<alpha+rhi<3，故r_j<3。实际source point由共同load support与真实increment有r_i≥rho>4，所以先于同一target j。新capture使真实中心距<[kappa+2/(a+1)]r_i<alpha beta r_i，因而原eligibility确实支付这些同一早期items。两权a_rad,b_soft≤1，实际weighted row≤a_rad*tau≤tau；不能乘theta^(-1)，也没有给实际mu假定随机law。

A/P不交，令ca=alpha^(-n)≥1。能量缩放ca²、度缩放ca、target mass缩放ca；difference identity给

Delta_a≥[ca² E−tau*ca*lambda(A)]/(ca*M*k)
 ≥[E−tau*lambda(A)]/(M*k)>gamma*n。

因此对所有a≥8n^4和所有原finite histories统一成立，不依赖selector或history极限。原cover factor、captured tau与terminal a^(-n)继续保留在finite inequality中。

对任意f∈L1取|f|，strict level set的每个点取自己的真实正共同radius witness，open小cube覆盖每个compact子集，有限子覆盖的每个selected dilate仍是真实witness。共同终端radius归一化保留上述所有guard。固定n,t,tau后仅在统一数值不等式中令a→∞，再以内正则性传到整个level set，得到Cn≤Phi66。没有compact support、radius lattice、scale cutoff或矩形operator替代；所有正共同radius均保留。

## 8. 精确log价格与完整root显示系数

读取已运行且完整前提现已满足的exact scalar证书，不重算log级数：price=log(1000/917)<4333/50000，

gamma−4333/50000−23/100=17/50000。

原精确价格B≤T+price*n+1给

Phi66≤T−(gamma−price)n+nR(L)，
R(L)=1/n+e^(−sigma L)+2sqrtL e^(−sigma L/2)。

保留A、cross、1/n后，以n≥L²/2和Taylor4上包R为

R4(L)=2/L²+4!/(sigma L)^4+2sqrtL*4!/(sigma L/2)^4。

三项按L^(-2)、L^(-4)、L^(-7/2)递减。L0处sigma L0>2^26、sigma L0/2>2^25、sqrtL0=2^15，故

R4(L0)<2^(-59)+4!/2^104+2^16*4!/2^100<17/50000。

于是整个Phi66<T−23n/100在L≥2^30成立。只比较gamma或只scalar基点均不足命名operator门槛；该结论使用前述全部统一guard。没有声称精确端点gamma+log(917/1000)本身成立。

## 9. 对旧Phi64的完整cost比较与valid域

旧R064完整cost为Phi64=F(A_n(t64),B_n(t64)−gamma64*n)，sigma64=91/1000、gamma64=61/200、x64=(909/1000)L。它的机制原结论仍eventual，没有从本轮移入显式域。

对原scalar B，下包保留正buffer：

B64≥T−log(909/1000)n−4n/L≥T+(91/1000)n−4n/L。

所以Phi64≥B64−gamma64*n≥T−[61/200−91/1000+4/L]n。准确Fraction在L=2^20已经给括号<9/40<23/100；括号随后递减。因此在本轮新有效域L≥2^30且旧R064自身机制有效的交集内，

Phi66<T−23n/100<Phi64。

这是不同x的完整root cost排序，不是两显示负系数的排序；不把旧未知N64变成2^20或2^30。其他已接受界仍按各自真实valid域取min，不宣称穷尽全部覆盖策略或公开记录。

## 10. 新验证与来源身份

verify_threshold.py仅做一批固定新L0的准确Fraction/integer gates，记录实际PID/时间、source/proof/code/旧价格receipt hashes及明确tails。它不调用作者R066的16gate程序，不重跑已有log/旧baseline/生产网格，不使用MATLAB/Git，不改其他作者文件。

所有q/r/face上的统一性来自原冻结全域界与本报告解析上包。L≥L0无限尾来自Taylor全域不等式、N与N/L单调及上文明确负powers，不能用有限采样外推。本轮门槛仍须soft独审，未ACCEPT前不进入Cn已接受账本，也不追溯改写原R066 eventual receipt。

# R063 · R059宽source band的显式operator门槛

## 0. 状态、目标与准确来源

本报告证明R059全部eventual条件在L=log n≥L0=2^33统一成立，继而给出

C_n≤Phi59:=F(A_n(t),B_n(t)−gamma n)<T−n/10，
T=nL−nlogL，x=(97/100)L，t=n/x，gamma=1637/12500。

本报告尚待独立threshold审计，不提前写入已接受Cn账本。冻结R059的已接受完整机制只有eventual原结论；这里的新门槛不追溯改变其历史run身份，也不声称门槛最小。F(A,B)=(sqrt A+sqrt B)^2仍保留全部cover成本。

实际读入R059作者源：
`/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/PROOF.md`，SHA256
`fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f`。

其full审报告为soft clone `audit/round059/wide_source_band/review.md`，SHA256
`c8aff39bf68749d25f9978cef11be82ccb7907ec7e72c5678ac66b06a65a33d4`。
局部coordinate引理为roots clone `audit/round056/one_coordinate/LEMMA.md`，SHA256
`97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec`。

## 1. 固定tuple与被接受预算

准确固定

d=391/1905，a_band=283/2032，b_band=63/64，sigma=3/100，chi=147/10000，
e=2^(-20)，delta=e，h_c=2^(-48)，h=h_c/(Ld)，
epsJ=epsN=epsS=epsT=epsD=e，epsRet=1/4096，epsA=1/65536。

baseline log tolerance为e。vmin=a_band*d，vmax=b_band*d，qmin=1−vmax，qmax=1−vmin，Rbase=qmin+d，rlo=Rbase+delta，rhi=qmax+d−delta。

源中已接受的准确预算是

P0=1−e−1280h_c，Ulo=P0(1−e)/(1+e)，Uhi=(1+e)/[P0(1−e)]，
M=2Uhi/qmin>2，b=(1−epsRet)2Ulo/qmax−2e，mL=159477/10^6，
Gamma=mL*b/M，gamma=1637/12500<(1−epsD)Gamma<Gamma<1。

(1−epsA)log(rhi/rlo)>mL、pin log<Qbar=26979/250000<1/9、meanrate<293/20000、varrate<1/900均依冻结full proof及其原准确log证书。本轮不重跑这些log认证。保留其精确scalar价格

gamma−sigma/(1−sigma)=121289/1212500=1/10+39/1212500。

## 2. 全部几何guard及局部coordinate接入

记N(L)=L^8/8!。Taylor全域给n=e^L≥N(L)及n≥L²/2。L≥L0时N递增，N/L递增；x∈[L/2,L]，x>1，

t=n/x≥N/L，t<n，1/2<qmin<qmax<1，1<Rbase<rlo<rhi<5/4，0<d<1/4<1/3。

准确基点检查给

h<min(delta/2,d/8,(Rbase−1)/2)，
1/[xN²]<h/2，N/L>16，qmin N/L>4rhi，N≥200000，Nd≥1。

第一项随L^(-1)递减，第二项除h后的ratio随L^(-16)递减，其余lower递增，所以对全尾成立。于是zeta t=1/(xn²)<h/2，kappa=1+1/t−zeta>1，beta=1+1/t<2。

整个source band的参考半径R=q+d−h>1；存在matched target的active source更有q≥qmin+delta−h、R≥Rbase+delta−2h>Rbase。固定target r的整个matching区间[q_-,q_+]=[r−d−h,r−d+h]都在原source band内，长度准确2h，log宽在[2h/qmax,2h/qmin]。baseline q_b=q_+使R=r；整个bridge上R∈[r−2h,r]，a0=(1−zeta t)q>0，

g=d−h+zeta t q∈(0,d)，g<R/2，rho=tq≥qmin N/L>4rhi≥4r。

另zeta t<d，m=(|w|−a0)_+≤g+2h≤5d/4。以上逐项满足被冻结coordinate lemma的d≤1/3、h≤d/8、q≥1/2、R≥1、rho≥4r、g≤min(d,R/2)及epsilon=zeta t≤d条件；不套用旧small-d机制。Target半径<2<t，source半径rho>4且rho≤qmax t<t，A/P不交。Inactive source的G/b仍有正定义，但column为零，不能对其误用R≥Rbase的active moment界。

## 3. 五个坏事件的统一定量概率

所有概率仅在每个固定full target面/半径上计算，一个pinned坐标加n−1个独立uniform坐标。实际history measure无需满足这种随机law。

### 3.1 Baseline

冻结局部矩给S=ΣlogJ(w_i;q_b)满足

|ES|≤[6d³L²/qmin²+dL/qmin]/n≤Mlog(L):=8![6d³/(qmin² L^6)+d/(qmin L^7)]，
Var S≤12d³L²/(qmin² n)≤Vlog(L):=12·8!d³/(qmin² L^6)。

pinned坐标只贡献mean。准确基点给Mlog(L0)<e/2及Vlog(L0)<(e/2)² epsA/8。所有powers递减，因此Chebyshev给P(|S|>e)≤epsA/8，uniform于每面和每r。

### 3.2 必须保留的K_w事件

d<1/4，故K_w≤4nd不自动成立。其law为1+Bin(n−1,2d/r)，mean≤1+2nd，variance≤n/4。Nd≥1时4nd−mean≥2nd−1≥nd，故

P(K_w>4nd)≤1/(4nd²)≤1/(4N d²)<epsA/8。

这正是匹配derivative求和所需额外事件，不能删除。

### 3.3 两个endpoint source-face count

每个q_±的Nsum均有准确mean n+(n−1)a0(q_±)/r及variance≤n/4。阈值epsN*n/4的两项失败率各≤4/(e²n)≤4/(e²N)<epsA/8。

基点另有2h+1/N<3e/4、2hd<3e/4，随后递减。冻结的cross-endpoint校准因此在整个matching interval给

1−e≤Nsum(q)/(2n ell_*(q))≤1+e。

这里保留n−1产生的1/n drift，没有把lost faces当作o(n)。

### 3.4 Target原cap

每个正target面的cap补集free fraction为1+1/(2t)−1/(2r)。基点1/rhi−L0/N(L0)>1/2，此后N/L递增，所以它<3/4；负面补集为零。用N≥65及准确(3/4)^64<epsA/8，得到逐面补集概率≤(3/4)^(n−1)<epsA/8。没有只在signed-face平均成立的因子1/2。

五项并集≤5epsA/8<epsA。因此lambda(A)≥(1−epsA)nlog(rhi/rlo)>mL n。无需五事件独立。

## 4. Product桥、individual J与inverse retention

在K_w≤4nd上，冻结a.e./AC derivative给Σ|∂q logJ_i|≤640Ld。整个matching长度≤2h，故log bridge≤1280Ldh=1280h_c。结合baseline tolerance e，准确P0=1−e−1280h_c给

P0≤∏J_i(q)≤1/P0。

独立逐坐标上包为|J−1|≤5d/(8rho)≤5dL/(8qmin N)。该式在L0<e，此后L^(-7)递减，所以J_i∈[1−epsJ,1+epsJ]全q统一成立，而不对n坐标作union bound。与上一节计数合并后，完整signed-face inverse row系数在[Ulo n,Uhi n]。

仅对active captured source face使用R≥Rbase。冻结准确moment常数给每个这样的face

E_inv V≤1/9+(293/20000)x，Var_inv V≤x/900，V=−logG。

chi−293/20000=1/20000；当x≥40000/9时，1/9≤x/40000，所以threshold chi*x距每face uppermean至少x/40000。于是

P_inv(V>chi*x | face)≤40000²/(900x)。

准确基点x0=(97/100)2^33满足40000²/(900x0)<1/4096=epsRet；这是本门槛中较紧的条件。上包按L^(-1)递减，故全尾成立。先逐face作Chebyshev再混合其概率，不把独立variance错误地赋给face mixture。Soft weight min(1,theta/G)遂保持至少1−epsRet的theta-scaled inverse row。

## 5. 宽source cap与tau的全域Taylor8尾

源cap必须使用宽band质量Csrc*n，Csrc=(qmax−qmin)/qmin<1，不能借用旧nd。准确eta=vmin/(2qmax)<1且eta−chi=499871/37603070000>1/100000。N≥200000使

(1−1/n)eta−chi>1/200000。

因此cap删减除k后的量≤(Csrc/h)exp[−x/200000]。令cS=97/20000000，它有全域Taylor8上包

S8(L)=8! Csrc*d/[h_c*cS^8 L^7]。

L0处Csrc/h<2^81、cS L0>2^15，故

S8(L0)<2^81·8!/2^120<2^(-20)=epsS。

S8按L^(-7)递减，覆盖全尾，不由有限source采样推断。

选择tau=sqrt[A_n(t)(B_n(t)−gamma n)]。原scalar积分界在x∈[L/2,L]、L≥64给

T−4n/L≤B_n(t)≤nL，A_n(t)≤ne^(−sigma L)。

这里B−gamma n>0：对L≥16有logL<L/4（在16由log2<1验证，随后导数1/4−1/L>0），所以L−logL−4/L−gamma>L/2，基点L0已进入此域。因此tau≤n sqrtL exp(−sigma L/2)。扣除theta的准确exponent是

nu=sigma/2−chi(1−sigma)=741/10^6，
tau/k≤(sqrtL/h)e^(−nu L)。

Taylor8全尾上包为T8(L)=8!d/[h_c nu^8 L^(13/2)]。L0处sqrtL0<2^17、sqrtL0/h0<2^98、nu L0>2^22，故

T8(L0)<2^98·8!/2^176<epsT。

T8按L^(-13/2)递减。最终ambient degrees≤Mk、row excess在支付cap与tau后≥b*k，target质量>mL*n，因此完整缺额≥Gamma*n>gamma*n。epsD留下准确正margin，并非取极限端点。

## 6. 所有finite dilation与一般L1/all-scale接口

本节把R059 §7–8中剩余“large n”逐项消去。对所有algorithmic dilation a≥8n^4，已有t<n且n>1，所以a>2t，c_cap=(a beta−1)/(a+1)>1。又beta<2给

alpha beta−kappa−2/(a+1)=zeta−(beta+2)/(a+1)>n^(-3)−4/(8n^4)>0。

在A/P上lambda_a=alpha^(-n)lambda准确成立；D⊂H_a的代数差为[((beta+1)m−1)/(a+1)]>0，对m≥1无额外n门槛。实际target increment遇A时由c_cap*r_j≤||u−w||<alpha+rhi<3得r_j<3；source increment有r_i≥rho>4，所以它在同一nonincreasing-radius history中先于该target。Center displacement小于[kappa+2/(a+1)]r_i<alpha beta r_i，故每个新capture仍是该同一target真实eligibility所支付的早期item。由于soft weight≤1，actual weighted row≤tau，没有乘theta^(-1)。

A/P不交、energy/degrees/mass按alpha^(-2n)、alpha^(-n)、alpha^(-n)缩放，以及alpha^(-n)≥1，正好保持gamma*n缺额下界。上述guard对每个a≥8n^4同时成立，不要求选择历史或selector收敛。

原finite cover/witness界保留cover factor、captured tau和terminal a^(-n)。对一般f∈L1取|f|；每个strict maximal-level点取自己的真实正radius witness，小cube构成每个compact子集的开覆盖，取有限子覆盖后，所有selected dilates仍是原f的真实witness。终端半径共同归一化不改变上述guard，所以不需compact support、radius lattice或尺度截断。固定n,t,tau后仅在uniform数字不等式中令a→∞，再由内正则性到整个level set，得到完整Phi59。这里仍是原共同半径cube operator，不是辅助矩形算子。

## 7. 完整scalar根式与严格一成扣减

保留准确B价格而非仅比较gamma：

B_n(t)≤T+[sigma/(1−sigma)]n+1=T+(3/97)n+1，
Phi59≤T−(1/10+39/1212500)n+nR(L)，
R(L)=1/n+e^(−sigma L)+2sqrtL e^(−sigma L/2)。

用n≥L²/2和Taylor8，三项被全域powers上包：2/L²、8!/(sigma L)^8、2sqrtL*8!/(sigma L/2)^8，分别递减为L^(-2)、L^(-8)、L^(-15/2)。在L0，sigma L0>2^27、sigma L0/2>2^26、sqrtL0<2^17，故

R(L)≤Rpower(L)≤Rpower(L0)
 <2^(-65)+8!/2^216+2^18·8!/2^208
 <39/1212500。

因此完整Phi59<T−n/10对所有L≥2^33成立；该门槛来自前述全部geometry/probability/history条件的交集，而非这个scalar基点单独命名。

## 8. 与R056完整root cost的比较

本节不是只比较两条显示负项。R056已接受的完整cost为

Phi56=F(A_n(t56),B_n(t56)−gamma56*n)，x56=(93/100)L，
gamma56=172531397335235056462877787662549380239/
1643616601687301779182642153860479582208。

其机制有效域为L≥2^20。原积分界I_n(t)≤2t/(n−1)、n/(n−1)≤2以及非负core项给，对x∈[L/2,L]，

B_n(t)≥nlogt−2n/x。

在x56处于是B56≥T−log(93/100)n−4n/L≥T+(7/100)n−4n/L，因为−log(1−u)≥u。根式正性给

Phi56≥B56−gamma56*n≥T−[gamma56−7/100+4/L]n。

准确Fraction在L=2^20已给gamma56−7/100+4/2^20<1/10；括号随后递减。因此在各自真实valid域的交集L≥2^33有

Phi59<T−n/10<Phi56。

这比较的是不同x的完整函数，保留R056正Bbuffer。不能把Phi56≥T−gamma56*n这个过粗下包当作有效比较，也不能把两显示上界排序直接等同于完整cost排序。结论只针对明确列出的这两个内部cost，不声称公开文献记录或穷尽所有策略。

## 9. 新有限证书与归档身份

verify_threshold.py只用Python标准库Fraction和integer，重建本报告新L0的finite guards；不调用作者certify.py，不重新认证64项log证书，不运行旧passed gates、生产扫描、MATLAB或Git。冻结source、full审、局部lemma、历史certificate及R056来源均记录实际SHA256。本次门槛基点验证与原R059历史证书身份分开。

所有统一半径/face界来自解析上包；L≥L0无限尾由明确负powers、N/NL单调及全域Taylor证明，不靠有限盒或数值网格。验证通过仍须另一次独立threshold审计后，才可把新operator门槛写入已接受账本。

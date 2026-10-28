# R070 LOCAL：二阶target profile及逐source-face中心化方差

**局部结论ACCEPT，以下逐项给出新证明。** 本轮不计入C_n定理，不审核完整target计数、径向列、actual-history/L1转移或Cantelli/BE应用。旧冻结c60及其他代理文件不改。后续完整机制的源码接入另写新binding文件，发布后的本证明、代码与运行记录保持冻结。

## 1. 定义与真实参数域

固定s=3/2、0<h≤2^(−20)、q≥2/3、1≤R≤r<3/2、ρ=tq≥4r。给定target r，q在

[(r−h)/s,(r+h)/s]，q_b=(r+h)/s，R=sq−h。

令ε=ζt≥0、a_0=(1−ε)q，g=R−a_0=(s−1+ε)q−h，要求εq≤h/s。于是ε<1/2、0<g/R≤1/3。所有q微分固定n,t,h,ε,r及target坐标w。定义

ell_R(z)=min(1,(R+κρ−|z|)/(2R))，κρ=ρ+a_0，
phi(z)=−logell_R(z)，Q(R,u)=−log(1−u/(2R))，
m=(|w|−a_0)_+，v=(g−m)_+，|w|≤r。

三个一维积分统一记为

J=(2ρ)^(−1)∫capture ell_R^(−1)dz，
A=J_phi=(2ρ)^(−1)∫capture phi/ell_R dz，
B=J_phi2=(2ρ)^(−1)∫capture phi²/ell_R dz，

积分source均在[-ρ,ρ]。本轮B不同于原始source-face V²：V是多个phi之和，故其V²还有交叉项。

baseline时R(q_b)=r。target概率law是固定full target face上的n−1个独立Uniform[-r,r]坐标加一个确定±r坐标。source inverse-face law则在给定target w、source radius及signed face后，对free Uniform[-ρ,ρ]坐标乘capture/ell再除以J。二者不互换，且均在D-cap/soft条件化之前建立。

## 2. 新闭式与Fubini baseline均值

一个source尾的profile原函数为

∫_0^u phi(z)/ell_R(z) dz=R Q(R,u)²，
**∫_0^u phi(z)²/ell_R(z) dz=(2R/3)Q(R,u)³。**

capture保留一条完整tail及另一条长度v的tail，故

**B(w;q)=R[Q(R,g)³+Q(R,v)³]/(3ρ)。**

m>g时额外删除的中央phi=0，所以没有−m项。相对guard及ρ≥4r保证只有一个source端点被删，此式也覆盖|w|>R。所有Q参数均在[0,g]，无负长度或非法log。

在baseline R=r，给定source z被target W接受的比例恰是ell_R(z)。Fubini取消一次inverse ell，得到

**E_W B(W;q_b)=E_fullsource phi²
=(2R/ρ)K_2(a)，**

a=ell_pin=1−g/(2R)，Q_a=−loga，
K_2(a)=∫_a^1(log y)²dy=2−a(Q_a²+2Q_a+2)。

这里不是E_fullsource(phi²/ell)，后者更大且有不同闭式。稳定的非负表达为K_2(a)=∫_0^(Q_a)u²e^(−u)du。

K_2'(a)=−(loga)²。若a=5/6、R/q=s，乘n后的理论上端常数为

**μ_2=3K_2(5/6)=1−5log(6/5)−(5/2)log²(6/5)。**

本机新Fraction用z=1/11的四项positive atanh和及严格几何tail包住log(6/5)，精确得

0.005289<μ_2<0.005290。

完整准确分数上下界保存在`verification.json`；它们不是浮点近似或数值积分。

## 3. 真实multiplicative结构给有限h的target均值下界

仅有a≥5/6只能给mean上界，不能给正的下界。本轮还使用精确关系

a−5/6=(2h−3εq)/(6R)∈[0,h/(3R)]，
R/q=s−h/q≥s−3h/2。

因Q_a²≤log²(6/5)<1/25，

K_2(a)≥K_2(5/6)−h/75，
2(R/q)K_2(a)≥μ_2−h(μ_2+1/25)。

后式用h≤2^(−20)使相乘的下界保持正；它包括ε正号的全部影响。

令S_b=Σ_(i=1)^n B(w_i;q_b)。n−1个free坐标提供随机项；baseline pinned target坐标精确为

B(±r;q_b)=R Q_pin³/(3ρ_b)<1/(250ρ_b)。

因此

ES_b/x≥(1−1/n)[μ_2−h(μ_2+1/25)]，
ES_b/x≤μ_2+1/(250nq_b)。

由本轮准确log区间，h≤2^(−20)、n≥2000、q_b≥2/3足以严格给

**0.00528x<ES_b<0.005295x。**

这些阈值均留有余量，不把pinned target坐标当成另一条随机Uniform law。

## 4. 新uniform导数、AC与wide-h匹配桥

对0≤u≤g≤R/3，有Q≤log(6/5)<1/5，

|Q_R|≤1/(5R)，Q_u≤3/(5R)。

R'=3/2、g'=1/2+ε<1。故完整tail的Q导数绝对值≤9/(10R)；剩余tail的Q导数绝对值≤6/(5R)，因为m>0时v=(R−|w|)_+，v'为0或3/2。

令S_3=Q(R,g)³+Q(R,v)³，则

S_3≤2/125，|S_3'|≤63/(250R)。

prefactor必须同时微分：

(R/(3ρ))'=(h/q)/(3ρ)。

所以

**|∂_qB|≤(21+2h)/(250ρ)<1/(10ρ)，**

特别满足所需safe 1/ρ，而不是只能用8/ρ。

|w|=a_0(q)处一般有真kink，左右B导数之差为(1−ε)Q_pin²/(2ρell_pin)>0。|w|=R处因Q(R,0)²=0而一阶导拼接。有限分支与上述界给闭匹配区间上的绝对连续性，故导数按a.e./AC积分。

若x=n/t≤L，则n/ρ≤3x/2，匹配q宽为4h/3。由强界得到

**|ΣB(w_i;q)−ΣB(w_i;q_b)|≤xh/5。**

h≤2^(−20)时h/5<10^(−6)。所以baseline事件

ΣB(q_b)∈[0.00527,0.00530]x

在所有matching q上保证

**ΣB(q)∈[0.005269,0.005301]x。**

这同时覆盖h=1/(256L)、L≥4096和h=1/(512L)、L≥4096（当然也覆盖8192）。若仅证明宽松1/ρ界，会得到2xh；对4096/256其证书数值超过10^(−6)x，不够闭合。这里新的<0.1/ρ界修复这一预算，并非靠有限差分猜测。

一阶profile的冻结引理在相同真实guards下给|∂_qA|<1/(2ρ)，故其全匹配shift≤xh。h≤2^(−20)<10^(−5)保留baseline一阶good条件到matching条件的原10^(−5)x预算。这只是宽h的guard接入，不重跑旧一阶baseline。

## 5. Target二阶profile概率：显式有限上界

在所有matching点，R<3/2且Q_pin<1/5，故

**0≤B(w;q)<1/(125ρ)。**

baseline pinned target半值为1/(250ρ_b)；该半值只属于真正pinned target坐标，不能拿去删任意source-face pin索引。

free target坐标独立，因而

Var(S_b)≤n/(125²ρ_b²)。

第3节mean区间到[0.00527,0.00530]两侧距离都至少x/200000。对n≥2000、h≤2^(−20)，Chebyshev给

**P{S_b∉[0.00527,0.00530]x}≤5760000/n。**

概率对每个full target face及compact r区间一致。不需要不同good事件独立；也不声称D-cap条件化后或actual μ坐标仍独立。完整机制须在ambient full-face概率建立后再与其他good条件及D相交。

## 6. 逐source-face真实中心化方差上下界

固定一个同时满足matching二阶profile区间的target w和source q。令

A_i=J_phi(w_i;q)，B_i=J_phi2(w_i;q)，J_i=J(w_i;q)。

在captured signed source face j上，pin phi是确定的Q_pin。free坐标在capture-and-inverse law下独立，其均值为A_i/J_i、raw second moment为B_i/J_i。因此**精确条件方差**为

**Var(V|face j)=Σ_(i≠j)[B_i/J_i−(A_i/J_i)²]，V=−logG。**

不能把它写成ΣB_i/J_i，也不能把完整V²的交叉项误加进独立方差和。

设逐点|J_i−1|≤η≤1/1024。相对坐标局部引理有|J−1|<1/(3ρ)，所以q≥2/3、x/n≤1/1000已足以给此η。另有所有坐标的确定界

A_i<3/(50ρ)=C_1/ρ，B_i<1/(125ρ)=C_2/ρ。

在S_2=ΣB_i∈[l_2,u_2]x，l_2=5269/10^6、u_2=5301/10^6时，去掉source pin和减去means²给

Var(V|j)/x
≥[l_2−C_2/(nq)]/(1+η)
 −C_1²[x/n]/[q²(1−η)²]，

Var(V|j)/x≤u_2/(1−η)。

去source pin必须用**任意**B_j≤C_2/ρ，不能误用pinned target的C_2/(2ρ)：source face j可对应一个central target coordinate。pin本身虽有phi²=Q_pin²，但其variance为零；正确做法是上式删去该free-coordinate profile，非加上Q_pin²。

对n≥2000、q≥2/3、x/n≤1/1000与η≤1/1024，本机准确Fraction给

下界≥9776676208/1862312890625>13/2500，
上界≤912/171875<531/100000。

故每个captured source face均严格满足

**(13/2500)x≤Var(V|face j)≤(531/100000)x。**

这个下界不需要一阶profile good条件；它用的是全体一阶profile的逐点O(1/ρ)上界来控制means²。Cantelli所需均值上界、BE第三矩或完整probability应用不是本定理的结论。

所有这些方差是人工ambient inverse-source概率law的条件方差，绝非actual increment measure μ的随机假设。混合face时total variance还有Var_face(E[V|face])，不能直接继承以上upper区间。完整机制可以在每face分别应用适当概率定理，再混合统一概率界；不能把mixture当成独立坐标和。D-cap或soft权条件化后的独立性也不在此处自动成立。

## 7. L≥4096的两个新宽度：局部finite guards

以下只证明本局部引理所需guards，不代替完整机制的finite gate表。取n=e^L、L≥4096、L/2≤x≤L，h_c∈[1/512,1/256]、h=h_c/L，source band下端2/3、上端q_max≤1。若target band按δ=2^(−19)从s q_min及s q_max向内trim，则

* h≤2^(−20)=δ/2，匹配区间位于source band，active R≥1；r上端<3/2须由该源tuple保持。
* e^L≥L²/2给n≥2000、x/n≤2/L<1/1000及ρ≥L/3>6≥4r。
* ε=1/(xn²)≤8/L^5。最窄h=1/(512L)仍有εq_max≤h/s，因为L^4≥6144；于是g/R≤1/3、κ>1及ε<1/2。
* 若要更强J归一eps_J=2^(−20)，由e^L≥L³/6得ρ≥L²/9，故|J−1|<3/L²<2^(−20)。
* Target二阶profile bad概率≤5760000·24/L^4，因e^L≥L^4/24；在4096处为16875/34359738368<2^(−20)，此后递减。

以上适用于R071的8192/256与R072拟用的4096/512，也适用于4096/256这一更宽局部候选；没有盲目转用旧轮维数阈值。若完整机制使用不同trim、r范围、x范围、eps_A分配或其他kernel，仍须单独核相应全链guards。

## 8. 不允许的推广与实际核验身份

仅有R≥1、g/R≤1/3时，可令g/R→0，target二阶mean也趋0。因此**正的双侧profile下端不由relative upper guard单独推出**；本证明必须用精确multiplicative结构及h小。另一项常见失效是用宽松1/ρ导数却在4096/256下声称已付10^(−6) transport，本轮已给更强导数修复。逐face方差界也不是mixture或actual μ的方差界。

本机SymPy 1.14.0一次新核验了phi²/ell原函数、Fubini二阶mean primitive、onset jump、R=|w|处C1及q prefactor。新Fraction准确核验mu_2区间、finite-h/n mean、H2 transport、真实去pin/中心化variance区间及4096新局部guards。未运行作者43个gates或旧passed基线，没有生产grid、共享MATLAB、n=7、下界或Git。

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round070/second_profile/verify.py
```

实际结果见`verification.json`。本LOCAL证明不自行登记为C_n全定理；其源码、代码及结果在发布hash后冻结，后续完整机制的source binding必须另写文件。

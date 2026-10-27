# R064 相对gap一维引理：新导数链、baseline及inverse-face矩

**局部结论ACCEPT，须保持本文件的active域和概率law。** 本任务不审target radial weight、column常数、target计数或完整C_n结论。它不援引旧absolute d≤1/3接口：本轮g可超过1/3，使用的是新相对条件g/R≤1/3。

## 1. 相对参数域与重要负导数

固定s=3/2、0<h≤1/1000、0≤ε=ζt，source q≥2/3、target 1≤r<3/2，且

q∈[(r−h)/s,(r+h)/s]，q_b=(r+h)/s，
R=sq−h∈[r−2h,r]，active R≥1，ρ=tq≥4r，
a_0=(1−ε)q，g=R−a_0=(s−1+ε)q−h。

要求εq≤h/s。它给ε≤h≤1/1000，以及

g−R/3=εq−2h/3≤0，0<g≤R/3。

g>0也由q≥2/3、h≤1/1000直接成立。闭式中g的+ζρ项不可丢掉。q_b是右匹配端点，R(q_b)=r；q=r/s时R=r−h，不能用作mean-one baseline。

ell_pin=(R+a_0)/(2R)=1−g/(2R)≥5/6，而

**∂_q ell_pin=−(1−ε)h/(2R²)<0。**

这与旧additive参数中的正号不同，不能照搬旧计数endpoint方向。本文只给该导数与guard，不负责后续count校准。

候选full source band可取q∈[2/3,909/1000]，但其中R<1的inactive点不在本引理内。匹配到r≥1+δ且h<δ/2的active点满足R>1；inactive band上定义source weight或证明column为零属于完整机制的另外责任。

## 2. 精确J与稳定Psi表达

令κρ=ρ+a_0，

ell_R(z)=min(1,(R+κρ−|z|)/(2R))，
J(w;q)=(2ρ)^(−1)∫_[-ρ,ρ]1_{|z−w|<κρ}/ell_R(z)dz，|w|≤r。

定义

F(R,u)=−2Rlog(1−u/(2R))−u，0≤u≤g，
m=(|w|−a_0)_+，v=(g−m)_+。

则精确闭式为

**J=1+[F(R,g)+F(R,v)−m]/(2ρ)。**

a_0>0使capture至多删掉一侧source端点；m≤g+2h≤R/3+2h≪ρ。因此公式也覆盖|w|>R、m>g时的中央段删减，后一情形F(R,v)=0。所有F参数都在[0,g]，不对负数或越界参数取log。

数值稳定的非负表达仍是

F(R,u)=∫_0^u z/(2R−z)dz
=2RΣ_(k≥2)(u/(2R))^k/k。

因为u/(2R)≤1/6，尾部可用几何级数准确向上包络；小gap时无需用两个近等数直接相减。

## 3. 新multiplicative导数链的显式证明

现在R'=3/2、g'=1/2+ε<1、a_0'=1−ε，而非旧链R'=1、g'=ε。令η=u/(2R)≤1/6，则

0≤F_u=η/(1−η)≤1/5，
|F_R|≤2η²/(1−η)≤1/15，
0≤F≤u²/[2(2R−u)]≤R/30。

故F_0=F(R,g)满足|F_0'|≤(3/2)(1/15)+1/5=3/10。

写B=F_0+F(R,v)−m。m=0时|B'|≤3/5；m>0时v=(R−|w|)_+，v'为0或3/2，而−m'=1−ε≤1。因此a.e.

|B'|≤3/10+(3/2)(1/15+1/5)+1=17/10。

同时R/q≤s、q≥2/3给

|B|/q≤R/(15q)+R/(3q)+2h/q
≤3/5+3h≤603/1000<1。

J≥1−m/(2ρ)≥1−1/24−h/(4r)≥11497/12000>1/2。由于ρ'=ρ/q，

J'=[B'−B/q]/(2ρ)，
**|∂_qlogJ|<(17/10+603/1000)/ρ<3/ρ≤256/ρ。**

这是无K_w条件的新unconditional导数界，不是将旧d域换名。

在|w|=a_0(q)处左右J导数之差为(1−ε)/(2ρell_pin)>0，一般存在真kink；在|w|=R处因F_u(R,0)=0，一阶导拼接。有限分支及上述界给局部Lipschitz与闭匹配区间上的绝对连续性。导数只能按a.e./AC使用，不能声称处处classical可微。

若n/ρ≤(3/2)L，匹配区间长度≤2h/s=4h/3，逐坐标积分得

**|Σ_i logJ(w_i;q)−Σ_i logJ(w_i;q_b)|
≤256·(3/2)L·(4h/3)=512Lh。**

h=2^(−64)/L时该量精确为2^(−55)。本轮不使用K_w。x=n/t=909L/1000、q≥2/3可给n/ρ≤2727L/2000<(3/2)L；其他t选择必须另保证同一guard。

## 4. Baseline target law与3/6/12矩

**本节W是target Uniform[-R,R]，不是source Uniform[-ρ,ρ]。** 仅在q=q_b、R=r时适用。固定source z时，其接受的target坐标比例恰为ell_R(z)，Fubini因此给E_WJ=1。

令T=(|W|−(R−g))_+、B_0=F(R,g)。闭式等价J−1=(B_0−A(T)/2)/ρ，其中g≤R/2给

B_0≤g²/(3R)，A(T)≤(4/3)T，ET²=g³/(3R)。

于是

E(J−1)²≤2g⁴/(9R²ρ²)+8g³/(27Rρ²)
≤(11/(27R))g³/ρ²≤3g³/ρ²。

此外|J−1|≤g/ρ≤1/8，故|logJ|≤2|J−1|、logJ≥(J−1)−2(J−1)²。结合E(J−1)=0与Jensen，严格保留

**−6g³/ρ²≤ElogJ≤0，Var(logJ)≤12g³/ρ²。**

这里仅用R≥1、g≤R/2、ρ≥4R，从未要求g≤1/3。本轮更强g/R≤1/3当然满足。pinned target坐标W=±R的确定项仍在[−g/ρ,0]。

对每个固定target r，N=n−1个free坐标独立，故全face log和S_b有

|ES_b|≤6N g³/ρ_b²+g/ρ_b，Var S_b≤12N g³/ρ_b²。

因为r<3/2、g<R/3<1/2且ρ_b≥2t/3，这些量一致为O(L²/n)+O(L/n)，趋于0；可在固定full target face上构造baseline typical set，再用第3节确定性桥控制所有matching q。此处不负责从该集合计算完整target annulus质量。

全匹配还有逐点界

**|J−1|≤max(F_0,m/2)/ρ≤(r/6+h)/ρ<1/(3ρ)。**

所以任意固定epsJ>0，足够大ρ后J≥1−epsJ，且无需对n个坐标作union bound。

## 5. Fullsource phi矩及新精确率

**本节Z是fullsource Uniform[-ρ,ρ]。** 令phi=−logell_R(Z)、Q_pin=−logell_pin≤log(6/5)。两个尾合计密度du/ρ，换元ell给

E_full(phi/ell)=R Q_pin²/ρ，
E_full(phi²/ell)=2R Q_pin³/(3ρ)。

中央ell=1的部分为零贡献。对captured free coordinate，inverse density相对于fullsource是1_capture/(ell J_i)。若J_i≥1−epsJ>0，非负性允许去indicator，故

μ_i≤B_1=R Q_pin²/[ρ(1−epsJ)]，
E_inv,i phi²≤B_2=2R Q_pin³/[3ρ(1−epsJ)]。

在固定captured signed source face，pinned phi为同一个确定Q_pin，其余n−1个坐标独立。因此V=−logG满足

E(V|face j)=Q_pin+Σ_(i≠j)μ_i≤Q_pin+(n−1)B_1，
Var(V|face j)≤(n−1)B_2。

epsJ=1/1024，R/q≤s及n/ρ=x/q给

E(V|face)≤Q_pin+x·s Q_pin²/(1−epsJ)，
Var(V|face)≤x·(2s/3)Q_pin³/(1−epsJ)。

不需要浮点log：对z=1/5，五项交错Taylor上界

Q_pin≤log(6/5)≤Q_5=z−z²/2+z³/3−z⁴/4+z⁵/5
=34187/187500。

准确Fraction核验得到

s Q_5²/(1−epsJ)=37400031008/749267578125<999/20000，
(2s/3)Q_5³/(1−epsJ)=639297430035248/105365753173828125<1/150。

故可统一使用每face mean≤Q_5+(999/20000)x、variance≤x/150。

## 6. Perface Chebyshev再混合，不能省略law区别

每face均值可能不同。对M=Q_5+(999/20000)x及任意u>0，先分别应用

P(V>M+u|face)≤(x/150)/u²，

再按inverse-row face概率混合，得到相同总概率上界。这个过程不要求所有face均值相同，也没有声称mixture variance等于within-face variance之和。

若一定要给无条件variance，total-variance须另加Var_face(mean)。特殊共同μ_i结构只允许估计该项≤B_1²/4（face j均值=Q_pin+Σ_iμ_i−μ_j），不能无依据删去。

上述独立law只在**fullsource、固定rho、capture和inverse G权**下建立；原corner cap条件化或soft weight改变law后不能继续默认为独立。完整机制应先在此law证明retention，再另外支付cap删除。target radial factor在固定target上的作用及列估计不属于本引理的审核结论。

## 7. 实际新核验与责任范围

本机SymPy 1.14.0新核验R'=3/2、g'=1/2+ε、relative guard等价式、ell_pin负导数、F完整新chain、onset jump及R=|w|处C1。Fraction检查新relative常数与512Lh=2^(−55)，并准确核验Q_5的两个新率预算。没有重跑旧baseline积分、Terra数值批或旧已通过案例；一般结论来自上述解析证明，非数值点外推。

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/audit/round064/relative_coordinate/verify.py
```

结果见`verification.json`。只写本任务目录，不用共享MATLAB、n=7、下界、网络或Git。完整新算子定理留给其负责审计者，源稿身份对照将另作局部增量锁。

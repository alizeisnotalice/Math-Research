# R066 first-profile：独立新局部引理

**Local ACCEPT。** 本文件只证明新的J_phi闭式、target profile矩与概率、全匹配导数桥，以及profile对每个inverse source face均值的控制；不审核target计数、径向column、actual-history转移或完整C_n结论。源定稿身份另写`SOURCE_BINDING.md`，不再改动此冻结数学引理。

## 1. 新区间的完整guards与两类概率law

固定s=3/2、0<h≤1/1000、ε=ζt≥0。source q∈[2/3,917/1000]，target r∈[1+δ,2751/2000−δ]，δ=2^(−30)。要求h<δ/2、εq_max≤h/s、ρ=tq≥4r，并在全匹配区间

q∈[(r−h)/s,(r+h)/s]，q_b=(r+h)/s

定义R=sq−h、a_0=(1−ε)q、g=R−a_0=(s−1+ε)q−h。于是active R∈[r−2h,r]且R>1、r<3/2，ε≤h<1/2，并且

g−R/3=εq−2h/3≤0，0<g≤R/3。

这里保留+ζρ项，不使用旧absolute g≤1/3域。新q_max=.917与旧.909不同，但真正使用的局部条件只有q≥2/3、R≥1、r<3/2、ρ≥4r、g/R≤1/3和上述匹配链，故此处重新核对后可接入。inactive source可能R<1，不能调用需要R≥1的baseline矩。

baseline q_b满足R=r。**target profile概率空间**为给定target r/face后，n−1个W_i独立Uniform[-r,r]，另一个pinned坐标固定±r。**source inverse-face概率空间**则给定target w及source rho/face后，在Uniform[-ρ,ρ]上乘capture/ell并除以J。二者不是同一概率law，均在原corner-cap条件化之前建立。

## 2. 精确闭式与baseline target均值

令

ell_R(z)=min(1,(R+κρ−|z|)/(2R))，κρ=ρ+a_0，
Q(R,u)=−log(1−u/(2R))，H(R,u)=R Q(R,u)²，
m=(|w|−a_0)_+，u_0=(g−m)_+。

新profile定义为

J_phi(w;q)=(2ρ)^(−1)∫_[-ρ,ρ]1_{|z−w|<κρ}
                   [−logell_R(z)]/ell_R(z)dz。

tail长度u的精确profile积分是

∫_0^u Q(R,v)/[1−v/(2R)]dv=R Q(R,u)²。

capture保留一条完整source tail、另一条剩余长度u_0，所以

**J_phi=[H(R,g)+H(R,u_0)]/(2ρ)。**

m>g时被删的额外中央段phi=0，故没有J中出现的−m项。此式也覆盖|w|>R及源一侧tail被完全删去的情形。

在baseline R=r，固定source z被接受的target W比例为ell_R(z)。Fubini因此严格给

**E_W J_phi(W;q_b)=E_Z phi(Z)
=(2R/ρ)[1−a+a loga]，a=ell_pin=1−g/(2R)≥5/6，**

其中W~Uniform[-r,r]是target law，Z~Uniform[-ρ,ρ]是fullsource law。它不等于source inverse均值，也不涉及未知actual increment measure。

K(a)=1−a+a loga的导数为loga<0。因此用R/q≤s与n/ρ=x/q，

n E_WJ_phi≤2sK(5/6)x
=[1/2−(5/2)log(6/5)]x=:μ_0x。

此处n只是把n−1个free坐标上界放大；pinned target项另付。由于μ_0含负的log系数，证明其上界必须使用log的**下界**。新的准确Fraction检查取八项even alternating下界

log(6/5)≥3988283/21875000，
μ_0≤386717/8750000
=221/5000−33/8750000<221/5000。

不能在这个方向错误地使用odd Taylor上界。

## 3. 新a.e.导数链与全匹配profile transport

对0≤u≤g≤R/3，

0≤Q≤log(6/5)<1/5，
|Q_R|=u/[R(2R−u)]≤1/(5R)，
Q_u=1/(2R−u)≤3/(5R)。

R'=3/2、g'=1/2+ε<1，所以|[Q(R,g)]'|≤9/(10R)。对于Q(R,u_0)，m=0时同此界；m>0时u_0=(R−|w|)_+，其导数为0或3/2，故统一有|[Q(R,u_0)]'|≤6/(5R)。

写S=Q(R,g)²+Q(R,u_0)²，则S≤2/25、|S'|≤21/(25R)。重要的是q依赖的prefactor也要微分：

(R/(2ρ))'=(R'−R/q)/(2ρ)=(h/q)/(2ρ)。

q≥2/3给

**|∂_qJ_phi|≤(21+3h)/(50ρ)<1/(2ρ)，**

因而源机制需要的8/ρ更保守且成立。另一种等价的H偏导界为H≤R/25、H_u≤6/25、|H_R|≤3/25；它也直接支持作者的宽松8/ρ链。

在|w|=a_0(q)处，profile左右q导数之差为

(1−ε)Q_pin/(2ρell_pin)>0，Q_pin=−logell_pin，

所以一般有真kink。|w|=R处Q(R,0)=0使Q²的导数平滑拼接。J_phi连续、分段C1且统一有界导数，故闭匹配区间上绝对连续；导数界必须按a.e./AC解释。

若0<x=n/t≤L，q≥2/3给n/ρ≤3x/2≤3L/2，匹配q宽度为4h/3。因此

**|Σ_iJ_phi(w_i;q)−Σ_iJ_phi(w_i;q_b)|≤16xh≤16Lh。**

h=2^(−80)/L时16Lh=2^(−76)。它可付额外x/100000：例如用更强16xh，只需h≤1/1600000；或者用固定2^(−76)，在x≥1时也足够。此为全连续q确定性桥，无须网格、K_w条件或实际μ浓度假设。

## 4. Target profile的uniform矩与good-event概率

所有target坐标|w|≤r、matching q都有

**0≤J_phi(w;q)≤R Q_pin²/ρ<3/(50ρ)。**

baseline pinned target w=±r时m=g、u_0=0，故

J_phi(±r;q_b)=R Q_pin²/(2ρ)<3/(100ρ)，

是确定项，不另赋Uniform law。free target坐标的profile方差至多9/(2500ρ_b²)。独立性只在固定full target face上使用，因此全profile和的方差≤9n/(2500ρ_b²)=O(n/ρ_b²)，均值≤μ_0x+3/(100ρ_b)。这些估计在整个compact r区间中一致。

以下给一个充分而非最优的明确概率界。令c=221/5000、η=1/100000。因ρ_b x=nq_b且q_b≥2/3，n≥9000保证pinned项≤ηx/2。于是

EΣJ_phi≤c x+ηx/2，
P{ΣJ_phi(w_i;q_b)>(c+η)x}
≤324000000/n。

该上界可能在小n平凡，但n→∞时一致趋零，足以构造baseline target profile event。它不是对已经D-cap-conditioned target或actual μ的独立性结论。profile与其他target good事件可取交，概率不需要彼此独立。

结合第3节transport，在good target和所有matching q上，

ΣJ_phi(w_i;q)≤(c+2η)x=(2211/50000)x。

## 5. 每个source inverse face的均值与安全方差

固定一个good target w及source q/rho。对每个free source坐标，inverse-capture density相对于fullsource Uniform[-ρ,ρ]是

1_capture/(ell_R J_i)，J_i=J(w_i;q)。

因此其phi均值**恰好**为J_phi(w_i;q)/J_i；此步骤没有丢capture。若统一J_i≥1−eps_J、eps_J=2^(−20)，任一captured signed source face j有

E_inv(V|face j)=Q_pin+Σ_(i≠j)J_phi(w_i;q)/J_i
≤Q_pin+[2211/(50000(1−eps_J))]x
<Q_pin+(223/5000)x。

最后一式已准确核验：2211/[50000(1−2^(−20))]=4390912/99296875<223/5000。省略的坐标是哪一个均无妨，因为所有profile非负，故全n项profile同一上包控制每个source face。target good条件只固定w，并不破坏给定w下source free坐标的独立性。

方差仍可使用fullsource二阶矩上包

E_full(phi²/ell_R)=2R Q_pin³/(3ρ)。

去掉capture indicator、除以J_i下界并对固定face独立free坐标加方差，给

Var_inv(V|face j)≤x[(2s/3)Q_5³/(1−eps_J)]<x/150，

Q_5=34187/187500是log(6/5)的odd-five上界。这里upper方向应使用log上界，与target μ_0处的lower方向不同；新Fraction同时核验了该二阶率。

χ=9/200且Q_pin/x≤1/5000（例如x≥1000即可），有

χ−223/5000−1/5000=1/5000。

所以对**每个**captured source face分别应用Chebyshev，

**P_inv(V>χx|face j)≤5000²/(150x)=500000/(3x)。**

再按inverse-row face权重混合，保持同一概率上界。不同face均值不必相同，不能把mixture variance直接等同于独立variance和。以上全部在fullsource、固定face、capture-and-inverse law下完成；source D-cap或soft b条件化后不再自动独立，应由全机制另付删除/retention。

## 6. 旧相对J接口的真实扩域说明

冻结R064相对引理c60中的抽象域为q≥2/3、active R≥1、r<3/2、ρ≥4r、εq≤h/s；它不以qmax=.909作数学上界。R066新source band .917与target上端1.3755仍满足这套真实guards，所以J基线和logJ导数可合法接入。这里核对的是guard，而非只凭旧轮名称或数值样本继承。R064冻结文件不作任何修改。

R066的J_phi本身则是本文件的新证明；target一阶均值身份、概率与source mean接口并未被旧J引理覆盖。

## 7. 实际新核验及冻结身份原则

本机SymPy 1.14.0一次新核验H=R Q²的tail原函数、target fullsource均值原函数、q prefactor导数、capture-onset profile jump及|w|=R处C1。Fraction新核验even-eight lower-log对μ_0的严格gap、target/bridge系数、source J归一预算、eps_J下的二阶率、源阈值gap、Chebyshev系数及target pinned记账。结果见`verification.json`。

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/audit/round066/first_profile/verify.py
```

没有重跑旧passed J/baseline、生产grid或Terra数值批；有限计算不是全域证明。只写本新目录，不用共享MATLAB、n=7、下界或Git。此数学引理、代码及运行记录在发布hash后冻结；正式R066源的接入身份单独写`SOURCE_BINDING.md`，不向此处追加版本史而改变既有依赖。

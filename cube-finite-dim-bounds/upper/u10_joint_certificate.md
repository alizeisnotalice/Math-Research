# U10：两个空间bump的全域联合修正充分接口

2026-10-25，固定三维体对角frame V={(1,±1,±1)}，L=−Σ_v|v·D|，各方向未单位化，jump率ds/(pi s²)。只给可接区间证书的充分接口，不给可行权重或新上界。已审U09的真实G、H_c、μ及全L1/尺度定理作为明确输入，不重证026。所有积分为Lebesgue。

## 1. 固定量、基与支持位于Ω

G=∫p_tdt，T0=4log2−1，c=T0/(24pi³)=G(1,1,0)/3，H=(G−c)+，LH=μ−δ0，μ为正概率、μ0=0。此c是固定edge值，不预设该edge是全cube最小点。

b0(x)=(1−|x|²)+²，b1(x)=Σ_(i=1..3,±)(1−4|x∓(3/4)e_i|²)+²。

令w0,w1≥0，B=w0b0+w1b1，W=w0+w1，a=5/4，F=H−B，K=F/D（D>0）。两基even且cube对称。

六小球radius1/2，中心orthogonal距离平方9/8、opposite距离平方9/4，均大于(2radius)²=1。因此闭支持互不相交，0≤b1≤1、0≤B≤W；所有支持S_B⊂closed B(0,a)。

输入主round009_bump_source_certificate.md给正Feynman表示G=3/(32pi³)∫_simplex 1/V及EV=|x|²/20、simplex体积1/6。Jensen给全域（允许∞）

G(x)≥5/(16pi³|x|²)。

这条已定位输入用于支持证明，不由采样G推出。log2<1，故c<1/(8pi³)。当|x|≤3/2时，G−c>5/(36pi³)−1/(8pi³)=1/(72pi³)>0。所以closed B(0,3/2)⊂Ω={G>c}；S_B⊂⊂Ω且dist(S_B,Ω^c)≥1/4。支持上更有G−c≥1/(5pi³)−c>3/(40pi³)。

## 2. C1,1、真实全空间源与质量

unit bump在球内gradient=−4x(1−|x|²)，Hessian=−4(1−|x|²)I+8xx^T；值和gradient在球面均连续为0，Hessian a.e. bounded（operator norm≤8），故b0∈C_c^{1,1}。平移/缩放和有限求和给b1同类（gradient Lipschitz可取32）。没有额外球面atom。

对任何这样的b，二差分以min(C|v|²s²,4||b||∞)控制sup，以min(|v|²s²||D²b||1,4||b||1)控制L1。故Lb为普通连续有界L1函数，并与分布L一致；∫Lb=0（逐方向平移Fubini）。因此

LF=(μ−LB dx)−δ0，ν:=μ−LB dx为有限signed测度、ν(R³)=1、ν({0})=0。

这里LB=w0Lb0+w1Lb1是带符号的总和；不同基的负部分可以支付其他基的正部分。不得替换成逐基独立支付。

准确质量及缩放为

M0=∫b0=32pi/105，M1=∫b1=6(1/2)³M0=8pi/35=(3/4)M0，

L[b0((x-z)/R)]=R^-1(Lb0)((x-z)/R)。

所以∫F=∫H−w0M0−w1M1；这只在全条件成立后才成为有效上界收益，不重复付六球的质量。

## 3. 完整充分条件与Ω外的自动支付

以下任意D>0、w0,w1≥0满足(A)(B)(C)时，K是合法cube核。

(A) 全域非负：F≥0 a.e.。只需在S_B核G−c−B≥0；其外B=0、H≥0自动。粗充分条件W≤1/(5pi³)−c可用，但不强制用它取代更精细联合区间。

(B) 全cube支配：F≥D/8 a.e.于[-1,1]³。S_B内需G−c−B≥D/8；cube\S_B需H≥D/8。D不能未经全域认证取8(Gedge−c)。如果已证g=essinf_cube G，可用D≤8(g−c)再检验减项，但本命题不假定g=edge。

(C) 只需在Ω内核联合源密度

q(x):=m(x)−w0Lb0(x)−w1Lb1(x)≥0，

其中m为μ在Ω的密度（下节）。点态核或a.e.核均可，因为q连续于Ω。

**为何不需另外手付Ω外的源。** U09给G连续有限于四frame线之外，所有frame线属于Ω；故Z_-={G<c}开放，Z0={G=c}由正齐次degree−2及每ray唯一level半径给零Lebesgue体积。S_B离Ω^c有正距离，B在其邻域恒0。在Z_-中F=0；若(A)已证，则对每个紧支撑于Z_-的非负测试，小jump项为0、大jump项是非负F的有限率卷积，由F L1和Fubini得LF≥0于Z_-。而Z0上LB dx不赋质量，所以ν|Z0=μ|Z0≥0。于是ν在整个Ω^c非负；不必另加Ω外源约束。

这自动结算必须使用(A)、支持buffer、Z0零体积及原源μ≥0，不能仅说“b有紧支撑所以Lb无外尾”。Ω内即使B(x)=0，LB(x)仍可能为正，仍须(C)。

由(C)和Ω外结算，ν为正概率。于是LK=(ν−δ0)/D≥0原点外，K≥0、even、L1且(B)给K≥1_cube/8。U09已审finite-frame比较立即给一般f∈L1及所有共同r>0的弱型常数

C(K)=[∫H−w0M0−w1M1]/D。

同一u/sigma/D_bad的全r结论沿用U09，不因本分母D重选坏集。四维的抽象C1,1/Ω外论证可依法另设基后使用；本文件具体支持、常数及有限尾只登记n3。

## 4. Ω内密度、近line及积分区间接口

h=min(G,c)为全局bounded continuous函数：G在line外连续，line上G=∞且lower semicontinuity使h在每个line点邻域恒c。μ=−Lh。Ω中h局部恒c，故μ有连续、locally bounded非负密度

m(x)=(1/pi)Σ_v∫_R (c−G(x+s v))+/s² ds，x∈Ω。

这包括line点甚至原点，不把G的∞值直接喂入不稳定公式；(c−∞)+定义0。若紧集E⊂Ω已证dist(E,Ω^c)≥δ，则|s|<δ/√3的各direction项严格0，无PV；截|s|>S的总尾位于[0,8c/(pi S)]。在剩余有限s段可接bridge的G/μ区间。未证明E⊂Ω的cell不能套m公式。

**可显式使用的近line屏障。** 对任意frame v，写x=t v+z、z⊥v、|t|≤A、ρ=|z|∈(0,1)。包含v的三列basis det4、inverse每行Euclidean norm≤1，剩余第四方向的basis坐标全±1。该剩余Cauchy变量落|C_u|≤u的概率恰1/2。对ρ≤u≤1，两个transverse basis密度各≥1/(5pi u)，v因子≥u/[pi(1+(A+2)²)]，卷积得

G(x)≥C_A log(1/ρ)，C_A=1/[200pi³(1+(A+2)²)]。

ρ=0按G=∞。因此ρ≤exp(−M/C_A)明确保证G≥M，可给line邻域的H−B≥D/8屏障（取M≥c+W+D/8），或给h恒c及small-jump零区。常数很粗，但是真实有限接口；可用更锐已有G表示替换，不可删去近line cells。

## 5. 明确有限的无穷尾归约

记S=∪_v span(v)，T_a=∪_v{dist(x,span(v))≤a}。若x∉T_a，则每条x+Rv均不碰B的支持，LB(x)=0，源自动μ≥0。只需考虑Ω∩T_a。

先准确给G的线外上界。对dist(y,S)≥d，取δ=d/(8√3)。如U09已审pair-law论证，X_t落该区域必至少两个Cauchy参数>δ。六pair分别补成det4 basis，其受限law不归一、density各≤t/(4pi³δ⁴)，因此

p_t(y)≤A_d t，A_d=6/(4pi³δ⁴)，同时p_t≤C t^-3，C=1/(4pi³)。

∫_0∞min(A_d t,Ct^-3)dt=√(A_d C)，所以

G(y)≤48√6/[pi³ dist(y,S)²]。

对x∈T_a取x=t v+z、z⊥v、|z|≤a；可由even/cube对称化取t>0、v=(1,1,1)。选w=(1,1,−1)，τ∈[1,2]时dist(v+τw,S)≥√(8/3)。t≥5/2给dist(v+τw+z/t,S)>1，因而

G(x+s w)≤48√6/(pi³t²)，s∈[t,2t]。

由log2>2/3（对1/x在[1,2]严格Jensen）、√6<5/2，可核t≥60时上述上界≤c/2。故对Ω内这样的x，m(x)至少该单direction贡献c/(4pi t)。

另一方面x在支持外，LB是非负jump-in。每direction打到ball B(0,a)的参数区间长度≤2a/√3，且|s|≥t−a/√3≥t/2（t≥5/2）。0≤B≤W，四direction求和给

0≤LB(x)≤40W/(pi√3 t²)。

定义显式T=max(60,160W/(√3 c))，R_*=√3T+a。若|x|≥R_*且x∈T_a，则|t|≥T；所以在Ω∩{|x|≥R_*}有m−LB≥0。Ω外已经第3节自动，T_a外LB=0，故**源条件(C)只剩有界集合Ω∩T_a∩B(0,R_*)需要区间认证**。

这是所有远tail/nearline tube的完整解析支付，不是有限源点替代无穷域；不要求所查有限点呈现正slack。T依权重可很大，但有限且公式可准确核。不得将这个充分尾半径误称候选通过；有限区域条件仍可能失败。

## 6. 可接证书的输入输出与未闭合项

一个候选认证包须锁定w0,w1,D（exact rational或严格interval）、c及所有来源hash，并提供：(i) S_B上的联合F非负；(ii)全cube F≥D/8，包括不在S_B部分；(iii)上述有界源域中q≥0的联合区间覆盖。phase-uncertain cell需先证Ω归属或用更一般分布/正核结算，不能套错m式；nearline cell可用屏障和zero small-jump；无穷s积分尾用8c/(pi S)或更锐已证界。源域外由本定理自动，不额外采样。

满足(i)–(iii)时给可核验一般L1/all-scale上界及精确质量式。未给这样的权重/覆盖即没有新上界；主LP的样本/容差warnings均不能补缺。本接口不是存在低于8修正的证明，也不排除两个基或更广联合修正。恒定Lebesgue实际用于密度/卷积/Fubini、源总质量、普通球质量、level-set零体积及cube归一；不能移植非恒定坏背景。

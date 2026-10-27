# R042 · Second-moment cut的显式门槛：s≥2^36

## 结论与冻结依赖

**显式门槛通过。** 令L=log n、s=L^(1/3)。对每个满足s≥2^36（即L≥2^108）的整数n，R042 second_moment_cut的全部eventual条件均已具体化，并严格得到

C_n≤F(A_n(t),B_n(t)−2^(−31)n/s)
 <nlogn−nloglogn−2^(−32)n/s，                                  (1)

t=n/(L−s)。本报告不主张该门槛最小。

冻结full proof：
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round042/second_moment_cut/PROOF.md`

SHA-256：`722764654cd15a8bf7fb18a799c298b54ced2f0188f5554ccae173240ddf90a2`。

冻结的一维moment审计：
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/second_moment_audit/PROOF.md`

SHA-256：`425a383952f341f42ff0dee9651924b00fff61e958c2957f355e8cc9fdc72d00`。

以下保持full proof的全部常数不变：cE=2^(−10)、CR=36、C2=144、H=2^20。较强moment常数48只用作已验证的充分依据，不替换其144。没有Arb搜索、生产计算、七维或Git/网络操作，也不改上述冻结文件。

## 1. 统一参数及基础域条件

置

S=2^36，s≥S，L=s³，n=e^L，d=1/s，
x=s³−s，t=n/x，β=1+1/t，ζ=n^(−3)，κ=β−ζ，
h=1/(100s²)，p=e^(−s/4)，M0=nhp。

目标和source半径为

r=1+u，q=ρ/t=1−v，u,v∈I=[d/4,3d/4]，|r−q−d|<h。

所有基础条件已由s≥8保证：x∈[L/2,L]、x>1，d≤1/8，q≥1/2、r<2。指数正项给

n≥L²/2=s^6/2，t≥n/L≥L/2=s³/2，ρ≥L/4>4。

所以target/source支撑不交，ρ>r；2d<r及2d<ρ使固定面矩公式与K统计正确。

真实κ shift满足

ζt=1/(xn²)≤8/s^15≤h/2，                                     (2)

最后一步只需s^13≥1600，在s≥2已成立。因此h+ζt≤3h/2<d/2，κ>1（因ζt<1），且κ<β<2。

反对角带的最短径向截面为h，而不是2h：其中心d−u或d−v均在I；h/d=1/(100s)<1/2使端点也至少有h长度。source log长度介于h和4h，target log长度至多2h。

## 2. Target统计和质量的显式界

固定一张target r面，固定坐标贡献d，其余坐标独立Uniform[−r,r]。S_w准确均值m_w=d+(n−1)d²/(2r)，方差≤nd³/(3r)。容忍t/100的Chebyshev坏率满足

P_bad,S≤10000nd³/(3r t²)≤20000/(3L)<1/4，                       (3)

因为x²/n≤L²/n≤2。最后不等式由S³>80000/3保证，且随s递减。

nd≥s^5/2>2，故E K_w≤1+2nd≤3nd；在12nd阈值的Markov坏率≤1/4。每张target面至少半数保留。

完整target annulus质量≥nd/4。对s≥4，

1/r−1/t≥1/(1+3/(4s))−2/s³>1/2。

右侧随s递增，其s=4基点为16/19−1/32>1/2。这使原corner交长ell_t(r)<3/4。n≥6给(3/4)^(n−1)≤(3/4)^5<1/4，所以被原corner捕获的target质量≤nd/16。于是显式有

nd/16≤a*=λ(A)≤nd。                                            (4)

所有target点都在r<t的原cap比较范围内，因t≥s³/2>2。

## 3. 全部matched row产品的统一e^(±2)包围

### 3.1 均值中心

在整个I×I上|rq−1|≤d且rq≥1−d。因此

|nd²/(4rρ)−s/4|
 =|(s−d)/(4rq)−s/4|
 ≤(1+d)/[4(1−d)]<1/3                                         (5)

对s≥8成立（基点上界9/28<1/3）。

### 3.2 Good-set容忍及阈值shift

写e=r−q−d、b=e+ζρ，则|b|≤h+ζt<d/2。所有受阈值移动影响的坐标都被K_w=#{|w_i|>r−2d}记录。exact row缺失和T_w因此满足

|T_w−nd²/(4rρ)|
 ≤[t/100+12nd(h+ζt)+2d]/(2ρ)
 ≤1/100+18Ldh+4d/L
 =19/100+4/s^4<1/5。                                         (6)

这里用q≥1/2、h+ζt≤3h/2、x≤L及t≥L/2。s≥8时最后一项最多1/1024，故严格小于1/5。

### 3.3 坐标平方误差与任意deleted product

每坐标归一化损失v_i=1−L_i非负，且

v_i≤(d+h+ζt)/(2ρ)≤4/s^4=:V。

由(5)–(6)，|T_w−s/4|<8/15，故s≥8时T_w≤s/2。V<1/2给

Σ_i v_i²/(1−v_i)≤[V/(1−V)]T_w≤4/s³。

删除任意一个坐标后，upper使用log∏_{i≠j}(1−v_i)≤−T_w+v_j≤−T_w+V；这里不声称相对完整log产品的增量−log(1−v_j)≤v_j。lower仍可使用完整乘积的lower。因此所有完整或删除一个坐标后的log产品距−s/4的额外误差至多

V+[V/(1−V)]T_w≤4/s^4+4/s³<1/10                               (7)

（基点s=8为9/1024，尾部是递减powers）。合并(5)–(7)，总误差小于1/3+1/5+1/10=19/30<2。故冻结full proof中的e^(±2)p产品包围对所有good target、全部matched source半径显式成立。

local source box每坐标至少一个面可用，是因为ρ>r且κ>1。normalized交长与面捕获的公式均在上述正域中，不需要另加eventual前提。

## 4. 原corner source的删除量

n≥4、x≥L/2、v≥d/4给

(n−1)(1/ρ−1/t)/2
 =(1−1/n)xv/(2q)≥3Ld/64≥Ld/32=s²/32。

完整source annulus的log长度≤d，因此其原corner捕获质量至多(nd/2)exp(−s²/32)。除M0可用稍宽的

Corner/M0≤100s exp(−s²/32+s/4)。

s≥16时−s²/32+s/4≤−s²/64。Taylor的单个正项e^u≥u²/2于是给

Corner/M0≤100s e^(−s²/64)≤819200/s³<1/36。                      (8)

最后不等式在s=S由Fraction严格核验，之后右侧随s递减。没有把一个固定阈值的数值指数样本外推为全尾。

完整source row下界为(e^(−2)/2)M0；扣去(8)，使用e²<9得1/36<e^(−2)/4，所以raw row至少r*=e^(−2)M0/4。上界4e²M0<36M0。由(4)，

E_raw≥a*r*≥e^(−2)ndM0/64>cE ndM0。                            (9)

cE=1/1024可由e²<9及9·64=576<1024直接比较。

## 5. Column二阶矩与energy截断

冻结moment审计已经对全部L≥64证明同一参数与同一Y具有

∫Y²dν_source≤48ndM0²。

本门槛L≥2^108当然包含L≥64。其域条件pinned source face、r_-=q+d−h≥1、Dgap=d−h+ζρ∈(0,d)、ρ>4以及κ>1均已由第1节和(2)明确保证。因此原full proof的C2=144也成立，未省略O(1)补偿。

原source限制到D∞只减小moment。按Y>HM0删列时删除的是energy而不是source总质量，安全方向为

E_removed≤∫_{Y>HM0}Ydλ≤(HM0)^(−1)∫Y²dλ≤(C2/H)ndM0。

精确整数比较

H=2^20≥4·144·1024=4C2/cE，H≥36

给E_good≥(3cE/4)ndM0，两侧degree均≤HM0。Y及截断集合都是联合Borel核的Borel限制，不需要可测选择器或μ浓度假设。

## 6. τ、预算正性与真实有限a接口

旧积分上界及x∈[L/2,L]给

B0≥n[L−logL−4/L]>8n，

在L≥64成立（基点用log64<9/2，随后L−logL−4/L递增）。δ=2^(−31)nd<n，根号B0−δ严格正。

又x≥L/2>e，使B0≤1+n(L−logx)≤nL；A0≤e^x=ne^(−s)。因此选择

τ=√[A0(B0−δ)]≤n√L e^(−s/2)

是合法的，而且

τ/M0≤100s^(7/2)e^(−s/4)
 ≤100·4^8·8!/s^(9/2)
 ≤100·4^8·8!/s^4<2^38/s^4≤2^(−106)<cE/4。                     (10)

Taylor8给e^(s/4)≥(s/4)^8/8!；中间粗化只用s≥1。所有尾界都是递减powers，常数100·4^8·8!<2^38由整数exact核验。

由a*≤nd，energy减去τa*后至少(cE/2)ndM0。故原抽象energy引理给

Δ≥[(cE/2)ndM0]/(HM0)=2^(−31)nd=δ。

对每个a≥8n^4，均有a>2t、c>1及

αβ−κ−2/(a+1)=ζ−(β+2)/(a+1)>0。

因为β<2，右侧至少ζ−4/(8n^4)>0。实际target半径<α+r<3，实际source半径≥ρ>4，故同一真实target的一次eligibility支付全部截断source；中心误差guard与冻结full proof相同。

两侧λ_a=α^(−n)λ准确成立，ambient energy乘α^(−2n)，degree和target质量乘α^(−n)。因此finite-a总缺额不小于上述δ，完全不需要TV或选择历史收敛。所有source/target增量仍来自同一实际历史。

有限输入的witness仍是原共同半径aQ_i，任意f∈L1、λ>0和所有正半径均允许。先对任意紧水平集取真实见证的有限覆盖，应用已量化的同一算法；固定n,t,τ后只对已证明的数值界令a→∞，再用内正则性，得到C_n≤F(A0,B0−δ)。不存在新的输入类别限制。

## 7. 新scalar余项及最终显式负修正

−log(1−d²)≤2d²（d²≤1/2）给

B0≤1+nlogt≤T+2nd²+1，T=nL−nlogL。

结合A0≤ne^(−s)、B0≤nL，准确上包为

F(A0,B0−δ)
 ≤T−δ+2nd²+1+ne^(−s)+2ns^(3/2)e^(−s/2)。

除以nd=n/s，使用n≥s^6/2及两个Taylor8，正余项比例至多

R(s)=2/s+2/s^5+8!/s^7+2^9·8!/s^5。                            (11)

其中最后一项从更强的2^9·8!/s^(11/2)粗化而来。每项在s>0单调递减，且s=S的纯Fraction核验给

R(S)<2^(−32)。                                                 (12)

于是F(A0,B0−δ)<T−2^(−32)nd，闭合(1)。这里没有复用旧x=L−3logL的W_n公式；当前scalar scale是x=L−s。

## 8. 计算证据、全部尾部与状态

`verify_threshold.py`只运行最小Fraction/integer检查：上述几何基点s=8或16、S=2^36处(3),(8),(10),(12)、truncation常数及冻结源hash。没有重跑生产参数、Arb搜索或MATLAB。

所有指数尾已被全域Taylor不等式化为显式powers；均值/row误差的相关rational函数在给定s域单调，或已使用统一基点的上包。不存在未覆盖的s>S无穷尾。Moment沿用本轮另一个已显式覆盖L≥64的冻结审计，未引入新的eventual假设。

本报告量化的是已经独审的operator证明；没有改变其常数，也没有声称s≥2^36最优。冻结full proof及moment的hash保持不变。结果见`verification.json`。

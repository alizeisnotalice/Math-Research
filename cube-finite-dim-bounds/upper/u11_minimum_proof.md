# 11K：四维purebody全cube Green下界准确证书

2026-10-25。固定V={(1,±1,±1,±1)}、L=−Σ|v·D|、真实G=∫p_tdt。结论：

**G(x)≥21/25000 对整个[-1,1]^4成立，G在active方向线上按∞解释。**

此处只证明g下界；没有I上界或新的弱型常数主张。独立审计状态另见report。准确算术用Python Fraction、python-flint fmpq_mpoly及Arb128；共享U11谱/真实Green表示按来源hash锁定，不把共享算术当独立审核。

## 1. 从完整1152锥公式到40个真vertex residues

U11已给384 chamber×3simplicial cones和q=ψ(r)，绝对Fourier逆正常数(2pi)^-4。generic x的Green公式是

G(x)=−1/(32pi³)Σ_cones α_R DD3(sign;b1,b2,b3,b4)，α_R=|detR|/Πq_j，b_j=x·r_j/q_j。

对每一ray r定义其simple-pole star

S_r(x)=Σ_(R含r) α_R / Π_(s∈R,s≠r)(x·r/q_r−x·s/q_s)。

合并有限和给G=−1/(32pi³)Σ_r S_r sign(x·r)。每个star是准确rational function，不从若干point值推断恒等式。

六ray orbits为A=(1,0,0,0)、B=(1,1,0,0)、C=(1,1,1,0)、D=(1,1,1,1)、E=(2,1,1,0)、F=(3,1,1,1)。checker准确证A,C,E,F stars恒0：各common-denominator numerator分别由48/12/12/12项组成，完全为零多项式（degree13/4/7/5）；signed permutation等变性覆盖其全部orbits。只剩B orbit24和D orbit16。此处不依赖额外Brion/顶点消除定理。

B,D各48项star用common degree26的线性因子分母合并，再exact gcd/商与原多项式cross-check。明确地，置u=a+b、v=ab、s=c²+d²、t=c²d²，

P_B=v u²(u²+v)−v(5u²+4v)s+4v s²+(13u²−32v−4s)t，

S_B(a,b,c,d)=u³P_B/[2(a²−c²)(a²−d²)(b²−c²)(b²−d²)(u²−4(c+d)²)(u²−4(c−d)²)]。

S_D=−(9/2)(a+b+c+d)³ divided by

(a−2b−2c+d)(a−2b+c−2d)(a+b−2c−2d)(2a−b−c+2d)(2a−b+2c−d)(2a+2b−c−d)。

这些individual star可有artificial projection pole，不能单独在pole点代入。

## 2. 20 line reps与三个完整物理区域

按cube群与even性令a≥b≥c≥d≥0。B的12条line reps为e_i±e_j (i<j)，其projection在严格排序内部全正；D的8条line reps为(1,±1,±1,±1)。对称给S_(-r)=−S_r，故两方向的S_r sign合并为2倍。利用signed-permutation协变，得到

C(a,b,c,d):=pi³G=−(Σ_12 B-line S_B(T^Tx)+Σ_8 D-line sign(x·r)S_D(diag(r)x))/16。

物理排序chamber由q=a−b−c分为三case：I q>d（8个D投影正），II q<−d（6正2负），III −d<q<d（7正1负）。边界随后按真G连续处理，不能自行给zero sign。20项按这三个已证sign patterns准确合并，common degree60，再gcd/cross-multiply化为各case的N17/D20。精确N,D字典在u11_minimum_global_results.json，生成源u11_minimum_global_revision.py；不是rational插值或采样拟合。

全物理cone覆盖是U11同一个非负参数表（现在用于x）：

|case|x=R u的rays，u_j≥0|覆盖条件|
|---|---|---|
|I|A,B,E,F|a≥b+c+d|
|II|E,B,C,D|b+c≥a+d|
|III|E,B,D,F|其他，即|a−b−c|≤d|

增量y=(a−b,b−c,c−d,d)下三参数式见U11；它们唯一、内部不交、并为全部排序chamber。因此没有未处理的全cube分层。分母D20的20个linear factors经相应R u代换，定向后均具有非负系数、至少一个正系数，整体常数正；故D>0于u_j>0。checker逐因子准确核，不靠一个sample判全域符号。

## 3. 固定目标的全局非负多项式证书

令g0=21/25000，p_up=22/7，k=g0 p_up³。Arb128严格核pi<p_up。对每一case准确构造

P(u)=[a³N17−kD20](R u)。

三个P都是齐次degree20，多项式展开的非零项数分别1739、1702、1702；所有系数为非负有理数，没有负系数，遗漏项为0。完整稀疏系数与指数在u11_minimum_certificate_results.json。分母的因子符号证书、四null star及source hashes也在该文件。checker没有调用任何I或cost数字。

因此对generic严格物理cone内部、a>0，

a³ pi³G=a³N/D≥k=g0 p_up³，

G≥g0 (p_up/pi)³/a³≥g0 当a≤1。

此系数证书是对全部实非负u的多项式证明，不是有限point证据。

## 4. D=0、zero/repeated节点与方向线的延拓

rational恒等式仅先在避开有限projection超平面的generic cone内部使用。该集合在每个case内稠密。026/U09给G在八条body方向线之外有限连续、在这些线上为∞；原点也∞。所以在case边界、排序边界、D=0的可消点以及任意zero/repeated节点点，用同case generic序列延拓下界；在线上不等式直接成立。没有把D=0设成0、没有把所有projection-zero平面错当奇性、没有遗漏零Lebesgue层。

对任意非零cube点，signed permutation后max坐标a∈(0,1]，属于上述cone闭包；degree−3齐次及第3节已给下界。x=0按∞处理。整个cube由此闭合。准确单点17/(640pi³)既未作为全域最小假设，也未进入多项式目标选择。

## 5. 复现与边界

最终固定checker u11_minimum_certificate_check.py SHA35489f9fa4c52a09154105ad949a8d3410f61bf226e8a24dd8d57359f0a3e060，运行会重建同一两star/三个global式并核四null identities、所有D因子和P系数。globaldata SHA43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f，certificate SHA23a8ce3a07c6c5ad94b29dd0b57aba4dd4e3f1fd4a9098374ec700b73b0feb3d；producer已去纯timing字段，数学主体现在锁字节不再改。

初始351 affine sign cells的失败（248单cone pole、103 finite但lower不足）仅定位了相消问题；唯一表示修订是上述全和rational化，没有机械扩网格。最终证明不依赖初始分区或其removed cap。初始的one-large-jump cap lower只作诊断账，不作为本证书缺省遗漏。

本结论限固定四维purebody真实G及Lebesgue cube。一般L1/全部共同r的maximal比较仍需合法核与质量/I认证，经026/U09传递；本g证明不独自给U4<16，不调整账本，也不推广tau或其他背景。

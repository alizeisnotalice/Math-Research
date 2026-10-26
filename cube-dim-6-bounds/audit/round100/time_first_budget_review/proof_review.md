# A100 独立证明审查

结论：接受100的新middle应用、joint u/rho局部预算、同域表达式比较、
条件构造预算的非竞争性及真实删列矩系数比较。接受范围是以下明确的解析接口，
不是global g、竞争性完整I、正S或Cn新界。原生成器的运行政策资格另见report。

## 1. 全时间接口及去重

正三resolvent积分J在062已证明和使用，含正参数重根极限。100的自含partial-fraction/
dominated-convergence复述正确；新贡献是把它用于043/074的middle补事件，
而非重新发现J或重新执行062公式基线。本独审没有数值重算该旧基线。

对actual split P+B=p，P是两个指定参数同时大的事件，B是其真实互斥补集。
估计B时可把w-small分支扩大，但不将扩大后的事件误认为互斥。
其他30个参数及所有其他大参数、抵消都仍在源积分内，没有丢掉source events。
043/074的31/30/29列正删除系数992/960/928保留。

补事件的第三距离沿保留参数fiber可以为0，不能先用虚假uniform shifted gap。
这里先保留真实Cauchy积分，使用原x到其他pair plane的g，再用已付的平均式；
保留参数在三维quotient中仍删除，所以原峰必须乘(10/9)^3=1000/729。
100的Krem仍为此损失乘完整Crem及1+sqrt(6/928)，没有借062的恢复顺序静默消掉该损失。
全时间积分后根归一为h/sqrt(A)，除以928*960*992，常数和π/2因子正确。

本100所用098解析接口也独立核对：q4的正删除上包在t=1给C/(|x|²+960)，
因而沿任何单位方向的线trace具有统一可积支配；连续marginal相等从a.e.延到每点。
q_t(re)=t^-4 q_1((r/t)e)，正变换s=r/t得到
`int t²q_t(re)dt=Q3(e)/(2r)`。全32个符号的二/四矩给ψ>=32|ξ|/sqrt3，
三维Fourier径向积分给Q3<=1/(π²c³)。实际两参数coarea J2=6及其事件caps再给
`int P_tdt<=sqrt3/(128π^4 r)`。这接受的是100使用的098解析normal-potential/K接口，
不替代主对整轮098历史运行及其他结果的独立验收。它是积分估计，不是各向同性的逐时间q_t估计。

## 2. 真正的联合gap及所有区域系数

R固定u∈[3/8,5/8]、0<rho<=1/64。真实face坐标x_A=(1-e1)1_A+e、x_B=u1_B+z，
sum e=sum z=0，rho²=|e|²+|z|²。由e1²<=2|e|²/3，
`|x-x(u)|²=rho²+3e1²<=3rho²`，每个distance可从中心下包减sqrt3*rho。

v-small原center residual的p>=1/2，nearest-body-line距离下界为
sqrt(3/2)min(u,1-u)。w-small中p=1-a、q=u+a为正，
`d[3(p²+q²)]/da=6(2a+u-1)<0`直到a=1/16，
所以origin gap至少sqrt(3[(15/16)²+(u+1/16)²])；
不能用较弱旧origin估计后硬塞cw=7/128。line gap至少sqrt(3/2)(7/8-u)。
其他pair-plane中心gap是min((1-u)/sqrt2,sqrt2u)，u>1/3选择第一项。
这些覆盖补事件实际用到的全部line/pair。

五个物理center gap在该u区间均>sqrt3/7，因此减sqrt3*rho后都>=
`(1-7rho)*h_center`。1-7/64>0。独立准确平方比较核了
1/sqrt2>7/10、sqrt(3/2)>6/5、sqrt928和sqrt960<31、sqrt992<32、
sqrt3/2>6/7、sqrt822/16>7/4。因此三个normalized roots的全部有理下包
`a=7(1-u)/310`、`bv=6min(u,1-u)/155`、`bw=6(7/8-u)/155`、
`cv=3/112`、`cw=7/128`方向正确，都是lower minorants。

J分别随三个正根递减、齐次degree−5，故
`B_G<=C(1-7rho)^-5 Z(u)`，C=πKrem/[2*928*960*992]。
这在取q次幂之前保留两补事件的合并Z，混合项没有遗漏。

A-face线性映射的Jacobian为sqrt3，B的u/zero-sum映射也是sqrt3，合计3。
e1为A中最大坐标的sector恰1/3径向测度；积分的是径向majorant，不假设补事件密度具有分支反射对称。
其他face clipping只缩小非负积分。因此Jac3与sector1/3抵消，四维径向面元为4v4*rho³drho。
u长度1/4，给体积v4*delta⁴/4和singular积分(5v4/14)*delta^(14/5)。

十个k3区域是A含1且|A|=3的10个选择。对R_A点，其到任何其他pair plane的距离
>3/(8sqrt2)-sqrt3/64>1/64。它们确实互不相交，不是把全部31个tube假作disjoint。
只用固定坐标1的置换对称复制相同预算。

## 3. 新有理/对数证书与径向积分

T=(a+b+c)/[abc(a+b)(a+c)(b+c)]等于三个正reciprocal-product的和；
每项f=prod l_i(u)^-1在各半区间所有l_i>0，
`f''=f[(sum l_i'/l_i)^2+sum(l_i'/l_i)^2]>=0`。
所以两piece的最大值只需三个解析端点3/8、1/2、5/8；这不是quadrature grid。
独立最小代码直接核三个完整有理值及M，不套用保存status。

对三个新有理integrands，独审通过每个有理pole的Laurent residue derivative生成局部系数，
逐式准确恢复原函数；没有调用作者sp.apart或作者程序。
逐项积分后所得有理部分及所有正有理log argument的**带符号**系数恰等于保存对象。
再在512-bit Arb重包一次该新signed-log和，结果完全落入原384-bit保存包围。
没有把负系数错误乘上同方向端点；所有log项相加的区间算术由Arb完成，
此外对关键S减法独立使用有理端点核上/下方向。

Z^q<=M^(1/5)Z因为q=6/5且0<=Z<=M。
rho因子的q次幂是(1-7rho)^-6，非−5或其他指数。
独立核F'(w)=(1-w)^3/w^6及w=1-7rho换元，
`Jr=[F(1)-F(57/64)]/7⁴=313/12033841140>0`。
因此`B_R=C[4v4 Jr M^(1/5) Ju]^(5/6)`。
先完成每个真实source event的全time积分，再做spatial Minkowski加P_R；
不使用逐时间F/L relaxation或未付的tangent/normal独立性。
rho=0为空间零测集，可让点wise majorant为infinity，q<4保证周围积分可积。

## 4. 同R比较及条件预算失败的准确含义

旧098的512/B0只是把已付uniform式转移到**同一个本轮新R**的预算，
不冒充098历史上已运行此R。新N_R约2.7980735e−5，旧B0同域式约3.9011791e−4，
新/旧约0.07172379；比的是同一个完整G范数上包表达式，不是实际G值或实际局部范数。
`B_R/N_R≈0.999514`是补事件在此范数预算的占比，不是实际源质量或真实积分的占比。

I_ten,U=10N_R^(6/5)约3.4374069e−5，g0=1/65000仍全域未证。
条件whole-face gate B=4/48750，I_target=64(g0/6)^(6/5)约1.2499782e−5。
独立核gate⁶=64⁵(g0/6)⁶、target端点及ten倍关系；
`S_plan`上端使用`target_upper-Iten_lower`仍严格负。
这只说明将此**构造上包**分配给十域、把余域乐观置0的方案失败。
不能从upper allocation超过门槛推出实际I超过门槛，不能下界真实模型cost或Cn。

仅在P_R和其他结构均固定、只把补系数统一乘γ时，十域要达到正预算需
`γ<[B_gate/10^(5/6)-P_R]/B_R≈0.430143565`。
这是这套固定预算下必要的系数目标；余域正预算会更严，global g变化也会改变门槛。
它不是对任何机制、真实补事件或真实Crem最优值的全局必要下界。

## 5. 删除一列的真实coercivity及原峰匹配

对实际遗漏b的31个forms，二矩恰为32n²-s²；四矩恰为
96n⁴-64sum ξ_i⁴-s⁴<=96n⁴-s⁴。
有限和Holder给`ψ_omit>=n sqrt[(32-x)^3/(96-x²)]`，x=s²/n²∈[0,6]。
使用四矩upper作为分母给lower coercivity，方向正确。
分母>=60，导数号由x²+64x−288给出，唯一minimum x*=−32+4sqrt82∈(0,6)。
独立整数sqrt82有理括号核出f(x*)包围完全位于保存c31²内；c31约16.55839748。
全R6 coercivity限制到任何三维quotient仍成立，故真实遗漏law峰<=1/(π²c31³)。

独审只读033五个原R_exact及normal rows，逐一核r=dimension=3、Gram determinant
(64,96,128,160,192)，并以自己的3×3整型det核coarea。
保存N=(435,1740,435,580,1305)与full类计数的29/32关系及ordered26970一致；
`24576 sum N_j R_j=44743079/2970`准确成立。
旧系数用的是原033 peak/coarea：q3/J3=R/(2π)³，再带遗漏损失1000/729。
新actual31-law峰已经删除b，不再额外乘一次该损失；coarea需保留1/sqrt(D)。
新系数`196608/(π⁵c31³) sum N/sqrt(D)`正确，其中196608=32³×6。

独立512-bit定向后处理的合法新/旧系数比落在原保存2.702150085包围内，严格>1。
因此这个uniform moment replacement比现有exact peaks更差，不能宣称改善。
共同Cauchy-average因子在比较中消掉；既未混入其他frame的峰，也未重跑033的fan/trace枚举。

## 6. 审核边界

独审运行两个新saved/exact组；未执行作者4配置、062旧J基线、043/074 F/L、
任意grid、native/G/H/p kernel或MATLAB。新的数值重包仅针对本轮signed-log primitive和
本轮actual-law/coarea相对系数；其他数值是保存dyadic关系/严格端点核查及精确源码检查，
不是独立重放完整四配置所有forward数字。解析证明和有限检查共同支持上述限定预算。
global g、完整竞争性I及正S均未闭合，实际上下界保持原记录。

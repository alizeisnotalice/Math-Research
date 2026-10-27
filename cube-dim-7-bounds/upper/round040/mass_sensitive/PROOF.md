# R040 · f-dependent预算与共同质量：精确flow障碍及cheap-hotspot planting

## 结果与范围

本报告给出一个具体coupling的严格反例，而不是再证明“所有margins可以是1+ε”。结论是：即使预算在整个允许区间|S_i|≤h_i≤M_i/λ中依据同一个真实f最优选择，R032的联合cover-flow证书若最终仍取combined load的essential supremum，其leading nlogn系数仍被R037真实历史饱和。

同一个反例历史具有共同star交点，故R039共同质量LP及top witness支付的真实selected-union比率却近1。这个差距明确定位了损失：精确f等式之后再次取pointwise maximum，而不是margins的独立不一致。

另给endpoint h_i=M_i/λ的cheap-hotspot planting引理。它说明“存在一个接近nlogn的load点”可用任意小f质量植入任意背景实例，因此仅以这一sup条件作全实例dichotomy的前提，并不能提供独立于目标uniform估计的捷径。

**没有新C_n上界。** 没有七维、低阶壳常数优化、生产计算、MATLAB或Git/网络操作。原算子始终是共同半径中心立方体，全部正半径与一般L1桥接见第5节。

依赖：主R032 `joint_cover_theory/proof.md` 的有限预算/CF定义；geometry R037 `general_cover/PROOF.md` 的真实全selected高load历史；geometry R039 `joint_mass_dual/PROOF.md` 的共同star支付。下文重新给出用到的接口和方向。

## 1. 固定原有限实例的精确f模型

固定a>1、b=a(1+1/t)、τ>0。输入有限小立方体Q_i，S_i=aQ_i，V_i=|S_i|。同一个非负f∈L1在水平λ>0严格见证全部S_i。记

M_i=∫_{S_i}f，V_i≤h_i≤M_i/λ。

预算在选择前固定；选中后保持真实增量I_i、v_i=|I_i|及p_i=v_i/h_i不变，仍按最大合格半径选取。资格为先前p_i在bQ_i中心上的和≤τ。

记输入并集U_all、selected并集U_sel、D=U_all\U_sel。R032的CF变量将D完整分配到selected项，每项接收c_i≥0且Σc_i=|D|。其combined load为

G_h(y)=Σ_i(v_i+c_i)/h_i ·1_{S_i}(y)，
L_CF(h)=min_合法flow ||G_h||∞。

每个合法flow都由同一个f给

∫fG_h=Σ_i(v_i+c_i)M_i/h_i≥λ|U_all|。                            (2)

在endpoint h_i=M_i/λ，等式恰好成立，不论如何分配flow：

∫fG_h=λ|U_all|。                                               (3)

因此endpoint保留共同质量，但任何只“优化f加权平均”的flow目标都是固定的未知比率，不能凭重新分流降低它。最后取essential maximum才是CF证书中的额外松弛。

更明确地令F=||f||1、π=f/F。在endpoint，

p_i=(λ/F)v_i/π(S_i)，
∫Σ_i p_i1_{S_i}dπ=λ|U_sel|/F，
∫G_h dπ=λ|U_all|/F。                                          (4)

所有π(S_i)来自一个共同概率测度；不能把它们当独立自由参数。

## 2. 核心命题：整个允许预算区间仍饱和CF

### 2.1 精确有限陈述

取R037定理1给出的任一个真实有限input历史Γ：所有小Q_i两两不交，所有input皆被选中，且0在每个S_i内部。原几何weights为

w_i=v_i/V_i=a^(−n)，W_old(0)=Σ_iw_i。

给定η>0，使用唯一共同输入

f_η=(1+η)λ·1_{∪_iS_i}。

它紧支撑且属于L1，所有真实质量准确满足M_i=(1+η)λV_i。

**命题1。** 对任何允许的预算V_i≤h_i≤(1+η)V_i，原算法仍选中全部input，真实增量仍为Q_i，而且

L_CF(h)=Σ_i v_i/h_i，
inf_h L_CF(h)=W_old(0)/(1+η)。                                  (5)

同时，对同一个f_η及同一个实际input，

λ|U_all|/||f_η||1≤c_a，c_a=((a+1)/a)^n；
|U_all|/mathfrak m(Γ)≤c_a，mathfrak m(Γ)=max_i V_i。              (6)

其中mathfrak m是共同witness质量LP，所有界均为原共同半径立方体实例上的界。

**证明。** 允许区间逐项给w_i/(1+η)≤p_i≤w_i。R037构造的全部input在原weights下始终eligible；将所有项降低不会破坏任何资格。small cubes两两不交，所以actual v_i不依赖内部同半径选择顺序。于是全部input仍被选中，D=∅，不存在任何非零cover-flow变量。

CF因此恰为||Σ p_i1_{S_i}||∞。所有stars有限且共同包含0的开邻域，故其essential maximum准确为Σp_i。endpoint h_i=(1+η)V_i同时最小化每项，给(5)。这不仅测试一个预算选择，而是测试R032允许的整个f-dependent预算区间。

共同star定理直接给(6)：若R=max r_i，任一y∈Q_i满足||y||∞<(a+1)R，所以互不相交increment总量≤[2(a+1)R]^n=c_aV_max；最大star的真实witness给λV_max≤∫_{S_max}f_η≤||f_η||1。共同质量LP至少V_max，并可在共同开交集中集中质量V_max，故其最优值恰好为V_max。证毕。

### 2.2 Nearlinear regime中的严格leading障碍

固定C≥0，n/(logn)^C≤t≤n(logn)^C且τ≥1。R037实际历史给

W_old(0)≥nlogn−O_C(nloglogn)−o(1)，

其中a可取任意充分大。取η=1/logn，并进一步选a≥n(logn)^2。则

inf_h L_CF(h)≥nlogn−O_C(nloglogn)，                              (7)

而同一真实f的支付满足

λ|U_all|/||f_η||1≤c_a≤exp[1/(logn)^2]=1+O((logn)^(−2))。        (8)

因此这个具体joint CF架构即使最优选h，也不能把其uniform pointwise certificate的leading系数降到1以下。真正共同质量LP和真实f支付却没有这种饱和。

还可保证U_all真的属于原{M_nf_η>λ}：若x∈Q_i，则S_i包含于以x为中心、共同半径(a+1)r_i的立方体，其平均至少

(1+η)λ/c_a>λ

（进一步加大a使c_a<1+η即可）。所以反例不依赖“只在中心有形式witness”的放松。它是一个有效的真实level-set finite-cover输入；式(8)只给有限cover（此处也即selected union）比率λ|U_all|/||f_η||1的上包。虽然U_all包含于整个level set，这个有限子集的比率上包不能反向控制整个level set的弱比率。本构造不给C_n下界。

## 3. Cheap-hotspot planting：sup前提无法控制共同f质量

本节仅针对endpoint h_i=M_i/λ，或其他逐项由局部M_i,V_i决定且在共同homothety下协变的预算规则。它不声称任意依赖全输入的预算oracle会保持背景选择不变。

### 3.1 原有限实例的精确直和

取任意有限背景实例B，同一a,b,τ,λ，所有witness严格且f_B支持在一个有界集合。另取命题1的probe H及f_H。对ε>0，把H的小立方体及其stars同时缩放ε，f_H变成f_H((x−T)/ε)，再平移到足够远T，使两个分量的全部stars、全部(b+1)Q及小Q彼此分离，且两个f支撑也分离。

所有witness质量、预算、increment体积在H内同时乘ε^n，因此p_i不变；B内真实M_i完全不变。跨分量任何中心都不落入另一分量资格区，所有资格和increment都精确解耦。全局最大半径规则只是两个原选择历史的半径归并，不改变任何内部资格决定。

于是联合input是一个真实有限实例，且

F_ε=F_B+ε^nF_H，
|U_all,ε|=|U_all,B|+ε^n|U_all,H|，
maxload_ε≥maxload_H。

准确比率为

R_ε=[F_BR_B+ε^nF_HR_H]/[F_B+ε^nF_H]→R_B。                      (9)

同样等式适用于selected-union支付。probe的真实f概率质量仅ε^nF_H/F_ε→0，而其高load完全不衰减。

如果还要求CF几何cover区与witness incidence完全解耦，前述分离已经足够：CF流多面体是两个分量直积，目标的essential maximum给

L_CF,total=max(L_CF,B,L_CF,H)。                                 (10)

所以一个任意便宜的probe也能保持CF的高pointwise瓶颈。这不是共享mass被重复计数；两个f支持完全分离，缺陷是最终sup本身不看到分量质量。

### 3.2 对局部sup dichotomy的精确限制

假定在固定参数下有一个“只要endpoint的sup load≥H0，就有全input ratio≤K”的命题，且H0可由probe达到。对任意背景B植入ε-probe，再令ε→0，(9)立刻给R_B≤K。

因此这种conditional命题对有限实例实际上已经蕴含不带高load前提的uniform比率界。它不能仅靠一个局部热点的几何证明来替代全局共同mass问题。这个结论并不否认带有额外f-weighted high-load质量假设的dichotomy；相反，它准确指出必须增加什么信息。

## 4. 一个有效但尚未闭合的mass-sensitive方向

由(4)，endpoint任何flow都有E_πG_h=R_all。若有已知G_h≤M并且R_all≥(1−ε)M，则对任意0<θ<1，

π{G_h<(1−θ)M}≤ε/θ。                                          (11)

证明只需E_π(M−G_h)≤εM。这是一个真正受共同f约束的必要条件；它不能由sup大而代替。

下一步可研究如下量化问题：当一个合法原历史在π质量的大部分位置都有接近nlogn的combined load时，是否能用各位置共同star cluster的几何，把同一π质量重复使用的representative支付统一压低？R039的cluster LP保留了这种信息，但本报告没有证明其uniform改善。

(11)只是必要条件，(3)也只是精确身份；它们不是新C_n上界。有限生产模型若由主协调请求，应保留同一incidence-cell概率π、原分量资格与actual increments，而不能把π(S_i)自由优化成不一致margins。

### 4.1 一个真正mass-sensitive的有限目标

对同一个实际selected历史Γ，使用R039的精确有限cluster列J_{C,j}：pivot arrangement cell C中的所有stars中，保留半径不大于代表j的项。每列有共同star交点，代表j为该列最大半径者。定义

S_f(Γ)=min Σ_{C,j}θ_{C,j}M_j，
θ≥0，Σ_{(C,j):i∈J_{C,j}}θ_{C,j}≥1 对每个selected i。              (12)

这里M_j不是自由价格，仍准确为原同一个f的star积分。每个primal可行θ都给

λU_sel≤c_a Σθ_{C,j}M_j，故λU_sel/F≤c_a S_f/F。

再接原有限cover界，可得

λ|U_all|/F≤(1+A/τ)c_a S_f/F，A=((b+1)/a)^n。

这是R039几何cluster定理的f加权有限版本；它不用对representative load先取pointwise maximum，也没有免费假设不同representative可只付一次。ΣθM_j本身已完整计入重复f质量。

该目标对命题1的共同star probe使用一列，立即有S_f≤M_max≤F；并且对第3节远离的两个分量，任何cluster不能跨分量，所以S_f精确加性：S_f,total=S_f,B+ε^nS_f,H。它因此不会被便宜热点的pointwise maximum单独主导。

其普通LP对偶是max Σ_i s_i，s_i≥0，且Σ_{i∈J_{C,j}}s_i≤M_j。一个dual可行点只给S_f下界，不能用于新的C_n上界；正向改善需要primal覆盖系数或uniform控制S_f/F。本报告未证明该uniform控制，也未运行生产LP。这个接口可供主协调请求后由Terra测试几何/历史一致的有限实例。

## 5. 一般L1及全部尺度mapping

任意一般L1原输入的有限compact-level-set子覆盖，只涉及有限个有界stars。删除f在这些stars并集之外的部分不改变任何witness质量，而只降低||f||1。若h_i已固定，或预算规则仅由每项局部M_i,V_i及固定参数决定并在共同homothety下协变，则这也保持预算与选择；本段的有界支撑归约及对原f支付的传回仅在这个规则范围内使用。任意依赖全输入或全局f信息的budget oracle可能因截断而改变h_i，因此不声称其预算或选择也保持。

反过来，命题1及planting使用的分段常数f都是非负L1真实输入，每个star有严格平均。所有变换是共同平移与共同homothety，原小立方体和stars仍是共同半径的立方体，没有矩形算子、尺度截断或独立分配的mass。τ与a,b保持不变，homothety使p保持不变。

这些反例足以测试一个声称对全部原有限实例成立的CF/conditional证书。它们不反推任何C_n下界。要从一个未来的正向uniform finite-family质量界推到C_n，仍必须按主R032的compact-level-set有限覆盖与Lebesgue内正则性步骤；本报告没有完成这个正向uniform步骤。

## 6. 最小exact核验

`verify.py`仅检查一个1维有理嵌套star/分离small-cube实例：所有witness来自一个真实f、预算区间端点、D为空时CF精确值、真实f平均身份、共同质量与planting的有理混合式。它不重算R037高维构造、不进行生产参数搜索。

任意维leading障碍来自R037真实全selected历史加本报告第2节的逐项解析证明，而非该小实例。结果见`verification.json`。所有原source及R032–R039保持不变。

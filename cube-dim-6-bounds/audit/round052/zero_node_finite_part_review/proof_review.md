# F052 零节点有限部分：独立数学审查 PASS

B052 最终源锁与25个输入哈希核验通过；以下给独立证明及实际精确参考结论。只讨论dimension6、正同时间Cauchy frame的scalar Green，没有gradient或Hessian推广。

## 1. 一般frame假设与真S2

令a1,…,aN为有限个非零实向量，span=R6。每列乘独立标准Cauchy(t)变量；正权重lambda_j可吸收到a_j=lambda_j v_j，不能把不同参数时间混在推导中。若有零权重应先删除，并保持剩余列张成R6。真实symbol为psi=Σ|a_j·xi|，正齐次且coercive：psi≥sqrt(xi^T AA^T xi)≥c|xi|。

S2定义为全部列子集rank≤2之span的并；因为frame张成R6，它等价于所有独立两列plane的有限并，也包含每条单列line和原点。不得用旧frame的S2，也不得假定任意三个不同方向独立。

正概率推前给p_t，exp(−t psi)可积，故其标准密度版本连续且非负。Fourier界给p_t(x)≤C t^−6，适用于所有x和t>0。

## 2. Rank-greedy 给 offS2 的一致小时间 O(1)

取紧集K⊂R6\S2，令delta=dist(K,S2)>0，M=Σ||a_j||，h=delta/(2M)>0。对参数向量z满足x=Az∈K，称|z_j|≥h为大参数。若其列的span rank≤2，则大参数部分和位于S2，剩下小参数部分的范数≤hM=delta/2，矛盾。因此每个这类参数向量至少含一个三列独立且都大的triple。

于是事件{Az∈K}被有限个真实事件E_I={|z_i|≥h for i∈I}覆盖，I遍历独立triples。每个I固定扩成六列basis B_I，det不为0。对其余N−6参数先积分；basis变量作真正六维线性换元，Jacobian为1/|det B_I|。三个受限密度各≤t/(πh²)，另三个basis密度各≤1/(πt)，所以该事件的推前密度≤1/(π^6 h^6 |det B_I|)，余参数概率积分为1。对有限triples求和得p_t≤C_K，独立于t。

证明应先对支持于K附近开集的非负test functions积分，不在零概率fiber上条件化。取略大的offS2开邻域仍有同样正delta下界，所得a.e.密度不等式由固定t的p_t连续性推广到每点。这里的large-set秩论证只选择存在的独立triple，不把dependent triples强行除以零det。

结合大时间Ct^−6，Green G(x)=∫_0^∞p_t(x)dt在offS2有限，并由局部共同可积支配和固定t连续性在offS2连续。该结论不提供任何空间高阶导数界。

## 3. 真S2上G无穷，必须拒绝

给x∈S2，选择独立pair a1,a2使x在其span内，再扩成六列basis B。其余列随机和Y_t=tY_1。选固定L使P(||Y_1||≤L)>0；N=6时Y=0，概率为1。

在此正概率事件上，B^−1(x−Y_t)的前两个坐标有界，后四个为O(t)。对0<t≤1，前两Cauchy密度各≥c t，后四各≥c/t；真实线性换元与余参数卷积给p_t(x)≥c_x t^−2。这个积分身份本身由basis密度与有限概率卷积给出，对每点有效。故G(x)=∞，包括line和交叠点。

即使某个有限部分值是负数、正数或零，它都不是这里的G。生产入口必须对真实新frame执行S2 gate；不能把computed finite part当作有限Green，也不能因节点零而不检查真实几何。

## 4. 完整finitefan的正epsilon恒等式

额外假设已证明一个完整有限simplicial fan：full-dimensional interiors只覆盖一次R6（或群multiplicity严格正确）、每锥R可逆、psi在该锥上线性且q_j=psi(r_j)>0、正确Jacobian alpha=|detR|/∏q_j。这些是独立几何前提，旧21cone或新全fan数据不能仅凭格式替代。

令b_j=x·r_j/q_j。对epsilon>0，Fourier锥积分绝对合法。由DD5[(t−ib)^−1]=i^5∏(t−ib_j)^−1，实部时间积分得到

    H_epsilon(x)=∫_epsilon^∞p_t(x)dt
      =(2π)^−6 Σ_fullfan alpha·DD5[atan(b/epsilon);b0,…,b5]。

有限epsilon时atan对所有实b解析，重复节点按normalized Hermite derivatives，零节点不得drop。完整fan求和后才讨论epsilon↓0。没有对单cone先积分到0，也没有把sign在0当作可微函数。

## 5. Laurent结构与负幂准确抵消

固定x和有限fan，非零节点b处

    atan(b/epsilon)=(π/2)sgn(b)−atan(epsilon/b)。

后项在epsilon=0解析且只有正奇幂；对b的所有正阶normalized jets也为O(epsilon)，常数项为0。节点b=0的value为0，normalized jets为

    k奇：(-1)^((k−1)/2)/(k epsilon^k)，k偶：0。

DD5只用最高五阶jets，所以每个cone的展开只有有限负幂epsilon^−1、epsilon^−3、epsilon^−5，加常数及O(epsilon)。重复非零节点完全由其真实Hermite coefficients处理，不能按distinct公式硬除零。

完整有限fan求和依然是这一Laurent形式。off真S2的G有限，且H_epsilon↑G；因此负幂多项式的三个系数必须依次准确为0：先最高非零负幂会导致H无有限极限，逐次排除。这里的取消是完整sum的数学事实，不是单cone正则化，也不是浮点近零测试。

定义纯linear Hermite functional：非零node的0阶数据赋sgn(b)，其全部高阶数据赋0；zero node的全部数据（含0阶）赋0。设Gamma为用正确完整fan alpha加权的这些DD5 functional之和。由于atan常数部分有π/2，故

    G(x)=Gamma/(128π^5)，x∉真实S2。

分母来自(2π)^6与π/2，不能误成64π^5或再次乘global负号2。Gamma是常数项functional，不是任何含零sign函数的经典导数值。结论仅scalar；x变化时zero/重复结构变化，不能直接微分此赋值规则以获得gradient/Hessian。

## 6. Tiny axis6 完整analytic reference及实际结果

另取6坐标轴正frame。它有64个orthants，每个q_j=1、|detR|=1，不使用任何body32/152frame证书。对x=(1,2,3,0,0,0)，其真S2为支撑≤2的向量并，此x在外。直接

    p_t=1/[π^6(t²+1)(t²+4)(t²+9)]，
    G=1/(120π^5)。

实际完整64term finite-part sum给Gamma=16/15，三个负幂系数全部0。这是一项真实tiny完整reference，不是0所有kernel工作。

拒绝gate对照x=(1,1,0,0,0,0)在真S2。此时

    p_t=1/[π^6 t²(t²+1)²]，
    H_epsilon=π^−6[epsilon^−1−3π/4+(3/2)atan(epsilon)+epsilon/(2(1+epsilon²))]。

G=∞，而H的常数有限部分为−3/(4π^5)，相应Gamma=−96。完整负幂epsilon^−1不抵消。若算法把该负数报作G即错误；本参考只能标S2_REJECT、finite-part diagnostic。

实际独立4配置一次通过（PID53463，0.230572250秒），完整结果见results.json。两项tiny frame各保留64个真实orthant项，合计128项；重复节点与零节点从未跳过。两个点的全部15个axis-pair距离分别最小为1和0，前者offS2，后者硬拒绝。没有production body32/152frame fullkernel或native调用。

## 7. 最终源文与代码审查

B052 theorem最终SHA256为f4ead1ef2785595fa71603f6a2ec651f5a31d16dd19296b963bfc4cd584233b3。其正交商证明和上文六列basis证明一致：对独立triple J，商symbol就是psi限制于U_J正交补，仍≥c|eta|，三维峰≤1/(pi²c³t³)；Jacobian为sqrt(det(A_J^T A_J))。三个受限参数尾带t³/pi³R6，相乘得1/(pi5c3R6sqrtGram)。先对正测度推前、再在开邻域用固定t连续性升为每点，没有零概率fiber条件化。

源文的大时间系数15/(8pi³c6)正确，S5面积pi³及Gamma(6)=120配(2pi)^−6；用C=||B^−T||F²、c=(1+C)^−1确给coercivity。rational列公共分母d_A令正Gram determinant≥d_A^−6，粗overestimate和全部pair残差规范均合法。parallel/repeated列无需删除或假设通用位置，只需非零且rank6；正权重已吸入真实列。

源check_scalars.py只实现正式scalar检查，无production oracle。已只读审其Hermite normalized jets、精确重复节点、resolvent相位和positive cutoff归一；没有执行其12case writer/checker。独立verify_reference.py使用Fraction负二项级数卷积形成局部pole系数，不导入作者代码，不做作者矩阵逆/导数残数/三角DD运行。另由直接Taylor remainder核新tuple[0,0,0,0,0,3]，得到(Gamma,L1,L3,L5)=(1/243,−1/81,1/27,0)，且时间导数严格为−1/[epsilon^4(epsilon²+9)]。一次symbolic jet family核k=0,…,5以及归一，负节点由奇性确认。

## 8. 审计裁决与后续硬边界

PASS：上述一般正spanning六维Cauchy frame标量定理、真实S2有限性判别、完整fan负幂取消、constant functional与1/(128pi5)归一。实际tiny reference通过并明确拒绝S2反例。

这不等于已审核任何新production engine。任何未来点必须绑定真实frame、rank6、完整fan、multiplicity及全部pair严格正距离；相同数值格式不够。零点相等性/符号必须精确决定，不可把近零当零、把无法区分的非零节点合并、跳过zero/repeated项，或重乘已合并的global-negation权2。finite part不能用来跳过S2 gate。

没有gradient/Hessian扩展、全域floor、box remainder、C6新界，亦未放行053 actual。主代理保留最终证明与实现放行责任。

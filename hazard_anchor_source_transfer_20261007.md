# 连续来源的有限 hazard 锚点转移：显式对数费用

2026-10-07，根代理。使用用户 J05（真实跳率/原source合同）和 N04（先确定全费用再求和）工作流。下列证明自含，不把 skill 当定理引用。上一轮“条件弱界不能直接平均”的禁令仍有效；这里通过正核比较及一个已证明的有限对数求和代替非法平均。

## 1. 有限弱 L1 求和，含质量熵版本

令 f_j≥0 可测，λ|{f_j>λ}|≤a_j，j=1..N；A=Σa_j>0。零a_j项为零a.e.。设 p_j=a_j/A。则
\[
\|\sum_j f_j\|_{1,\infty}
\le A[3+2\log2+2\sum_jp_j\log(1/p_j)]
\le [3+2\log(2N)]A. \tag{1}
\]
证明对给定λ，将 {max_jf_j>λ} 以 A/λ 支付。在其补集，低部分 f_j≤b_j=λp_j/2 的和≤λ/2。剩余若总和>λ，必有Σ f_j1_{b_j<f_j≤λ}>λ/2。由层饼
\[
 \int f_j1_{b_j<f_j≤λ}
 \le b_j|\{f_j>b_j\}|+\int_{b_j}^{λ}|\{f_j>u\}|du
 \le a_j[1+\log(λ/b_j)].
\]
Chebyshev即得(1)。所有积分均有上下截断，不使用无限测度空间上的无截断L1。熵≤logN由Jensen或相对熵非负给出。这个对数并非免费Minkowski。

对数依赖确有必要：在[0,1)的N个等长格上，令 f_j(k)=1/[1+((k+j) mod N)]。每个弱准范数为1/N，总a=1；和恒为H_N=Σ_{r≤N}1/r，故总弱范数为H_N。它是一般测度论求和的压力样本，不是原传播或cube反例。

## 2. 正 cocycle 和来源一次预算

设 U(s,t) 是正演化族，U(s,a)=U(s,t)U(t,a)，s≥t≥a；来源 F_t≥0 为一份固定 Bochner L1 时空密度，支持 t∈[0,T0]，总质量 M=∫||F_t||1dt。可同时允许有限多个L1密度时间脉冲。预设N个互不相交时间格I_j及左端a_j，覆盖来源支持；若
\[
 U(t,a_j)\ge\theta I\quad(t\in I_j),\qquad 0<\theta≤1,
\tag{2}
\]
则对每个s及t∈I_j,t≤s，正cocycle给 U(s,t)≤θ^{-1}U(s,a_j)。令 ν_j=∫_{I_j}F_tdt（包括所分配脉冲），则Σ||ν_j||1=M，而且
\[
 \sup_s\int_{t≤s}U(s,t)F_tdt
 \le\theta^{-1}\sum_j\sup_{s≥a_j}U(s,a_j)ν_j. \tag{3}
\]
未来尚未输入的那部分ν_j只在正上界中加入，没有重启实际过程或改变来源。若各锚点初始密度最大算子有同一弱费A_*，应用(1)得
\[
 \|\sup_s\int_{t≤s}U(s,t)F_tdt\|_{1,\infty}
 \le\theta^{-1}[3+2\log(2N)]A_*M. \tag{4}
\]
可用真实质量p_j=||ν_j||1/M替换logN为其熵。对有符号F先取|F|，其总变差进入M；不继承来源的旧winner/父组语义。

初始密度最大界和可测版本是本定理明示假设，不由正性自动推出。对所用有界jump模型，正幂级数和局部共同支配给逐点右连续版本；时间脉冲取右连续包含端点。亦可先取有限/有理输出时间集，界不依赖该集合，再以共同版本取极限。

若U另有参数L，(2)对L一致，且F_t同一、不依赖L，则(3)右边可换sup_{L,s}U_L(s,a_j)ν_j，使用相应初始密度联合弱费。若F_{t,L}随L改变，不允许把它们当同一ν_j；需要另有共同正来源支配及其质量证明。

## 3. 原有序噪声：early来源、任意future终点

固定物理尺度ell。原K_{s,t}=[(1-r)I+rG_{1-t,ell}]^{tensor n}，r=(s-t)/(1-t)，构成真实正cocycle。来源仅early t≤1/2，输出s可一直到1（或先s<1再极限）。取
\[
 N=2n,\quad a_j=j/(4n),\quad I_j=[a_j,a_{j+1}),\quad j=0,..,2n-1,
\]
末格含1/2。这里a_{2n}=1/2。对t∈I_j，未跳跃权
\[
 p_0(t,a_j)=\left(\frac{1-t}{1-a_j}\right)^n
 \ge(1-1/(2n))^n\ge1/2.
\tag{5}
\]
最后用Bernoulli不等式，完全不需数值指数拟合。因此θ=1/2。锚点 c_j=1-a_j∈(1/2,1]，每个固定锚点的完整future K最大界若为A_ord(n)，则
\[
 \|\sup_{0≤s≤1}\int_0^{\min(s,1/2)}K_{s,t}F_tdt\|_{1,\infty}
 \le2[3+2\log(4n)]A_{ord}(n)M. \tag{6}
\]
这没有对c取免费sup，也没有每窗重新领W；输入分割的是一次来源F_t，不是反复存活的η_{a_j}。分格宽度1/n仅导致logn求和费，不导致nW。

**尚未证明 A_ord(n)=polylog。** 固定冻结H的Ologn不能代入A_ord；窄窗正支配只覆盖r≤1/(2ceil√n)，而(6)要求每个锚点全部future。旧D_n强包络仍给A_ord≤D_n，但不能把这个旧界报作新增一般上界。新结果是连续来源参数的严格有偿消除。

## 4. 已完成的冻结时空来源最大界

对于一个固定 c∈[1/2,1]、固定ell，U(s,t)=H_{s-t}^{c,ell}，s≥t，来源t≤1/2，而输出s≥0不设上界。相同2n时间格上总跳率q=n/c≤2n，未跳概率≥exp(−q/(4n))≥exp(−1/2)>1/2。故(4)与已证 C_fr(n)=Ologn 给
\[
 \left\|\sup_{s≥0}\int_0^{\min(s,1/2)}H_{s-t}^{c,ell}F_tdt\right\|_{1,\infty}
 \le2[3+2\log(4n)]C_{fr}(n)M
 =O(\log^2(n+2))M. \tag{7}
\]
此处F为任意固定共同正L1时空输入，不需产品分布、有限原子数或输入深度；有符号输入以总变差支付。C_fr的完整holding/faces保留。

若还对ell∈[a,b],b≤2a取sup，且同一F不依赖ell，上一轮已付尺度网给
\[
 O(\sqrt n\log^2(n+2))M. \tag{8}
\]
这是已完成的冻结时空来源工具，不是原R_dagger。若M本身为O√nW，(8)会变成O(nlog²n)W，不能把两个√n隐藏。若来源是原ρ_{t,L}，其moving-L共同支配仍未证。两项限制均不得省略。

## 5. 数值范围与下界目录

预登记 hazard_anchor_source_transfer_registration_20261007.json，三轮n8/32/128。核原B/G/log cocycle和exact未跳预算、early右端点及future端点1；另以精确Lebesgue阶梯函数测有限弱求和、循环调和的必要对数及非均匀来源质量。原核数值用80位Decimal，仅为实现诊断；时间/holding/阶梯弱范数用Fraction。正cocycle与任意输入正支配由解析证明负责。

重新读取下界总表A/B：任意非负联合权、增长XOR/Cantor/宽峰不受本工具输入限制；不以峰数或深度收费。本轮不是这些cube源的实际模拟，也不从噪声符号检查推断FIRST/CPGP/history全部合格。

运行终态：三轮n=8/32/128，共68262项谓词通过（3266/13010/51986）；全部2n格、五个格内点及future端点1均保留。原符号80位Decimal的最大cocycle残差分别9e−80、2.6e−79、6.3e−79；仅诊断、不作interval或渐近证明。精确Fraction核holding和Lebesgue阶梯弱范数。循环调和例的真弱比为H_(2n)，对数求和损失不能删。脚本/注册hash保存在结果中，主cube账不变。

独立解析审计：gated_hazard_anchor_source_transfer_audit_20261007.md 全部§1–4通过；tensor代理亦只读独核相同常数、cocycle方向和源时间范围。二者未重跑数值。

右连续版本补全：原ordered模型对固定source t展开有限mask，v_A(t,x)=G_(1−t)^A*F_t(x)，Σ_A∫v_A(t,x)dt有L1范数2^nM，故a.e.有限。原系数关于s连续且有界，变动上限的支配收敛给连续版本；有限L1时间脉冲产生右连续跳，s=1还可从左极限恢复。Frozen在每个紧输出窗以∫Σ_k(qT)^kP^k|F_t|/k!dt作共同支配，定义性L1范数e^(qT)M有限，不进入最终费用。再取尺度sup时，正一维Gc,L≤2Gc,b提供共同参考和相同支配方式。

实际正退出源的更强旧接口：若ρ_t与η_t来自同一固定物理尺度、同一初始ν的真实killed演化，则Duhamel给∫_0^s K_(s,t)ρ_t=N_sν−η_s≤N_sν；early正子源亦成立。因此这类真正同源正ρ可直接归约到c=1的完整N_s最大界，不必额外领本锚点对数费。本稿主要新增用于一般时空源和有符号forcing的合法有偿转移。原initial-data弱A_ord和moving hard平均仍未解。

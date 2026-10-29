# 第38轮高维 driver 独审

审计结论：在真实均值/PB-CDF 上下包有效、普通浮点计算按已声明前提解释的条件下，本职责范围内通过，未发现分层方差、条件 Cantelli、counts 对称性、proposal RN、贡献帽或336项统计CI的致命错误。均值积分和BE/continuity-correction常数由另一代理单独审查，本文件不替其背书。没有重建cache、没有抽新counts、没有导入执行driver。

输入是 `round038_highdim/probe_highdim.py`、`HIGH_DIM_PLAN.json`、`HIGH_DIM_REGISTRATION.json`、`HIGH_DIM_RESULTS.json` 和现有来源/缓存源码。读取结果的SHA256为 `391c502d32366d7f441af7dbc2e4dd20ec18b37696b565e26f5546d54e68f56d`。只读脚本 `read_saved_stats.py` 与实际输出 `HIGH_DIM_SAVED_STATS_AUDIT.json` 保存核对链。

## 1. 来源质量、对称条件化和缓存索引

各表box留在coarse cell内、每格峰质量相同，背景D对所有表相同，所以各表在第j格的总质量是共同m_j，通常不是1/K。联合cell按∏m_j赋质量，表由源格class-count S的固定整数段决定，则gated joint质量仍是Σ_cells∏m_j=1。χ使用前127格，其实际α=pchi≈.4980478009。driver的cut以 nominal nα+[-1,0,1]σ舍入后冻结；它明确是固定整数切点，不是精确分位点。

joint4使用四段及表0,1,2,3，joint2用中间cut的两段及表0,2。joint平均不再乘π；prior2/4才使用真实 Bin(n,α) 段概率的上下包，未把π设为1/R或公平bit概率。空段如果发生，仍须保留零系数；此次三个cut严格递增。

基准查询law是各坐标独立uniform。给定ordered bin labels后隐藏坐标独立；只给counts时是坐标排列混合，但所有坐标使用同一表核，gating只依赖源class-count，共同cube平均对坐标置换对称。因此每种ordered标签排列的目标条件分布相同，选代表排列计算条件均值/方差是合法的。这里没有假设表、class或长度之间独立。

feature缓存轴为(4,4J,Q)，目标row=bJ+j。矩阵乘后reshape(batch,4,4,J)保持这个顺序：z=ΣN_q log total；PB mean=ΣN_q p_q；PB variance=ΣN_q p_q(1−p_q)；最后轴给ΣN_q log(1+r_qb)。索引转换与builder一致。

## 2. 分层方差与 Cantelli

对固定b，每坐标class核的relative variance由 `relative_variance_class_max` 控制，故表内任一正monomial及段和的relative variance≤

`r_b=expm1(Σ_q N_q log1p(r_qb))`。

若段均值上包为u_b，则跨表中心化L2 Minkowski给绝对variance上包

`Vhat=(Σ_b u_b sqrt(r_b))²`。

driver `log_mixture` 用 `2*logsumexp(logu+.5*logr)`，恰是该式；没有额外R²，且没有把各表方差直接相加。fixed表也可以用同一class-max帽：其total核是两个正class核之和。prior的表权重用区间上端，得到同样安全的绝对方差上包。

`tail_enclosures` 的lower对每L使用均值下包m_l和Vhat：m_l>λ时给 `(m_l−λ)²/[Vhat+(m_l−λ)²]`，否则0。upper在均值上包m_u<λ时给 `Vhat/[Vhat+(λ−m_u)²]`，否则1，并与Markov `min(1,m_u/λ)`再取min。均值上下包有效时这些方向均正确。对L的lower取max，是counts层先选择一个共同L；upper取各L上包之和再截到1，是union bound。未替换为avg表max、逐坐标选L或独立尺度模型。

严格阈值边界处理保守：m_l=λ时lower=0；m_u=λ时upper仍可取1；零variance且确定值大于λ时lower=1。driver微检查覆盖这些代表情形。含−inf lower logmean的数值分支被回退为0；upper的非有限保守回退为1。普通浮点下这些保护并不构成interval认证，但方向没有反转。

## 3. 实际 counts proposal 与 RN

proposal是50%基准Multinomial(n,1/Q)，另50%均匀选择四表与J尺度，再按nominal `exp(feature0)` 显式归一后的bin概率抽Multinomial。无需把这些nominal均值视为精确真实倾斜：只要将所用proposal明确作为一个正概率law，其counts RN由该law计算即可。

令A_bj=exp(z_bj)、Z_bj=mean(exp(feature0_bj))，则在实数算术下相应component RN为 A_bj/Z_bj^n，multinomial coefficient与Q^{-n}恰好消去。因此

`D(N)=1/2+(2·4J)^{-1}Σ_bj A_bj/Z_bj^n`，`W=1/D≤2`。

driver的logD正是对z−nlog(norms)作logsumexp，使用marginal mixture RN，没有错误地条件于抽中的latent标签作单分量权重。norms由exp(feature0).mean()计算，与nominal tilt定义一致；未以原模块mean(total)的相近值代替。抽样的浮点exp/归一化与log RN往返的舍入仍属于明确未认证的FP层，不能把公式正确改写为逐样本严格数值证书。

summary中uniform65行、tilt63行，tilted b计数[10,21,15,17]相加正确。权重最大1.94210≤2，均值1.00985、描述SE .05106。保存summary没有逐样本counts/权重，故本审计不能独立重算每个RN或paired covariance；count-stream hash只提供原次运行身份，不代替原始数据复核。

## 4. C_Z C_E 与两侧贡献帽

记C_Z≥max(1,max Z_bj^n)，C_E≥max exp(loginflate_bj)。若均值inflation包有效，各目标模型在某L的条件均值上包m_u,L≤C_E Σ_b A_bL：joint的段概率上端≤1，prior的每个π上端≤1，fixed是其中一项。故

`D≥Σ_bj A_bj/(2·4J C_Z)`。

lower envelope满足λB≤目标均值≤m_u；upper envelope含Markov/union截断，故 `λU≤Σ_L m_u,L`。因此两个统计贡献都不超过

`2·4J C_Z C_E`，

同时不超过2λ。upper也有此帽，不是误用了仅lower active时成立的推理。

只读重算C_Z=1.0000000014887647、C_E=1.0001050575742725，统一帽1864.1958300937927，与注册一致。每个已保存max contribution均≤其注册 `min(2λ,统一帽)`。

## 5. 336项CI与独立重算

21个冻结h=0,.5,…,10，八模型，每项lower mean用单侧下CI、upper mean用单侧上CI，总336项。每项α=.05/336。cache、cut、N与H均在RNG抽样前冻结；N128只由确定矩阵cost gate选出，没有依尾收益增样本。

对独立0≤X_i≤cap，归一后指数函数的弦上界给 `E exp(tX_i/cap)≤1−μ+μ exp(t)`。Chernoff优化给Bernoulli KL偏差界，反演sample mean即可得到driver所用CI，无须假定贡献为Bernoulli。跨模型/阈值有共同样本依赖不妨碍union bound。于是对有效的条件上下envelope，336项同时将对应真实固定阈值tail系数夹在lower-CI与upper-CI之间，名义失败概率≤.05。

`read_saved_stats.py` 不导入driver，用Decimal55位、120步二分独立反演336项CI。最大绝对差2.8422e−13；由保存mean与SE重构ESS的最大差1.8475e−13。全部模型经验winner、最大lower-CI的h和数值一致。代码/plan/注册/结果/source/builder的hash链全部通过。只读重算耗时1.663秒，零新样本、零cache构建。

`best_lower_CI=max(1,max_h CI_lower)` 的1锚来自L=1。结果带有 `lower_CI_anchor_includes_L1=True`，因此该字段对应 `{1}∪自然窗口grid`，不能将1锚误称为裸I=[10/9,11/9]网格的下界。此次各模型max CI实际都>1，数值winner未受锚影响。

## 6. 赢家与未分离结论

以下是按经验lower mean选h后的同h上下包；CI按336项family调整。它们是固定有限grid的tail数值，不是完整weak排序。

|模型|h|lower积分估计|upper积分估计|lower CI|upper CI|
|:---|---:|---:|---:|---:|---:|
|joint2|3.5|7.106987|10.881034|1.839421|21.744445|
|prior2|3.5|7.103487|10.762993|1.837713|21.599718|
|joint4|3.0|5.448822|8.645798|1.709115|15.616174|
|prior4|3.0|5.477128|8.721494|1.724762|15.703277|

joint2−prior2的两条lower差+.0034996；它不是尾概率差，更不是joint胜出证据。同h joint lower−prior upper为−3.65601，联合CI端点差−19.76030。R4对应−3.27267与−13.99416，亦未分离。summary没有保存逐样本paired贡献，不能从两列SE推paired差的SE。

这些数据完成了真实高维gated来源、共同尺度、条件矩与counts抽样的可算接口。上下包和统计区间仍宽，没有证得joint优于实际prior，也不能从R4 lower略低就证明joint来源没有潜在收益。缺少完整阈值及连续尺度上控制，更不能宣称完整weak范数排序。原匹配下界参数资格仍失败，普通FP尚未interval认证；本轮不产生渐近C_n或log增益结论。

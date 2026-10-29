# 第44轮：旧来源的非贪心MILP count-selector pilot

2026-10-29。固定第43轮的K2/R4、n32正背景内部盒source和18个共同尺度，不调source。先实际执行N64成本试验，再对N128的训练样本运行三档MILP，每档time_limit=30秒。取实际incumbent、round为onehot count映射后，按原始真实DP重算训练weak值，选择并冻结一个输入与λ，最后独立8192点测试。本文只报告浮点探索；没有向外区间或连续query格认证，不构成新的C_n下界证明。

## 固定输入及真实平均

只读引用43轮 `round043_binary/probe.py`；source源码SHA256写入PLAN并由verify_saved重新核验。T=R/Z分两格I_0=[0,1/2)、I_1=[1/2,1)。ε=2/25、δ=1/25，四表相位(.15,.35),(.35,.85),(.65,.15),(.85,.65)，中心c_bj=(j+phase_bj)/2，峰盒B_bj=[c_bj−δ/2,c_bj+δ/2]完全位于相应格内。

ρ_b(t)=ε+(1−ε)/(2δ)∑_j1_Bbj(t)。每格质量恰为1/2。

对全部query与尺度共用一个onehot count映射s:{0,…,32}→{0,…,3}，令k(y)=#{i:y_i∈I_1}，输入为

f_s(y)=∑_bρ_b^⊗32(y)1_{s(k(y))=b}。

由于所有table都给k相同的Bin(32,1/2)分布，∫_T^32 f_s=Σ_k C(32,k)2^−32=1。实际prior g_s=Σ_bπ_bρ_b^⊗32，π_b=Σ_{s(k)=b}C(32,k)/2^32，保存了精确整数分子。

真实window/class质量z_bij=L^−1∫_[x_i−L/2,x_i+L/2]ρ_b1_Ij，用实际周期区间交叠计算。D_bk=[X^k]∏_i(z_bi0+z_bi1 X)，按非负Poisson-binomial系数DP递推。平均A_(L/2)f_s=Σ_kD_{s(k),k}。18个L固定为.25,.4,.6,.8,1,1.04,1.08,1.15,1.25,1.4,1.6,1.8,2,2.1,2.3,2.5,2.7,2.9；全部32坐标共用同L。整数窗L1、2的class质量分支精确设1/2。无Gaussian代理。

## 有限训练MILP及cap等价

训练threshold λ∈{1.5,2,3}，margin=10^−7，T=λ+margin。w_bk∈{0,1}，Σ_bw_bk=1。对query i及尺度ℓ设二元witness y_iℓ，Σ_ℓy_iℓ≤1，最大化Σ_iℓy_iℓ。每个witness约束为

Σ_bk min(1,D_bk(x_i,L_ℓ)/T)w_bk≥y_iℓ。

没有big-M。其精确实数等价性如下：对任意非负a_j和binary w_j，令S为选中指标。若Σ_{j∈S}a_j<1，则Σmin(1,a_j)≤Σa_j<1；若原和≥1，则要么有选中a_j≥1，其cap项本身为1，要么所有选中a_j<1、cap前后和相同。故Σa_jw_j≥1等价于Σmin(1,a_j)w_j≥1。y=0时约束因非负性自动满足，y=1时由该等价性，witness恰对应原始平均≥λ+margin>λ。binary假设必需；不能对fractional relaxation作相同等价主张。

对每个固定onehot映射，在数学模型中有witness当且仅当至少一个共同尺度达到阈值。但有限time_limit下，solver的y incumbent未必把其w incumbent已覆盖的全部query选出来。因此训练挑选依据是rounded w的原始D平均重算，而不是solver目标值。实际线性系数仍为双精度，solver也会丢弃小系数、容忍可行误差，数学模型等价不构成数值incumbent正确性证书；所有incumbent均重算验证并记录分歧。

## 实际执行顺序与预算

成本试验使用独立seed、N64、λ2、time_limit3秒。实际耗时约3.012秒，status1/time limit，gap约.5143，incumbent objective为35个witness；rounded原始DP重算也覆盖35/64，weak1.09375。该成本观察后才写正式PLAN，正式训练使用新seed202610294410、N128，三阈值各30秒，time_limit均保留为实际status而不称全局最优。独立测试seed202610294411、N8192，在FROZEN写出后才生成测试query。

正式模型每档2436 binary变量、2465约束，稀疏矩阵存储308868项。PLAN、COST、逐档TRAINING、最终FROZEN、SUMMARY和原始简短execution均保留；未存全feature缓存或冗余完整MIP日志。

## 训练incumbent、gap和独立测试结果

三档均status=1，message均为time limit，`success=false`；没有一档证明全局最优。下表“最优性gap”只按HiGHS返回值转录，不能当source弱指标误差；solver在最小化负witness数。

|λ|实际秒|solver witness|rounded重算覆盖|gap|dual bound|重算训练weak|训练bestfixed|
|---:|---:|---:|---:|---:|---:|---:|---:|
|1.5|30.022|113|114/128|.132743|−128|1.335938|1.324219|
|2|30.025|71|71/128|.633803|−116|1.109375|1.218750|
|3|30.021|43|43/128|.209302|−52|1.007813|1.007813|

每档都把w按每个k的argmax转换为严格onehot选择器，并完全丢弃solver的近整数数值作为输入权。原incumbent最大整性误差分别2.78e−15、3.19e−10、7.03e−8；round后重算所有真实D。三档solver所声称的query均在round后的真实平均超过λ+guard；所声称尺度的最小真实margin slack分别.0113055、.0002972、.0009521。没有用MILP可行容差代替这一步。λ1.5的w覆盖114个query，solver仅记录113个y，是尚未优化完其witness变量的现象。

按重算训练weak选择λ=1.5。其π精确分子为(4164793799,130168008,5489,0)，共同分母2^32；实际质量约(.96969162,.03030710,.000001278,0)。不是按33个count的table使用次数平均。该训练gate略高于训练bestfixed，但训练actual prior=1.347656反而更高。

冻结后独立8192点的结果如下；weak=λ×成功数/8192。

|模型|成功数|经验weak|探索性CP weak区间|
|---|---:|---:|---|
|MILP gate|6739|1.233948|[1.216861,1.250442]|
|actual prior|6749|1.235779|[1.218737,1.252227]|
|fixed0|6821|1.248962|[1.232246,1.265068]|
|fixed1|6165|1.128845|[1.109658,1.147569]|
|fixed2|6118|1.120239|[1.100913,1.139113]|
|fixed3|6803|1.245667|[1.228867,1.261859]|

MILP gate对bestfixed经验差−.01501465，对actual prior差−.00183105。训练至测试weak下降.10198975；训练中对fixed的微小优势没有在独立测试保留。没有新的C_n改善信号。

另以完全相同的训练query、阈值和来源运行43轮原平滑坐标搜索，作为优化器控制。该控制的PLAN在读取MILP heldout结果前写出、单独FROZEN后才评价相同独立测试query，不改变MILP选解。贪心控制选λ1.5，训练weak1.171875，测试weak1.123169，明显小于本次MILP incumbent的经验值1.233948；两者仍都未超过fixed0的1.248962。它说明本组某个贪心随机初始化确实会漏掉更好的count映射，不能从“贪心塌缩”直接推断所有count映射无效。MILP也有显著未闭gap，本轮同样没有排除更好的映射。控制只作为描述性比较，未追加CI或新显著性主张。

## guard、复核与证据边界

test对全部三阈值、六模型均保存raw/guard/near计数、最小绝对距离及全部near样本。10^−9 guard范围内零near，最小距阈值约5.106e−6，所选阈值raw与guard计数逐模型相同。该诊断没有变成浮点误差证书。

`verify_saved.py` 实际执行1728个有理数cap/binary等价小例，验证全部incumbent的count范围与精确prior质量，再独立重新生成相同训练query、原始DP并重算三档训练覆盖，最后从冻结selector与test seed重新stream全部8192点；六模型成功数完全复现。训练guard成功query的实际共同L witness逐个留在VERIFICATION中，可检查它们没有逐坐标独立选窗。

报告区间是12个单侧方向合计α=.05的Clopper–Pearson区间，仅对理想iid uniform采样下的“冻结双精度guard分类器”成功概率解释。尚未以向外区间包住source积分和DP，也没有把离散dyadic query扩成整个连续query格，PRNG的iid性不由程序证明。因此这些区间不构成连续Lebesgue超水平集测度的认证区间。没有重要性采样或RN权重；uniform点测试仍可能看不到更稀有的高阈值事件。

## 周期到L1的归约与本轮结论

固定本轮输入、n32与有限L。对整数N，F_N=f_s1_[−N,N)^32具有质量(2N)^32。m=2大于全部窗口半边长。内部盒[−N+m,N−m)^32内所有窗不受截断影响，周期真实超水平集E的每周期测度为|E|，故

C_32≥λ(1−m/N)^32|E|，令N→∞得C_32≥λ|E|。

这是对真实连续E的严格归约；本轮浮点随机query本身没有给其测度严格下包。与已有source一样，每个窗平均积分1，λ>1时整数窗恒1不可触发，16个非整数窗给source-only union/Markov帽λ|E|≤min(λ,16)，未采样尺度不在该帽内。

本轮实数cap模型避免了big-M，并实际找到优于同数据一轮贪心控制的incumbent；独立测试仍未胜过最好单表或实际prior。它提供优化器局限及训练过拟合的实证，不提供C_n增长下界，也不提供该source/selector族的全局否定。核心正式计算约91.61秒，成本试验约3秒，控制及复核只需数秒。

复现：将最终脚本复制到新输出目录，保持相对旧43源路径可达，按 `python3 probe.py --cost`、`python3 probe.py`、`python3 greedy_control.py`、`python3 verify_saved.py` 顺序执行。主脚本防覆盖SUMMARY；勿覆盖本次冻结证据。训练MILP的time limit使不同硬件/求解器版本可能得到不同incumbent，status/gap须以各次实际保存值为准。

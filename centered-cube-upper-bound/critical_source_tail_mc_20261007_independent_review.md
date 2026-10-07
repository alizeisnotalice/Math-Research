# Critical source-tail MC：独立只读概率及实现审计

2026-10-07。审计者 gated_radial_energy_audit；指定 6.1-sol high。已读取新 registration、原 critical_source_tail_arrival 工具/md 和新 run.py 全文。只读审查，不导入或执行新/旧 main、旧 probes、guard 或 sampler 函数，不改实现者文件。

**结论：固定256 source、每 replica8 inner draws、正常完整执行的理想连续概率接口及新实现映射通过，未发现阻塞数学错误。** 两个限制必须随结果保留：(i) 浮点 proposal/log/weight 不经 interval 认证，exact dyadic winner 不补其连续抽样误差，所报统计 CI 只属 ideal law；(ii) 外部中断后的 prefix 恢复/outer CI 尚非 run.py 自动完成的行为，须独立读出完整记录、后处理并撤销 stopped-mean unbiasedness。本稿不验执行结果、不认证 actual FIRST 或渐近阶数。

## 1. 被审文件和冻结身份

- [registration](critical_source_tail_mc_20261007/registration.json)：SHA-256 cc6aabb10f51ef2e961405056641c1667884da214169e1f2ea47b962ed194446。
- [run.py](critical_source_tail_mc_20261007/run.py)：此次全读版本 SHA-256 45aa740d552d330502f7d7f82cd71bdec7296df89a3b600f474d2a5a5ebfa5a2。
- [原 arrival module](critical_source_tail_arrival_20261007.py)：SHA-256 58fb95ab7dc77c0d1413471068d602f414d6bd57a3cedfeffdb19967454b58dc，和新 registration 登记一致。

三个既有 gzip 原输入 profile 的 hash 及参数在新 registration 中封存，run.py 启动时逐项核 hash。审计只读结构/数学，不重新计算旧24 probe 或9 RN config。新脚本的 import 使用 pure definitions，原 main 在 __name__ guard 下，不随 import 启动。

本次引用行号属上述 run.py 哈希：整数 oracle31、partition106、sampler guards114、receiver144、Hoeffding216、main224、source law246、profiles260、conditional CI294、outer summaries308–329。

## 2. 完整 source law 与 coincident contributions

完整输入四个 tensor grids k=1,2,4,8，每坐标 m/k、-Ak≤m≤Ak，qk=2Ak+1。原 μ=Σ_l w_l μ_l，w=(1/2,1/6,1/6,1/6)，W=1。

run.py 用均匀整数0,…,5 的前三值选择 k1，其余各一个选择 k2,k4,k8；每个坐标独立均匀整数 [-Ak,Ak]。在理想随机整数法则下，source 实际物理点 y 的概率正好
\[
 \Pr(Y=y)=\sum_{l:y\in{\rm grid}_l}w_lq_l^{-n}=\mu\{y\}.
\]
重合位置贡献真实相加，component label 仅记录抽样来源，不会在 full response 或 source first moment 中以标签重领质量。source_grid_membership 的整数整除判别也正确。没有 outer tilt/RN、没有丢 support edges、没有先筛超阈来源。

固定 seed 给可复现伪随机序列，不独立证明概率随机性。下文 iid/conditional independence 是预登记 ideal experiment 的法则；实际 numpy generator 与有限浮点计算没有被声明为严格连续法则证书。

## 3. cone/core prior 的 Lebesgue 密度

记 Cn=1+n log2、d∞=||x-y||∞。原 positive hard envelope 为1于 d∞≤1/2、(2d∞)^{-n}于1/2<d∞≤1。其积分为 Cn。因此
\[
 q_y(x)=C_n^{-1}\sup_{1\le R\le2}h_R(x-y)
\]
是概率密度。

core 概率1/Cn，原 offsets uniform[-1/2,1/2]^n。shell 概率 n log2/Cn，log r uniform[0,log2]；选2n个 signed face，free offsets uniform[-r/2,r/2]。对一个 face，
dx=(1/2)dr dx_free=(1/2)r^{n-1}dr du_free。
shell prior 在 (r,x_free) 中的密度为
\[
 \frac{n\log2}{C_n}\frac1{r\log2}
     \frac1{2n}\frac1{r^{n-1}}=\frac1{2C_nr^n}.
\]
除以 face Jacobian1/2 得 dx 密度1/(Cn r^n)，和 q_y 完全一致。face ties 是理想连续法则的零集；源原子并不把 receiver 的 Lebesgue proposal 变为 source probability。

## 4. conditional global-mixture RN：固定 face 必须抵消

冻结 source y、core/shell、r、signed face。free 坐标的原 interval 用真实 finite k1 reference 捕获数=2的区域划分，p_i为 exact interval 长度比例，finite source 边缘/0-count 均进入原 partition，不假设 p_i=1/2。

T 为登记 tilt odds，Z_i=1-p_i+p_iT。整个 tilted product 相对原 free product 的密度为
\[
 L=\frac{T^{\sum_{\rm free}H_i}}{\prod_{\rm free}Z_i}.
\]
代码在坐标循环外只选一次 tilted_branch，实际 proposal 是全 product mixture
Q=(1/4)P+(3/4)P_tilt，不是逐坐标独立混合。所以恢复原 P 的 RN
\[
 w=\frac1{1/4+(3/4)L}\le4
\]
正确。固定 face 的参考 H 为常数；即使形式加入 T^{H_face}，归一化中同一项相消，所以它不得进入 logL。run.py 正确只让 free hs/logZ 进入 logL，face_high 仅作诊断。core 全部 n 坐标 free。

tilt 高/低区选取、按实际 segment 长度选区间后均匀采点，是上述 product law 的理想实现。p=0或1时只采非空分支；新 finite-coordinate guard 使用真实 finite-body 边缘的 p_i 而非自由 Bernoulli。

run.py 对每个实际浮点 offset 重新使用原 exact axis_count 分类 actual_high；tilted branch 若因 rounding 产生 status discrepancy 会 assert，不静默沿用所选 flag。它还精确检查 offset 在 dyadic width 区间内。p_i/T/logZ/RN 的后续浮点计算仍不属于 exact probability certificate，范围见§7。

## 5. 同一点的 exact dyadic full winner/capture

IntegerArrival 将每个 realized float offset 化为确切 binary rational。公共 D 至少为8，且所有 binary denominator 和 source k 都整除 D。整数 x_i=codes_i·D/k+offset_i·D，正是 exact y+dyadic offset；没有先把 y+offset 以 float 相加而改变输入点。

各 coordinate arrival key=2|x_i-mD/k_l|；初始 side1闭 interval 用 ceil/floor count，随后仅登记 (1,2] 的 actual finite grid arrivals。相同 key 的所有 layers/axes/multiplicities 先更新完再评 full mass。零 count/product 状态更新正确，不除以0；所有 source layers 按原正权求和。

共同 mass 分母6Π_lq_l^n、layer 系数 (3,1,1,1)Π_{j≠l}q_j^n准确实现完整混合。候选含1、2和全部 actual coordinate arrivals；相邻候选间 mass 恒定且 response 随 R下降，所以 exact integer cross multiplication 给 complete continuous [1,2] winner，没有 radius grid。> 才替换赢家，保持最小 R tie；strict E=M>tau 也是 integer 比较。

source capture 用2max|exact offset_i|≤winner key，和同一个 x-y 完全一致。shell fixed face offset=sign·width/2，binary除2精确；free offset的 exact bound保证 max∞=width/2。因此 shell geometric factor(width/R)^n 使用的 prior r 就是这个 dyadic点的实际来源中心半径，未出现“oracle点变了却沿用旧 face”的错配。core 的 envelope density为1/Cn。

new crosschecks 仅调用旧 pure arrival_oracle 于3个新实现点/n，并对 winner/M精确相等；不会调用旧 main 或旧 probes。本稿未执行这些 crosschecks，只有实现只读审。

log gap/near threshold/near max 用浮点仅记录诊断，没有进入 winner 或 E 决策。g=tau/M及 score 再化为浮点，仍需保留未认证 rounding 限制。

## 6. score、two replicas 与完整来源 I1/I2

理想 proposal 下
\[
 Z=1_E1_{\rm capture}\frac{\tau}{M(x)}
         \frac{R(x)^{-n}}{q_y(x)}w,\qquad0\le Z\le4C_n,
 \quad \mathbb E[Z\mid y]=S(y).
\]
因为完整 capture 保证 h_R≤sup h、M>tau，score 上界方向正确。S(y)≤Cn 是原 positive envelope；没有删 source、重选 winner 或按局部受限输出重新归一化。

每个 source 的两个独立 conditional receiver replicas 分别取8个 score平均 A,B；给定 y 后，各自均值 S(y)、相互独立。因此
\[
 \mathbb E[(A+B)/2]=\int S\,d\mu=I1/W,\qquad
 \mathbb E[AB]=\int S^2\,d\mu=I2/W.
\]
W=1；I1=τ|E|，正 Tonelli 精确恢复 receiver Lebesgue 量，因为 source 对 selected Q 的积分为 m(x)，不是把 μ 当 dx。q_y 支撑覆盖每个可能 capture 的 x，support外 kernel 为0。μ有限、receiver有界窗口及 source有限体支撑使相关 first/second moments有限；所有最初换序均非负 Tonelli，σ-finite 无隐藏可积假设。

代码为两个 replicas使用登记的独立 RNG streams，source RNG独立且先抽全部 source；没有用 replica1选择 replica2 proposal。独立性是 ideal experiment 的条件法则，不是由观察两seed不同推出的数学随机性。

I2 cross product 在 fixed sample ideal model 中无偏；I2/I1 的点估计作为比值一般不无偏，代码正确叫 point_estimate。plugin max((A+B)/2-K,0) 为有偏 nonlinear diagnostic，明确不称 unbiased source tail。任何 ratio/order置信论断还需传播两条区间并要求 I1下界正，脚本没有假称 ratio CI。

## 7. 同时 CI 和浮点声明边界

outer 每个 source 的 I1 score范围[0,4Cn]、I2 score范围[0,16Cn²]。3 rounds×2 targets×256 allowed prefixes=1536；在 ideal independent source experiment 下，Hoeffding 两侧概率及 union δouter=.025 正确。正常最终固定256调用这个较保守 unionCount 仍合法。目标的已知上界 I1/W≤Cn、I2/W≤Cn²可与置信区间取交。

profile 为同source16 conditional draws，union3×256=768、δprofiles=.025；逐 source条件Hoeffding再对随机 y取期望/union合法。outer+profiles总 ideal coverage≥.95，随机profile与outer估计有关联不影响union；不要求二者独立。

这些预定区间很宽：profile半宽为
4Cn√[log(2·768/.025)/32]，约2.35Cn；
final I2半宽为16Cn²√[log(2·1536/.025)/512]，约2.42Cn²。
上界裁剪后通常不提供source-tail/order分辨率。数据点估计、ESS、replica差别都只是诊断，不允许在缺显著区间时外推阶数。

有限浮点 RNG、exp/log radius、p_i转float、T转float、segment插值、logZ/RN及score都有未界定的连续法则/rounding误差。integer oracle只保证对每个已实现dyadic点的M、R、E、capture精确，不证明其采样点按理想q/Q分布。注册和run.py输出明确用 ideal exact continuous model 统计区间及 floating_proposal_law_not_interval_certified；本次通过审计严格保持此边界，不把CI提升为实际浮点实现的certified95%区间。

## 8. 外部中断：注册合同与已实现行为的区别

正常固定N完整执行没有 optional stopping问题。登记若外部中断，uniform finite-prefix事件仍可在恢复一个已完成 source prefix 后使用；停下的样本均值/AB均值不能继续称 unbiased。

当前 run.py 写逐source jsonl gzip，每32 source flush；只有整个 round完成后生成outer summary。它没有自动保存每个prefix outer CI、没有中断异常收据或gzip partial recovery/续算实现。故目前不能声称“代码已自动提供任意中断prefix置信收据”。如果被中断，需核 gzip中哪些完整source行真正可恢复、只取一个完整prefix，不把半个replica/source计入，然后按预登记 union1536单独后处理outer CI并明确撤销stopped-mean unbiased。不得重新启动已完成source或假称正常fixed_N_complete。

此项是中断时的额外操作范围提醒，不阻塞正常预登记fixedN run；已及时告知root及实现者，不修改其代码。

## 9. 最终审计范围

被审snapshot正常固定预算的source sampler/coincident components、core/cone Lebesgue prior、global free-coordinate RN、两replicas条件独立接口、完整continuous winner/capture和I1/I2校准均通过。sampler guard、3新oracle crosschecks和main是否实际通过，须看作者最终执行收据；本审计不执行也不冒称已验结果。

这些原finite atomic模型统计不认证 actual FIRST/fullfuture/LCA/history，若出现疑似一般L1反例仍需正宽化及winner/阈值稳定性。当前所有统计只能提供理想模型有限探路，没有新一般sqrt(n) endpoint结论。脚本若后续更改hash，需补读变更而不能沿用冻结身份自动宣布审过。

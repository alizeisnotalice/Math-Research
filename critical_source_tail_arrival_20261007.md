# Critical source-tail：continuous arrival oracle 与合法 cone tilt 登记

2026-10-07；指定 6.1-sol high；只编辑新 prefix。承接 [source_occupation_tail_geometry](source_occupation_tail_geometry_20261007.md) §7，旧33记录封存。本轮先建立准确 arrival/RN 接口，只登记小型三轮原空间交叉校验；不登记大样本 I1/I2/source-tail 估计。

完整来源固定四个 n 维正 tensor grids：k=1,2,4,8，半宽 A=16n，每轴 m/k、-Ak<=m<=Ak，q_k=2Ak+1；正混合权1/2,1/6,1/6,1/6，W=1。所有重合物理 atom 的质量真实相加；源 layer 数在 n4/16/64 中保持4，没有以标签数代替维数。它实现 B41 型共同分辨率正混合机制，尚未实现 A13/A14 完整下界输入。

在参考 k1 层、R*=3/2 的 interior uniform Lebesgue lattice phase，坐标捕获数是1或2，二者概率1/2。先以精确二项尾选择 t_n=min{t>n/2:P(Bin(n,1/2)>=t)<=n^-2}，得到 t4=4、t16=14、t64=47，预设 tau=(1/2)(32n+1)^(-n)(3/2)^(-n)2^t_n。该规则只映射 reference 的真实计数；不是断言完整 mixture 或 source-cone 的超阈概率等于二项尾。完整 max 全部重新计算，不 componentwise 选 winner。

Receiver x 的 continuous [1,2] candidate 是1、2和所有坐标到轴 grid atom 的 arrival sidelength 2|x_i-m/k|。这些值精确 Fraction 去重；同一半径上全部坐标、层、左右同距离 atom 同时更新 capture counts，再用完整正 mixture mass/R^n 比较。任意两个 arrival 之间每个坐标 count 恒定，故完整 mass 恒定，response 随 R 单调下降；所以候选给真正 continuous winner，无 radius grid。候选可以包含没有完整 tensor atom 入场的冗余坐标 arrival，但不会漏完整 mass 的变化。

合法 tilt 保持 source y、shell r、face/sign 的原 prior；只改变 free coordinates 的条件 law。令 H_i 为参考有限 grid 在 R*=3/2 捕获两个坐标 atom 的实际 high-set indicator，p_i 是该集合在 free-coordinate 原 interval 的精确长度比例，包含 source finite-body 边缘。T=e^beta 是登记的正有限有理数，Z_i=1-p_i+p_i T。条件 product tilt 的密度相对原 cone prior为 T^(sum H)/prod Z_i；再混入 eta=1/4 的原 base branch，真实 RN 为

\[
 w_{RN}=\frac1{\eta+(1-\eta)T^{\sum H_i}/\prod_iZ_i}\le4.
\]

所以 base 支撑完整保留。face 的固定坐标不当随机 free coordinate；core 分支则在其原 cube1 里改变全部 n coordinates，使用同一长度/密度计算。这个 conditional RN 与 cone Jacobian 一同恢复 receiver Lebesgue 积分，不把 source probability 当 receiver volume。给每个来源仍须使用整份 mu；未来 outer source tilt 若有需要，必须另给完整 mu 的 RN。

三轮参数、精确 critical尾、seed、8个 bounded rational receiver oracle probes 和3个 actual-cone high-set RN configs 已在 registration 保存。复杂度为 O(n sum k) arrivals、O(H log H)排序及 rational integer 位长成本；不会枚举 q_k^n atom 来处理主 n64 样本。独立 scalar all-atom 实现只用于 n1/2/3 小源交叉校验。主三轮 rational probes不构成随机 volume/source-tail 样本，故无概率CI或阶数结论。

## 1. 完整连续 oracle 的精确处理

设 c_li(R) 为第 l 层、第 i 坐标在闭 interval [x_i-R/2,x_i+R/2] 捕获的实际 grid 点数，q_l=2Ak_l+1。完整 response 是

\[
 U_R(x)=R^{-n}\sum_l w_l q_l^{-n}\prod_i c_{li}(R).
\]

初始 R=1 用闭 interval 的精确 ceil/floor 计数。之后把每个实际坐标 atom 的 2|x_i-m/k_l| in (1,2] 登记到 Fraction radius key；同一 key 的 (layer,coordinate) multiplicities 用 Counter 聚合。全部同距离到达先更新各层 count/product，再评价**完整**上式。来源重合 atom 不是拆标签合作：每个 layer 的质量按其正权贡献，故共同物理 atom 的质量相加正好是原 mixture。

产品更新保存每层 zero-count 数和 nonzero integer product；不能遇零时除以旧零。响应用 Fraction，各候选比较 > 才替换赢家，所以真正 equality 保留最小 R。R=2 一并含所有边界 arrivals；没有 tolerance grouping、半径 grid 或浮点 winner。主 probes 的输入坐标为明确有理数，它们的全部 continuous R 都被 oracle 覆盖；这不表示将来浮点 receiver 近似的数学误差已被认证。

Independent all-atom scalar 方法在小维数直接构造全部实际 atom，重合位置质量合并，按 l-infinity distance排序、同距离全部质量累计，并在1、2及实际 atom arrivals 比较。它不使用 tensor coordinate arrival counts。主大源无法枚举全部 atom；在每个 probe 的若干实际候选处另以闭 count product 从头重算，核 incremental 结果。

## 2. RN 证明的真实条件 law 与实际几何

给定来源 y、shell r、face i/sign之后，free receiver 坐标 x_j 的原 interval 是 [y_j-r/2,y_j+r/2]，密度1/r。参考层 k1 的 count在中心 x_j 越过 m±3/4 时改变；用这些**真实有限 grid**边界分割该 interval，midpoint精确计数为2的 segments 组成 high set。令 p_j=这些 segments总长/r。finite-body edge、来源 fractional phase、r 的变化全部进入 p_j，而不是把它设成1/2。face 坐标已固定，单独计其参考 high count；不会被 tilt 改动。

选择 tilt density相对uniform interval为 T^{H_j}/Z_j，每个 interval分支的精确归一化为1。可执行抽样是：先以 p_j T/Z_j 选 high，否则选 low，再按各 segment长度选一个 segment，最后在其内部均匀采点；p_j=0或1时仅采非空分支。给定 source/r/face，各 free coordinates独立采用它。全局再以 eta=1/4选择 base 或整个 product-tilt分支。它的 mixture RN 正是上式，而非逐坐标独立混合的 RN。

Core分支的原坐标 interval为[y_j-1/2,y_j+1/2]，所有n坐标都 free；同一证明成立。原 core/shell概率、cone face/sign、log radius prior均未改变，conditional RN因此也是对整个 receiver prior的合法RN。原 shell积分的Jacobian为 n r^{n-1}dr dnu；它与前稿 density sup h/C_n合用，任何后续 S estimate仍恢复receiver Lebesgue，不以 source law充体积。

小守卫的 reference status law是由这些实际 p_j 得到的 Poisson-binomial，而非自由toy coins。原law系数P_h和tilt law系数P'_h由精确polynomial乘法计算；proposal Q_h=eta P_h+(1-eta)P'_h。逐 h 的 w_h=1/[eta+(1-eta)T^h/prod Z_j]检查：sum Q_h w_h=1，sum h Q_hw_h=sum p_j；加入实际face高计数后critical tail也完全相同。RN<=4给E_proposal w²<=4，另以精确Fraction直接核验。

## 3. 冻结阈值与有限边界资格

精确 binomial tail 分别为1/16、137/65536、2092620625940629/2^64。相对旧baseline(32n)^(-n)，预设 tau 的 alpha_n约1.5318、12.0885、366.4130；它确实随维数改变，不再使用固定alpha .1/.5/2冒称覆盖critical。

然而binomial规则只适用于参考single grid的interior uniform Lebesgue phase；完整 mixture、sourcecone prior、finite边界和 continuous winner不是该law。选择tau时未使用测试probes，且每个决策仍由 full response给出。在count=t_n的reference response贡献等于tau，其它正层可让full M严格大于tau；所以不能把full mixture levelset概率当作登记binomial tail。

RN guard数据中的 prior_critical_tail/weighted_critical_tail 指实际参考 **high-count status** 的 tail，而非完整M的tail。finite-body边缘还可能出现count0；此时low不一定为1，sum H>=t_n不能自动保证参考component response达到tau。两者仅在所有参考坐标count属于1/2的interior条件下按前述公式对应；原闭边捕获和full M阈值始终由完整oracle另算。

将来可使用这个oracle和RN探索full source的I1/I2/截断尾，但仍须抽真实outer mu、对所有inner样本重算continuous winner，并对两独立conditional replicates分别加完整RN。当前只使用有限有理probe和conditional status积分；没有用它们估计总体矩、volume或概率。

## 4. 三轮完成状态与可复算记录

全部Fraction守卫通过：

|n|oracle probes|实际候选数量范围|同时到达groups|最大同距multiplicity|critical超阈probes|
|---|---:|---:|---:|---:|---:|
|4|8|9–29|116|9|4|
|16|8|17–125|352|26|1|
|64|8|17–509|1112|85|1|

超阈probe数是设计的有理位置统计，包含一个刻意在reference high-set内部选取的prime-denominator receiver，**不是**Monte Carlo事件频率或总体压力。主24个连续oracle probes有143个实际候选direct tensor重算全部相同；n1/2/3的18个全部atom independent scalar oracle cases全同。22/24主probes的full winner至少与某一component winner不同，直接记录componentwise选择为什么不能替代完整mixture。

三轮各3个真正conditional core/shell RN configurations全通过。n64的三种prior reference critical tail分别约.00011344、.01796825、.00015670；中间值与uniformphasebinomial尾差两个数量级。n4某一face条件下reference count达到4是解析不可能，conditional tail准确为0；这个**固定参考count事件**的空性不代表full winner E或source-tail为空。各配置RN second moment在2.0824–3.7513之间，均准确<=4。

首run已保存n4/n16 profiles，在n64输出一个大整数精确Fraction时触发Python默认4300十进制位的序列化限制，未发生数学assertion失败。单独保存failure receipt，修复仅为解除该输出限制及读取已保存profiles续算；n4/n16 profiles没有重算或覆盖。恢复run正常exit0，约0.4679秒；18个微型scalar cases作为入口guard复核，n64未保存部分重新计算。结果记录当前script hash与resume flags；failure receipt记录初始script hash和前两profile hashes。

三个gzip JSON profiles包含完整source参数、输入hash、所有真实有理source/receiver位置、完整候选响应trace、原winner、full/component差别、high-set实际segments、原/proposal状态laws、RN-by-count以及tail identities。无浮点response或tolerance判定；np RNG仅选择有限有理实现probe，不被解释为continuous MC抽样。

## 5. 范围与下一缺口

已完成continuous arrival oracle、确切critical映射、相对真实cone prior的合法RN，以及小型原空间交叉校验。没有运行大样本I1/I2/source-tail测试，没有新的维数阶证据，也没有refutation或付款一般source-tail合同。

下一轮若启用rare-event估计，需要独立登记sample预算、tilt/base mixture、outer来源proposal、原source近边界率、完整continuous winner和confidence；本登记不自动开启它。完整source是有限positive atomic模型，若候选违反一般L1合同，还需正宽实现及严格阈值/winner稳定迁移，不能仅凭atomic oracle宣布actual geom反例。主一般sqrt(n)预算仍未闭合。

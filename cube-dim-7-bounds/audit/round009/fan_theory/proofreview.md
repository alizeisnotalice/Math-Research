# R009 fan theory独立Sol审计

2026-10-27。新独立clone起点 `75eedb49b0e95dff748081d658d5a9f53a593974`，仅写 `audit/round009/fan_theory/`。全文读取作者fan_proof.md、fan_normals.json、registration及source/code-review/receipt，另只读Terra generator/checker/report和独立cofactor审计报告。本审没有执行数学脚本、科学枚举、MC、Green或数值积分。下列有限数值计数来自解析组合恒等式或明确归属的Terra实际运行，不是本审新算术输出。

**结论：接受35/29 restricted arrangement、36行active-rank完整性和有条件的full-fan H_epsilon/Fourier/DD6接口。** 作者证明没有把ray catalog当full fan或三角剖分，也没有从conditional公式宣称g/I/C7。46条ray的实际计数由Terra及另一个完整cofactor审计负责，本审不冒称重枚举。

一项使用约定须明示：fan_proof §7的k阶directional-derivative常数按单位方向成立；任意方向h须乘|h|^k。本报告按此正确版本接受该估计。未发现要求改动核心证明的实质缺口。完整chambers、triangulation、连续算术和g/I均尚未给出。

## 1. 坐标、lambda截面和整数归一

在fundamental ordered cone ξ1≥…≥ξ7≥0，w_i=ξ_i−ξ_(i+1)、ξ8=0，逆式ξ_i=Σ_(k≥i)w_k。矩阵M上三角、diagonal1，detM1；M和逆差分矩阵均整数。因此gcd(w)=gcd(ξ)，primitive w对应primitive ξ，频率Lebesgue Jacobian为1。

lambda=Σ_iξ_i=Σ_k k w_k，而不是Σw_k。非零w≥0保证lambda>0；截面lambda1是compact6-simplex，归一点为ξ/lambda。u_k=k w_k给标准simplex坐标；若用u做完整七维频率积分，det(dw/du)=1/7!，不得沿用w的Jacobian1。该论断涉及频率体积，不把六维section Euclidean measure误当同一对象。

对v=(1,±1,…,±1)，交换求和得到v·ξ=Σ_k s_v(k)w_k，s_v(k)为prefix sum。u式系数s_v(k)/k可能有分母；选整数w normals是正确exact-kernel输入。作者JSON的lambda向量1..7和七个coordinate walls符合这一约定。

## 2. 35非cut、29cut、first-hit16+8+5及边界

所有walk第一步+1。relative simplex interior有每个w_k>0；若prefix normal全非负，其第一项1保证form严格正，不切interior。若有负prefix，则也有第一项正，可在两个相应axis附近取interior points使form异号，再沿内部segment取zero，所以确切cut interior。

非负length7 meanders按终点7−2j分类，j0..3，reflection count是binom(7,j)−binom(7,j−1)。telescoping得binom(7,3)=35；它们自动第一步+。余29条必须第一次到−1在step3,5,7：

- step3仅+--，四步任意tail，16条；
- step5先走length4 Dyck path ++--或+-+-再−，两步tail，8条；
- step7由length6五个Dyck paths后接−，即+++----、++-+---、++--+--、+-++---、+-+-+--。

三类不交且穷尽，数16+8+5=29。已按文本文义逐行核JSON的29 prefix rows及固定lex tails顺序，含两个step5 prefix和五个step7列式；没有在本审运行64-path生成器。

两个internal normals首项均1；若其lambda1截面hyperplane相同，restricted affine forms成比例。差linear form在整个lambda1 affine plane为0，因该plane的linear span是R7，差必为0；首项1固定比例。prefix consecutive differences还原v，所以无global duplicate section。interior cut使其trace含相对开片，故也不能只由bounded-simplex同trace产生别的hyperplane。boundary上trace重合可以发生，不能据此删掉global distinct rows。

discarded nonnegative-prefix form的orthant零集合正好是所有positive-prefix坐标的w_k同时0的coordinate face；其normal是这些active coordinate-wall normals的线性组合。它不会带来遗漏的boundary subdivision或active rank。故七walls加29cuts足够描述完整restricted arrangement，但body-zero label必须仍看全部64 original forms。

## 3. 六active normals与射线无增无漏

对任意36行中的六行矩阵A，signed6-minors z_j=(-1)^j detA_without_j满足Az=0。非全零等价rank6，其kernel一维；全零是rank<6，不足以certify unique ray。若kernel能定向为nonnegative且非零，则保留零分量、按正gcd primitive；mixed signs不在orthant，必须reject。不同active subsets可能指向同一退化ray，正确dedup key是primitive w。

无extraneous：对accepted w，取正orthant中避开全部29hyperplanes的generic点趋近w；有限sign patterns给一个full-dimensional chamber子列，closure包含w。所选六行都是该chamber的有效weak inequalities且在w active，不必每行都为irredundant facet。六行rank6、lambda(w)>0使加lambda后的rank7；它们在section中孤立该点，故w/lambda为chamber-section vertex。

无missing：每个full-dimensional chamber只可能由这些36supporting hyperplanes界定。其section vertex的active normals加lambda rank7；所有homogeneous active rows又annihilate非零ray，rank≤6，因而恰rank6。选其中六个独立者，枚举必含，cofactor即该ray。即使有超过六个active rows或w多坐标0，仍存在这样的subset。

因此遍历binom(36,6)=1,947,792的逻辑完整。按选c墙行分层的count为binom(7,c)binom(29,6−c)，c0..6给475020、831285、498771、127890、14210、609、7；这是组合覆盖恒等式。是否实际执行每一个子集仍需要独立运行证据，不能由该恒等式或source的assert语句替代。

## 4. 46计数的证据归属与source读取边界

Terra候选run.py/check.py已只读审源：rank/nullspace、首尾pilot cofactors、nonnegative orientation/gcd、ξ cumulative sums、lambda、全36incidences、全64body signs及nonzero6minor witness与上述数学规格一致。prefix iterator消耗前4096，separate suffix iterator取末4096，中间只再取TOTAL−8192，避免pilot重复；1903 terminal spans和wall strata是coverage合同。source review不是本审runtime证据。

Terra候选实际receipt/report宣告全1,947,792组合执行，rank<6为552979、rank6 mixed为1066827、feasible hits327986，dedup为46 primitive rays。另一个独立Terra审计 `audit/round009/fan_rays/sourceaudit.md` 记载以七个signed6×6 cofactors独立重算全量集合，46 rays、normalized ξ、incidence、body signs和witness逐项一致。这是46计数的实际数值来源，本审未再执行它。

该独立审计有必须保留的执行故障：原程序完成所有计算/比较/hash并输出完整PASS JSON后，末尾except BaseException误捕自身SystemExit(0)，又输出FAIL并exit1。原stdout/receipt和修复portable checker均被其封存；修复版未重跑全科学枚举。可说实际完整计算及比较已完成，不可说那个原始命令正常exit0。本审阅读其sourceaudit并锁hash，未替其重新运行或消除故障。

46ray catalog仍只列极射线，没有给每个full-dimensional chamber的H/V描述、interior witness、coverage/disjointness、positive triangulation或全部cone列表。ray数量、ray validity或section volume总数不能单独替代这些缺失输入。

## 5. Spectral catalog与group multiplicity

w=e1给ξ=(1,0,…,0)，是六coordinate walls的arrangement ray，但64个body dot均为1，没有zero body normal，ψ在附近只有一个linear branch。ξ/ψ(ξ)在谱单位球facet relative interior，不是其vertex。这一例已经证明arrangement catalog不能等同旧21-orbit spectral-vertex catalog。

一般r_body=rank{v:v·ξ=0}必须用全部64forms。ψ单位球的active supporting gradients由一个与ξ点积ψ(ξ)>0的gradient，加zero-body-normal span生成，rank=1+r_body。因此真正谱vertex当且仅当r_body6；排序walls给的rank6不能代替它。discarded forms可在wall上zero，尤其不能仅看29cuts判断body rank。

signed-permutation group order2^7·7!=645120，ordered nonnegative ξ给ray orbit canonical representative；boundary stabilizers须计数。未来可将完整fundamental-cone triangulation的真实group images用于覆盖R7至measure-zero boundaries。对一般x，单cone项经g变换是原项在g^-1x处的值，不能把一个cone项乘群阶。若已保留全部opposite cones，不再额外乘2；若只保留每个opposite pair一个，才用2Re。只有完整process after sum享有所有对称。

## 6. 有条件fullcone Fourier乘积与absolute Fubini

必须先有真正完整、interior不交的positive simplicial fan。设C=R[0,∞)^7，七ray独立且同一sign chamber，q_i=ψ(r_i)>0。weak signs在全部ray一致使ψ(Ra)=Σa_i q_i。频率change of variables带|detR|，逐个正a_i积分给

`p_t(x)=(2π)^−7 Σ_C |detR| Π_i(tq_i−i x·r_i)^−1`。

phase的−i号由原Fourier inverse e^(+ix·ξ)固定，t>0实部tq_i正。nonunimodular determinant不可丢；ray positive rescaling同时缩det和对应denominator，恰抵消。w ray用ξ=Mw进phase，detM1；u ray要记1/7!频率Jacobian。

时间从epsilon>0积分，absolute double integral为 `∫ e^(−epsilon ψ(ξ))/ψ(ξ)dξ`。由ψ≥8|ξ|，近0的radial行为≤常数r5，远处exponential decay，所以绝对有限；所有finite cone/time/Fourier交换合法。单项real part可负，必须保留complex phases及signed cancellation。该绝对性也覆盖下节partial-fraction后再合并，不能把尚未给的fan本身视为已免费取得。

## 7. 七维DD6正号、重节点和branch

设alpha_C=|detR|/Πq_i，b_i=x·r_i/q_i。要积分的是Π_i(t−i b_i)^−1。distinct knots的partial fraction coefficient为

`A_j=i^−6/Π_(k≠j)(b_j−b_k)=−1/Π_(k≠j)(b_j−b_k)`。

∫epsilon^T每项为A_j[log(T−ib_j)−log(epsilon−ib_j)]。upper端logT系数ΣA_j=0，余下log(1−ib_j/T)→0，所以负的lower boundary再乘A_j的负号，得到**正**DD6：

`∫epsilon^∞ Π_i(t−ib_i)^−1dt = DD6[log(epsilon−i·)](b1,…,b7)`。

n7恰order6且i^−6=−1，是此正号的关键。右half-plane epsilon>0选log principal analytic branch，无real b的pole或branch cut。

重复knots使用confluent/Hermite divided differences，因为解析节点函数的连续极限存在。直接微分给

`f_epsilon^(k)(b)/k! = −i^k/[k(epsilon−ib)^k]`，k≥1。

all seven knots0时DD6=f^(6)(0)/6!=1/(6epsilon^6)，与∫epsilon∞t^-7dt一致；all knots同c时同样给1/[6(epsilon−ic)^6]。零/重复real knots因此分析上安全，不允许在distict-node显式分母中直接代0。rational x/rays使b相等可准确决定；overlapping interval balls不证明exact equality，未来稳定confluent/near-confluent算术仍需专门实现。

于是conditional H_epsilon公式为 `(2π)^−7 Σ alpha_C DD6`，没有负号、没有额外factor6!，且branch与cone乘积相符。本审这些是解析推导，没有运行DD数值小例子。

## 8. smoothness、单位方向界及原G transfer

真实H_epsilon Fourier multiplier e^-epsilonψ/ψ给C∞；|Imx|<8epsilon时e^{−(8epsilon−|Imx|)|ξ|}/(8|ξ|)仍可积，故holomorphic tube。在此tube每cone denominator t q_i−ix·r_i实部≥(8t−|Imx|)|r_i|>0，DD log branch兼容。

对单位direction h，|(h·ξ)^k|≤|ξ|^k，sphere areaS6=16π3/15，radial Gamma积分得

`|∂_h^k H_epsilon| ≤ Γ(6+k)/[120π4·8^(7+k)epsilon^(6+k)]`。

这正是作者§7常数；任意非单位h的正确界乘|h|^k。未来interval AD若沿未归一化方向必须支付这一因子。该界并不替代full-face arithmetic cover本身。

若未来在closed face证明连续H_epsilon≥g0，则G≥H_epsilon在boundary，G的−6 homogeneity把lower带到cube interior；这里使用G homogeneity，不能误称fixed-epsilon H本身homogeneous。positive-part全空间delta-source/all-common-r/generalL1原接口仍属G；LH_epsilon=−p_epsilon不是原delta source，不能把截断核直接塞入原finite-source theorem。

若未来B≥∫face H_epsilon^(7/6)，原R001全face short-time预算给 `I_G≤[B^(6/7)+E_epsilon]^(7/6)`，E=(7/8)2^(9/7)π^(-13/7)epsilon^(8/7)。也可独立复用已有full-G I budget配H-derived g。没有g或I连续certificate时，ray count、谱体积、单节点和conditional cone公式都不能产生C7 bound。

当epsilon↓0，先对完整cone sum取消后才取limit。某些单cone zero-knot terms可各自diverge，即便offS3真实G有限。原局部p_t=O(t)与连续性证明支撑offS3整体limit，S3上真实G∞。不能把individual denominator奇点当原kernel奇点。

## 9. 状态与停止

本审接受上述理论和conditional接口；没有执行cofactor、fan/chamber枚举、数学脚本、quadrature、MC或MATLAB。只运行文件SHA256/manifest等元数据核验。46计数严格归属Terra实际执行及独立cofactor审计，不归属作者Sol source review或本审运行。

尚无full fan/triangulation/H值、g、I或Cn新结果。原root repository除本新audit目录无修改，无commit/push或子代理。此批结束停止，不自动启动R010或其他策略；主集成人保留采用及后续证明责任。

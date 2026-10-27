# R040 共同f、CF障碍及weighted cluster接口独立审计

**ACCEPT，按以下精确scope。** 初审`upper/round040/mass_sensitive/PROOF.md`身份为`4d6a22fcfb343cbdb851cea6d017c799e8b04ece5d9b9f0813cec76afba530b1`。两处scope文字现已修正，最终正文`8bdd6b862daf9dda2a021ee85fb8ad7a719a2f541d55915a70067e5559356910`及新增树附记均已增量核验并ACCEPT，详见第5节。

核心CF no-go是固定R032架构、h在原允许区间、near-linear t窗口且τ≥1中的pointwise证书障碍；不是所有t、所有选择算法、所有共同质量方法的障碍。证明保留同一个原f、所有有限witness及真实选择，基本支付/cluster接口适用于任意维数。leading渐近继承R037的充分大n结论，而不是由最小一维检查推出。

两项解释须始终保留：

* 近1上界是λ|U_all|/||f||以及selected-union比率。即使另证U_all⊂{M_nf>λ}，也不能由此上界整个λ|{M_nf>λ}|/||f||；正文§2.2“其弱比率很小”仅可按有限cover支付比率理解。
* §5截断f保持质量后，须保持原h作为固定参数，或使用正文§3限定的局部质量/体积预算规则。任意依赖全f的oracle可能在截断后改变预算；本报告不把这种oracle纳入不变性断言。原允许预算区间本身不变，因此在整个区间上取inf的no-go仍成立。

## 1. 整个允许h区间的CF no-go

f_η=(1+η)λ1_(∪S_i)给每个M_i=(1+η)λV_i，所有witness严格一致。V_i≤h_i≤(1+η)V_i给w_i/(1+η)≤p_i≤w_i。

仅“某条旧路径eligible”一般不足以保证改预算后同一selected结果；但被调用R037构造有更强性质：在每个source center，**全部source**原权重之和已低于τ/2，最后最小cube由全部source零捕获。故所有remaining source始终eligible，降低每项p保留该性质；任意半径合规tie都选全族。small Q两两不交，actual increment不依赖tie。于是D为空，CF无分流自由度。

有限开stars共同包含0，故含一个共同正体积邻域，真正esssup恰为Σp_i。对h逐项取上端点同时最小化每一项，严格得到inf_h L_CF=Σw_i/(1+η)，不是仅测试一个h方案。

R037给W≥n log n−O_C(n log log n)−o(1)。η=1/log n只扣除O(n)，且可在R037 a≥a_0的量词后再取a≥n(log n)²及a>2t，因而leading系数1障碍成立。有限actual histories足以反驳声称对所有尺度所有原历史成立的该pointwise改善；此事实不把原算子变成尺度截断版本。

共同star几何把全部disjoint increments放入半径(a+1)R的cube，故|U_all|≤c_aV_max。原f给F≥λV_max，mass LP又由单约束下界V_max及共同开交集上集中质量的可行解得恰好V_max。因而有限支付比率≤c_a，a大时近1。使c_a<1+η后，以每个x∈Q_i为中心的扩大cube还严格见证x属于真实level set，但这个包含关系不反向控制完整level-set测度。

## 2. endpoint与cheap-hotspot planting

endpoint h_i=M_i/λ时，每项(v_i+c_i)M_i/h_i=λ(v_i+c_i)，故任何合法flow都满足∫fG=λ|U_all|。以π=f/F取加权平均只是精确的未知原比率；改变flow不能改变该目标值。取esssup才再次松弛。

probe共同缩放ε时，M、V、h、increment、F和U都乘ε^n，而p不变。背景与probe有限且支撑有界，可平移使跨分量stars、资格区、cover E及f支撑均不交。这样选取的全局半径顺序可按两条内部半径历史归并；若存在同半径tie，应保留各分量原内部tie顺序。所有局部eligibility与实际increment精确解耦。

endpoint或局部homothety-covariant预算规则中，真实比率严格为正文加权混合式，ε→0回到背景；高sup热点则不缩小。若cover区及witness区都分离，CF可行流集合是直积，esssup是两个分量load的max，因此min也精确为max(两分量min)。没有把独立质量重复计数。

由此“只要高sup就有uniform原比率界”的全实例命题已蕴含去掉高sup条件的同一uniform界。这个no-shortcut推论只针对sup前提，不反驳额外要求高load拥有正π质量的dichotomy，也不保证任意global budget oracle保持背景不变。

## 3. f-weighted cluster LP与方向

θ列是同一几何族的radius-truncated common-star clusters，并保留最大半径代表j。每列中actual increments仍相互不交；覆盖每个selected i至少1次给λU_sel≤c_aΣθM_j。所有M_j仍为同一个原f的实际积分，故不同代表重复覆盖的f质量已在ΣθM_j中重复付费，没有免费“一次付清”假设。

有限列覆盖可行（每个i至少有含自己的列）、M_j严格正且有限，故min S_f存在。每个primal可行θ给上界，取最小值仍给λU_sel≤c_aS_f。普通cover LP dual是s_i≥0、每列Σs_i≤M_j，目标maxΣs_i；**dual可行仅给S_f下界**，不能据此推C_n上界。

原cover估计使用h_i≥V_i，给|D|≤(A/τ)Σv_iV_i/h_i≤(A/τ)|U_sel|，所以接上(1+A/τ)c_aS_f/F方向正确。若放宽h到V_i以下，不能无检查复用该cover系数；本文固定的原区间没有这个问题。

共同star probe可一列覆盖全部，S_f≤M_max≤F。远离分量没有跨分量common-star列，列/行分块且代表原质量随缩放乘ε^n，所以S_f精确加性。这里没有uniform S_f/F估计；任一树/贪心构造仍须支付其全部代表成本，不能从层级结构自行删去overcount。

(11)由0≤G≤M及E_πG≥(1−ε)M对非负缺额M−G直接用Markov得到。它是共同f加权的必要条件，sup大不推出该条件。

## 4. 一般L1、全尺度及实际核验范围

在固定有限input中，应截断到**所有input witness**的并集，而非仅selected stars，以保持拒绝项的原见证。该并集有界，原质量不变、F只减；固定原h仍合法，或局部预算规则保持不变。一般有符号输入可先取|f|。未来若对全部这样的有限原实例给uniform质量上界，再经原compact-level-set有限覆盖及Lebesgue内正则性处理全尺度算子。本文尚未完成该正向uniform步骤。

本审计没有重跑作者已有一维toy、R037构造或生产LP。实际只用本机Fraction核查：任意固定标量预算比的区间bracket、endpoint之和及混合恒等式；另给一个新微型mass-consistent cluster成本检查。后者取stars(-1,1)⊂(-2,2)与远离的(10,11)，f=2于其并集、λ=1，实际M=(4,8,2)；一列共同cluster与远离singleton给primal成本10，显式dual s=(0,8,2)也为10，检验LP方向及分离加性。未使用自由不一致M，未调用solver。它只是接口核验，不是uniform界。

复现：

```sh
PYTHONDONTWRITEBYTECODE=1 python3 cube-dim-7-bounds/audit/round040/mass_sensitive/verify.py cube-dim-7-bounds/upper/round040/mass_sensitive/PROOF.md
```

关键只读依赖身份：

* geometry R037 `general_cover/PROOF.md`：`6261341afbdafe868616e3e365f3484e0dab7e9e8e6ca09f83c56c93751df90e`。
* geometry R039 `joint_mass_dual/PROOF.md`：`3fa6117e186343466f3a1cbeb18ed4084a12f4b303e966f14480687a1edd24f6`。
* 作者R040 `verify.py`：`a0016d33e740bdaa78f6b7bf289802fb1f4ef391159f6f92be9bf6d97dd545ac`；原receipt文件：`cc43381dbac8a24599cc9a900900c57697e2629218a658115f739720f3be9734`。

本审计只新增本目录三个文件，未做Git修改、生产搜索、旧案例重算、MATLAB或n=7工作。

## 5. 最终文字锁定与新增树附记

最终实际源身份：

* `PROOF.md`：`8bdd6b862daf9dda2a021ee85fb8ad7a719a2f541d55915a70067e5559356910`。
* `editorial_adaptation.json`：`6e1748499ae498fe584c530752c57b2fa58b38279f450f94cb3c97eb6c80d452`。
* 新`CLARIFICATIONS.md`：`ace4081168847f1bc5ef2ff976be643193ed1c95c8c5a58c5c9e0c05956c0a97`。

实际按adaptation的两段new→old逆替换，重建正文SHA恰回到初审`4d6a22fc…530b1`，确认只有有限cover/full-level-set与fixed/local budget两处scope文字变化。作者`verify.py`和原`verification.json`的实际字节哈希仍分别为`a0016d33…545ac`、`cc43381d…9734`，未重跑历史toy。本审计自己的原`verification.json`也保持原字节，仍标识初审时的小型计算。新增`incremental_lock.json`单独记录此次版本锁定及新附记审核，不冒充原run。

**ACCEPT树附记的解析结论，范围如下。** 树来源`upper/round040/star_lp/PROOF.md`的实际hash为`6aece373926c7dfc65e7373b7245d117ba3fe4ae8f17fcd09492d1acc55a2a95`。节点半径8^(−k)，子中心偏移父半径/2，children严格包含于parent而跨分支stars分离。共同star列恰是祖先链，代表为最高节点。

在a>8/3、1≤β≤2时，实际small Q两两分离；祖先资格区包含后代中心，非祖先资格区不捕获后代中心。深k候选恰受k个祖先原a^(−n)项，因此D a^(−n)≤τ保证全选，actual v_i=a^(−n)V_i>0。任意h_i≥V_i只减小这些项，保留全选；若同时作为原witness预算，还需保持h_i≤M_i/λ。

单个真实f把总质量F均分至2^D个leaf内部。层状包含及跨分支分离给M_i=F2^(−k)、V_i=V_root8^(−nk)，所以F>λV_root足以严格见证每个节点。这里M_i来自同一f，而非自由指定价格。

unit-index共同star LP下界由非leaf s_i=2^(−k−1)、leaf s_i=2^(−D)给出：从深k代表向下任意链的s和至多2^(−k)=M_rep/F，而全树s和恰为1+D/2。上界用完整v-to-leaf链列，root代表总权重1、其余每节点总权重1/2；深j节点覆盖量2^(−j)+(1/2)Σ_(k=1)^j2^(−(j−k))=1，成本也恰1+D/2。因此S_f/F的精确值适用于全部fractional链列，非仅某个greedy规则。arrangement列可保留代表扩大链，故同值。

对于改收真实volume的共同star列，代表j仍在J中，故c_J=Σ_(i∈J)v_i/V_j≥a^(−n)。于是同一unit-index cover的cost至少a^(−n)(1+D/2)。这个下界与真实支付的有界几何级数相比保留深度overcount；它本身不否定所有参数范围内的uniform改善。固定a、τ时D受限，不能声称无限深；允许a随D增长才可使volume-free S_f/F无界。

containment-root列可收全部actual increments，因为它们都包含于root star。其cost为U/V_root=a^(−n)Σ_(k=0)^D(2·8^(−n))^k。反向对任意节点anchor j，M_j/V_j=(F/V_root)(8^n/2)^depth(j)≥F/V_root；所以任何unit-index覆盖按真实Σv_i/V_j收费的总f-cost都至少(F/V_root)U。root列达等号，严格证明containment接口的精确值。F↓λV_root仅为严格witness极限下所述有限selected支付逼近证书值，未产生operator下界。

此次仅另外以Fraction检查D=0,1,4的三条纯标量恒等式（链dual telescoping、全树dual总量、每节点fractional覆盖量）；没有重建树几何、调用新toy/生产LP或重算历史案例。任意D与任意n的结论来自上述解析证明。增量文件的加入使本审计目录共四个必要文件，未改变他人源文件。

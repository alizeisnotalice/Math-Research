**122：118双父锚点零H复用的独立整数审计（2026-10-26）**

全部有限门通过：6492个child、12984个单锚点整盒候选、207744个闭角点独立重算一致。新增1132个通过叶，分别为H64的1080与H128的52；其中996个成功者使用child外的H64锚点。剩余5360个child及原877个未消费叶保持待证，最终98007 passed、51884 split、6237 queue，所有52360根精确闭覆盖。整体g110仍未证明，本轮全局S记0，没有更新I、Cn或C5。

只写新audit_reuse122/，在116完整冻结之后另行登记，没有扩大116登记或重算其成功数值。旧116、108、100、090及既有250件锁保持原SHA。未提交、推送、使用MATLAB、改生产或派877/child H。本程序只加载只读081 common.py的JSON与整数几何工具，没有导入arith.py、任何kernel、Arb、log或作者函数；实际新H、梯度、kernel、log调用均为0。

118 manifest SHA `1d280a59665564b0f50b04f298f0fbd9f941d7711705b5ef673e72d7a163684d`，report SHA `8d14f20cd374e6f9226c81938ce7642ad9cf06d57f6bfc4c187aa625e382e789`。proof_sources.json SHA `aa179f0ca6244c5974c81f751d2e7d5f208c95f117821c4805b5a3235e451d54`，state SHA `10c3f96617c5e6573779b284925cff50438f63651b7f3ed2423bdfcd4854c737`。这些输入与全部源文件在计算前后按字节核对；实际路径及完整SHA见登记和manifest。

源107首批已经由116完成4096个新H128 jet的20480分量独立值审、fresh Delta/物理梯度比/全角点与partial树门，其manifest SHA `111840a3858a2a8b616d6fe731c16a28f8c04a18515dc9ac8fa39e879a484f2d`，结构结果SHA `33567c24b16532d08542fe74bac5b1e4e299d20cb85dc946550323d7e073d38a`。源H64新110目标的Delta与梯度资格已由108全量独审，manifest SHA `e4c6ee6f1a4e9b26bd8cc8d93d09203536a2d4b66c3533a942080a74fdfac813`。本轮消费这两条既有资格，不重新调用log或旧H。

对每个107失败父，逐字典绑定两个不同数学source。H128 source的原certificate和record SHA必须等于116独审回执；H64 source的完整原105新目标失败行SHA、certificate SHA与source index必须等于108的成功审计块。优先级文件新增的gap_rank在比对原证书行时仅作为排序元数据剔除。point、epsilon、typed pair key、H原球、Delta、四个物理log-gradient区间与1/(16ε²)罚项完全保留。125个父的H64/H128 point不同，也分别用各自point计算位移，没有混用球或梯度。

proof_sources中的作者qualification文字没有作为通过依据；当前H128来源资格由已冻结116独审回执实际补齐，H64由108回执实际锁定。6492个proof source恰为3246个失败父各两条；12984个候选是每个source分别在两个不同child上的整盒检查。全部6492个child及其source id、闭盒、父节点、路径，与107的3246次准确闭二分完全对应，没有重测父盒或删去子叶。

075/099的全空间不等式为

\[
\log H_\epsilon(c+u)\ge\log H_\epsilon(c)
+\nabla\log H_\epsilon(c)\cdot u-\frac{|u|^2}{16\epsilon^2}.
\]

它不要求c属于待证child。若已知 (0<L\le H_\epsilon(c))、\(\Delta\le\log(L/g)\) 及同epsilon梯度比区间，每个固定source的安全函数

\[
M(v)=\Delta+\sum_{j=1}^4\min(\gamma_{j,L}v_j,\gamma_{j,U}v_j)
-\frac{1}{16\epsilon^2}\sum_{j=1}^4v_j^2
\]

全局凹。只要同一个source在child全部16个闭角点的最小值严格大于0，就给该整盒 (H_\epsilon>g\)，进而给同一完整G下界。H下端只需严格正，不要求L>g；Delta允许负。这里是g证书资格，没有将anchor外盒规则迁到积分均值代码，也没有给混合H使用同epsilon齐次性。

check.py以root坐标与每四位交织path位直接生成每个child的有理闭盒。对每个source分别计算box−anchor位移，允许各坐标偏移两端同号；将Delta、梯度区间与全部位移化为公共整数分母，穷举16个整数残差，准确罚项为H64的256或H128的1024。独立最小有理margin、最坏bits、inside标志、HL>g标志和Delta负号与冻结候选结果逐条一致。没有导入作者Fraction角点或其separable verifier函数。

选择先检查H128整盒，再检查H64整盒，选第一个完整16角点通过的source。没有按不同角点切换锚点，没有对多个anchor的pointwise max套凹性。实际没有child同时由两source通过；H64候选1080通过、5412未证，H128候选52通过、6440未证。每个epsilon各有3244个盒外anchor候选；成功盒外anchor恰为H64的996个，H128为0。H64的204个负Delta候选得到正常有限检查，全部未过；算法没有因负Delta或anchor HL≤g提前排除它们。

1132个new_passed证书逐字典重建，准确复制所选source的Delta、梯度、罚项、record/certificate SHA及typed key，只新增child余量、最坏bits及075盒外规则标记。每叶只能由同一所选source的完整16角点支付。所有新叶epsilon、target、对象G及推论一致；原96875个已通过叶和全部189086个typed缓存对象保持不可变来源引用，旧历史H64障碍没有覆写或重贴标签。

未证5360个child的best_single_anchor_margin只是两条完整整盒最小值中的较大者，逐条核对其来源和排序；它只用于优先级，不用于跨角点混合证明。new_passed、remaining_children、queue与source索引一一重建。原877个未消费叶的node、冻结point、source leaf index、H64失败证书SHA及其他mapper字段逐项保留，只新增pending118_index。所有6237条pending mapper和精确gap顺序核对一致，没有新split、inflight或H请求。

最终把source105的96025个通过叶、107的850个通过叶、新1132个通过叶及6237个pending一起做完整整数树审计。所有52360根的前沿prefix-free，dyadic整数区间排序后精确铺满单位区间，Kraft严格为1；51884个split恰为前沿全部严格祖先，全部闭二分不变。活叶104244，来源和未证区域完整，没有用计数等式代替覆盖证明。

实际程序PID23029、exit0、5.504211500秒，单进程nice5，无pool或wall alarm。每256个child持久保存成功前缀，未出现fault、未重跑成功候选。audit_registration SHA `25e784b7a1758b03675f6184ea5067dc72fcf90d99af76186ff1dd48d48fdf8a`。result.json SHA `62380f0dfb98c405af7df299d9a0bb11f6b17869a15bf82f8222478e6bef5a95`。锁前的实际PID检查及250件旧锁结果见final_lock_check.json。

独立性范围明确为：新child物理盒、全部单source整数角点、source链/选择/metadata及整树独审；共享已锁081 common整数几何与075/099理论，H64数值/Delta依赖108及其先前值资格，H128数值/Delta依赖116。没有重新计算任何source数值或log，也没有把作者label当成来源资格。本轮提供局部g110叶证据，6237未证仍在，全局S=0；不能由局部叶数提高全域g110或直接改变I/C5验收门，主保留最终接纳责任。

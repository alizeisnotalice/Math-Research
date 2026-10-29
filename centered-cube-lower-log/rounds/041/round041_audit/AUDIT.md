# 第41轮后处理、统计预算与新upper独审

2026-10-29。只写本目录，40轮source/pilot/payload原件未修改。lower阶段756个CI、MCD与uniform-Lipschitz正文、同样本replay验收及FINAL全部1008个CI均已完成只读独审。本人未启动MC或replay；父线的重建是旧协议再次执行并逐位验收，无新增独立样本。

## 1. variance/3与保存贡献恢复

复核40轮已证接口：当前h=δ/6，class range帽恰为整bin的Lip×h，N_partial覆盖所有端点斜率；独立uniform差平方恒等式给单核variance帽÷3。正乘积展开及正系数Minkowski给整体Vnew≤Vold/3。对旧conditional lower B，φ(B)=3B/(1+2B)合法、单调，与共同L的max交换；upper不能盲变。

`postprocess_checks.py/json`先实际运行24个Fraction恢复/φ例、64个max交换例及60个FP极端恢复观察。x=λWB可通过B=x exp(logD−h)恢复，但x=0只保守留0；已丢失的underflow信息不能复原。等价稳定表达为 `new_x=3x/[1+2x exp(logD−h)]`。一般range、两点隐藏law均不能使用uniform-bin的÷3。

当前 `round041_post/postprocess.py`恢复的是conditional B，φ未作用于贡献x或旧CI。它重新计算旧lower、新lower、旧upper的样本CI；旧upper不改变。实际raw B均在[0,.966]内，本次clipping没有改变任何值，新lower逐样本≥旧lower且≤旧upper，原cap仍成立。原40原件及其504个CI均保持不变。

## 2. 1008预留预算与事后分析口径

当前252个冻结目标=4个R×3个模型×21个h。预留四种单侧函数：旧lower、新lower、旧upper、新upper，共1008项，每项α=.05/1008=1/20160。当前只计算756项，union费用≤.0375，余.0125留252个合法新upper；不应把分母临时改为756。

若new upper是对同一个真实条件tail的refined-Cantelli、MCD、Markov等**已证逐点upper的min**，再按同一冻结共同L聚合，它是一个合法upper函数。对只有各自边际coverage、没有共同pointwise支配的一般统计CI，不能任意取min免除选择费用；本轮满足下述方法统一条件的包络则另论。预留252这一数字不能补救换目标、换阈值或未证代理。

必须保留的事实：40结果及一次c=1/3的描述性预览先于本轮PLAN，41是事后确定处理，不是盲测或新独立确认。c完全由源的Lipschitz解析证明固定，R、模型和h族未按41收益改变。PLAN已明确披露这些事项。

在以下充分条件下，可以给出固定族的**无条件**名义coverage陈述：各R样本行按固定提议iid；同一个真实条件tail贡献在所有样本上被对应lower/upper夹住；共同source-only cap也上界真实贡献；所用算术精确或误差已认证；252个目标本身固定。然后有限族Chernoff–KL反演和union给总体失败概率≤.05，跨模型、阈值及旧/新界可相关。

披露事后提出的时间并不自动推翻这个特定固定族的数学无条件陈述；但不应声称“在本研究为何被选中、为何继续报告”的选择事件条件下仍有95%覆盖，也不为其他轮的模型/目标搜索提供总错误率保证。当前完整FP未认证，因此报告宜继续称数值诊断及名义CI，不当作有限真概率/Cn证书。

方法统一的正式证明见 `METHOD_UNIFORM_CI.md`：设固定真实贡献为Z_i∈[0,C]，任意合法lower l_i≤Z_i、upper u_i≥Z_i。KL两端点均随样本均值递增，所以lower端点(l)≤lower端点(Z)，upper端点(u)≥upper端点(Z)。同目标上即使方法或每行bound依赖整批数据，只要所有实现均逐点合法、同cap，仍同时被同一个真实Z的target-direction事件罩住，不需按方法数union。它覆盖EZ，不声称selected-bound自身期望的iid CI。252个固定target的504个direction事件足够形成这个数学陈述；本轮仍保守使用1008预算，不改冻结端点。1152个单调邻点、512个全样本自适应夹逼检查及data-tuned cap/target反例均先实际计算，再给自含证明。

## 3. Lower阶段保存后处理：756 CI复核通过（历史阶段）

已只读核对 `round041_post/postprocess.py`、PLAN/RESULTS及其全部40输入hash。程序实际约.026秒，没有新样本、source重算或seed变化。

自有 `post_result_checks.py/json`以60位Decimal独立计算exp、贡献恢复和直接φ表达式，并以100次高精度bisection重算全部756个当前CI（α精确1/20160）。最大CI差1.82×10^-12，均值差7.11×10^-15，paired均值差8.89×10^-16。所有已保存输入hash匹配；raw B无需裁剪，原零贡献保留为零。

84个同阈值sharp-joint-lower减old-prior-upper经验差，以及对应的同时CI端点差，均未正分离：

|R|最大经验差（h）|最大端点差（h）|
|---:|---:|---:|
|4|−.08477411（0）|−.62503869（0）|
|8|−.07400318（0）|−.61453368（0）|
|16|−.04113418（10）|−.61379493（0）|
|32|−.05415056（0）|−.59506018（0）|

这仅是没有正比较witness，负的lower−upper不证明真实joint没有增益或劣于prior。CI覆盖冻结h族；不能作为所有阈值或连续尺度weak上包。理论锚1仍明确来自额外L=1、λ↑1。

## 4. MCD/replay审查前提与当前状态

已读 `round041_mcd/probe.py/RESULTS.json`，并先在自有 `mcd_checks.py/json`独立运行n3、8个hidden状态的正多项式：48个Fraction坐标替换比、10个数值log-MGF及4个sum-tail检查通过。此处允许独立两点law，未将它与uniform-bin的variance/3混用。

现已读完整MCD正文及observer代码并通过预审，详见 `MCD_REVIEW.md`：正多项式ratio、Doob/Hoeffding常数、Jensen方向、u_b/Σu阈值分配、S0与空项均正确；固定C=4.05覆盖精确C0。将S写为C log(1+R_old)时确实使用**旧quarter/range relative cap**，未把r_old/3误塞入原C。source r0预算在原真实均值误差前提下成立，仍非完整FP认证。

replay只能恢复同一40样本。必须逐R核对原int64 count stream hash、sampled nominal anchor hash及原四特征统计stream hash；seed相同本身不够。保持n、N、R、grid/cuts、rng调用顺序、source/feature算术及组件遍历顺序；新upper所需额外统计不能替换原特征后再声称原stat hash匹配。若hash不匹配，应停止将其称作同数据重放并记录差异，而非追加seed或把重放当新独立确认。

实际重建已在230.1765秒后accepted，峰值RSS561938432 bytes。自有 `replay_checks.py/json`重新逐R验证三stream hash，24个旧payload数组全部逐位相同，side NPZ及计划/证明/代码hash链匹配。`sharp lower≤new upper≤variance-only upper≤old upper`及原cap的保存数据链也通过。这里只验收旧数据身份和保存bound-chain；完整逐尺度upper公式的独立数值复算由另一独审负责。

## 5. Uniform-Lipschitz factor6及FINAL1008复核

`UNIFORM_LIP_REVIEW.md`审查通过：1-Lip函数的非降quantile仍1-Lip，中心convex-order支配uniform，精确MGF≤sinh(t/2)/(t/2)≤exp(t²/24)；Doob未来积分保留Lip，给factor6尾。当前source的Lip(log F_b)h由q/(1−q0)控制，仍被相同C×Σlog(1+r_old)覆盖。这里确实使用Lip surrogate，不能由仅oscillation的旧McD直接换6。

独立穷举81个四段斜率−1/0/1函数，405个Fraction stoploss、567个分片指数积分数值检查通过。连续uniform前提与target隐藏law相符，不适用于旧两点law。

已读 `round041_post/aggregate.py`、`uniform_upper.py`、FINAL_PLAN/FINAL_RESULTS。primary upper采用factor6及variance/3-Cantelli/Markov的逐尺度min，再作233-scale union；factor2与variance-only仅描述means，没有新增CI族。exponent低于−700时向上floor为−700，避免仅因exp下溢把正upper写成0；这支付该局部方向，不代表完整FP认证。

`final_ci_checks.py/json`以60位Decimal、精确α=1/20160、100次bisection核对全部1008端点；最大差2.17×10^-12，均值差7.11×10^-15，paired均值差8.89×10^-16。CI反演、exp权重及sharp-lower恢复独立计算；upper核心调用已审冻结formula，逐尺度完整独立复算另有专审。FINAL hash链全匹配，**最终检查不依赖phase1 RESULTS.json**，只使用FINAL、40原payload、已验收per-scale输入及必要PLAN。

最终84个same-h joint-sharp-lower减prior-sharp-upper的经验差及对应同时端点仍全负。各R最大端点差均在h0，约−.605285、−.605776、−.604268、−.585493。R8高h10经验差接近0（−2.82×10^-21）但端点差仍约−141.59，不能解读为尾已精确确认。

事后提出、已见40描述性预览、同样本重放/确定处理及非独立确认声明均保留。METHOD统一lemma解决同fixed target合法包络的方法选择计费，不认证尚未证明的pointwise FP夹逼，也不覆盖data-tuned cap或新的source/model/h。最终无正分离witness，不证明真实joint无增益；完整FP、全阈值/连续尺度weak及Cn目标均未闭合。

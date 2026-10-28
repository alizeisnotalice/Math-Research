# 第37轮分位数来源独审与实际小模型

结论：真实分位数来源的质量、共同尺度PB公式、sampler与全pool边际RN通过审查。实际小模型证明了与mod分类实质不同的现象：在同一个L>1及正体积query区域上，joint/prior点态倍率可达R阶；但固定L、全尺度M和weak范数必须分开。该构型没有给出weak增益，实际MC也未确认weak排序。

本人先运行自己的Fraction/真实box模型，再写THEORY.md；未重复代理尾MC，未修改他人来源、冻结文件或Git，未启动代理。全部自有文件仅在本目录。

## 1. 独立实际计算及最强结构发现

quantile_probe.py以K4、n63/255/1023、R2..32作12组exact Binomial/box检验。分位段按整数CDF首次达到b/R，pi使用准确格点质量。每表每坐标的box积分在三个有理偏移上准确核对，PB top段用准确Bin(n,3/5)有理尾，不用Gaussian。

取delta=1/(64KR)、L*=1+2delta，query在top表plus箱中心减半周期处。此时top表posterior plus概率3/5，其余表为1/2，likelihood比a=(5/4)^n。准确joint/prior比为(a u+1−pi)/(a pi+1−pi)。n1023、R4/8/16/32的实际结果约4.06856/8.51839/16.95670/32.99566。THEORY.md进一步用Chebyshev和准确binomial递推证明n=2R^4+1族的R阶点态增长，而不止有限数值。

该增长在一个准确描述的plateau Q上成立。额外窗口全长2delta、源箱全宽delta，完整包含的center区间宽delta，任意相交区间才宽3delta。K4有两个plus原子，所以每坐标union长度2delta，总体积(2delta)^n=(1/(128R))^n；没有将单个区间误写成2delta。

固定L*在Q上的尾贡献最多是THEORY(6)的极小值，n1023/R32的log帽约−8280。**这个帽不能用于full M。** Q同时是table0、全minus source cells的箱支持，已用Fraction准确核对半周期映射。故full M在Q内准确等于密度峰H，限制到Q的weak贡献为H|Q|=2^-n，log约−709.09。它仍只描述Q，不是整个来源weak上界。full Mjoint/Mprior在Q内介于1/pi_max与1/pi_0之间；n1023/R32的下端约20.0971，固定L的32.99566不能替代这两个max的比值。

此外joint源密度自身为高度H、体积1/H的平顶，weak系数准确1；prior密度weak为max_k k pi_(k)≤1，L1又已给Mprior weak≥1。峰高与支持体积存在真实归一化抵消。需另研究Q之外的整体共同尺度尾，不能从点态R倍率跳到weak或log增益。

## 2. 代理来源pi与共同尺度PB实现

只读round037_quantile_cell/PLAN.md、完整probe.py、registration/calibration及RESULTS。quantiles使用Python整数comb/CDF，重复cut保留为明确空段，未把pi伪装成1/R。n15/R8与n31/R8各有7个active段，固定的8张表仍保留；joint/prior抽样不会选择空段，fixed对照仍是单位质量来源。

每joint粗格先均匀抽取、据plus计数判表，再放该表箱；theta与1−theta及delta保证箱严格留格。于是每cell质量K^-n、joint总质量1，表标签真实先验恰为记录的Binomial段质量。prior则先独立按该pi抽表，源坐标格仍均匀；两种输入确实不同。

PB实现先形成每query/table的单坐标总权重，A_b的乘积用log保存；对active pair以plus/total和minus/total正递推形成PB概率，直接加该表的整数段系数，最后以logsumexp加A_b×概率。没有用两个接近1的CDF相减，未给joint再乘pi，prior才按实际pi加权；zero total与空段给严格零分量。同一调用全表/全坐标共用L。保存的K4/n3全cell直接枚举及Fraction控制覆盖空段情形，与DP一致。

真实后验/Bayes表达及精确积分守恒见THEORY(1)–(3)：ell_b=PB_b(I_b)/pi_b，joint/prior是按pi_b A_b归一化加权的ell_b，且每个固定table-anchor律下E ell_b=1。可出现1/pi_b级倍率，但必须付query律下的概率；该恒等式不提供maximal weak上界。

## 3. sampler、实际RN与数值稳定性

draw_proposal先选uniform或全部来源/全部L；joint先cells后classify，fixed固定b，prior用独立n个fair bits再classify，准确生成格点pi，避免用Gaussian或浮点categorical替代。source box噪声和cube位移逐坐标独立，同一行共享L。prior抽中的b与本次source cells独立，正确。

evaluate_queries在log域累计全部R+2来源×J62尺度的真实平均，RN为D=.5+.5 sum/((R+2)J)，不是抽中标签项或max。W≤2，lambda-tail贡献帽min(2lambda,2(R+2)J)有效；L1保证log max≥0，所有模型都是同一有限尺度集合。

PB总质量残差最多约8.22e−15，仅为一致性诊断。subnormal-loss预算约e−274至e−133，明确只支付理想概率递推的subnormal量级，不涵盖其他roundoff、logsumexp丢项、精确pi转float、实际RN误差或尾指标边界。不能把这一小数值当完整浮点证书。有限box正和的解析合法性也不自动认证所有FP输入。

## 4. 保存统计及范围

自有audit_saved.py仅读保存摘要，以独立整数算法重算9组quantile cut/pi/空段，再反演全部holdout CI。probe、PLAN、calibration、registration SHA全部匹配。每n的模型数为4+6+10=20，三个n总60目标、双方向120 CI，alpha_each=.05/120；端点独立复算最大差0，alpha总和.05。

各caseN2048、前后1024分半。discovery用ESS≥32选log阈值，再nextafter向下，后半只验证该阈值；适用于该单尾系数。全样本unrestricted与ESS64扫描为描述，未给完整weak上带，因此不能排序两模型的最优weak。不同来源的holdout通常选择不同阈值，更不能以两个CI区间当全weak排序。

稀疏候选必须保留：n15/R2未过滤joint约1.38517但ESS15.0；n15/R8约1.28486但ESS4.89；n31/R8约1.44600但ESS6.35，ESS64后分别退到约1.10494、1.03650、1.02541。实际joint九个holdout下端都小于1，未新增超过L1确定锚的已确认下界，也未确认相对prior/fixed的weak改善。

独立macro诊断在固定64个query、L10/9、7/6、11/9上的残差未按mod分类指数消失；例如n63/R8最大约.1661，与两机制不同相容。它仍只是有限query探路，不能替代理论全x/L比较或weak结论。

所有工作仍限于有限来源与有限MC尺度；THEORY的明确plateau结论单独注明何时涉及全部L。一般Cn、匹配原来源、增长阶和完整FP认证均未完成。

## 5. ESS下溢修订的最终独审

后续只读ESS_FIX.patch、ESS_revision.json、ESS_REVISION.md、replay_stable_ess.py及RESULTS_stable_ESS.json。最小补丁仅替换weak_scan的累计权重/平方/score/ESS计算，改用logaddexp累计logW及2logW；来源、RNG调用、query grid和RN没有改变。1e−13的log-ESS筛选容差只作用于描述性ESS门槛，不影响holdout有效性；它不是额外样本或独立确认。

首次没有保存raw logM/logD，所以研究代理实际将主数据生成执行了两遍：首次26.1449329853秒，按原9seed、原版本和batch协议重放25.8362388611秒，合计51.9811718464秒。它们是同一协议的纠错重放，不能写成只作了一遍计算，也不能当两份独立验证数据。宏观诊断直接复用首次输出，未重跑。

自有ess_audit_checks.py只读这次重放保存的9个NPZ：文件hash及逐batch logM/logD流hash全通过；两版metadata、质量权重统计、source-b计数、贡献上限、PB诊断和subnormal记录逐位相等。首次raw未保存，故这些一致性不能被夸大成首次与重放完整逐query hash比较。NPZ的原始prefix公式在n63/R8确有joint13、prior57个非有限ESS前缀，各fixed表1–12个；其他8case没有该问题。确定性小例logW=[−400,−401,−402]也复现平方下溢，log域ESS与独立scaled-fsum结果约1.95869865一致。

数值主表应以RESULTS_stable_ESS.json为准，原文件保留历史，不静默覆盖。自有audit_saved.py对修订主表重新核对9组准确pi及120个CI端点，差0，家族费用仍.05。实际所有holdout字段与首次逐位相同；经验赢家系数仅舍入级差，没有新胜出候选。n63/R8稳定joint最佳约1.0014778442、ESS1750.0673，prior约1.0002581997、ESS1745.1303。旧版“忽略非有限前缀”不可视作数值方法，但修复后结论未变。

父线另从保存NPZ以独立scaled math.fsum重算120个经验赢家并完整聚合ties，系数最大差约2.25e−14、ESS最大差约6.49e−11，未见遗漏ties；这是零新MC的后处理核验。本独审没有重复生成样本。普通FP积分、指标边界、实际RN等其余认证缺口仍保留，稳定ESS不等于整场计算已认证。

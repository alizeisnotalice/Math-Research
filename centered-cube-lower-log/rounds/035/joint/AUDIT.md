# 第35轮joint-cell来源短独审

已只读PLAN.md、完整probe.py、registration.json、calibration.json及RESULTS，sampler、共同尺度DP、全pool边际RN和完整finite-grid weak带的主要接口通过。六点均未确认joint的最优weak值超过固定来源或实际先验混合。本人未重复大模拟、修改代理脚本或冻结文件，未启动其他代理；自有audit_saved_checks.py只复算保存摘要。

## 来源质量与共同尺度DP

每joint cell j向量分配K^-n，选表b=majority(chi(j))后在该cell内放一份边长delta的均匀box。两表的原子到原格边界均至少1/(2K²)，而delta/2=1/(128K²)，因此box严格留在原joint cell，来源总质量准确1。fixed0/fixed1也各质量1，prior_mix使用实际先验。平衡chi和奇n使majority两侧cell数准确相等，π0=π1=1/2；这不是声称给定b后各source格仍独立均匀。

coordinate_weights逐坐标计算真实周期box落入同一个长度L的窗口之质量，再除以KL。x在圆周、L≤1.25时三份周期平移−1/0/1足够。L1的快捷返回全1准确。正权重多项式prod_i(z^-_ib+t z^+_ib)的系数计数source cell的plus个数；b0取少数、b1取多数，再直接相加，即joint真实cube平均。此处若再乘π或除π都会改来源；当前实现没有该错误。fixed表用全部系数和，prior_mix为两者各半，角色分开。每次DP中所有坐标和分支共用同一个L，未误把逐坐标/分支最优尺度放进乘积。

保存的小控制K4/n3枚举64个cell、majority各32个，12个查询的DP与直接枚举最大相对差约1.91e−16。保存Fraction在L3/8、5/8、1下准确相等，分别joint=40/27、128/125、1；prior_mix=8/9、24/25、1。本人用本机Python Fraction重读这些保存值，核对joint≤fixed0+fixed1均通过，未重新运行代理控制。小例不认证所有浮点查询，但与上述代数一致。

## 提议与weak比较接口

计划的真实提议若按如下顺序实施则合法：先选uniform，或选一个来源及共同L；joint先均匀独立选source cells，再据这些cells算b，随后在相应joint box中独立抽平滑噪声。prior_mix才先按π选b，再独立抽各cell。令T为真实来源样本，各坐标独立抽U_i∈[−L/2,L/2]，X=T−U mod1，即抽中该来源的A_LF密度。共同L不等于所有坐标共用同一个位移U。

边际RN必须为D=.5+.5/(4J)sum_(model,L) A_LF_model(X)，包括全四来源与全部J尺度，不能只用实际抽中latent项或只取max。由D≥.5得W≤2；当target最大平均>lambda时，匹配项给lambda W≤8J，故指标贡献帽min(2lambda,8J)有效。该帽依赖所有目标尺度都在proposal中。

所有来源因L1而有有限grid max≥1，密度峰值均≤(64K)^n，故完整weak筛选可把lambda≤1部分确定上界设1、高于峰值部分设0；内部阈值段用lambda_(k+1)P(M>lambda_k)上包，方向正确。严格大于号在峰值处尾为0；weak下界至少1可由lambda向1的左极限得到，不要将P(M>1)误设为1。逐阈值、模型、case及方向需联合分配CI。固定查询点的平均值增大或两条点估计弱曲线交叉均不构成完整weak排名，更不构成log或Cn阶。

## 两表机制的范围

在每joint cell只选择两表中的一份质量，正性立即给F_joint≤F_fixed0+F_fixed1=2F_prior。cube平均及共同尺度最大值也保留此支配，所以W(M_joint)≤2W(M_prior)。保存的有理查询与该支配相容。它不排除有限常数弱收益，但只靠这个二表选择不能相对该prior产生随n发散的log因子；不能将点态joint/prior比解释为weak比。

浮点box积分、DP、权重和置信求值尚未区间认证；直接有限和的解析精确性不自动支付普通浮点误差。

## 完整实现与保存统计复核

draw_proposal确实先均匀抽source cell向量，joint才由该cell算majority b；fixed0/1分别固定b，prior_mix另抽独立公平b。每坐标的平滑box噪声、窗口位移均独立；同一行共享的只是L和来源类别，并未误用单个标量U平移全部坐标。样本T属于对应真实来源，X=(T−U_L)mod1的密度就是其真实A_LF。

evaluate_queries将所有L上四来源的DP均值求和，形成D=.5+.5 sum/(4J)，同时独立保存各来源的共同L最大值。分母不是max或抽中标签项。W≤2与W·max≤8J均有assert，指标贡献帽正确。L1快捷值保证max≥1；完整weak下CI显式取max(1,各阈值下CI)，没有强设严格尾P(max>1)=1。

成本定标先对固定确定query点计时，按冻结120/180秒规则选N；六点均N16384，seed202610290351至356各一次，未按收益增样本。注册冻结了全部L和阈值。J为K4的68、K6的70；K6保留了两个浮点上极接近2/3的不同长度，属于冗余但合法的两个query，proposal和帽均使用实际J，没有少算。

自有audit_saved_checks.py仅读保存摘要，用标准库独立反演全部CI和整段包络；probe、PLAN、calibration及registration的SHA均通过。六case×四来源×双方向共有11112个CI，alpha_each乘此数准确为.05，家族分配覆盖数据选择的所有最大下端及最坏上段。CI复算与保存值最大绝对差0，整条weak上下端复算精确一致。上包使用next lambda（最后段裁到密度峰值）乘前一阈值尾概率上界；峰值及以上的严格尾置0，方向正确。

| n,K | joint完整weak带 | prior_mix完整weak带 |
|---|---|---|
| 3,4 | [1.13324,2.19156] | [1.12333,2.39657] |
| 3,6 | [1.10298,2.40232] | [1.10955,2.29054] |
| 5,4 | [1.13651,2.27184] | [1.12705,2.27004] |
| 5,6 | [1.10336,2.39994] | [1.13087,2.35790] |
| 7,4 | [1.13755,2.24469] | [1.13606,2.18566] |
| 7,6 | [1.10924,2.49510] | [1.13412,2.33527] |

六点joint下端均未超过任一对应baseline的完整上端，不能确认joint weak改善；也不能据点估计宣告其真实weak更差。固定网格最佳点的非负尾ESS在joint约9352–13463，并非本轮结论受限于个位数尾样本；主要是完整阈值上围栏仍宽。自选empirical_weak及ESS过滤结果属于同数据的描述量，不带额外置信；遇重复max值时所报阈值对应从左逼近的经验弱函数值，不应解读成该精确阈值的严格尾概率。

上述带只围住列出的有限L集合，不能以上界部分控制所有连续尺度。所有普通浮点积分/DP、实际采样分布与计算RN的一致性、阈值分类及CI舍入仍未区间认证，hash和摘要复算不补这层误差。本轮通过的是数学与实现方向审查及保存统计一致性，没有新的Cn或log阶结论。

# 第40轮 streaming driver 独审

结论：本次预pilot只读静态审查通过，未发现阻塞执行的实质错误。均值/BE上下包的数学前提沿用38轮已审链；普通FP、BLAS、erfc和随机数实现仍未向外认证。没有运行新MC、没有构建cache、没有改stream/source冻结文件。本审不声称已有第40轮结果或统计收益。

审查对象为 `../round040_stream/probe.py` 和 `helpers.py` 的落盘版本；范围包括source量词、make_cuts、stats轴、RN、Minkowski/Cantelli、两侧帽、504项CI与paired口径。版本SHA256分别为 `9ec993e2e25e4415eb5341f90b9dec604efc0cc63deca547c254ecbae94e52b4` 与 `73c8ae324a2e9dffddc382314e73ec31a852b39b32389323491dd7d3cb98d1d3`。

## 来源与真实class law

R4/8/16/32各表只改变cell内峰位置，common D与每格峰质量相同。R32的box允许触格边但不越格，端点零质量不改变共同cell law，故所有表使用同一η≈.498047800888639。gated表由固定源class-count段决定，各产品律下源S都为Bin(n,η)，共同分段使总质量1。prior使用该实际Bin段概率，而不是1/R。

`make_cuts` 采用rounded Gaussian-quantile候选整数，然后将这些整数作为定义冻结；结果字符串明确不声称精确Bin分位点。CDF真实上下包用actual η及其guard。pi_l的首项是CDF lower−0，中间项是右CDF lower−左CDF upper，尾项1−CDF upper；pi_u方向相反，均正确，且没有重新归一化破坏包含。当前n和R≤32的cut严格递增；代码若不满足会报错，不会静默删除标签。

## Stats轴、query conditioning与PB段包

counts每行是各query bin的整数数目。来源、固定产品和prior均对坐标置换对称，所以给定counts的目标分布与任意代表ordered bin排列一致；该代表排列下隐藏query位置独立。此处没有宣称histogram条件下原坐标本身独立。

`feature` 按列提供log total、p、p(1−p)、log1p(class相对variance帽)。counts@columns给N×4，转置赋值stats[:,:,b,j]恰是shape(4,N,R,J)的对应轴。z、PB mean、PB variance和log product-variance帽的索引均正确。

`bounds` 的低cut为[-1,c0,…]、高cut为[c0,…,n]，以(1,R,1)广播，与各table自己的(N,R,J)PB law逐元素配对；没有混用别表的PB CDF。首低端强制0、末高端强制1，真实支持精确处理。段概率lower=右lower−左upper、upper=右upper−左lower，方向正确。

joint按表乘对应段posterior，不再乘π。prior才乘π lower/upper。每表r_b由ΣN log1p(r_qb)构造；跨表variance采用 `(Σ_b u_b sqrt(r_b))²`，不是variance相加或显式R²。fixed0也可用同一class-max表内帽，因为total是两个正class核之和。

## Tilt与marginal RN

抽中的anchor和后续全表stream都调用同一个 `feature(src,length,b)`。nominal向量明确为exp(log(total))，抽样概率nominal/sum nominal，normalizer为mean nominal。后续component RN为 `exp(z−nlog norm)`，与该nominal概率law的理想定义一致；不是拿source原ρ当Y倾斜，也不漏实际归一化。

每个sample只选择一次uniform或(b,j)标签抽counts，之后RN使用全部R×J anchor的混合和。stats[0]为(N,R,J)，reshape(N,−1)正好枚举全部b,j；`logD=log(.5+.5/(RJ)Σ A_nom/norm^n)`。没有用抽中单一component的RN代替marginal RN，W≤2正确。

counts在int64状态下hash后再转float64。每项0≤N_q≤n=6704800<2^53，整数到float64精确；不涉及大计数丢位。np.multinomial的数学总数固定n。stream统计hash及保存的logD/component_sum/payload有助后续复核，但不是完整FP或原counts证书。

counts层importance抽样无需把nominal倾斜看成真实连续query tilt；它只需是明确正的counts proposal，RN按所用nominal分布计算，真实target则由条件矩/误差包单独承载。

## Source-only帽与条件尾

error inflation与PB参数误差规则写在冻结源码中。C_Z、C_E在stream后才汇总，但只依source/evaluate及固定grid，不依counts、标签、贡献或观测max，因此不是自适应统计cap。C_Z≥max norm^n，C_E≥各A真实/nominal inflation。

各模型每尺度均值上包≤C_EΣ_b A_nom,bj；对prior也成立，因为各π upper≤1。D≥Σ_bj A_nom/(2RJ C_Z)。lower采用active共同尺度时λB≤目标均值；upper采用union及Markov时λU≤Σ_j mean_upper,j。因此两个方向都满足

`0≤λWB,λWU≤min(2λ,2RJ C_Z C_E)`。

这一步upper不依赖仅lower成立的active条件，也没有从观测max设置cap。`tail_enclosures` 的lower/upper Cantelli方向与38轮一致，lower取L上的max、upper取L上的union sum截1；同一个L在全部表形成平均以后才选，不假设尺度独立。

## 统计family与pairs

四个R、三模型、21个预冻阈值、两种单侧CI，共504项，每项α=.05/504。N256及四seed在plan前定好，plan保存全部cuts、π区间及source-only cap规则后，才进入各R的RNG。CIs用预定有界贡献的Bernoulli-KL反演；不同模型/阈值共享样本不影响union bound。

pairs按相同R、相同h、同一行counts构造joint lower−prior upper。保存paired SE只作描述，不提供未经预算的paired CI；`simultaneous_CI_endpoint_gap`用joint的已有lower-CI减prior已有upper-CI，在504项覆盖事件下仍是合法固定阈值差的下界，不需额外α。它不能换成两条lower相减或完整weak范数差。

保存的lower/upper数组、logD、component_log_sum_normalized和tags足够后续只读重算摘要/paired统计；没有保存巨大的feature cache。full自然R32 cache会超过11GB，而该driver只持counts与4×N×R×J充分统计，结构上符合streaming方案。

`best_lower_CI_with_L1_anchor`明确带`anchor_scale_is_extra_L1=True`；1锚属于额外L=1，不属于裸I=[10/9,11/9]的233点grid。实际tail阈值λ≥1时L1不贡献严格超水平事件，这与旧窗口profile分开。

## 边界

本审无须修正的公式或量词。执行后仍须核对source/driver/helper hash、504项CI实际重算、实际帽、质量/分支计数及固定阈值比较。未审尚未存在的数值结果；未给continuous-L上界、fullweak排序或C_n/log改善结论。普通FP前提继续保留。

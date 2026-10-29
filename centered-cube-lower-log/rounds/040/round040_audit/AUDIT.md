# 第40轮 general-R来源与streaming pilot独审

2026-10-29。预运行版审查通过后父线完成一次冻结pilot，本报告已追加实际四组RESULTS/NPZ的只读核查。没有运行或重跑父线tail pilot、未改其源码、未抽新MC。

## 来源资格：通过

只读审查 `../round040_source/source_tables.py`、`build_cache.py`、`validate_source.py`、`SOURCE_VALIDATION.json`。Source扩到R4/8/16/32，n6704800、K255、v8、w4、B4/M3及common-D保持固定。R是表数，**不是Sidon标签字母B**；field_labels只认证原B4标签资格，不认证原匹配概率链或新B32来源。

plus offset `(label/4+(b+.5)/R) mod1`、minus offset `(label/4+(.5−b)/R) mod1`由奇整数numerator/(2R)实现，modulo及符号正确。自有 `source_checks.py/json`实际做15300个Fraction offset/留格测试，通过全部理想有理盒质量检查。R32有510个触格边端点，非严格整数≥判据正确；uniform盒端点无质量，不导致跨class漏质量。

各表每格峰质量p/K，共同D背景不变，故完整class质量共同为

`η=127/255+(1/1020)Σ_{m≤3}c_m sin(2πm·127/255)`，

而非1/2。共同背景严格下界437/918≈.4760348584>0；D有限余弦积分0，峰质量加背景质量为1。κ正性及common η是解析来源性质，浮点η的具体包围仍需FP认证。

同一class内所有offset仍处于平移的四标签quarter格，因此最小内部相邻atom gap≥1/(4K)=1/1020，各R的精确检查均达到该值。它大于δ+h，保证class峰均值算法每端最多一个partial候选；不能仅凭generic R32半格宽推出错误的双候选结论。跨class两族分别处理。

R4的atoms/labels与38来源相同；生产验证还在16个L/table全Q evaluate上逐位一致。R8/16/32的midpoint phase网格与较小R不嵌套，不能把本试验解释成同一源逐次加表的单调实验。新的R4 normal-quantile cuts也不同于38轮±1σ cuts，因此不是38轮完整复现。

## Kernel均值与误差接口：通过数学方向

trapezoid-CDF均值、full/partial开闭端点、common-D class-window分片中点积分与38轮相同。R32触格边在零测端点没有新增质量。来源验证保存独立point-kernel Gauss参考与有理peak检查；这些只认证所测有理峰部分，Gauss与完整FP不是区间证据。

每class名义均值a,b各误差≤e且a,b>e时，p的极端在(a+e,b−e)及(a−e,b+e)，所以 `|p_true−p_nom|≤e/(a+b)`。PB均值和variance各≤n e/minS；总kernel log乘积双方向用 `n[-log(1−2e/minS)]`，不缺额外2倍。这一传播以e真包含误差为前提；经验guard10^-13未认证，之后BLAS也未自动被e覆盖。

## 冻结streaming脚本预运行审查：通过

对象 `../round040_stream/probe.py` 与 `helpers.py`，审查哈希保存于自有 `pilot_checks.json`。最初17共同尺度意图已明确升级为全233 natural∪old网格；当前源码用 `exact_grid`，各R共用同一精确长度列表，含I两端点。不是连续尺度sup上包。

每个R固定N256，seed=202610290400+R；当前模型fixed0/joint/prior，h0至10步.5。来源、cuts、N、阈值、cap规则在PLAN写入后才抽样，source-only normalizer/cap虽流式算在后面，但不依赖观察尾结果。

先按每行一次选择uniform分量，或均匀选(b,j)后从该正nominal bin概率抽n坐标计数，即multinomial(n,q_bj)。统计只存N×Q counts和4×N×R×J矩，不建4×R×J×Q全cache。RN对完整n样本的混合密度为

`D=.5+.5/(RJ)Σ_bj M_bj/Z_bj^n`，

代码使用全部R×J的z减n lognorm后logsumexp，没有误用逐坐标混合再取n次方。norm严格按同一nominal=exp(logtotal)向量定义；理想实数算法里multinomial系数在比值中抵消。实际概率归一化和BLAS/log/exp舍入尚未认证。

自有 `pilot_checks.py/json`对n3、Q3、R2、J2做独立Fraction count枚举：proposal总质量1、E_proposal(1/D)=1，RN与直接多项抽样概率比逐项相等；包含非1的Z，60个双方向cap小检查通过。没有运行高维tail。

## 共同L、Minkowski与双方向贡献帽

各b先在同一个j计算PB段概率，再对表正求和，得到每j真实条件均值包与分层Minkowski方差帽。此后lower在j取max、upper在j作union，维持一个共同cube尺度；未对每坐标或每表分别取最优L。

joint variance使用 `(Σ_b upper_mean_b sqrt(r_b))²`；prior也用实际π区间的上端均值组合，方向正确。fixed0沿同一个表的kernel乘积。π上包不必和为1，不能再强行归一；只需每项≤1即可支付Markov贡献帽。

记C_E为最大kernel均值膨胀，C_Z=max(1,max Z_bj^n)。任一目标局部均值上包≤C_EΣ_b M_bj，upper-tail union还被Σ_j mean_upper/λ控制。因为ΣM≤2RJ C_Z D，故

`λ W U ≤2RJ C_Z C_E`，同时W=1/D≤2给≤2λ。合法lower也≤真实条件tail≤U，因此两方向均可用

`cap=min(2λ,2RJ C_Z C_E)`。

实际assert容忍及FP未向外认证不能转成严格数值证书；上述是理想算法的数学cap推导。

## 实际π、CDF与统计族

normal-quantile公式仅选定整数cuts，实际prior仍按共同真实η的Bin(n,η)做BE包围。空段逻辑支持，但当前n/R下cuts全不同。coarse quantile不声称精确Bin分位数。

Source的宏观范围和η∈[.49,.51]保证所有真实p∈[.3,.7]；C=.5606的有理检查给CDF BE≤1/√n。CDF的μ/variance矩形角点取min/max正确。low cuts [-1]+cuts与high cuts cuts+[n]逐表对应；首尾直接置0/1，内部上下包作相反CDF端点相减，没有重复count或漏掉第二端BE。

自有动态R mock核对，在四个R下PB参数与prior一致时，joint/prior的全部mean/variance区间逐位一致；prior lower/upper sums允许分别小于/大于1，R32最小π下端仍约.030318。

统计族504=2方向×4个R×3模型×21个h，α_each=.05/504。所有固定h/模型/R的端点同时成立后，选择经验赢家或报告joint lower CI减prior upper CI无需另付多重选择；paired SE仅作描述，未作额外显著性声明。504个区间**不覆盖所有阈值的fullweak上包**。理论锚1明确来自额外L=1及λ↑1，允许完整cube算子下界使用，不属于严格233点有限算子的无说明附加值。

## 实际四组保存结果：通过一致性检查，无正分离

父线实际四组R4/8/16/32、每组N256、每组233共同尺度已完整落盘。RUN_RECEIPT总耗时243.401秒，进程峰值RSS558514176 bytes；没有4×R×J×Q全feature cache。与最早17尺度意图不同，最终冻结PLAN明确采用233尺度。所有source/builder/driver/helper、PLAN、四RESULTS及四paired NPZ的hash匹配，预运行独审的源码未变。

自有 `result_checks.py/json`已实际只读运行：

- 从payload中的 `component_log_sum_normalized` 重构 `logD=logaddexp(log .5,log .5−log(RJ)+component_sum)`，四组最大误差均为0。没有重生未保存的raw counts，也不声称由payload独立重建每个component统计；RN的组件构成还依赖已审源码。
- 四组所有lower/upper contributions非负有限，逐样本lower≤upper，每个h的双方向cap均满足；最大contribution/cap约.97005。
- C_E与exp(max_logmean_inflation)一致；C_Z直接由保存的normalizer_max重算为 `exp(max(0,n log(normalizer_max)))`，四组均为1.0000000014887647。normalizer_min/max也准确重现max_abs_residual；R16的最大绝对残差来自负侧，不应直接把它当正残差指数化。四个RESULTS均已在冻结pilot前保存这些signed min/max字段。
- 每个保存方向均值、SE及ESS以独立math.fsum/缩放平方和重算，最大差分别7.11×10^-15、3.56×10^-15、3.98×10^-13。paired差和paired SE也直接从同一NPZ两方向相减重算，最大差8.89×10^-16。未重复反演全部504个CI，该项由另一独审负责；保存endpoint gap与两端CI相减一致。
- 各R的21个same-h `joint lower − prior upper`经验差全部为负；84个从504项同时CI产生的对应endpoint gap也全部为负。

|R|最大的same-h经验差（h）|最大的同时CI端点差（h）|
|---:|---:|---:|
|4|−.13214101（0）|−.65024871（0）|
|8|−.09166156（10）|−.63499109（0）|
|16|−.04257966（10）|−.64247116（0）|
|32|−.10280369（0）|−.62117876（0）|

这说明冻结候选上没有正分离证据。**负的lower−upper并不证明真实joint劣于prior，也不证明不存在分类增益。**不能用joint最高lower与prior最高lower的差替代本次应有的同阈值lower-versus-upper比较。

own结果校验首次在R16的metadata assertion停止，原因是检查器漏读signed min/max并把绝对残差误当正残差；最终版本直接使用已保存normalizer_max，简短留痕于result_checks.json。该过程仅重读已有数据，没有source重算、tail重采样或修改父线结果。

## 最终结论与残留误差

新source、冻结stream数学接口及保存payload一致性审查通过，未发现必须修改的实质数学/RN/端点错误。所有阈值仍是严格 `>λ`；zero variance在等号处lower取0、upper可保守取1，没有把等号当正下界。不同R的phase/source与cut规则不嵌套，不能解释为同源加表的单调证据。

当前结论仍属于普通FP有限pilot：三角函数、erfc、全部名义source矩、BLAS点积、proposal概率的实际机器归一化及统计端点未获完整向外认证。504区间只覆盖冻结有限族，不是连续尺度/全阈值weak上包。没有新Cn阶、没有已证log收益，也没有排除真实共同尺度改进。

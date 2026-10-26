# 157：155 新曲率与公共 retarget 的独立预审

理论与已读实现接口未发现阻断缺口。157 只拥有本新目录；154、152、145、141 等冻结件不改。当前只读生产源码、真实 schema 和 byte hashes，并做 AST 解析；没有运行157的新 q、log、Delta、角点或候选 floor 算术。生产155已完成的结果仍仅是待独审候选。

## 新曲率及尾质量资格

153 的公式已由独立154证明：对同模型、同 ε 和同 c，令 A=1/(16ε²)，有 `Qε≤min(qold,2A/3,A(7/10-8L/(15U)))`，其中 L≤H2ε(c)、Hε(c)≤U，旧 q 必须是同一 Qε 的合格上包。157 不导入153 penalty API，而以整数共同分母直接构造 split 分子 `A(21 U_n L_d−16 L_n U_d)`、分母 `30 U_n L_d`，再与 uniform 和原 q 取 min。

155 固定策略为 EPS64 无 H32 证据取 L=0、q=512/3；EPS128 原075对象无 H64 mass 证据取 L=0、q=2048/3；8135 个原129对象仅复用原 H64g105 支撑 L=21/200000 和同点原 H128 上端 U。这里 H64 的积分下端大于 H128，确实是 H2ε；反向用 H128 下界作 H64 下界无效。新 q 可以小于256，不能沿用旧129下限或 clip 回去。所有原 Hε Delta/整组梯度保持不变，H64 只进入质量比。

预审实际读取的生产准备计数是107530 canonical、111201闭叶、15 shards、3812原H64支撑叶；canonical rule 为97727 EPS64、1668 EPS128原075及8135 EPS128原129。未来独审将在任何新数值之前，逐项绑定149 canonical/aliases/closed leaves与151的原jet projection；对3812支撑回读真实098完整 state，核原 leaf/certificate SHA、EPS64、g105 target、原严格正 margin、原 accepted105 source SHA，再核每个129锚点在所用原闭叶内。不重算这些旧支撑角点，也不从 global H128/gstar 或新raw查找补 L。

## 一次16角点同时支付两目标

固定一 canonical 的原 g0=11/100000 下包为

`Φold(u)=Delta110+Σ min(γL_j u_j,γU_j u_j)−qold|u|²`。

因为 qnew≤qold，同一叶的每个 u 都有 `Φnew110(u)≥Φold(u)`，因此新最小 margin 不小于原保存正 margin。无需重算旧成功角点，就能核新角点结果与原正 margin 的次序，并保留完整旧域。

g1=23/200000，g1/g0=23/22。若公共 Ulog≥log(23/22)，则 `Delta115=Delta110−Ulog` 是合法下端，故 `Φ115=Φnew110−Ulog`。常数与角点无关，最小值和最坏角点都可由一次新16角点计算同时得到；同一最坏角点也给g115 margin。负 Delta、盒外 anchor、穿越奇异子空间均合法，仍是一 anchor/一 ε/整组梯度/全16闭角点，不拼接不同 anchor。

157 的角点实现先为整证建立一个共同正整数分母，以整数分子枚举16角点，区别于155 Fraction逐项求和。只核新曲率的角点，并检查 `new_margin110≥saved_old_margin110>0`。g115使用严格正门；零或负只进入该目标的pending，不损害已保留的完整g110覆盖。

## 公共 log 的独立正有理尾证明

令 z=1/45，准确 `log(23/22)=2∫0^z(1-s²)^-1 ds`。前64项是积分正多项式，余项

`R=2∫0^z s^128/(1-s²)ds`，满足 `0<R≤2z^129/[129(1-z²)]`。

157 为64个正有理项先形成统一整数分母，再直接累加整数分子，核生产 lower 和 positive tail 上端；不调用生产 shared_log、通用 log、exp 或Arb。此为一个新公共常数的独立证书，全部 canonical 共享其同一上端和hash，不重新求旧 log。

## 完整覆盖与后续接纳接口

源149的所有111201叶和52360根保持不变。157 将独立核原理上的 dyadic root box、闭边界、antichain、全部前缀 split 数58841和逐根 Kraft=1；同时核 canonical_index%15唯一 owner、每对象一次新q、每叶一次16角点和全部journal hash-chain/checkpoint。52360个有序指标root及坐标置换对称性覆盖整个 F，半群模型与对称性未变。各叶 cutoff 仅EPS64/EPS128，所以全叶的正新 margin 给整个 F 的 H128 下界。

完整最小值 m=min所有叶的新margin110>0，可给 `H128≥g0 exp(m)≥g0(1+m)` 在F上。生产保存的 linear floor 只作此次匹配对象；某些g115叶pending与这条完整域结论严格分开，不能从部分g115通过叶推出完整g115覆盖。

主后续拟采用 `P8(m)=Σ(k=0..8)m^k/k!`，只要m>0，正Taylor余项给 `exp(m)≥P8(m)>0`。因此主可另登记计算

`g192=2^-192 floor(2^192 g0 P8(m)) ≤ g0P8(m) ≤ g0 exp(m)`。

若主选择先下量化P8再乘g0，同样保持下界方向，但需明确最终网格口径。157本轮只核该接口合法，不计算P8或新g192，也不执行根/范数合成。这个H128下界仅在F上；完整G≥H128后才可借G的−4齐次性与cube对称扩展到Q，不把固定ε的H当齐次函数。后续成本区间下端只夹住所选证书表达式，不能据此声称C5的新下界。

## 有限实际审计合同

准备代码为 `sources157.py`（来源与整树门）、`integer157.py`（独立新q/log/整数角点）、`audit157.py`（15个显式 source-shard worker 与联合结果）、原字节复制的 rational152 gcd核心及jsonreader。源门先于任何新log/q/角点。总范围固定为107530个新q、111201×16=1779216个新角点、一个公共log证书；0旧H/gradient/角点/log基线重放，无新raw L、H、split或cutoff。

真实155输出须先以完整只读SHA形成binding，主审准确plan与binding后才启动全新owned输出目录。15worker各只加载自己的源码分片；已完成或失败输出不自动覆盖、重启或追加case。所有源与实际候选文件在运行前后检查hash。若任一字段/公式/结构不符则保留失败，不能以partial结果推全域。正式g/C5接纳仍由主负责。

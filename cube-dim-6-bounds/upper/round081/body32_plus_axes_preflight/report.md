# 081：body32 加六轴精确同模型去重封存

状态：**DUPLICATE_EXACT_MODEL — STOP，actual 0 新科学配置**。

这不是新的上界机制。012 已定义并实际计算本次指定的精确模型，
015 已写出匹配的 band 曲率接口，018 已完成其中一个固定短时势候选的
严格障碍证书。依主任务“先查旧记录，精确模型已有则去重停”，081 不再
构造 beta 表、不重新枚举 fan 或 pair planes，也不申请 native。

## 1. 同一模型的原件定位

旧 C012 定理的第一式就是

\[
A_\beta=[V_{\rm body32},\beta e_1,\ldots,\beta e_6],\qquad
\psi_\beta(\xi)=\sum_{v\in V_{\rm body32}}|v\cdot\xi|
                    +\beta\sum_{i=1}^6|\xi_i|,\quad\beta>0.
\]

所有 38 个独立标量使用同一个标准 Cauchy 时间 t，轴的尺度由列 beta e_i
吸收。原件明确固定 Lebesgue 测度、一般 L1 输入及所有正的共同 cube
半径，未改原始 C6 算子。它不同于055的 body+sparse4。

原件与主仓归档逐字节对照均相同。完整路径、大小和 SHA-256 见
`source_binding.json`；其 49 个文件全相同，4 个原有 seal/receipt 中的
31 个文件哈希全匹配。C012 和 A012 没有名为 manifest.json 的旧封存件，
这里如实绑定其 input_manifest/kernel_manifest/delivery_receipt，而不虚构。

| 原有对象 | 本轮实际读取的关键原件 | SHA-256 |
|---|---|---|
| C012 精确模型与 fan 证明 | upper/round012/augmented_frame/theorem.md | `9c722245d3da04254537dcb9776c08d4d099c09cef05cc7e809740f625ecf136` |
| C012 完整交付 receipt | upper/round012/augmented_frame/delivery_receipt.json | `0256b1229e6bcc7757a2745740d8f3638e4bbf3839edeca3a4f01a75e06447d2` |
| A012 全 beta 正的精确距离分类 | upper/round012/closest_subspaces/theorem.md | `5e34d41d83331cae2900f6952f34f86bc8bc4edcec9eb0938aa67e467d4bf65f` |
| A012 输入登记 | upper/round012/closest_subspaces/input_manifest.json | `dd4d1f2eda4ad31007c9d500b0c7c0e20232af24a696446568a396e7bc1a55c6` |
| F012 一般正 Cauchy frame 独审 | audit/round012/augmented_frame_theory/manifest.json | `c45fed7739a1c1d64cea0647f511d883d686eebfbcdde46cab0261dcaad99375` |
| A015 匹配 band 理论 | upper/round015/band_taylor_cover/theorem.md | `922f4132fec891fafa30cfe9aa8d89f91409e7a6c1a50852d34a92e56aa2b901` |
| F015 独审 | audit/round015/band_taylor_review/manifest.json | `841cbfd78b6e5986786a859c09055691c685c155c4e96f75f966d7f8665ab192` |
| B018 固定候选障碍 | upper/round018/augmented_obstruction_probe/manifest.json | `43d1912b43bbbd012484039afc771a3420680fb3af00652fcd1285cb894f17b5` |

C012 原件在 curvature-astra clone，A012/A015 在 theory clone，F012/F015
在 literature clone，B018 在本 upper-astra clone。归档对应主仓
`Math-Research-cube-dim-6/cube-dim-6-bounds`。路径映射均见机器来源表。
本轮只核上述静态原件与对应 seal 的文件绑定，不宣称重新认证历史上所有
传递依赖或已经改变的 STATE/ledger 文件。

## 2. 已经付过的 fan 和软件接口

C012 已证明：旧 21 个精确锥位于正有序 dual chamber，其每个坐标形式
在整个锥上非负。因此添加轴项在每个旧锥上仍然线性，原几何是新符号的
refinement，不需重新三角剖分。每个 ray 的新 rate 与 cone weight 是

\[
q_{\beta,j}=q_{{\rm body},j}+\beta\|r_j\|_1,\qquad
\alpha_{\beta,R}=\frac{|\det(r_1,\ldots,r_6)|}{\prod_jq_{\beta,j}}.
\]

完整 signed-permutation group 仍是 64×720；H 使用更新后的
`b_j=x·r_j/q_beta,j` 和同一个 DD5 公式。**q 与权重都要更新**，仅改符号
名称不够。零点、重节点及完整 Laurent principal-channel 消去要求仍在。

旧 C012 实际证据为 beta=1：21 基锥、126 ray occurrences、18 unique rays、
46080 group elements、967680 complete contributions，38 forms 的弱常符号
与新 rate 已经实际核对。其 fan_beta1.txt 的 SHA-256 为
`450caf04d124334474b91f3d9fdabce4eadf3259ea47def9094ae5eb056b6c4e`。
旧 B018 已用该 fan 运行一次 C014 native kernel。081 没有重新跑这些量，
`saved_slots.json` 中相关字段是原 JSON 的逐字复制，保留 source hash 和 pointer。

故“旧完整 fan 能否用于这个新加轴符号”不是未付的新问题，至少 beta=1
的软件实例也已实现。任意新有理 beta 的整数编码、int64/int128 范围及其
完整 node/confluence 输入仍需针对新输入单独验收；这里没有新验证或授权。
不把012/018的 beta=1 运行说成任意 beta 软件批已通过。

## 3. 已经付过的真实 S2 分类与曲率接口

A012 给出了适用于所有 beta>0 的解析分类（beta=0 不属本模型）：

- body/body：496 个不同的几何 pair planes。
- body/axis：192 对，去重后96个平面，全部已在上述496个中。若 w 是 v
  在第 i 坐标的翻转，则 `span(v,e_i)=span(v,w)`。
- axis/axis：新增15个互异的坐标平面，均不等于 body/body 平面。
- 总方向对数703，真实几何平面总数511；不是703个互异 tube。

令 z 为 |x_i| 的递减重排，`W(a)=sum(a_i^2)-(sum a_i)^2/len(a)`，则旧定理为

\[
\operatorname{dist}(x,S_{2,\beta})^2
=\min\left\{W(z_1,\ldots,z_k)+W(z_{k+1},\ldots,z_6):1\le k\le5;
                    \ \sum_{i=3}^6 z_i^2\right\}.
\]

这包括 ties、零坐标及符号置换，不意味着固定五个物理平面，也不意味着
tubes 不相交。body/axis 有真实 dependent triples，短时事件选择必须按
actual rank 而不是“不同三列必独立”。C012/F012 已按有限正 frame 的
rank-greedy 事件覆盖证明：G 与 K_T 在真实 S2 外有限连续，在 S2 上发散。

A015 已给出 `A_beta A_beta^T=(32+beta^2)I`，对应正 band Hessian 下界系数
`2/(32+beta^2)`，以及 cell loss `D*U2_upper/(32+beta^2)`。其中所有 H、
gradient、U2 数据必须来自同一 augmented frame。旧 body32 数值与31-tube
预算均不随这个 Gram 恒等式自动转移。

## 4. 历史数值结果的准确范围

C012 的 beta=1、T=2/45 两个旧点曾通过必要点测试；那是2G+2H四次旧完整
pass，并非全 cube 下界。其后 B018 在另一个精确有理面点得到严格反例。
B018 的完整 dyadic/rational 原值、点 numerator/denominator、frame hash
及 source pointers 已逐字保存在 `saved_slots.json`，未从显示小数还原。
旧证书的保守有理摘要为

\[
\frac{10048510817}{10^{15}}<K_{1,2/45}(z)
 <\frac{100485108171}{10^{16}}
 <\frac1{92160}=\frac{T}{4096}.
\]

旧真实 augmented S2 距离平方是
`6020805515960314567115330493/40000000000000000000000000000 > 0`。
因此旧连续性定理将这个点障碍提升到有正体积的 cube 相对邻域；并非
可忽略的孤立点。这里只引用已封存严格证书，不重新计算这个比较。

**它只排除 beta=1、T=2/45 的短时势全域证书**，不排除其他 beta/T、混合
时间、全 Green 正部势，亦不证明原始 C6 的下界超过64。旧 point-cost
约69.109只是这个固定辅助候选的必要成本障碍，不是 C6 的数值。

## 5. 没有由去重产生的新全域结论

C012 已证明的源/成本接口是
`L_beta K_beta,T=p_beta,T-delta_0`、`mass K_beta,T=T`，以及在真全 cube
essential lower bound m 成立时 `C6 <= T/(64m)`，适用于一般 L1 与全部
共同半径。这不等于一份同模型全 Green 正部 `R_c=(G_beta-c)_+` 的已付
竞争性 g/I 证书。

本轮不把“off-S2 有限连续”替代成全 face `I_beta=integral G_beta^(6/5)`
的已付数值预算；也不把有限511平面分类替代成含全部重叠层的 I 估计。
新的正部源/停止接口如需使用，应明确引用一般有限 frame 定理并逐项匹配
其假设，随后提供同一个 beta 的严格 g/I；081 在触发去重后没有另开这个
证明或计算分支。原 body32 g、I、31-tube 常数均不得继承。

因此本轮唯一验收对象是“精确重复模型的来源绑定与准确范围归档”。没有
未报告的 native 候选、未闭合后仍继续扩算的计划，亦没有新的 C6 上界。
要重新开放此旧模型，应先由主提出一个区别于旧点测试的新具体竞争性
g/I 或 cover 机制；不是再次执行已经做过的几何和 beta=1 旧点。

## 6. 实际执行与停止

只执行 metadata-only `python3 bind_sources.py`，PID85879，tool chunk
`0bf08d`，exit0；wall0.016263s，CPU0.034210s。49份原件及其主仓副本相同，
31项内部封存文件匹配，48个旧 JSON 字段按原值复制，provenance gap=0。
脚本设置 CPU30s/alarm170s，无子进程、无数值库。

实际新科学配置0，新 beta0，kernel0，fan/pair 重枚举0，旧 author program
执行0，MATLAB0；没有 stage/commit/push，旧目录未写。封存后停止本分支。

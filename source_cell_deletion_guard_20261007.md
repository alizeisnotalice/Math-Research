# 来源 cell 一般 p / 删除 scalar 守卫

按新登记，一次完成对冻结 spatial_packet_collision 批 15 个完整来源/74 个 receiver 的三轮后处理，没有重跑极大值、旧 oracle 或旧实验。64 个非空 row 的 p=1/4、3/4、7/8，q=2、4、16 共 576 项 Bernoulli 期望比较，以及 192 项删除 scalar 比较，全部严格有理 log 区间通过；unresolved 数均为零。10 个 empty、14 个 top 分支全部保留。

这只核验本批原真赢家后验上的有限 scalar 合同。没有把它作为一般引理的证明，没有对 16W 全空间积分作数值认证，不提供 Lebesgue receiver 覆盖或原 geom/近优资格。

## 1. 冻结原空间输入与 cell 分配

完整 source、a,b、源 weights、receiver x、真实连续 winner R/M、保存 D 与 half-open cell 分配，全部取自冻结原 profiles。每个原 cell μ_C 包含该 cell 全部标签（含重合标签及未捕获标签）；π_C=μ_C(Q_R)/m_R，χ=Σπ_C²。Σπ=1。未捕获 cells 的 π=0，但原 full-cell labels/weights 仍保留；没有按捕获量替换 μ_C 或重新最大化。

三轮沿原分组：第一轮 n=1,2,3；第二轮 n=4,8 的四种来源；第三轮 n=16,64,256，同时保留原 N 和窗口变化。完整种子、坐标、weights、source hashes 在原 input profiles 与新 case 副本中，输入文件 hashes 在 registration/final_receipt 中冻结。没有新的随机采样；以下 coin law 是 DP 精确积分，不是 Monte Carlo。

原 top R+2d>b 分支仍可核以下 posterior scalar 身份，但未借此追认其原扩张几何预算资格；empty row 不定义 π/V，显式跳过而不剔除数据记录。

## 2. 同 cell Bernoulli(p) 的真实期望

固定实例的全部 cells 上定义独立 η_C∼Bernoulli(p)，同 cell 的全部原标签共用同一个 η_C，ν=(1/p)Ση_C μ_C。因此 Eν=μ。对原 Q_R，记

    U=Σπη, V=U/p−1, 1+V=ν(Q_R)/m_R,
    ell=log q−log_+(q(1+V)), log_+(0)=0.

未捕获 cells 的 η 对这一查询无影响，故只在 DP 中边缘积分掉，不从 source 删除。每个 row 的 π 是原完整 μ/原 winner 所生，不是自由后验或组件单独 maxima。

设 p=A/B，π_C=u_C/T 为保存的整数 share。在逐 cell DP 中，reject 对整数权重乘 B−A，accept 乘 A；每个 subset score s 的概率为 N_s/B^K。全部概率之和为一。精确核验

    EU=p, EV=0, EV²=((1−p)/p)χ.

192 个 weighted DP law 均通过。K 最大 31、DP states 最大 159，最大概率共同分母 8^31=2^93；这是压缩积分所有 coin outcomes，没有抽样或舍弃小概率零接受分支。192 个 DP law 均保存 score=0 的正概率。

核验合同为

    E ell ≤ (2+4log q) ((1−p)/p) χ.

p=1/4、3/4、7/8 的 variance factors 分别为 3、1/3、1/7。原逐 row τ=M/q 只作为原输入身份随记录保留，不把不同 row 的 τ 当共同阈值。

## 3. 删除式及 π=1

对每个 q 核

    Σ_C min(log q,−log(1−π_C)) ≤ 1+(1+4log q)χ.

π_C=0 的未捕获 cell 项精确为零。π_C=1 时 −log(1−π_C)=+∞，min 项显式定义为 log q，不向 log helper 传零分母。28 个 receiver 出现这一分支，对三个 q 共 84 个 term，全部保存。

其余 0<π_C<1 时，使用精确有理 branch criterion：π_C≤1−1/q 等价于 1/(1−π_C)≤q，选 untruncated log；否则取 log q。等号合法，实际 log_q 缓存包围一致，不使用浮点比较决定截断。全批 term branches 为 618 个 untruncated deletion log、18 个 truncated log q、84 个 π=1，合计 720=3×240 个原 captured-cell incidents。

## 4. 严格 log 证据与三轮结果

只从已审 spatial_packet_collision_probe_20261007.py AST 提取 outward、atanh_log_interval、log_interval、log_plus_interval 四个函数定义到新的 Fraction namespace。旧 module/main/geometry/response 函数均不执行，也不导入会写 bytecode 的旧模块。旧文件 SHA、所选四函数的 AST SHA 及 100-term/192-bit 参数均在 registration/results 中冻结。

原 helper 用 range reduction 和 atanh 的正项截断，严格余项上界 2z^201/[201(1−z²)]，再以 Fraction 向外取整到 192-bit dyadic 格。每项 expectation 用精确 DP 概率乘 log 区间；最终将 Eell upper 或 deletion sum upper 与 RHS lower 比较。所有数据区间为有理端点，没有 Decimal 正残差替代认证。

|轮|n 范围|非空/总 row|top/empty|weighted DP laws|coin q 区间通过|deletion q 区间通过|π=1 row|
|---|---|---|---|---:|---:|---:|---:|
|1|1,2,3|12/18|3/6|36|108/108|36/36|9|
|2|4,8|20/24|3/4|60|180/180|60/60|2|
|3|16,64,256|32/32|8/0|96|288/288|96/96|17|

三轮 coin gap 的最小严格下界约为 0.252648、0.0239112、0.0531314；deletion gap 约为 1.278515、0.128134、0.284424。按 p 分组的全批 coin 最小 gap 下界约 0.515642、0.0557527、0.0239112。近似小数只供阅读，判定实际使用并保存精确 Fraction gap 下界。

ell 是 signed scalar，程序未把 Eell 正部化。p=3/4,q=2 与 p=7/8,q=2 各有 33 个严格负 expectation interval，其余 510 项 expectation interval 严格正；完整 raw intervals 全部保存。这不影响上界合同，也不能把负值删除后声称核验了另一个更强的正损失预算。

原 session 70184 正常 exit 0，主后处理约 1.982 秒；没有 active handle，没有扩大样本或挑选输入。exact finite DP 与 log enclosures 不需要概率置信区间。

## 5. 数据、hash 与适用边界

15 个新 gzip profiles 保存完整原 case/source、原 input profile hash、每 row 身份 hash、x/D/M/R、原 full cells/captured cells/π/χ、所有 p 的概率 DP、mean/variance 值、所有 q 的期望与 RHS/gap 区间、删除 branch 与逐项区间。empty 和 top 原分支完整保留；输入来源未改变。旧代码仅作为已审 log helper 的只读数据。

registration、执行结果/receipt、本文与 final_receipt 全部采用 exclusive-create。复现应复制到空的新目录，并保留原冻结 profiles 与 helper 文件引用；不能覆盖本批或旧批。全体输入/输出文件 SHA256 和 helper AST hash 在 final_receipt/registration 中。

本批没有新 partial-cell 几何、receiver 积分、source-tail 阶数或一般 16W 支付实验。它检验的是同一原 source 与真实 winner 的已保存 π 上，规定一般 p 与单 cell 删除的两个新 scalar 合同。

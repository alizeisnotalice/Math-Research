# 固定 physical-cell 异质保留时钟：局部 scalar 守卫

新登记后对原 15 个完整 source/74 receiver 一次后处理，原 winner 与 μ 不变。三种固定 physical-cell patterns 的 192 个 exact DP law 均满足 EV=0、EV²=Σπ_C²(1/p_C−1)；576 个 q 检查全部通过，无 unresolved。37 个零方差 patterns（111 个 q）按 exact algebra equality 通过，不能说成有正 gap。其余 465 个 q 有严格正的有理 gap 下界。

这里仅检验 fixed patterns 的局部时钟损失 scalar 门，不是自适应政策模拟，不是全空间 22 质量-work 积分或已闭合 geom 上界。

## 1. 全源 cell p_C 的共同、预先分配

使用 spatial_packet_collision_probe 的冻结真实连续 winner 后验 π_C=μ_C(Q_R)/m_R。cell 是原完整半开格，不拆标签；原所有 sources、完整 cell weights、未捕获标签、x/D/M/R、source hashes 和原 seeds 均保留。

每个实例（原 a,b,d 不变）一次冻结完整 cell p_C 分配。对 shift s=0,1,2，

    p_C=[1,1/2,1/4][(first full cell index+s) mod 3].

Python floor/mod 给负坐标同样固定的 residue 0/1/2。登记从完整 full-cell 列表生成分配，不读取捕获 π 来选 p；运行逐 row 断言 full-cell 列表与登记相同，所有 receiver 查询使用同一份 map。窗口改变的旧实例仍保持其原独立 grid 定义，不把不同 grid 的 index 混为同一物理 partition。

每个完整 μ_C 用同一个 η_C∼Bernoulli(p_C)，ν=Ση_C μ_C/p_C；同 cell 的标签共同保留，不按 receiver 重分 η 或 p，Eν=μ。未捕获 coin 只在查询边缘 law 中积分掉，完整 p assignment 仍保存在新 profiles。没有实际抽 coins，没有 Monte Carlo 或新的随机 seeds。

## 2. 真实响应 DP 与 variance

对原 Q_R，V=Σπ_Cη_C/p_C−1。π_C=u_C/T 为保存的整数份额，p_C∈{1,1/2,1/4}，所以每个 accepted response 的整数位移 u_C/p_C∈Z。DP reject/accept 概率按每个 p_C 的有理分子/分母精确累乘，共同概率分母为所有 captured p_C 分母之积；p_C=1 的 reject 分支概率零，不保留虚假 outcome。

精确核验 probability sum=1、EV=0 与

    EV²=Σ_Cπ_C²(1/p_C−1).

所有 192 patterns 均通过，最大 compressed DP 状态数 181。完整概率 numerator、denominator、response shifts、captured p_C 和 full-source p maps 均保存。

## 3. log interval 与零方差

固定 q=2,4,16，ell=log q−log_+(q(1+V))，核

    Eell ≤ (2+4log q) EV².

只 AST 复用已审四个 log helper，100-term atanh 严格尾界、192-bit dyadic 外取整；旧 module/main/max 未执行。每个 outcome 若 q(1+V)=q，直接使用 ell=0 的精确恒等式，不把同一 log q 当作两个独立区间相减。

EV²=0 iff 全部被捕获 p_C=1，此时 DP 唯一 score=T，V=ell=0，RHS=0。37 个这样的 pattern、111 个 q 的 [0,0] equality 均保存。原完整 source 在未捕获 cells 仍可能改变，零行仅表示原查询没有 coin 方差，不推出全空间损失零。

对非退化 law，用 signed expectation upper 与 RHS lower 比较，全部 465 项严格正 gap。ell 未正部化；一个 expectation interval 严格负，其 raw 值保留。其他接近零的 expectation 可受 log 包围宽度影响，但不改变已认证的正 RHS gap，也未改输入避免它们。

|轮|n 范围|非空/总 row|top/empty|pattern laws|零方差 patterns|q 通过|nondeg 最小 gap 下界（约）|
|---|---|---|---|---:|---:|---:|---:|
|1|1,2,3|12/18|3/6|36|10|108/108|0.471396|
|2|4,8|20/24|3/4|60|5|180/180|0.160610|
|3|16,64,256|32/32|8/0|96|22|288/288|0.118399|

主 session 33712 正常 exit 0，约 1.490 秒；无旧 oracle、重新最大化或 sample expansion。empty/top 原分支全部保留；top 行只做 scalar，不追认其扩张几何资格。上表小数供阅读，判定使用保存的 Fraction interval endpoints。

## 4. 保存与未证范围

15 新 gzip profiles 含完整原 source/case、输入 profile hash，三个 full-cell assignments，每 row x/D/R/M/原分支、row hash，每 pattern 完整 DP 与所有 interval/equality 分支。registration 冻结 code/helper hashes、完整分配、q 和预算。代码与输出 exclusive-create；复现应复制至新空目录并保留 frozen input profiles/helper 引用。

适应过去的 p_C 政策仍需以真正 conditional coin law 和逐时 source-work 结算证明，本批只是固定控制参数的必要 scalar 守卫；不是该自适应一般结论的数值证明。没有量 receiver 体积、Φ、nearmax 或 √n 阶。最终所有 input/output/hash references 在 final_receipt。

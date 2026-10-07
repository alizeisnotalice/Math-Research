# 质量倾斜 source spine 的有限停止复杂度守卫

新登记三轮 N=3、5、8，完整正来源各为 3、5、8 个不同固定半开 cell。全部 Fraction：293 个非空 subset 的真实 Q 率/后验、879 个有限时间状态律、20,484 项 mass-tilt 与 direct-tag transition 比较，三种停止规则的 879 项 G=w_A H 与 879 项 F=w_A 比较、586 项 count harmonic 公式均通过。失败数零。一次主运行 exit 0，约 0.530 秒，无未结进程。

本批不计算 Φ、cube maximal、nearmax 或原 geom/history 资格。n=1、16、256 只确认 source 位于不同物理格，绝不从这三轮拟合维数阶。原 thinning 损失合同能否在允许预算内使任意近优来源进入已付类别，仍未由本守卫解决。

## 1. 工作流、完整来源与真实跳率

已读 E04/J01 的 SKILL、provenance、method、cube-interface，以及 source_cell_thinning_budget_20261007.md。这里使用实际有限状态首跳与停止边界流程；没有引用未核外部文献或把 maximal 算子假设成内生过程，也不使用无限时域 optional-sampling shortcut。

三轮 raw cell 质量固定为 (1,2,4)、(1,2,3,5,11)、(1,1,2,3,5,8,13,55)，除以各自总和 7、22、88 后完整 W=1。每 cell 一原子 y_j=(3j,0,...,0)，窗口 metadata [1,2]、d=1/(8n)，固定 cell index=(24nj,0,...,0)。三輪 cells 分别恰为 3、5、8 个不同物理 cell，源标签不拆格，不重整为被保留子类。整个 source、位置、权重、source hash 与预算均在响应/过程表计算前登记；无需随机 seeds，所有状态 exhaustive。

P 下每个 cell 独立寿命 Exp(1)，当前 alive set 为 A；来源为 ν_t=e^t Σ_{j∈A_t}μ_j，增长补偿后的质量 W_t=e^t w_A。原 P 各活 cell 死亡率 1、总率 r=|A|；质量漂移增长 w_A 与总删除质量 Σw_j=w_A 相消。非空 A 的 mass-tilted Q 为

    tag posterior(j|A)=w_j/w_A,
    lambda_j(A)=1−w_j/w_A, sum lambda=r−1.

所有 1,116 项死亡率属于 [0,1)，293 个 posterior 和为一、293 个 total-rate 及常数生成元残差为零。单标签 state 的 Q 总率为零，无 Q empty state。停止前的这些是真实 rates，非均匀死亡近似；停止边界另将过程冻结，不把 terminal r>1 的原 rates 当成仍在继续收费。

## 2. 固定时间的精确 change of measure 与 hidden tag

p=1/2、1/4、3/4 表示 e^(−t)，无需数值 log(t)。对每个非空 initial A 及所有 B⊆A（含空集），直接独立 survival 得

    P_A(B)=p^|B|(1−p)^(|A|−|B|),
    dQ/dP(B)=w_B/(p w_A).

独立校验不从这个密度公式反推：先以 w_j/w_A 抽 tag，tag 永远 alive，其余 cells 独立保留 p；对每个 tag 分别作乘积 survival 概率，再加权求和。全部 20,484 个 outcome 与 tilted P 完全相等。每个初始状态的 P、Q 总概率为一，empty 在 P 有正概率但 Q 密度与概率精确为零。879 个 law 均有 E_P[e^t w_B]=w_A。

这部分是未停止链的 deterministic-time law；没有擅自将固定时间密度当作任意停止后的概率身份。下一节以 finite 首跳递归直接核停止预算。

## 3. 三种 adapted 停止规则与补偿质量

规则是 count≤1，count≤ceil(N/2)，以及 count≤1 OR max_j w_j/w_A≥3/4（等号停止）。在当前 retained source 上判定，适应当前状态；不是使用 future/source-unknown 参数。每个规则在全部非空 initial A 上计算，terminal 可以有多于一个 cell。

记 H(A)=E_Q T_A，G(A)=E_P∫_0^T e^t w_{A_t}dt，F(A)=E_P[e^T w_{A_T}]。特别地，F 是补偿增长后真正来源的终端质量，不是未乘 e^T 的 raw survivor mass。

terminal 上 H=G=0、F=w_A。nonterminal r≥2，P first hold T_h∼Exp(r)、deleted label 均匀，故

    E exp(T_h)=r/(r−1),
    E integral_0^T_h exp(t)dt=1/(r−1).

其后 source 的增长继续带 e^T_h 补偿。这给精确递归

    H(A)=1/(r−1)+Σ_j lambda_j/(r−1) H(A−j),
    G(A)=w_A/(r−1)+Σ_j G(A−j)/(r−1),
    F(A)=Σ_j F(A−j)/(r−1).

每个 finite table 按 subset count 递增计算，保存所有 child contributions、holding 因子和 terminal 状态。最多 N−1 次删除，非终态 hold rate 至少 2，e^holding 的矩有限；全部递归是有限且正的，不将未核可积性藏入 optional stopping。

精确逐项得到 G(A)=w_A H(A)、F(A)=w_A。固定 count k 时，H=0 对 r≤k；其余 H=harmonic(r−1)−harmonic(k−1)，共 586 个状态核验。不能对 r<k 硬用负的 harmonic difference。

|轮|N/n|非空状态|finite-time laws|tilt/tag outcomes|停止 rule-state 比较（G/F 各）|count harmonic 比较|
|---|---|---:|---:|---:|---:|---:|
|1|3/1|7|21|78|21|14|
|2|5/16|31|93|726|93|62|
|3|8/256|255|765|19680|765|510|

完整 initial A、w_A=W=1 的 H（亦等于 G）为：

|N|count≤1|count≤ceil(N/2)|singleton OR share≥3/4|F|
|---:|---|---|---|---|
|3|3/2|1/2|8/7|1|
|5|25/12|7/12|193/132|1|
|8|363/140|319/420|31013/31680|1|

share 恰等于 3/4 的停止状态三轮为 0、1、2 个，均按 inclusive threshold 保留；未用容差移除。不同 masses 改变 share stopping 的时间，count 规则只由 r 的总率 r−1 决定。这是 input-count 复杂度，不是维数复杂度或任意 geom 路径费用。

## 4. 冻结输出和范围

3 个 gzip profiles 包含完整 μ/权重/源 cell index，全部非空 A 的 tag posterior/death rates，每个 A,p 的全 B 分布及逐 tag 概率，三种停止规则的完整 H/G/F/terminal 表与 child contributions。零来源 outcome 不隐藏；Q singleton 零率与 count/dominance 终态在数据中显式区分。

同前缀 .py 是 standalone 新 Fraction 检查；register/run exclusive-create，registration 冻结脚本及读入 workflow 快照的 hashes。旧 max、MC、thinning 脚本未执行。结果、本文、执行收据、所有 profiles hashes 在 final_receipt；复现应复制到空的新目录，保留显式 source 公式，无需旧数据或外部库。

一般来源 work 预算与随过去控制的 cell 时钟需要另外的解析合同；本批仅核有限 mass-tilted spine 和三个 prescribed 停止 rule 的来源一次结算，不能将 F=w_A 本身说成已支付 receiver 或原 Φ 损失。

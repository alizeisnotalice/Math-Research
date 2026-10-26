#105：有来源资格的 sharp Khintchine 与任意三维角平均迁移

## 结果与限定

**既定三组 exact 认证一次通过。** 从公开 Xie v1 Theorem1 的稳定性输入，
得到对任何三维线性子空间 S 的真实 full-body32 time-one peak 上界

    Q3(S)<=eta/[pi^2(16sqrt2)^3],
    eta=95408781/97336000
       =15/64+(49/64)(114/115)^3<1。

与098/104沿用的 moment 系数相比，sharp-only 与105 cap 的固定变化为：

|对象|准确因子|展示值，仅供阅读|
|---|---|---:|
|sharp-only 的 norm/线性系数|(2/3)^(3/2)|0.5443310539518174|
|105 cap 的 norm/线性系数|eta(2/3)^(3/2)|0.5335534881029437|
|sharp-only 的 I 预算|(2/3)^(9/5)|0.48198745386564384|
|105 cap 的 I 预算|eta^(6/5)(2/3)^(9/5)|0.4705584278108021|

这些是同一固定上预算表达式的变化因子，不是实际 G 或 I 的下界。
所需完整有理端点在 results.json，未从展示数恢复输入。

本轮作为公开定理迁移的独立备选来源。主另有108候选机制；本105仅对其
给定表达式作条件代数比较，没有证明或接受其机制。若108的统一界
Q3<=1/(4096pi^3)最终获接受，则其系数相对旧098约0.49007012926381516，
相对105 cap约0.9185023436099459，严格更小。

按运行前锁定的修订范围，比较通过后省去完整104预算重算。没有执行旧104/100、
kernel/native或旧网格；没有新 Cn，更没有 C6 下界推进。global g_L、完整
同模型 I_U、far/补域以及完整正 S 仍未支付。105理论仍须后续独立审查。

## 公开来源资格未随数值通过而改变

主取得的 Xinyuan Xie 原始8页 PDF、完整 TeX 和 tar 已按指定 SHA 绑定。
实际全文审读的是原始 TeX；本作者没有另作8页PDF视觉审查。该文是
[arXiv:2610.07526v1](https://arxiv.org/abs/2610.07526v1)，2026-10-05预印本，
不是期刊论文；使用的是 Theorem1 的不等式，不需要其“最佳常数”主张。

[FW26官方原HTML Theorem1及证明](https://arxiv.org/html/2608.27908v1) 已读并保存，
不是全文独审；其证明仍调用经典 sharp real Khintchine。经典 Szarek/Haagerup
原始 PDF 没有读全，主记录403/截断，不能说本项目独立读证了那些原始定理。
Kershaw原始页面亦403；PROOF.md AppendixA以 digamma-Laplace 正差与极限
给出本次所需 gamma-ratio 特例的独立解析证明，不冒称读过1983原文。
资料中的非数学对话/操作材料没有作为数学来源或工具操作指令。

|原始证据|SHA-256|
|---|---|
|Xie PDF|5731a47fcbd91d116778df9fa0033353a82f9e24db3672e6cdb48f6552bc31d4|
|Xie原TeX|794b7efaead1ea2288baac88ef0d6c7fd87812e9225abd001439f155bd41f074|
|Xie source tar|01029ebde7bde564bfb173ede10aa4789eb0c5c5f16e88d50937cbf9bbc87bb8|
|FW26原HTML|3317d13947b2b3e69f63ef7f5f5edc01476450eca2cc44cbc445123d7ca7845b|
|完整迁移证明与gamma特例|d33fb0e9cfd974c6691b6645528df363d64c20d41219058df5fbda37cc8a6118|

原 TeX 中采用范围的 case split、Fourier 正负方向、cosine tail、Gaussian
remainder、concavity endpoint 等已逐段分析，未见内部缺口。这个审读状态
不能抹掉外部 classical 输入的资格。完整路线与独立Gamma特例见 PROOF.md，
逐源阅读范围见 SOURCE_QUALIFICATION.md。

## 任意三维子空间的机制

full body32 保持原32个 unit-weight、alpha1、common-time Cauchy、Lebesgue
及 Fourier(2pi)^-6。psi=32 E|sum xi_j eps_j|。Theorem1 经排序、符号与置换
作用到60个 signed two-coordinate unit equality vectors，给

    psi(theta)>=16sqrt2+32c_*d(theta),
    c_*=(2-sqrt2)^(3/2)/8。

对于任意3D S 和任一 ambient unit equality vector u，令 r=|P_Su|<=1。
chordal 距离||theta-u||<delta 等价于 theta·P_Su>1-delta^2/2；若投影太小
cap为空，否则其面积比是[1-(1-delta^2/2)/r]/2<=delta^2/4。
所以即使 u 不在 S，cap 也不会更大。delta=1/8 时，单 cap<=1/256，
60 cap 的并至多15/64，无互斥假设；至少49/64的球面有d>=1/8。

严格有理平方比较给 c_*>1/20 和相对增益>1/114。于是余下球面上的
psi^-3 至多原sharp upper的(114/115)^3。正球面积分再接真正三维
Fourier peak，得到eta。这里没有把平均改进误写成全空间点态 coercivity
超过16sqrt2；equality directions 上不可能有这种提升。

相应系数精确简化为

    K_sharp=sqrt2/(192pi^4), K_105=eta*K_sharp,
    D104_sharp=Sigma/(16pi^4), D104_105=eta*Sigma/(16pi^4)。

它们只适用于实际 full-symbol 3D restriction，不能移给 omitted-frame law。
104恢复真实参数后的 Q3 正是这种对象。其所有固定 P/B norm 对该共同系数
线性，其 I allocations 为6/5次齐次，故表中的新旧预算关系严格成立。
105没有替104认证其几何或恢复证明，也没有重算其 source/time/spatial预算。

## 条件108比较与停止

本轮第三组只测试条件表达式，不测试108的新证明：

    C_BL/C_105=2sqrt2/(pi*eta)<1。

有理数eta>31/32，pi>3，而93^2>8*32^2，所以pi*eta>93/32>2sqrt2。
正有理根包围另给严格分离的系数区间。由同一固定104预算的正齐次性，
该条件系数将严格优于105 sharp-only和cap的所有原计划 norm/I预算。
因此按主在运行前的修订，不再作重复完整104重算，也无参数扫描或自动fallback。

108作者另有三组通过的状态不改变这里的资格：其独立理论审查尚待，本105
比较没有把该前提当成已审定事实。105自身的任意3D cap迁移保留为有来源
资格的备选，供独审和主最终采用。

## 实际执行与最小保存检查源码

主写入 parent_release.json，SHA440192f04575f76db9329b91b072eb2cc638b0d5d13d3634885ffcbc0412efb7，
先静读完整证明、Gamma特例、完整脚本，并核24个源与prepared pins，才放行
唯一既定3组。原prepared材料保持原hash，不能把它们的“未运行”历史状态
当成本报告的当前执行状态；本报告及最终manifest记实际结果。

唯一数学程序 exit0，PID51537，3组、91断言。结果保存前数学CPU0.002888秒，
receipt前CPU0.003249秒；完整观察到的tool roundtrip0.059秒，包括最终写入、
stdout、exit/reap和工具开销。24源及prepared pins在程序前后相同。
没有失败或自动重试，0kernel/native、0旧104/100/网格、0MATLAB、0递归代理。

用户当前已取消程序运行时间上限。旧177秒/20CPU只保留为本极小批的自选
短反馈设置，绝不写成当前用户要求。没有因runtime门槛增加研究批次或刷计数。
数学用本机标准库Fraction、整数幂和4个已保存正根包围，精度固定60位。

check_saved_minimal.py 是新准备、**本作者未执行**的最小保存检查源码。
它只核文件hash、保存根的整数幂关系和有限系数恒等式，不搜索根、不执行
旧程序、也不证明外部定理或108前提。实际运行证据只有上述唯一3组。
复核可另记账运行该保存检查器，不能把未运行源码写成额外PASS收据。

必要材料为来源资格/原HTML、限定scope、解析证明、固定输入和源清单、原
prepared代码及绑定release、完整有限结果/收据、最小保存检查源码与本报告。
旧104快照只保留为原计划未执行输入，明确不计重算。所有写入仅在本人105
目录，未stage/commit/push。封存后停止，不自行启动随后108审查任务。

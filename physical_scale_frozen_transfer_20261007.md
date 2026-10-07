# 物理尺度移动与 frozen maximal：可付 Gaussian 尺度网，尚缺实际硬平均转接

2026-10-07。本稿只编辑独立文件及同前缀守卫。当前 actual 余项仍为保留原全部门的 R_dagger；本文不更改总账。新得到的是固定 c 的 frozen semigroup 在完整物理窗口上 O(√n log(n+2)) 的弱 maximal 工具，不声称已经支付实际 FIRST 流。

## 1. 读取范围与实际字段

应用 J05，实际读取 SKILL、provenance、method、cube-interface；它只规定核对真实生成元和补偿合同，不提供本稿定理。重新读取原 `概率接口阶段证明.tex` 4873–5070（NB）、5430–5585（EI/FJL），`short_shell_jump_face_audit_20261007.md` §3、§5，`actual_gate_contract_audit_20261007.md` §2，以及新 all-visited 稿。原下界总表 A/B 保留任意正联合来源与增长 label 复杂度，本稿没有 atom 数、产品输入或固定占用假设。

实际目录仍是同一完整 μ、原高来源分配、原 y,z、共享 θ/唯一 distinct-child LCA、出生 t_C、失败币、重捕获、严格 CP/GP、owntrace/score/时钟，原可测选择 σ(x), L_s(x), R_h(x) 与 full future cap。原早首跳 t<T_n，c_t=1−t，t≤min(σ,1/2)，跳点 w=y+r e_i 必须退出原格。实际 continuation 为
\[
h_{L_s}^{\otimes n}*K_{\sigma,t,L_s}
=\big[r_t h_{L_s}+(1-r_t)(h*G_{c_t})_{L_s}\big]^{\otimes n},
\qquad r_t=(1-\sigma)/(1-t).
\tag{1}
\]
原门在 y,z 上，不移动至 w；实际 history 不能由新参考先验重新判断。

NB--L 只对同一个固定 η_s 比较一维轴向退出，给 E_{s,L(s,x)}η_s≤2E_{s,b}η_s；不比较各尺度自己的 killed/live 源，不支付选时图。NB--S 的等待时间加权列不能倒向原响应。EI/FJL 的一份小退出源已有 N_h D_n 付款，不重加。新 frozen 定理每个固定 (c,L) 有各自障碍域，不能免费把不同 L 的障碍源合并。

## 2. 一个全维 Gaussian 的两端点尺度界

对任意 1≤m≤n、固定正对角协方差与位置，Gaussian dilation 密度（无关常数省略）为
\[
g_u=\exp[-mu-Qe^{-2u}],\qquad Q≥0,
\tag{2}
\]
其中 u=log(L/L_0)。在区间 [u_0,u_0+δ]，若最大值在端点则无需误差。若内峰 u_* 存在，Qe^{-2u_*}=m/2。取离峰最近的端点 v，d=v−u_*，|d|≤δ/2，精确地
\[
\log[g_{u_*}/g_v]
=m\{d+(e^{-2d}-1)/2\}
≤m d^2e^{2|d|}≤mδ^2e^δ/4.
\tag{3}
\]
Taylor integral remainder 给中间不等式，正负 d 都成立。若 δ≤1/(2√n)，则 mδ²≤1/4、δ≤1；e^δ<3、log(4/3)≥1/4 给
\[
\boxed{\sup_{u\in[u_0,u_0+δ]}g_u
≤\tfrac43\max(g_{u_0},g_{u_0+δ}).}
\tag{4}
\]
使用 e<3 的严格余量，(3) 上界≤3/16<1/4。m=0 的共同原子分量不随 L 变化，亦满足相应正测度不等式。这里是 centered Gaussian 的共同 dilation，不允许把不同均值随 L 移动的 Gaussian 直接套入 (2)。

## 3. 固定 c、完整 holding/face 的移动尺度 maximal

固定 c∈(0,1]，令 H_t^{c,L} 是新 all-visited 稿已证 maximal 的**完整**冻结 semigroup，包含 holding 与 proper faces。原 NB--A 的正热钟表示给一维 G_c 的 centered Gaussian 概率混合。固定 t 的 compound-Poisson 展开，在各坐标 jump 标签与所有热钟条件下，终核是某个 m≤n 维坐标面上的 centered Gaussian；其热钟、标签权和 holding 权都不随 L 改变。多次同坐标跳跃仅使 Gaussian 方差相加。零热钟可并入未活动面。

因此 (4) 对每个潜在标签可用，**先以同一潜在律正积分**，得到对每个同一 t 和 L∈[L_j,L_{j+1}] 的正测度比较
\[
H_t^{c,L}≤\tfrac43(H_t^{c,L_j}+H_t^{c,L_{j+1}}).
\tag{5}
\]
面维数及协方差不改变常数；没有将 partial faces 丢掉、没有作 posterior independence 假设。可对任意非负 L¹ 密度 ν 卷积。可测性及完整 supremum 不单凭 L¹ 强连续推出：正潜在混合定义一个共同 Borel 版本，固定 t 时每个内部 L 的 Gaussian 密度由两个端点的和支配，端点卷积对 a.e. x 有限，故 dominated convergence 给该格内连续。要同时覆盖 t，可取一份固定参考轴向正算子 \(\overline{\mathcal J}=2c^{-1}\sum_iG_{c,b}^{(i)}*\)；一维 Gaussian 混合比较给 \(\mathcal J_L≤\overline{\mathcal J}\)。任意固定紧时间窗 [0,T] 的正幂级数由 \(\sum_{k≥0}T^k\overline{\mathcal J}^{,k}ν/k!\) 支配，该函数 L¹ 范数≤exp(2nT/c)||ν||₁。逐有限项用上述尺度 dominated convergence，随后尾的一致支配，得到 a.e. x 的联合紧时间连续版本。可先在有理端点及可数 T 上去除共同零集。这项有限但可能很大的定义性范数不作为付款常数；(5)–(7) 的真实常数不变。

若 a<b，预先取 M=ceil[2√n log(b/a)]≥1、L_j=a exp[j log(b/a)/M]（j=0..M）；单尺度 a=b 只需一个节点。每格 log 宽≤1/(2√n)。于是
\[
\sup_{L∈[a,b],t>0}H_t^{c,L}ν
≤\tfrac83\max_{0≤j≤M}\sup_{t>0}H_t^{c,L_j}ν.
\tag{6}
\]
令 C_fr(n)=[(2+8C_0)log(e1024n²)+5/4]，沿用新 frozen theorem；正水平集 union 给
\[
\boxed{\alpha|\{\sup_{L∈[a,b],t>0}H_t^{c,L}ν>\alpha\}|
≤\tfrac83(M+1)C_{fr}(n)\|ν\|_1
=O(\sqrt n\log(n+2))\|ν\|_1.}
\tag{7}
\]
这是**同一输入**的有限并集付款，来源质量不按尺度重定义；O(√n) 个端点的费用已经显式计入，并未宣称一次 dilation 就免费消除 supremum。常数对固定 c 一致，尚未得到 sup_c。另若预先固定来源拆分 Σ_jν_j≤ν，且每个 j 只用一个固定物理尺度，则对共同阈值 α 的集合并，仅有 α|∪_j{H_j^*ν_j>α}|≤Σ_j C_fr(n)||ν_j||₁≤C_fr(n)W；预先固定的互斥 receiver 分区取对应 j 也可用此界。这不估计 sum_j H_j^*ν_j 的弱 L¹ 范数：弱 L¹ 不可直接相加。actual 输入没有上述源先验分配合同，不能只由 receiver 分组生成它。

## 4. 原实际转接仍有两个不同缺口

本稿 (7) 只处理空间 dilation；非齐次 softness 的同参数冻结转移由另一任务审计。本稿不把 K_{σ,t,L} 假定等于某 H_τ^{c,L}，也不重新付旧 near-full 或低 β 支。

即使未来有人给出 K≤C H 的正支配，原 (1) 还有随 L 移动的 h_L^{⊗n}。若把 h_L*ν 当每个端点的新输入，质量均是 W，却已不是同一个来源函数。直接用 h_L≤(L_{j+1}/L)^n h_{L_{j+1}}，在 (7) 的格宽约1/√n 内损失 exp(O(√n))；要常数控制则需1/n格和 O(n)节点。这就是 NB 末尾仍留下 N_h 的具体物理接口，而非新的未付行条件。

Gaussian 论证亦不能在条件 uniform seed u∈[−1/2,1/2]^n 后免费用于 h_L*H：给定 u 后，Gaussian 中心是 Lu，随尺度移动，(2) 不成立。均匀硬平均有实际面跳；作为纯尺度工具，h_L 单个原子在窗内到达时，取到达点正上跳的硬核与相邻格端点比为 (L_{j+1}/L)^n，足以否定 1/√n 网格的同常数硬核比较。此仅为**转接不等式障碍**，不是满足 FIRST/CP/GP/history 的 actual 反例，亦不否定目标弱费。

最小尚欠的可付款合同为：在不更换原 y,z 与共享 seed 的情况下，将当前 R_dagger 的 actual soft continuation/hard 接收交通正支配或水平集控制于 (7) 的一个共同固定密度 ν（或预先源分配 Σν_j≤Cμ），系数≤polylog；同时保留 actual FIRST/LCA/history，且不再以 h_L 改写为逐输出新源。若只有 λE 的 weak 控制，不能自动积分 actual traffic；但完整实际行上界 C_h q 与 paid 输出集合可给实际交通≤C_h q|E_paid|，前提是有阈值比较把该集合接到 (7)。本稿尚无该比较，故 (7) 不进入原 cube 主账。

## 5. 下界压力与新守卫预登记

下界总表 A/B 的 arbitrary positive 联合标签、增长 Cantor/XOR/Gibbs 及宽峰输入均满足 (7) 的任意非负密度范围；没有按 label 数收费。它们的真实 cube 响应含硬平均及原尺度选择，不能因新 Gaussian 网只用 O(√n)节点就套成 cube theorem。旧 lower 模型不重跑。

预登记仅新尺度组成守卫：三轮平方维数 n=576/1024/4096；纯 Fraction 区间 Taylor 检查 (2)–(4) 的 Gaussian log 比较、不同 face 维数、内峰/两端峰、共同尺度端点和闭端点，及固定正混合的两端上包代数。另核 mδ²≤1/4、3/16<log(4/3)，有限并集的系数和来源固定性。用 exp 的正 Taylor 余项及 log(1+x) 的交错余项作有理外包，不以浮点 PASS 认证端点。该守卫不认证 G_c 数值混合、原 softness 转移、h_L 硬面或 actual FIRST/history。执行后仅在本稿补终态，不重启已有实验。

执行终态：`physical_scale_frozen_transfer_exact_guard_20261007.py` 与 `_results.json` 已保存，无随机种子。三轮分别 1707/2603/8747 项，总计 13058 项精确谓词 PASS（含一项共同 log 常数检查）。每轮对所有 m=0..n 核连续误差预算，对 m=1,n/4,n/2,n、7 个固定 q、17 个含闭端点的 u 核 476 个有理 Gaussian log 外包；另核 holding、正混合及来源固定的有限并集系数。最小 log margin 的近似仅供阅读，原 PASS 使用 Fraction 正下端及 SHA，未由浮点决定。没有运行旧 all-visited/Beta/core/下界实验。

本轮可独立提交的严格增益仅为 (7)，它是从单个冻结物理尺度到整个 [a,b] 窗口的实付弱 maximal 推论；目前没有将 (7) 接入 R_dagger。因此原实际主账费用与未付项均保持。

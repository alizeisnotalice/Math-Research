# B055 postrun 独立审计：固定 frame/time obstruction PASS

**PASS_SAVED_ONLY_FIXED_FRAME_TIME_OBSTRUCTION。** 对真实 sparse4 frame beta=1/80、唯一 T=32/825 和 x0=(1,1,1,9/10,27/50,0)，已完成的 C055 区间满足

\[
 K_T(x_0)\le U<\frac1{105600}=\frac{T}{4096}.
\]

结合已审 off-true-S2 连续性，这给出 cube 上 essential infimum 的严格上界障碍：此固定 frame 与固定时间不能通过 m_T>T/4096 的 strict-64 必要条件。它不排除其他 T、其他 beta 或其他 frame，不涉及原始 R_c 的排除，不是新的 C_n 上界。

本次只读取冻结输出和回执，实际运行一次 6 配置 saved-only 审计；production native、旧 6 项预检重跑、group/pair 重新枚举、矩阵计算、作者 writer、MATLAB 均为 0。原 precheck 目录及 producer 文件均未改动。

## 1. 源锁与运行证据

producer 为 `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round055/sparse4_first_scalar_probe`。核验开始与结束均检查同一冻结字节链，共 108 个 pin 条目：full manifest/report/11 个 full 白名单，C055 precheck manifest/18 个 local/22 个 external，B055 precheck manifest/11 个白名单/23 个原输入，13 个编译头文件、binary、3 个连续性/接口源文和 2 个本审计源文件。无 mismatch。

| 对象 | SHA-256 |
|---|---|
| C055 full_manifest.json | 01ee335a507d48c8292da462eee1da4b796621249703f0c52a3a61c5356b3a54 |
| C055 full_report.md | f1292e18905659d963564339f33be25b54ea9d306c8c38b5fd730e2a8cb201df |
| C055 precheck_manifest.json | e5a61613fa0af7d235c5bd8017e8a9b1171b3577d06509c31ccb75f74088892c |
| B055 precheck manifest.json | c1171119bbefd07f9120633791d13d5c76dbc43a43034a5481462329f392e1ab |
| sparse4_probe.cpp | d44781595a04e8be9a4292e724fa7f3b3d63703494c12a2790399c63859c5af2 |
| native_input.txt | 610b5e73a07e354ab644ef27c2653bdf39ec9341a04dc3a2577cb502c8034174 |
| binary | 13a4b2ee9067d101859b83afdae48df8b0e6028f9a8fd26cb1b0c0378263c98f |
| native_full.stdout.json | b0a9a0adc8702c58db72b079e6024109eb31315cbd8408298370612b06b9515f |

full_release 参数与 precheck 唯一点/时间一致，时间早于 native 启动。full_execution 与 tool_session_receipt 交叉一致：attempt 1，native PID 59613，tool session 70248，launch chunk f43ca9，terminal chunk ea43f8，exit 0，restarts 0。进程 wall time 为 91.59672154099098 秒，kernel time 为 91.5061 秒，低于 1,800 秒上限；stderr 为空。审计时 PID 59613 已不存在。

这些历史执行事实以冻结回执与源/输出绑定为证据；本审计没有重新观察过去的操作系统事件，也没有重新运行 native。tool session 身份来自原 tool receipt，主代理可与原实际工具调用记录核对。当前 PID 不存在只证明审计时未运行，不能单独证明历史 exit 0。

## 2. 全部保存 checkpoint 与字段计数

逐行检查 61 条记录：前 60 条依次完成 32、64、…、1920 个 images；每一行 paired terms=images×13229，两个字段各 DD5=paired terms，logical terms=24×paired terms，total DD5=2×paired terms。elapsed 严格递增，unique-node counts 非递减；最后一条 computed/complete=true 与最后 running 记录及 native stdout 一致。没有以抽样代替全 ledger 检查。

| 字段 | 完整计数 |
|---|---:|
| paired image/cone tuples | 25,399,680 |
| G logical contributions | 609,592,320 |
| H logical contributions | 609,592,320 |
| G DD5 evaluations | 25,399,680 |
| H DD5 evaluations | 25,399,680 |
| 两字段 DD5 总数 | 50,799,360 |
| zero-node paired terms | 259,937 |
| confluent paired terms | 84,987 |
| unique exact nodes | 303,744 |

zero-node 和 confluent 的 logical counts 分别为 6,238,488 与 2,039,688，与 point_comparison 的乘 24 记录一致。这两类可以相交，不将其相加当作互斥分区。完整 ledger 支持已审 source 的全遍历计数；ledger 本身不是序列化 partial sum 或独立重算每项值的证书。

## 3. 独立 dyadic 区间和精确阈值证书

审计未调用作者 verify_full_saved.py。使用本机 Python 的整数与 Fraction，把 G、H、K 的原生 exponent -511、-508、-512 统一到单位 2^-512。若三个原生 mantissa 区间分别为 [g_l,g_u]、[h_l,h_u]、[k_l,k_u]，则保存差区间在该共同单位下为

\[
 [2g_l-16h_u,\;2g_u-16h_l].
\]

实际整数核验得到 k_l 比其下端小恰好 1 unit，k_u 比其上端大恰好 1 unit。因此保存的 K 区间包含完整保存 G/H 区间的差，不只是与它相交。三个区间端点有序，G>H>0 且 K>0。

阈值核验直接使用整数不等式

\[
 105600\,k_u<2^{512}.
\]

其右减左的正整数差为

```
772551557857646814812454158317104870941172593219457297001296516067216531632909835027357328522271508133213982222517101144374078728178306870152271085920896
```

这与阈值间隙下界的精确分子相同，分母为 105600·2^512。完整共同单位整数端点、约分后 exact Fraction 端点和区间宽度见 results.json；无需依赖显示小数判断 strictness。

下列 50 位小数由整数 floor/ceiling 除法独立外向舍入，未使用作者 Decimal 实现：

| 量 | 下界 | 上界 |
|---|---:|---:|
| G | 0.00001681892544263902088994249411166539840647790147 | 0.00001681892544263902088994249411166539840647790148 |
| H_T | 0.00000789486793799146978863105303299382098826704766 | 0.00000789486793799146978863105303299382098826704767 |
| K_T | 0.00000892405750464755110131144107867157741821085380 | 0.00000892405750464755110131144107867157741821085381 |
| 1/105600−K_T | 0.00000054563946504941859565825589102539227875884316 | 0.00000054563946504941859565825589102539227875884317 |

这里的严格算术证书验证保存输出的区间关系；把区间解释成目标 scalar 的 enclosure 仍依赖已接受的完整 fan、F052 scalar 定理、冻结 engine 与 Arb 的区间语义。没有用未经证明的浮点近似替代这些前提。

## 4. 负代表、归一与真实 S2 的继承

B055 precheck 已独立检查完整 frame/fan 输入、46,080 个群作用、全部 11,476 个 pair，以及 native 正负代表映射，本次只继承其冻结 PASS，不重新枚举。真实最小 pair 距离平方为 363/4375>0；x0 严格 off 新 frame 的 S2，所以 G 的 finite-part functional 确为有限 Green 值。

native 的 `y<neg` 取首非零坐标负代表，独立 checker 内部取正 canonical representative。H 的 atan 数据为奇函数数据；G 的非零 sign·pi/2、零 value=0、高阶 constant jets=0 也满足相同反射规则。五阶 DD 同时取反全部节点时，阶数因子与奇性因子给 (-1)^(5+1)=1；normalized confluent jets 的反射同样成立。因此选正或负代表所得每 cone scalar 相同，无额外负号或 factor 2。已包含稳定子及 global-negative 配对的 weight 24 只乘一次。

G 已含 pi/2，故 G 和 H 都除以 (2pi)^6=64pi^6；纯 sign-functional 的 128pi^5 形式只是等价写法。共同固定 T、512-bit Arb、zero/repeated 节点处理均由未变的预检 source PASS 覆盖。本次不推广到梯度/Hessian。

## 5. off-S2 连续性到 essential 固定 T 障碍

令 Q=[-1,1]^6，K_T(x)=G(x)−H_T(x)=∫_0^T p_t(x)dt，并令 m_T=ess inf_Q K_T。

已审 F052 定理给出 G 在 R^6\S2 上连续：小时间密度在离 S2 正距离的邻域上一致有界，大时间有可积 Ct^-6 控制。对于固定 T>0，H_T=∫_T^∞p_t dt 由固定 t 连续性和同一可积尾界连续。因而 K_T 在 x0 的某个 off-S2 邻域连续。这仅使用已证明的 scalar 连续性，不需梯度界、Hessian 界或任何新数值点。

x0 的前三个坐标等于 1，确在 cube 边界，不能只靠孤立 boundary value 否定 essential floor。不过 Q∩B(x0,r) 对每个 r>0 都有正的六维 Lebesgue 测度：可向前三个坐标的 cube 内侧微移，取得仍在球内的 interior 点及其小球。

记已证 gap=1/105600−U>0。连续性给足够小邻域，使其中 K_T(x)<K_T(x0)+gap/2≤U+gap/2<1/105600；该邻域与 Q 的交集有正测度。进一步让连续性误差任意趋零，得到

\[
 m_T\le K_T(x_0)\le U<\frac{T}{4096}.
\]

这闭合了 essential 接口，不只是“找到一个 measure-zero 反例点”。没有宣称一个显式定量邻域半径，也没有建立全域下界。

B044 已给该正保守 semigroup 的共同尺度 time-average 接口，其条件成本为 T/(64m_T)（当 m_T>0）；strict cost<64 必须 m_T>T/4096。本次证书恰好否定 beta=1/80、T=32/825 的这一必要条件。若 m_T≤0，则该正 floor 接口本来就不可用于此时间；无论哪种情形都不产生 strict-64 证书。

不据此声称 all-T、all-beta、旧 body32/原始 R_c 的排除，也不改变 C_n 或主项目已证上界。主代理保留最终证明审查与集成责任。

## 6. 本审计实际执行与封存

registration 在审计前登记最多 6 saved-only 配置、600 秒和所有零新增核约束。输入与审计源码随后冻结。独立程序 review_saved.py 一次实际运行：PID 61114，工具 chunk 788cb7，exit 0，2026-10-26T09:08:24.958420+00:00 开始，09:08:24.995510+00:00 完成，checker 耗时 0.03637558300397359 秒；准备、源审和报告时间不冒充 checker 时间。运行即时退出，未产生持续 tool session。

results.json SHA-256 为 `196bf1fe44ed64afe96be9fb1e34780e41df29dc5abeb91de52d8ed0fe283236`，由 execution_receipt 引用。执行源码拒绝有回执时重跑。复核封存使用 `python3 verify_sealed.py`，它只校验现存字节哈希与原回执引用，不重做 6 配置或调用 native。

裁决无阻断性缺口。保留原 full/precheck 所有文件与回执；本目录为独立 postrun 附件。未 stage、commit 或 push。

# B055 sparse4 首次标量探针独立预检：最终 PASS

**裁决：PASS_PRECHECK_STRUCTURE_AND_SOURCE_ONLY，最终源绑定完整，无已识别的阻断性预检缺口。** 可交主代理决定是否对冻结 C055 放行唯一一次 native。B055 自身不放行、不创建 release/token、不启动 native，也没有计算 G、H、K 或给出任何新 C6 界。

本次接续只完成原证据核验、完整源文审读和封存。原六项科学配置没有重跑；本次新增科学配置、矩阵操作、作者 writer、编译、production native、MATLAB 均为 0。历史 registration 的等待状态和 input_manifest 的 pending 字段原样保留，它们是执行时快照；本报告及 final_binding.json 记录已完成的最终绑定，不能通过改写旧快照伪造先验注册。

## 1. 冻结对象与证据链

唯一模型为六维真实 152 列 Cauchy frame，beta=1/80，x0=(50,50,50,45,27,0)/50，T=32/825。producer 目录为 `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round055/sparse4_first_scalar_probe`。

| 对象 | SHA-256 |
|---|---|
| C055 最终 precheck_manifest.json | e5a61613fa0af7d235c5bd8017e8a9b1171b3577d06509c31ccb75f74088892c |
| sparse4_probe.cpp | d44781595a04e8be9a4292e724fa7f3b3d63703494c12a2790399c63859c5af2 |
| native_input.txt | 610b5e73a07e354ab644ef27c2653bdf39ec9341a04dc3a2577cb502c8034174 |
| 原 build_record 所绑定 binary | 13a4b2ee9067d101859b83afdae48df8b0e6028f9a8fd26cb1b0c0378263c98f |
| B055 原 input_manifest.json | 782e34a9c2cd04b6979648fb37abb1456b881d570eb9a11a6c0afcde9dcdb0c4 |
| B055 原 results.json | 57f880ad248c0ae94a7d1fabbdc5461c1ecbcf965e42362e21089ab5f4523f62 |

本机 Python 实际完成 96 项字节哈希核验：原输入 pins 23、原审计源 pins 4、receipt 对 input/result 的引用 2、最终 C055 manifest 本身 1、其本地白名单 18、外部只读依赖 22、编译头文件 13、相同 13 个头文件与锁定 archive 内文件内容的比较 13。全部通过。source、input、binary、runtime 的现值与原检查时完全相同；最终 manifest 未替换原科学输入。F050 manifest 的 PASS 范围为完整数学 fan；F052 manifest 的 PASS 范围为 scalar theory 与 tiny references。两者分别是本审计继承完整 fan、primitive determinants 与 off-true-S2 有限部分公式的依据。

本次核验只确认现存二进制与已记录 build 的指纹相符，没有再次编译或形式化验证编译器。原执行及 build 的历史事实按冻结回执记录；不能把本次哈希复核称为再次运行原计算。

## 2. 原一次执行及六项检查结果

execution_receipt.json 记录 PID 55570，于 2026-10-25T22:35:12.939617+00:00 开始，于 22:35:13.232037+00:00 完成，exit_code=0，实际 checker 计时 0.291640749986982 秒。六个注册配置一次 PASS，22,952 个矩阵操作，低于 30,000 上限。计时为原 checker 主体，不包括准备、等待、后续人工源审或本次封存。收尾时进程查询确认该 PID 已退出；未将“all_handles_closed”字段单独当作进程已退出的证明。

1. 独立生成 32 个六支撑方向和 120 个四支撑方向，逐列匹配真实权重 frame。Gram 为 (2561/80)I；六坐标 Cauchy rates 均为 33。五个相邻交换与一个坐标符号翻转的 912 个列像验证群生成元不变性。
2. 对全部 11,476 对方向以 Bareiss 的 detGram3/detGram2 精确计算平方距离，不采用作者 normal-equation 展开。最小 pair 距离平方为 363/4375>0，唯一指标对 [39,47]；最小 line 距离平方为 2991/10000。因此固定点严格位于新 frame 的真实 S2 之外，非旧 body32 资格的继承。
3. 全 46,080 个 signed-permutation 作用得到 3,840 个不同像、每像 multiplicity 12；global-negative 配对得 1,920 个代表、每代表权重 24。首坐标为零的 320 个代表被保留。原 checker 用首个非零坐标正作 canonical line；native input 用字典序较小的负代表，两种选择只差每对统一取反，集合比较已覆盖全部对。
4. 全 478 条 rays 的真实 rates 与同时缩放逐条匹配；全 13,229 个 simplex 的 determinant/rate 权重逐条一致。primitive determinant 继承 F050，未新增 13,229 次矩阵计算。整个 native stream 的整数记录、顺序、cone IDs、point images 和 canonical digest 均已匹配。
5. 核定 scalar constant-functional 的非零 value=(pi/2)sign(b)、零 value=0、所有高阶 constant jets=0。H 使用真实 atan(b/T) normalized Hermite jets。两个字段均除以 (2pi)^6；这是 G 已含 pi/2 时的正确归一，等价于纯 sign-functional 除以 128pi^5。
6. 绑定完整 driver/input/build/runtime；每字段逻辑贡献为 609,592,320，压缩后的共享 tuple traversal 为 25,399,680，G 与 H 各作同样多的 DD5，两个字段共 50,799,360。共同固定 T、512-bit Arb、必要点阈值及失败语义均纳入审查。

**代表方向的审计澄清。** `results.group.representative_rule` 中的 “first nonzero coordinate positive” 仅描述独立 checker 的 canonical representatives，不能当作 native stream 的方向规则。实际 stream 的第一条代表为 `(-50,-50,-50,-45,-27,0)`，第二条为 `(-50,-50,-50,-45,0,-27)`；native `load` 强制 `y<neg`，故它取首非零坐标负。原 `review.py` 第 73–74 行已先把 input 映到正向 canonical line 作集合比较，再逐条断言 input 取负向，因而没有漏核这个差异。冻结 C055 report 第 48–49 行只写 “lexicographic first NONZERO coordinate, not the first coordinate”，未明确正负；应按源码理解为字典序较小的负代表。没有修改这些冻结文件。

取反不会改变两个字段的值或权重：H 的 atan 数据为奇函数数据，G 的 constant-functional 数据也满足同样的反射规则；对五阶 DD 同时将全部节点取负，反射的阶数因子与奇性因子相乘为 `(-1)^(5+1)=1`，重复节点的 normalized Hermite 数据也遵守该规则。因此正负代表互换后每个 cone 的 DD5 值不变，alpha 不依赖点代表，现有 weight 24 仍只乘一次，无新增符号或因子 2。此处只解释已通过的全 stream 与群检查，没有重跑群或新增数学配置。

这些值来自原冻结 results.json 与 execution_receipt.json；收尾没有重新生成它们。原 results 保留全部 pair 精确值，结果哈希与回执相符。

## 3. 完整源码审读结论

本次完整阅读冻结 sparse4_probe.cpp、launch_once.py，并对照已锁定 F052 proof_review.md。原 review.py 中的字符串断言仅是绑定辅助，不是任意 C++ 语义的独立证明。现给出实际人工审读要点：

- 输入严格限定全部 rays、cones、images、点分母、群权重及唯一 T；canonical int64 parsing、exact reduced Rat、int128 点积/交叉乘积保证节点排序、相等与符号由精确整数决定。Ray/rate/determinant 范围及 cone 内六个不同 ray IDs 均受检查。
- Node 构造 atan(b/T)，jets 由 atan 导函数的局部 Taylor 系数递推并除以阶数取得 normalized jets；重复节点使用该 Hermite 数据。G0 按精确符号取正负 pi/2，零节点保留零；所有重复高阶 G 数据显式赋零。无跳过零节点或重复节点的分支。
- sorted DD 的相等节点分支与 distinct 分支分开；distinct 节点的 Arb 差若包含零则抛错，不容差合并、不假装分离。缓存只复用节点 jet，不减少 tuple 或任一字段的 DD 次数。
- 每 cone 的 alpha24 先乘权重 24，再逐个 Arb 除六个 rates，未以 int64 构造巨大 rate 乘积；没有额外乘稳定子 12 或 global-negative 2。DD5 对 odd value/相应 Hermite 数据具有全体节点取负的不变性，故可用完整配对代表；该性质与已审 scalar 公式一致，不引入新配置。
- 全部 tuple 完成、精确项数一致且 G/H/K 有限后才输出 COMPLETE。归一使用 (2pi)^6，K 用 interval subtraction。输出给出整数 mantissa 与共同二进制 exponent，后处理用 Fraction 作精确阈值比较；double 仅用于耗时记录。
- launcher 要求主代理写独立 release，并同时绑定 B055 manifest 与 C055 manifest；重新核全部冻结 local/external hashes、严格 S2 资格和 source/input/binary。full_execution 一旦存在就拒绝第二次，checkpoint/token 存在也拒绝重启。它只允许一个 1,800 秒 child，超时只终止该 child，必要时 kill 并 wait；失败、非零退出、stderr、count 错误或不完整 checkpoint 均不输出成功比较。

封存时 producer 尚无 full_release、full_execution、full_token、full_checkpoint 或 native_full 输出；B055 没有创建这些文件。没有实际运行表现、耗时可达性、数值区间宽度或数值结果的预先保证。

## 4. 主代理后续边界与未完成工作

**未识别需修复才能进入主代理放行判断的预检缺口。** 主代理仍须独立接受本报告并决定是否放行；本报告和 manifest 中 production_launch_authorized_by_this_audit=false 不是 release。

若主代理放行，只限当前冻结 source/input/binary/runtime 下的一个 native。后续必须审核实际 process receipt、完整 checkpoint、两个字段计数、有限区间、输出 hash 及必要点分类。源码、输入或依赖变化使本次绑定失效；本 PASS 不自动覆盖变更后的对象。

本点 K_lower>T/4096=1/105600 仅通过该点的必要条件，不证明全域 m_T；K_upper<T/4096 可排除此固定 frame/T 的对应必要条件。launcher 对等号和区间重叠均保守返回 unresolved，不自动重跑。62.5 对应 T/4000=1/103125，但当前 launcher 只分类 64 的阈值。本报告不宣称 gradient/Hessian、全域 floor、C6 改善或最终证明完成。

## 5. 复核与封存

在本目录执行 `python3 verify_saved.py --verify`，只读核验所有冻结 hash 和原回执一致性，不导入作者模块、不执行 review.py、不调用 kernel；manifest 存在时同时核本审计白名单。final_binding.json 留存此次实际核验摘要和进程/生产文件状态的时点观测。manifest.json 绑定全部必要原证据、复核脚本、本报告和最终绑定记录；manifest 自身 SHA 通过交付消息传递，避免自引用。

本轮只新增封存与核验文件，不覆盖原证据，不改 producer，不 stage、commit 或 push。

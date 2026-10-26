# A090：真实首批八次覆盖试点的 Sol 保存证据独审

结论：**PASS_SAVED_ONLY_ACTUAL_COVER_REVIEW_WITH_RETAINED_REVIEWER_HARNESS_ERROR**。
接受本次精确冻结的原八次 native、数学导入、局部整盒结果、状态生命周期和运行时收据之间的绑定。
这是保存证据独审：没有新 native、梯度参考算法、Arb 指数/角点或旧 fixture 重跑。
首次独审程序在完成一次 core replay 后发生归档读取错误，原错误没有抹除；主明确批准的有限续审完成后续绑定，三组总额未增加。

审查对象是主 clone 的 `cube-dim-6-bounds/upper/round090/cover_batch001_actual/science_manifest.json`，
SHA256 `5b9f004437d34fc0fd58be76ef762f1a7f975fa7cc4ce0e61dac2a89cc02f82d`。
其每一条文件的字节数和 SHA 均核对。最初只读 source，收到主的 final-freeze 通知后才读取并哈希实际整批证据。
专用 no-checkout clone 是 `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-cover090-audit`；
全部写入只在本审计目录，未 stage、commit、push、MATLAB、递归分派或修改原代码/结果。
模型是 GPT-6.1-sol，遵守当前禁止新增 Astra 的协调修订。

## 真实数学结果与剩余范围

八个保留输入均为同一 `body32_unit_cauchy_v1`，同时间带 `B_[1/32,1/2]`，目标 `1/65000`。
八个中心的原 Arb512 区间下端均严格高于目标；这没有被当成八个整盒通过。
原保存 Arb384 证书中，仅 `r000.a0.n000008` 的整个五维面矩形 `[0,1/8]^5` 通过，
其余七个整盒未决并按原规则细分。通过盒的严格有理下包保存在 `core_results.json`，
来源是原 `verification.json`，约 `5.275e-5` 的显示值只用作定位，验收采用完整有理数。

| 原 token | 原 tool session | 原 wrapper / native PID | 整盒结果 |
|---|---:|---|---|
| r.a0.n000001 | 2739 | 19118 / 19120 | 未决，细分 |
| r0.a0.n000002 | 13324 | 19344 / 19346 | 未决，细分 |
| r1.a0.n000003 | 20429 | 19695 / 19697 | 未决，细分 |
| r00.a0.n000004 | 70210 | 20773 / 20775 | 未决，细分 |
| r01.a0.n000005 | 23541 | 21438 / 21440 | 未决，细分 |
| r10.a0.n000006 | 89741 | 22284 / 22286 | 未决，细分 |
| r11.a0.n000007 | 28784 | 23251 / 23253 | 未决，细分 |
| r000.a0.n000008 | 27356 | 23709 / 23711 | 整盒通过 |

最终 `sequence=25`，head SHA256
`30ffc83e8b7643250429cb5267c57a11d8131fc31123fd88abd165ffa8d28cf2`。
`passed_leaves=[r000]`；待处理队列严格为
`[r001,r010,r011,r100,r101,r110,r111]`；`empty_leaves=[]`、`inflight=null`、`halted=null`，无 pending snapshot。
原八个 distinct started nodes、八个累计 reservation、八个已观察 native 和 serial9 一致。
主的最后一次 reserve 确实返回 `BATCH_BUDGET_EXHAUSTED_QUEUE_RETAINED`，
其 before/after head 都是上述 SHA；没有借耗尽预算清空队列、退款或隐式刷新批次。
恢复接口的一般性 action 文本不增加额度；真实 control 的预算分支已拒绝第九次 dispatch。

## 一次状态重放与原件数学绑定

本审检查真实 068 core、Context、math_contract 和冻结 063 几何转移源代码，
保留 063 持久化 API 原拒绝。068 的 genesis 固定为完整 `[0,1]^5`，
而非任意自洽 snapshot；每个 typed event 从已验前态生成完整后态后逐字段比较。
根、闭分割边界、排序域 tightening、兄弟和唯一 unresolved frontier 均受原确定性转移约束。
原 token/input、handle、已导入 observation、passed leaf 和累计预算不得删除或回退；
没有 reset/refund/requeue/reissue/new-batch 事件。

第一组只调用一次原 `V.recover`。唯一适配是 reviewer `SavedOnlyContext.verified_observation`：
替换原会重新调用 Arb 的 reconstruction 为保存原件的精确有理核查；
原 core/geometry/release/handle 校验不改。程序把 `M.reconstruct` 和 `checker.certify`
替换为一旦调用就报错的内存 guard，因此没有隐含的 oracle 或 Arb 重放。
原8个 import 触发了全部8份 saved verification 检查，共256个带标签角点。

每份检查绑定 reservation、唯一 attempt 目录、输入哈希、intent、精确 native 命令、
model/definition/window/center、原 release、source/binary/fan/runtime SHA、原 native bytes/stderr、
原 wrapper/native PID/worker PG 和 observed tool session 表。
六个物理梯度、band、两个端点和 native Euler 都由原三字段 dyadic 解码；
检查 `x·grad(F)+5F+T*p_T-a*p_a` 的独立有理区间含零，然后核对保存的完整六分量与 Euler 字段。
这是一条标量一致性关系，不能声称重新独立确定了六个梯度。

每份 native 都保留 `483840` computed sign-pairs / `967680` logical terms。
原零节点/重复节点计数真实非零，全部核对原输出、允许范围和偶性；源码的零节点与 Hermite
重复节点路径在完整循环中处理，没有把其跳过。本审没有重新枚举这些计数：

| 原 node | confluent_terms | zero_node_terms | unique_node_values |
|---|---:|---:|---:|
| r | 605520 | 367200 | 133 |
| r0 | 541920 | 301920 | 189 |
| r1 | 422208 | 226272 | 278 |
| r00 | 527760 | 282240 | 258 |
| r01 | 391488 | 197088 | 394 |
| r10 | 451152 | 254880 | 236 |
| r11 | 346824 | 236928 | 290 |
| r000 | 507840 | 282240 | 303 |

保存32角点证书表逐行核对全部标签、位移、梯度点积下端、共享中心值 clamp、
固定 frame32 二次项、指数表达式、有序区间、上下最小值、逐角点 passes_target 与整盒 status。
保存证书、checker hash、原 native provenance 和完整 verification 对象也核对，
不只读取 PASS 字符串。

本审没有重新认证指数函数数值本身；原 `check_saved_box` 也明确只认证精确几何/表达式，
不认证 exp。指数端点通过冻结原件、真实 verify/import 收据及已接受 A057/068 实现证明绑定。
八个新点没有重新运行 F061 的全梯度参考算法；F061 的有限已登记案例不扩张为本次八点全部独立 crosscheck。
source proof 和 provenance 继续明确 `independent_full_gradient_crosschecked_for_this_new_attempt=false`。
主090 report 中“065的两个有限crosscheck”是文字索引误写；真实 release/context/hash 绑定的是 F061
manifest `91b26762f71229fbaa4c011be12b5a678ce968afc298e23a2c5320403c51fbda`，按真实绑定解释，原封存报告未改。

## 真实运行时、所有权和恢复状态

068 main release `8102bf03…` 与082 runtime release `a316a087…` 分别绑定真实独审 acceptance，
082 runtime preparation 为 `a6894939…`；原071的真实拒绝没有改成PASS。
检查了 fixed soft170/group178/parent179/watchdog179.5/supervisor179.8 及清理/收据仍受 deadline 的源码。
本批没有触发 body/native timeout，不能把已有合成故障证据说成本批真实 timeout 实验。

42个原 production operation = 8 native + 34 control，原5文件的 UTF8 字节/hash与
冻结 `runtime_receipts.json` 逐项一致，原物理文件也逐字节一致；另外读取其原 stdout/stderr。
每个 operation 的原 supervisor/worker/watchdog handle、armed事件、task_result、release、
正常组清理ACK、direct worker/watchdog bounded reaps 和 completion 均匹配。
42个 task_result 的 logical_exit_code 都为0，无 timeout/error/receipt_failure 或 unresolved cleanup。
全部 completion 的清理含收据时间低于180秒，最大为 `8.46855154199875` 秒。
原 native_exit_code=0 且原 stderr 为空；worker末尾−9是正常清理身份锚点，不能改记成 native 错误。

原session表来自主的实际 exec yield 与同session write_stdin，含原 launch/completion tool chunk。
本审能交叉核对这些观察记录与落盘原件，不能冒充独审当时亲自观察了历史工具调用。
没有重发任何原session或native。

第三组只对原件明确记录的134个任务 PID 做了一次 `/bin/ps -p ...`，exit1且stdout/stderr为空，
观察时间 `2026-10-26T15:03:53.663814+00:00`。没有按名称扫描或发信号。
这表示所列原PID当时不在进程表；不推出 supervisor 回收了所有孙进程。
所有输出继续保持 `all_descendants_reaped_claim=false`。

## 失败与有限续审账本

主原先的两次setup错误分别保留：系统Python缺flint；预检对独审manifest错误使用
`artifact_sha256` 键而真实core audit用 `whitelist_sha256`。两者发生在原 registration/release/
state transition/native 之前，0 native、0 transition，无已记录实际耗时则不补造。
它们与本审以下归档读取错误分开，均不改称科学重跑。

本审第一组 `check_saved.py` 原样保留，SHA `dbf84d54…`，tool chunk `06194a`，exit1。
堆栈到达 main225 → ref57，将已验hash的Markdown报告JSON解析而失败。
根据该精确源码路径，191行一次 recover 已返回，192–222行的8份检查、256角点、8/8、
1passed/7queue、原工具/PID和source pins断言已通过；但本次没有写出最终 core_results。
不能把这次整体程序说成首跑全PASS，也不能补造当时未写出的输出。
原 failure receipt 记录 wall `0.23348362499382347` 秒、CPU `0.21158` 秒、PID27220；
工具包自己的 wall字段为 `0.148163` 秒，分别保留，不互相覆盖。

主明确批准剩余有限续审并保持原3组总额。第二组只完成剩余source/review字节hash绑定，
不再 V.recover、不重跑已通过的数学核查；`core_results.json` 的汇总字段从冻结原件读取，
清楚标注第一组证据范围和未产生最终输出的事实。第三组完成 runtime/预算/准确比较。

| 审核接口组 | 实际结果 | wall / CPU 秒 | native / Arb / recover |
|---|---|---|---|
| 1 原saved-only core | core已过后发生harness归档错误，exit1 | 0.233483625 / 0.21158 | 0 / 0 / 1 |
| 2 明确批准的source-binding续审 | PASS，exit0，chunk1bd3cc，PID27850 | 0.044889167 / 0.028897 | 0 / 0 / 0 |
| 3 runtime+准确比较 | PASS，exit0，chunk70a3a4，PID27887 | 0.067363334 / 0.04133 | 0 / 0 / 0 |

每组hard alarm179.5秒、CPU limit20秒，包含清理/输出阶段；累计3组、原core recover仅1次。
没有第四组、自动重试、native/gradient/Arb重算、旧基线/fixture重跑、源文件写入或MATLAB。
本报告和manifest仅为行政归档。主另报的科学封存脚本 parseSyntaxError 发生于任何程序运行前，
不属于本审接口组，也不与上述两次主setup或本审harness错误混账；本审未修改科学seal。

## 准确比较与未解决缺口

第三组读取并保存075当前完整准确有理区间，核对其证书 SHA
`dec62274abcf8204cf698cd7665dd6c754db3fee87a55beb5bfa1ece676e2817`，没有重算旧下界。
精确Fraction比较验证当前075完整下端高于冻结003完整区间上端，
当前完整区间下端严格高于显式有限L1阈值 `2073015624830/1000000000000`，
当前完整区间上端低于64。完整区间而非显示小数存于 `runtime_results.json`。
当前内部下界仍为075，当前/public已核上界仍为64；本次090不登记C6新界。
`dual_checkpoint=null`、`post_dual_completed_rounds=0`，没有把本次局部成功计为双端突破。

本批不是全ordered-face cover；七个队列区域仍缺足够整盒证书。
因此没有global g、竞争性完整I、严格正S，也没有general L1/全部正尺度的新C6上界。
090试点及068/082 release或局部盒本身不替代这些全域数学前提。
原批额度已经耗尽；未来工作需新有限登记与合法机制审查，不能清空原队列或静默续跑。
审计已停止计算，最终材料由主集成。

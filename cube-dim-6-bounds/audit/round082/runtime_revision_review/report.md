# A082 运行时最小修订独审

结论：**PASS_INDEPENDENT_RUNTIME_SUPPLEMENT_REVIEW**。接受本次封存源码的运行时
修订及下述证据范围，可以交由 main 绑定实际 release。独审本身不创建 release，
不启动 production；八新 node / 八 reservation 上限及原 C068 队列不变。
这不是实际 body/native 验证、全域覆盖、global g 或新的 C6 结果。

锁定作者目录：
`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round082/runtime_revision/`。
preparation_manifest SHA256：
`a68949397da1dc7825c1d4979ca1fcc096fe5e49b81884d36c4fb9170694af33`。
17 个被 seal 绑定的 regular artifacts 及全部 readonly pins 在本审计前后均一致。
原 A071 失败 seal、C069/071 源和 C068 core/math 均未修改。

## 补丁判断

原 A071 的真实满管道配置暴露了 envelope 重复 killpg 的未捕获 EPERM，导致收据
缺失；当时没有观察到 180 秒超限。本次不是撤销该历史发现，而是接受新目录中的修订。

新 supervisor 只从原 watchdog 事件管道、且 PID 与原 watchdog 相符的消息接受
终止 ACK。已有该 ACK 时省略重复 kill，因此修复了观察到的重复动作路径。
没有 ACK 时，kill OSError 被记录为 termination_confirmed=false；不向 watchdog
发送虚构 K，不凭 ps 观察升级为已终止。原 worker 在可能仍需要组操作时不被先行
reap，避免提前释放 PG 身份锚点。watchdog 确认或退出后才尝试相应有界回收。

receipt 写入异常产生 receipt_failure，明确 final_receipt_written=false；若写入
本身挂住，仍保持独立 watchdog 和默认总 alarm。若所有 checkpoint 写入都故障，
只能保留已经落盘的 start/handles/原 reservation，不能保证最终收据一定存在。
这是明确的部分结果边界，不应把最终收据的缺失误称为成功或自动重新发起的许可。

worker、watchdog、native_task、control_task 四个文件与原已审层逐字节相同。
production 固定 soft170 / group178 / parent179 / watchdog179.5 /
supervisor179.8；native 继承 worker 的独立 PG。native cleanup 和 supervisor
清理/回收均受原绝对期限限制，清理及收据阶段没有取消总 alarm。本文接受正常 OS
调度下的有界终止机制；不声称对整机暂停、内核不调度或故意 setsid 逃逸作数学保证。

新增 record_native_stop 的静态检查：它重读原 C068 state，要求当前 inflight token
匹配，使用该 attempt 的 input hash；已有 failure 不覆盖，新记录绑定原运行时
收据和已有原件。缺 original_process 时使用 UNKNOWN_ORPHAN_STOP 和
native_launch_observed=None，符合冻结 Context.verified_failure 的接口；有原件时
采用允许的 timeout/error reason。它不修改 replay state，不退款、不重新保留 token，
不 relaunch；main 仍须处理原 token 的 HALT。该生产 helper **没有实际执行**。

production_gate 的 release、source/receipt/独审文件哈希绑定仍保留，并再次验证
C068 main release。当前 runtime_release.json 不存在。native/control 的入口、
guard 与同 model/window 的既有数学检查未被修改；本轮未重跑数学或真实 kernel。

## 作者已有证据：只读核验，不重跑

merged_fixture_evidence.json 内 44 个原始 regular 文件的文本字节 SHA 全部匹配，
存在的原文件也逐字节一致；四个 observation 与 results 的相应对象完全一致。
除读取 status 外，还重新检查了关键字段：

- 真 cleanup 满管道留读端，原 watchdog 的 GROUP_SIGKILL_SENT 持久化，已 ACK，
  supervisor 未重复 kill，before/receipt/completion 完整。
- 真实 child/grandchild 同 worker PG；其后均不运行，无关自有 PG sentinel 存活，
  随后由作者仅清理原 sentinel。
- supervisor 的真 receipt 满管道由独立 watchdog SIGKILL，只有 before checkpoint，
  没有伪造 final receipt/completion。
- 注入未 ACK EPERM 被保留为未确认；真实 FileExistsError 写收据异常产生失败记录，
  没有误标 final_receipt_written；缺 release 和参数覆盖请求被拒绝。

作者实际四配置、一次 exit0、harness 2.453384375 秒，只作为封存的已有证据。
本审计没有重新调用作者程序或旧 normal/soft/math。

## 三项新独立组合故障

生产运行时代码按 seal 字节复制到独审所有权目录，只有 fixture_stubs.py 被替换为
新的独立注入模块，并更新副本对应 pin。原作者源码不改，策略/时限不改，没有虚构
runtime_release。测试组合与原作者四项不同。

| 组合 | 原 supervisor / worker / watchdog | 实际结果 |
|---|---|---|
| A1 错误管道伪 ACK + 真 cleanup 满管道 | 93738 / 93739 / 93740 | PASS；exit1，0.879154084 秒 |
| A2 watchdog 早死 + cleanup 满管道 + receipt 目录异常 | 93745 / 93746 / 93747 | PASS；exit1，0.969015541 秒 |
| A3 watchdog 早死 + supervisor receipt 满管道 | 93751 / 93752 / 93753 | PASS；默认 SIGALRM exit−14，1.619232792 秒 |

A1 从 worker result 管道发送与原 watchdog PID 相同的伪 GROUP_SIGKILL_SENT；消息
确实进入事件记录，但没有成为合法 ACK，也未导致提前清理。随后原 watchdog 在
0.85 秒策略截止后发出的真实 ACK 得到接受，supervisor 不重复 kill，收据完整。
这检查了管道来源，而非仅检查消息内 PID/status 字符串。

A2 在真实 cleanup 阻塞标记后，仅向原 watchdog 93747 注入 SIGKILL。supervisor
仍在固定截止附近完成自己的 worker 组清理，group_ack=false，没有伪造 ACK；写入
receipt.json 目录触发真实 FileExistsError，receipt_failure 明确 false。完成记录
中两个直接子句柄均被回收，watchdog returncode−9。

A3 在 supervisor receipt 满管道标记后，仅向原 watchdog 93753 注入 SIGKILL。
worker 已被原清理路径终止并回收，supervisor 继续真实阻塞，随后由自身固定
1.6 秒默认 SIGALRM 终止。外部计时包括 Python 启动，故为 1.6192 秒。未生成最终
收据或 completion，符合部分结果的契约。证明了外部 watchdog 失效时总 fuse
仍覆盖收据阶段；不是实际生产 native 的 timeout 测试。

三项新组合一次运行均 PASS；continuation harness PID 93737，exit0，
3.527681959 秒。没有额外自动重试或第四项组合。所有 12 个本审计已知原 PID
（包括下述 setup 失败的三 PID）在后续只读 ps 中均无行。这个观察不是“全部后代
由 supervisor 回收”的断言；all_descendants_reaped_claim 保持 false。

## 首次审计脚本故障和明确授权续跑

原 check_review.py 的 await_marker 在文件刚创建、尚未写完整 JSON 时即读取。
第一次调用在 handles.json 上遇到 JSONDecodeError，harness 0.086721583 秒 exit1。
已经启动一个 envelope，但没有 entered/forged/cleanup 标记：**setup started=1，
组合 fault entered=0**。这属于独审 harness 的 setup 故障，不能据此拒绝 C082。
初始 supervisor/worker/watchdog 为 92738/92739/92740，没有创建 sentinel。
原脚本、登记、输入 pins、失败 observation/results 全部保留。

主审明确授权一次 continuation：先只读确认原三 PID 均不存在，再于新文件/新输出
中仅修复 marker 读取，在原绝对 deadline 内等待完整 JSON，不延长期限、不改被测
源码或注入条件。supplement_registration.json 绑定原失败哈希及授权；然后一次
完成原登记三组合，首失败停止规则不变。原 handler 在失败时有仅针对其新建进程组
的条件清理分支；原 observation 未记录该紧急分支，不能把空 reviewer_signals
列表解释成“首次故障没有发送信号”。后续只读观察本身不发送信号。

总账：两次 harness 调用；首次一个 setup 中断、零实际组合注入；一次明确授权续跑，
三项实际组合全部通过；总 envelope launches=4。body/native/Arb=0，旧数学及
normal/soft 重跑=0，自动重试=0，无观察到的 180 秒超限。

## 验收与归档边界

本独审接受 **a6894939… 对应的精确封存版本**，不追认旧失败版本，不认证所有
未来输入已重算。主审须另行写入合法 acceptance/runtime release，并继续使用冻结
C068 的来源、数学、队列和八 node/八 reservation 约束，才可运行真实首 batch。
本轮未执行实际 body native launch/kill、完整 g cover 或新 Cn 计算。

manifest.json 是明确归档白名单。两次源码/登记和证据均独立保存；临时 runtime
副本排除，必要原始输出已合入对应 observation（完整文本及 SHA）。所有写入仅在
本 audit 目录。086 保持 WIP，与本 PASS 不混账。未 stage、commit、push 或 MATLAB。

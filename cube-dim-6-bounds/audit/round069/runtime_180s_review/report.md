# 运行时补充独审：首个真实阻塞配置失败后封存

目录保持 `audit/round069/runtime_180s_review/`；主审安排本对象归档为第 071 轮。

结论：**REJECT_RUNTIME_FOR_RELEASE**。本结论针对运行时错误处理及未完成的验收覆盖，
不撤销 A068 的 core/math ACCEPT，不是 global g 证书或新的 C6 上界。
没有观察到程序超过 180 秒，没有运行 body/native kernel 或 Arb，也没有建立真实覆盖。
生产 gate 继续关闭；只有主审能够作出后续 release。

## 固定输入和执行范围

作者源目录为
`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round069/runtime_180s/`。
锁定的 preparation seal 为
`ea1bdbbfe261b07f675c82da4f62e71eb06188496e79e5555b7ab747f8dba9df`。
作者 report 为
`c145579d1c2c8a3ca64eca47081ff175c9b93b7d01dafcd1c483bde14f19b39f`，
原始 execution receipt 为
`f4ea654c3dc9961e3a2e948c50e2786bb30d0c532a0136ae564391316b408004`。
54 个作者 regular artifacts 以及全部 readonly/source pins 在独审前后均匹配；
068 源、旧审计、原有登记和作者失败记录均未改变。

作者只实际启动三个配置：normal/soft 通过，第三个 cleanup FIFO 注入提前产生
BrokenPipeError，未到独立 watchdog 截止。作者第 4–6 配置未运行，不能算通过。
本审计没有重放这些 normal/soft、旧数学配置或作者 check_runtime.py。

独审先固定 registration.json 和 input_manifest.json，再写定审计代码并生成
pre_run_manifest.json，然后只启动一次 review_runtime.py。登记上限六个新配置，
程序上限 180 秒、审计 harness 默认 SIGALRM 170 秒、每个被测 envelope 外部等待
最多 4 秒、stub 独立默认 alarm 4 秒。第一项失败即停止，不重试、不扩展参数。

写入仅发生于本审计目录。运行时生产代码逐字节复制到本目录的临时副本，唯一有意
替换的是 fixture_stubs.py（采用独立注入 independent_stub.py）及其副本 input pin。
envelope、watchdog、worker、runtime_support、native_task、control_task 和
check_runtime.py 的副本 SHA 与封存源码一致。没有虚构 runtime_release。
副本只是确保测试输出留在独审所有权目录，未改变被测运行时控制逻辑。

## 阻断发现：已实际允许到达的清理路径抛出 EPERM

**R1_retained_reader_full_pipe 实际运行，但完整验收失败。**

新 stub 收到原 worker 的 soft TimeoutError 后建立匿名管道，始终保留读端而无人
读取；先用非阻塞写填满管道，再切换到阻塞写。实际记录为：

| 观测 | 实际值 |
|---|---|
| supervisor / worker / watchdog PID | 84002 / 84003 / 84004 |
| worker PGID | 84003 |
| 进入真实阻塞写 | supervisor 起点后 0.3559547500 秒 |
| 已填充且读端仍保留 | 65536 bytes；retained_reader=true |
| cleanup 返回标记 | 不存在 |
| supervisor 退出 | exit 1；外部测量 0.8953347910 秒 |
| before_receipt / receipt / completion | 均未生成 |
| 随后按三个原始 PID 查询 | 均无进程行 |
| harness | PID 84001；exit 1；0.9407727500 秒 |

封存 envelope.py 第 108 行无条件再次执行
`os.killpg(worker.pid, signal.SIGKILL)`。它在本次实际路径抛出：

```text
PermissionError: [Errno 1] Operation not permitted
```

该段只捕获 ProcessLookupError（第 110 行），不处理此次实际发生的 PermissionError。
因此代码在第 134 行 before_receipt 写入之前退出。后续 harness 因找不到
before_receipt.json 抛 KeyError；**KeyError 是缺失收据的后果，根因证据是被测
envelope 的原始 stderr**，两者均保存在 R1_retained_reader_full_pipe.json。

本次真实挂起得到终止，且早于 worker 的 1.1 秒 fallback 和 supervisor 的
1.6 秒 fallback；时序与 watchdog 0.85 秒组截止一致。不过 watchdog 事件留在
supervisor 内存而未写成收据，故不将 GROUP_SIGKILL_SENT 或 ACK 标为已持久化观测。
不对 EPERM 的操作系统内部原因作超出证据的断言；也不把本次失败表述为 180 秒
时间约束违规、真实 native 失败、数学证明失败或已观察到遗留活进程。

独审 observation 中 reviewer_emergency_cleanup=true 表示进入失败处理函数。
此时 proc.wait 已取得 supervisor exit 1，函数中的发信号条件均为假，只有重复读取
已完成句柄。R1 没有 reviewer 实际发送的 kill 信号；execution_receipt.json 单列
此区别。后续只读 ps 检查也未杀任何进程。不存在“审计先杀死后宣称运行时清理”的情况。

修订需在新目录进行，不能改原 seal。一个可审查的方向是：若已收到独立 watchdog
的有效终止 ACK，避免再次向同一已结束组发送冗余 kill；其他 kill 错误应生成明确的
未确认/失败记录，并保留独立时限及身份锚点。**不能简单吞掉 EPERM 然后宣称组已
清理**。具体实现和重新验收由主审组织。

## 静态结构核验及其限度

以下结论来自封存源码逐行检查，不把尚未运行的 R6 当作动态证据。

1. runtime_support.py:11–14 固定 production soft=170、group=178、parent=179、
   watchdog=179.5、supervisor=179.8、worker_fallback=178、reap=1。
   fixture 使用固定的短时策略，CLI 没有延长生产时限的参数。
2. envelope.py:4–6 在后续导入前启动默认 SIGALRM；第 35 行按同一个 START
   调整总 fuse。代码没有在清理或收据写入前取消它。watchdog 同样设独立总 alarm。
   这是正常 OS 调度下的执行机制，不是关于整机暂停或内核不调度的数学定理。
3. envelope.py:75–81 创建 own-session worker，而 watchdog 位于该 worker PG 之外。
   worker.py:17 核对原 PID=PGID 和父 PID；只有 watchdog armed ACK 后才向 worker
   发 gate。native_task.py:71 明确 start_new_session=False，native 继承 worker PG。
   原 worker 在逻辑结果后保留为未回收的组身份锚点，直到清理，避免先 reap 再 kill。
   不声称能清理故意另建 session 的恶意子进程；真实 body 未在本轮启动。
4. 原 native wait 使用距 soft 截止的显式 timeout；native cleanup 通过
   bounded_reap 使用不超过 1 秒且受 group 截止限制的等待。supervisor 的 ACK/reap
   也受绝对截止约束。worker 的 signal.pause 循环被 worker/watchdog fuse 限制；
   不把这种保留身份的等待误判成可无限延长整个程序的路径。
5. 生产入口先调用 production_gate，再创建 operations 目录或启动 worker。
   该 gate 先核验固定 pins，再读实际 runtime_release.json。当前没有该文件，
   因而按源码不能进入实际工作。worker 的 production 分支再次 gate；native/control
   的 main 首行 require_worker，要求生产 guard、原父/组身份、未过 soft 和相同 release。
   两个 standalone 入口会主动拒绝。本轮没有动态运行这四个拒绝 probe。
6. 后续有效 release 还须绑定本 preparation seal、main 的 runtime acceptance、
   独审工件哈希和 C068 main release；再调用 C068 Context.verified_release。
   native/control 仍使用冻结的 C068 replay/math contract 和八 node/八 reservation
   规则；本轮不重验这些已接受的数学配置，不把源级绑定当作未来每个实际输入的验证。
7. receipt 写入、异常路径本来仍受总 fuse 约束；但已观察到的未捕获 EPERM 会使
   正常硬停止路径提前失去审计收据。因而不能给出所要求的完整运行时 release PASS。

## 完整配置总账

| 配置 | 状态 |
|---|---|
| R1 真满管道，保留读端，不排空 | 实际启动；完整验收 FAIL；收据前 EPERM |
| R2 早返回留子孙进程，并保留无关 PG sentinel | NOT RUN |
| R3 supervisor receipt FIFO 无读端 | NOT RUN |
| R4 receipt.json 目录导致写入异常 | NOT RUN |
| R5 真满管道时故障注入杀死原 watchdog | NOT RUN |
| R6 缺 release 的四个拒绝 probe 与 AST 断言 | NOT RUN；相关源码已人工静态审查 |

实际 1/6，完整通过 0，失败 1，未运行 5。native=0、body kernel=0、Arb=0、旧
数学配置重跑=0、normal/soft 重放=0、自动重试=0。program/tool 运行收据与注册哈希
已封存。原作者部分结果仍保持 3/6 实际、2 通过、第三注入失败，不能与本次合并成六项通过。

未执行的后代 PG、无关 sentinel、FIFO 收据、收据异常和 watchdog 死亡路径仍是
后续必要验收范围。本报告不允许生产，也不启动新的补测。真实 body native 的
launch/kill/timeout 路径从未在本轮尝试。C068 core/math ACCEPT 的既有范围保持；
runtime 必须修订、重新封源并独审，再由 main 决定实际 release。

## 归档与复现边界

最终 manifest.json 给出唯一归档白名单，含登记、输入 pins、预运行代码 seal、
独立 stub/harness、单次开始标记、R1 合并原始证据、results、执行收据和本报告。
临时源副本为本审计可再生成数据，不属于归档工件；没有修改、stage、commit 或 push。
R1 的合并证据含每个实际 regular fixture 文件的完整内容和 SHA，另保存 supervisor
stderr/stdout。不得尝试读取 FIFO special file。

复现仅可在主审授权的新有限登记与新输出目录下进行。原 review_runtime.py 使用
排他 run_started 标记，拒绝原目录重放；本轮不通过覆盖标记来继续其余配置。

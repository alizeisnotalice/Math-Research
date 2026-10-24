# 第45轮准备：第44轮候选的两套有限全域核验

状态：第44轮solve已实际终态并导出候选，本文件据保存结果绑定输入。尚未执行本轮checker；由主侧另审wrapper后安排唯一运行。本轮不做LP、参数修复、模型续接或新点诊断。当前已证B1仍不变，不能先称新C2上界。

## 1. 实际发现结果与证据边界

第44轮build实际完成7092个新T AS模型、15336个方向区域记录、20448条新连续行；旧三父格16272行/12204模型原样复用。唯一solve外包装terminal exit0，42.64252837499953秒，内部42.373197167005856秒；1次MATLAB LP，exitflag1。

保存浮点LP成本3.590862617332495，原尺度最坏浮点slack约−8.6841e−12，不能将其误读为真实生成元负值。按事前固定λ=1023/1024收缩和唯一结构量化后，36720条准确连续行零失败，最坏项在父格(52,1)、node658、piece1、Bernstein(2,0)，准确slack为正。两批共103148个缓存行值的原/加强K和A目标均无超容差违反，仍只属数值筛查。

候选的保存成本约3.590862654958542；exact_cost.json给出严格上下包围、improves_baseline=true、below_user_target=false。仅准确质量/结构和有限条件通过，尚不能将该数字写为全域一般C2已证界，更没有达到17/5。

## 2. 固定输入及运行前准确gate

- plain候选：`../round25/candidate.json`，SHA
  `7081cfef6fb47e5d47c95d3740c7a6ba48d48bb6197928a84d8cb6340c4166d5`。
- structured核：`../round25/structured_kernel.json`，成本证书绑定的SHA
  `00c50260612ec3e8b15e810af2fe62269c10a09db0a7cd136c624df8a90bc02f`。
- 准确成本：`../round25/exact_cost.json`，SHA
  `0120e7a55961ae7858541fb95319a1e4454c3c4e953c603478e51e488b7f8b33`。
- solve源码执行前后及result记录的SHA均为
  `e9074068c755b7635f7fcdfb2acf262e6a1d48e49143b4490e5d581978c753cd`；solve.run.json已terminal，result的candidate_exported=true。新轮wrapper必须再核这些来源与输入，不仅依赖bool。

两个checker读取同一个不修改的plain候选，准确继承α=47/40、R=33/20、h=1/32、N=80。plain schema严格只有alpha/R/h/N/c/a，不能把beta或structure交给旧schema checker；plain六字段须与structured逐项完全一致，structured的beta=alpha。准确成本的certificate_sha256必须等于structured SHA。

预定second的claim为准确有理数 **359087/100000（十进制参数3.59087）**。每个checker启动前用标准库Fraction重核 `cost.upper < cost.baseline_lower`、improves_baseline=true，且 `cost.upper < 359087/100000`。若此准确claim gate不符，停止不运行，不用近似成本放行或临时改claim。该claim不是3.4，只是本次候选明确可认证的上界目标；最终仍以全域核验和一般输入接口为准。

## 3. 固定两checker及预算

每套只有一次全域运行，各 **180秒外部hard、6进程**，总数学预算最多两次360秒；可由主侧并行安排，仍各自严格计时。输出/日志/run文件已存在就拒绝重跑。只清理各自wrapper创建的独立进程组，不触共享GUI。预算包括进程启动、全部检查和报告写入；超时立即终止进程组并标incomplete，不能把局部PASS当全域通过。

### 修正first

使用第36轮固定修正副本，不能误用未修正上游hlcert：

- `../round17/fixed_checker/cert_run.py` SHA `079906db354febdd57d537f0d747428bc90bd37fa88297449cd1d46c7446310d`。
- 同目录hlcert.py SHA `4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05`。

固定参数 `--Q 8 --maxdepth 7 --Kdepth 6 --procs 6 --only MKA`，独立输出first_report.json。沿用其原精度与公式，不为失败修改阈值、深度或精度。该副本的有限积分log已传递analytic契约；源码before/after SHA都须一致。

### independent second

- `../upstream-v4/checker/independent/hl_check.py` SHA `0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8`。
- 同目录hl_core.py SHA `c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf`。

固定 `--claim 3.59087 --procs 6`，完整默认C1/C2/C3/C4，不指定cells缩小范围，不使用tight追加探索。默认深度明确为maxdepth8、kdepth10；如wrapper显式列出，应与此相同。独立输出second_report.json，不把两个共享Arb和解析族的实现称为全部数值原语完全独立。

## 4. 产物与判定

每套保存不可变候选SHA、命令、全部source/dependency SHA before/after、启动/结束/超时、PID/退出码、实际时间、日志和完整报告SHA。候选、structured和成本文件不改写。大报告可留本地并保存原hash、最小摘要和失败定位，但不能只保留PASS文字而丢失完整终态。

成功要求两套均在预算内完整终态exit0、first的M/K/A要求均通过、second的C1/C2/C3/C4均通过，且准确成本claim gate通过；再由主侧核一般L1/所有尺度接口适用性、独立审计后才报告新已证C2≤3.59087。局部连续LP行或solver optimal都不能替代该门槛。

若报告存在失败盒，负充分下界只表示该checker没有闭合，不等于严格负点或数学反例；若超时，则为预算内未完成，不把C4写为数学失败。保存父格/盒/下界与方法，下一轮可据实际范围另定固定点严格诊断或实质模型改进，但本轮不自动追加点、加深、换精度、重跑checker、修补参数或求第二候选。任一套未通过时B1继续保持，3.4目标明确未完成。

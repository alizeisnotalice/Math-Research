# 第34轮草案：固定收缩候选的两套完整核验

状态：仅静态接口及预算草案，未运行。依据第29轮实际run元数据及当前公共CLI源码提取；没有重查旧证书、重新评价核或调用checker。主侧可在repair终态通过后冻结本方案并唯一执行。

## 启动条件与schema

先完成本轮预定lambda=1023/1024准确仿射、全部16272条Q half行及51574缓存数值行门槛、准确整数根成本，并证明新成本严格低于已证B1及本轮固定claim=3.5909。失败即结束，不先启动全域checker，亦不换lambda或claim。局部安全行仅覆盖三个父格，不能由与已证基线混合自动认定新核全域可行。

两CLI输入同一份plain `candidate.json`：仅alpha、R、h、N、c、a六个算术字段，值保留准确有理字符串，alpha=beta=47/40。若repair原件还有lambda/structure/provenance等元数据，另导plain副本并逐项核六字段完全相同，保存schema_binding及双方hash；不要把元数据塞入严格的Spec，也不能把alpha改成其他beta而同时改变算子。不得二次舍入系数。

运行前固定候选SHA、4个公共源码SHA、命令、cwd及输出路径。核c非负、规范D4键及实际非零样条外包U>=R：第一程序的有限窗口U由样条索引推得，不自动把R纳入。固定内支撑T与正s应继承基线外部非零系数，但仍需在候选schema层核窗口。此为输入检查，不是新核采样。

## 唯一命令与预算

cwd固定为 `../upstream-v4`（相对本round15目录）。PY为本机数学venv Python。候选与输出使用相对此cwd的路径，避免错误地从round15 cwd解释旧命令：

```text
PY -I checker/cert_run.py ../round15/candidate.json --Q 8 --maxdepth 7 --Kdepth 6 --procs 6 --only MKA --out ../round15/first_report.json
PY -I checker/independent/hl_check.py ../round15/candidate.json --claim 3.5909 --procs 6 --maxdepth 8 --kdepth 10 --only c1,c2,c3,c4 --out ../round15/second_report.json
```

若plain文件另命名为candidate_public_schema.json，只在冻结前统一替换两命令输入路径并更新schema绑定；不得两checker检查不同候选。

建议每程序**一次180秒外进程组hard**，各6进程，可并行；无加深、重启、剩余格补跑、改only或tight。两个程序预算合计360进程墙钟秒，若并行不等于承诺总实际耗时。第29轮同配置曾分别88.852秒与47.447秒终止，因此180秒是有既存运行参照的有限预算，**不能保证本新候选足够完成**；若树扩大导致超时，结论unknown，不自动延长到旧600秒。本轮repair的90秒是另一个阶段，不含在上述各180秒内。

第二命令中的显式深度8/10与旧默认一致；主侧当前wrapper省略这两个参数及only，由已pin源码保证同一默认配置。不得传 `--cells`，否则只是局部C4；不拆成先K、后A两次以隐藏重复调用。程序内部按质量、C1/C2、C4顺序执行，完整条件包含在每套唯一调用中。主侧wrapper采用绝对脚本/输入路径、cwd=upstream-v4，和上述相对命令等价；输出名称应统一为first_report.json/second_report.json。

`--claim 3.5909` 是主侧本轮冻结的新门槛，低于已证B1，已不沿用第29轮的3.59091；只有准确成本上端严格低于此claim才运行。它仍不是3.4目标，最终报告更准确的新数值须用独立准确成本上界。不要直接把宽度约1e-60的准确成本上端当128位Arb claim而造成纯包围宽度误拒；第一程序M选项也只输出质量包围，不自动同目标比较。不得观察checker输出后临时改claim。

## 依赖及完整通过标准

当前静态核对的源码SHA256：

- cert_run.py：`079906db354febdd57d537f0d747428bc90bd37fa88297449cd1d46c7446310d`。
- hlcert.py：`4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5`。
- independent/hl_check.py：`0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8`。
- independent/hl_core.py：`c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf`。

两套有不同包围/程序实现，但共享python-flint/Arb后端及已审数学公式，不能称完全独立算术或形式化证明；本轮无Lean build。第二程序保存主/worker来源与spec hash；主侧也应记录第一套实际加载的输入/源及运行环境。

第一程序应正常exit0，完整JSON同时含mass_half、U、K_fail、A_fail、A_cells，K_fail/A_fail均空，日志末尾PASS且argv为完整MKA；不能只凭exit0或默认空列表判断，因为只运行M也可能exit0。M只报告质量，准确成本另验。

第二程序应正常exit0、完整JSON有c1/c2/c3/c4且四者ok均true，c4.prov_consistent=true，spec_sha256等于候选，实际完整cell数量/结果数吻合；日志OVERALL完整四条件PASS。CLI的allok只聚合实际存在项，故不能遗漏某条件后仍据PASS称全域成功。

任一超时、被杀、缺报告、输入hash变化、非有限或条件未闭合均不得升级上界。负的盒下包围失败不自动是真实点反例。记录真实状态/exit code、elapsed、源和候选前后hash、日志及report hash；终止后冻结最终日志快照，避免把shutdown追加日志混入先前hash。只清理本次专属进程组，不触碰共享GUI或其他任务。

完整大JSON可留本地，但归档摘要必须逐条件保留终态、fail计数/覆盖数量、source/input/reporthash及实际预算，明确没有包内全叶重验。两套全域条件及准确成本/解析桥全部通过后才可声明阶段改进；即使如此，若准确成本仍高于17/5，目标C2<3.4仍未达成。

## repair及阶段wrapper的最终静态审读

已完整读过 `repair_once.py`，最终实际SHA256 `038e19156bb3be6aeb54b839554dda2c29a025444cc867b62c45239e42ce54b6`；阶段wrapper实际为 `b3cfa2de8c88290268b0d7113f475030262eafc14f867d4bb46c0593a7bdc688`。无剩余阻断，已通知主侧可唯一执行。本审读只读源码及hash，没有生成新参数、行算术、缓存矩阵乘法、成本或checker。

固定lambda的c/s/t准确仿射及z=abs(t)=lambda*old_z正确，前后准确参数盒核验，不二次舍入。16272行复核baseline身份/正余量、非负signed penalty及新参数half目标；仅全部通过才进入一次51574缓存数值筛查，K/A原与加强目标分别保留1e-9/1e-8，非finite拒绝。没有AS/Psi模型、LP或第二lambda。实际T支撑再准确核<=3/4；c/s及a=s*a1+T身份保持，不能据此推AK外区。

成本只调用已pin的旧 `structure/enclose`，不运行其main或Arb；旧B1算术文件另pin。准确上端同时严格低于B1下端与35909/10000才导plain候选。plain只保留六字段且逐项同structured，candidate/structured/cost/parameters的hash均在result绑定。

最终wrapper在每个全域任务前重新核candidate/result、structured、cost及plain字段，核准确upper双门槛、四个公共源码pin；任务后核候选hash未变。两阶段专属组超时/终态只清理自己组，无重试或修改深度。repair90s与各fullchecker180s分阶段计；不能把未运行的静态放行写成通过。

旧solve_result/build_result/parameters/rejection四个小输入文件pin已实际核符；其余大行/缓存将在唯一repair中按源码核hash，本审计没有重复读取模型或重新核完旧行。后续实际结果另作终态审计。

# R009 fan rays：独立全量精确审计

结论：完整的 1,947,792 个六行子集均已实际计算，得到 rank<6 552,979，rank 6 混号核 1,066,827，rank 6 非负核 327,986；三类和等于总数。46 条本原非负射线和候选逐项一致，包括 w、未归一化 ξ、Σξ=1 的精确有理坐标、36 条法向量的关联位、64 条体对角形式的符号、字典序最早的六行 witness 与非零余子式。这只核验射线枚举数据；不声称计算完整胞腔、三角剖分、Green 核或上界。

输入锁定：实际引用的候选 manifest 全部散列在计算前、后核验；Sol 理论输入目录各文件散列同时锁定。input_lock.json 记录完整 SHA-256。由 (1,±1,…,±1) 字典序独立生成 64 条路径，逐项求前缀和，得到 35 条各项非负、29 条混号。将首达负值的长度 3、5、7 模板按冻结顺序生成，在 7 条坐标墙之后排列；得到 36 条法向量，与候选输入和 Sol 的冻结法向量逐项相等。

方法：对固定字典序每个 6×7 子矩阵，使用 python-flint 的 fmpz_mat.det 分别精确计算七个带符号 6×6 余子式。七个全零等价于 rank<6；否则矩阵 rank=6，余子式向量给出其核。若核两符号兼有则排除，否则统一方向为非负并除全部坐标的正 gcd。保留每个本原向量的最早字典序 witness。该枚举不调用作者 rank/nullspace 或引用其 hits 作为答案，不做浮点筛选。按 2048 项分块、最多 30 个在途任务，15 个单线程进程；按顺序验证 952 块的首尾衔接及完整长度，并按含墙行数比较所有七个组合数分层。验算 46 条射线的其余数据后与候选逐项比较。使用该方法得到的全部射线及实际每进程 CPU/RSS 在 exactresults.json。

**真实执行状态须区分：**原执行程序先完成全部计算、全部比较及源散列后输出一个 PASS JSON，但末尾 `except BaseException` 错误捕获自身的 `SystemExit(0)`，接着输出 `FAIL: SystemExit(0)` 并以状态 1 退出。run.stdout 原样保留两条输出，receipt.txt 保留 time 输出和这一退出状态；exactresults.json 是从 stdout 首个完整 PASS 对象原样解析，没有改动算得的内容。该进程的计算结论可以从已保存的完整 PASS 对象与源代码检视，但不能声称该原始命令正常退出。为避免重复 1,947,792 次科学计算，未再次运行完整枚举。

portable_checker.py 只将原执行器结尾的 `except BaseException` 改为 `except Exception`，以允许正常的 SystemExit(0) 通过；该修复版已做语法解析与单行差异检查，但尚未完整运行，不能冒充运行验证。需要 Python 支持 ProcessPoolExecutor.map 的 buffersize、python-flint；运行方式：`python -B portable_checker.py --candidate CANDIDATE_DIR --theory THEORY_DIR`。脚本只读源，完整结果写 stdout；路径须指向与 input_lock.json 一致的文件。原计算脚本 recheck.py、原 stdout、两版散列及其它材料由 manifest.json 锁定。

全部子集为本轮新候选验证，没有重算旧基线或进行 lower、cells、Green、MATLAB、提交或推送。实际独立 clone 从本地源 commit 75eedb49 建立，配置 SSH origin、promisor/blob:none、sparse checkout；不宣称网络 partial clone 已执行。理论上射线足以刻画相关截面的证明由独立理论审计负责。

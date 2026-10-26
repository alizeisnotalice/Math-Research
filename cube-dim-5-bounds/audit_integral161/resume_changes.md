# 161 v2：仅兼容缩减 packet 并续接未完成后缀

run001 的主进程 80086 实际退出码为 1，工具报告墙钟 3.136044792 秒。
15 个 worker（80094–80108）均记录 `KeyError('level')`，已查无存活 PID。
旧 driver 没有保存每个子进程的 returncode；本报告不把推断的退出码称作实测。
完整 UH、DFS 和最终接纳均尚未发生。

根因是 range092 的已锁原 packet 使用四字段缩减 schema：
`endpoints/key/original_raw_sha256/point`。旧独审把适用于 adaptive095/103
完整 packet 的冗余 `level/indices/weight` 检查无条件用于它。
原 112 row/certificate/tree 中的这些字段及真实点偏移、orbit、16-child 分割
已由冻结准备步骤独立核查；故障不是一条均值不等式失败。

`audit161_v2.py:packet_check` 只给锁定的 range092 来源允许上述精确四字段 schema；
其路径、全部 packet SHA、key/point、73 个时刻仍检查。其他来源仍须实际核查
原字段 `level/indices/weight`，没有 `.get` 默认值或放宽到任意缺失字段。

`prepare_resume161.py` 只作元数据绑定，实际 PID 81062、退出码 0、2.707985541987 秒。
它逐成功行核对原 source 顺序、当前输入 packet SHA、原执行代码身份以及
producer 的实际字段；15 个 gzip 都完整可读。它不重算 factor、band、mean、gain、根和。

冻结成功前缀数为
`2790,2790,2789,2789,2789,2789,2789,2789,2789,2789,2789,2788,2788,2788,2788`。
合计 41833 个 mean／3053809 个 band products；292 个 factor 证明也已经完成。
旧 run001、旧代码、plan、binding、manifest 全部保留原字节。

v2 实际执行须另获授权，使用排他 run002。它逐字节复制旧 factor_checks、coefficients
和每个完整前缀 gzip，核对 SHA。前缀只读元数据和收集已得 self，不调用 weighted_mean
或 gain 算术。仅余 56074 个 mean／4093402 个 band products 在 15 个 worker 中计算；
均值、ceil160、divmod gain 算术保持原独立路径。15 个 worker 全 ready 后才放行新乘法。
旧几何/orbit/16-child 准备门继承，不重放。

全部新后缀成功后，第一次完成整个原 97972 节点／52360 根的独立 DFS，
复用旧150的73 factor／64 root值，保留1个原 self，不重放旧150算术。
仍比较所有 producer effective_changes 和准确 UH；一项失败就保存并停止。
v2 同时保存实际 worker returncodes，修复旧 driver 的退出证据遗漏。
资格仍受原 theory.md 的同一模型、112严格活动弦／范围、正时间质量、完整 F 和029限制。
这里只能交 UH 给主代理组合，不自行合成或宣称新的 C5。

# R022 有限参数批：首错停止

本轮预登记 a=50,…,80 除64共30个新参数，逐一用512位 Arb 求候选公式的外向区间，并与封存 R021 完整准确区间比较。首个参数 a=50 的区间已计算并写入 terminal.jsonl（约77.8256807772983，严格高于旧 a=64 的完整上端）。随后写入持久化环节调用 `os.fsync` 时，程序遗漏 `import os`，抛出 NameError；错误处理分支再次调用同一未定义名称。进程 exit 1。原 terminal 保留先 completed 后 failed 两条 a=50 记录：前者表示算术对象已生成并 flush，后者记录持久化步骤失败，不代表两个独立参数。

依首错停止规则，a=51,…,63、65,…,80 共29个新参数未执行，无本有限集合的胜者，更无连续参数最优性结论。未重试、未续跑、未计算旧参数、phase、MC或 lower；未提交推送或运行 MATLAB。源代码、原 stdout 和 stderr/time 均原样保存。唯一执行的终端 a=50 与R021原区间比较布尔值可从所存 exact Fraction 离线复核，但这不构成完成整批。源目录本地 clone 来自 commit 71570d6469db24e037fa763e32536bbe7852585e，配置 SSH origin、promisor/blob:none、sparse checkout。

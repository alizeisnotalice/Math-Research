# C149实际运行结果补充：100 NEW密度通过

这是冻结准备包之后新增的实际结果说明。原source、proof、packet/manifest原件均不改。执行者为Terra，原controller PID68873；实际原session已exit0并reap，operator_terminal_receipt记录automatic_restarts0、exact_density_generations_started100、unobserved_case_ids空。这里不重新运行任何密度或旧functional。

## 实际科学终态

Terra result/exit状态为`COMPLETE_100_NEW_CANONICAL_DENSITIES`：100/100 NEW canonical density cases通过，380完整精确stage检查，3634 literal endpoint substitutions，per_case_counts_complete=true，全部原worker已reap。15个实际worker PID被记录；没有失败/未观察case。源runner与core仍是已批准的50e292d2…0727e、6e5c41aa…bfba6。

结果SHA256：`e2012d1a0d70912a3e7bd5bf154246f3e2a73da645ba5e2358d588f8a0feae81`。

完整operator wall为6862.296279625007秒，约两小时。operator所收原child/descendants user CPU为24692.25454秒、system CPU403.004964秒，合计约25095秒。controller exit另记录其完整reaped child CPU25095.066199999997秒；两者scope不同，不能互换当同一个精确计数。无wall/CPU ceiling或自动重启。

100 NEW canonical laws对应106个原cone/body pair来源，准确聚合weight为2698731413/2247297614217216000。3个C147已通过keys仅引用，未重复生成。每个case/stage的source与证书引用保留在Terra输出内。

## 最后完成的实际硬案例

最后completion是case_050，key `0adf3f872414493dc9cf81c867b8a1ba698ac234a87d05a10f67f7c61172b97e`，worker68883，与已有采样对应。它的完整wall6821.838279332995秒、CPU5903.582865秒；6 active列、6方向、5个完整stage。记录的最大非零项19、最大可见expression节点878、最大多项式degree3、最大系数38bits；累计Laplace monomial summands78、literal substitutions172。

这些显式计数不是CAS内部操作数。较小的可见输入/输出与长耗时同时出现，给后续独立checker/producer优化提供了冻结参考，但不能据此宣称已定位某函数百分比，不能重编此已PASS key。后续性能静审在本clone `upper/round151/performance_diagnostic/`；当前原生产源不改。

## 独立saved审计与范围

155v2已经实际核固定5个案例/15 stages，receipt引用而不重放。C158将按独立新登记核其余95案例/365 stages；本补充不把C158计划写成已终态通过，也不重新跑它。原数据和source hashes供独审直接消费。

本批只证明这些source-bound二维density law的完整LT/归一化接口。0 k_h integral、0特殊函数积分、0旧6D网格/原生kernel重算；没有新的beta2数值上界、M/globalg或C6结论。全局数值上界仍须所有keys和全部积分/尾/误差齐备。

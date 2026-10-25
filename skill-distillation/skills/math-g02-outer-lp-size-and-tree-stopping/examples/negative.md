# 大剩余条件与旧全局范数目标的反例

在 16 点集上令 ω 计数、μ(nonempty)=1、F≡1、r=2。条件 μ(A\B)≤μ(A)/2 迫使 B=A，真量为 4；误换成 μ(B)≥μ(A)/2 允许只保留一个点，量为 1。

另取 f=(1,2,…,1000)。单点窗口的平均值为 1000，而全局 RMS 为 √(Σ f_i²/1000)<1000，故无条件的 `sup_window average ≤ global RMS` 为假。原 R13 脚本实际检验较松的 `supavg ≤ L2·√(N/10)`；该通过不能支持 plain-L2 目标，已在 `audit_current/g/history/` 原样封存并从当前结果撤回。

以上反例及局部 Cauchy 正目标均由当前 G02 案例脚本独立计算，结果见 `audit_current/g/case_results_20261007.json`。

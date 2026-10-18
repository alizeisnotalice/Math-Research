# 第348轮：真实行权容量扣减与已到达上尾的原时钟预算

**一般根阶一致上界仍未闭合。** 本轮保留完整原时钟行权a=qω，对全部relative行精确拆分a=min(a,τ/U)+(a−τ/U)_+，τ=2λ(1+√n)/q_-。第一部分以与347完全相同的C_hit,n Wex支付，覆盖所有旧Top-hit，并在每个高score行支付真实部分；余量保留soft扣减，不放大回整个高score交通。

再将soft余量的固定行mask乘回原Ω，用真实径向逆质量二阶尾和原ell（含b原子）证明：当前prefix上尾捕获>1/q_-的分支只需γ_n Xq，其中γ_n=2(λ/q_-)e^-√n。highpast可非空维数中严格有γ_n<ε/16；必要维度条件失败时原分支为空，令γ_eff=0；否则γ_eff=γ_n。

完整槽为

`B_H ≤ G_late[g_ex] +[(13/3)ε+γ_eff]Xq +C_ε,n W +(C_top,n+C_hit,n)Wex`。

旧347付款不重复添加。剩余G_late保留真实ηΩ且当前prefix上尾≤1/q_-，一般未来联合交通预算尚未证明。历史CZ同一外账条件下接收系数841/2304+(256/119)γ_eff<853/2304。软层、质量可比重尾、极端尺度、初始/删除/未来被删来源及全部外积分继续独立留账。

实际本机复现、独立预审/全文后审及具体证据范围见包内main、numeric、geometry、core、audit。精确有理结果与Decimal诊断分别标识；有限原子测试不认证原一般立方体实现。没有编译PDF。

[研究包](research.zip)包含详细LaTeX源、完整总账及必要脚本。解压后用Python标准库运行r348_numeric.py、r348_core_checks.py、r348_geometry_checks.py。承接[第347轮](../347/README.md)。仅新增两个仓库文件并更新索引。

归档SHA256：`0f99f982c094da0566aa6c216ef4a9cee3efee62a2d091f68e17cf7eae54b1a6`。

# 第46轮诊断源码独立静态gate

全文审冻结diagnose_once.py SHA `53195f55dc40dad9da6749e2bea7ad55d748cd9d3597fdea53d38f36948dad3a` 及wrapper SHA `5bb9140fc42aac834c134209f4a20e5d8e13602fdf82ce471d1ae643e11fb6a6`；结论GO，可按固定4调用/192bit/60秒预算唯一执行。审查未运行新数学。

- 候选、两终态报告、round04适配器及修正first/independent源码逐文件pin；P1从first完整A_fail最低下界/原序并列规则的physical盒重建，P2只取second唯一失败父(11,5)首保存local盒，经h*(11+x,5+y)转换。与两个预定Fraction中心及父格断言吻合，不改最低样本或追加点。
- 导入helper（不main）、修正hlcert和hl_core后设置ctx.prec=192，再构造Kern/Spec，避免模块默认精度覆盖。输入全为准确有理数。
- 复用已审OriginalAdapter.AK的精确B/D样条求和及内区Calpha*r^-2alpha+4I0+2I0+2I0，c仅乘一次，样条sc=h^-alpha*κ完整。IndependentAdapter用Spec.AS_point/APsi_point；两正坐标/r<R断言排除外区积分分支。不调用盒模型或K点评价。
- scalar返回后同次保存已算AS/cAψ分项；man_exp只对有限exact上下端转换Fraction，不从展示mid取端点。严格符号、差球含零和区间交叠分别记账，非finite/异常停止不重试。
- started/completed/failed/interrupted及not_started、append-only事件、源/输入前后hash和traceSHA齐全；最多两点×两方法=4完整诊断，外wrapper拒覆盖并清理自有进程组，60秒hard不扩。
- AK局部连续性和严格负→开邻域负→非负测试函数负配对的结论适用于这两内区点；若值正/跨零只记录点事实，不能认证父格或全域。共享Arb及解析族的两路径，不称独立积分原语。
- 保留first/second失败下界与本轮真实点值的证据区别；完成点评价之前不预判符号，不更新B1或声称达到3.4。

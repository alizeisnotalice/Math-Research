# 第34轮：唯一固定向已证基线收缩

状态：只静态准备，不运行；主侧负责90秒hard包装及后续全域核验。第33轮仍按原协议拒绝，不能追改成已通过。其失败是额外half-baseline加强余量的一行小负残差，独立核验原安全行仍为正，并非真实AK负点。

固定λ=1023/1024，旧准确参数q来自第33轮parameters.json，基线x0=(c1,1,0)。只作一次准确Fraction仿射组合：c′=(1−λ)c1+λc、s′=(1−λ)+λs、t′=λt，z′=|t′|=λ|t|。不再次舍入，不选择λ，不作第二修复或LP；任一门槛失败即停。

读取pin的第33轮build_result/directions/rows/parameters/solve_result及rejection_verification SHA，核完整12方向/16272行、原顺序与准确旧z=|t|。参数盒保持c1/2≤c≤3c1/2、1/2≤s≤3/2、各粗/细|t|≤1/2或1/4，z上界相同。新参数应仍在此凸盒，实际准确检查。

逐16272行准确核新安全值≥baseline_slack/2，无容差；同时核已存baseline_slack=c1*c_column+s_column>0。收缩的解析恒等式为h′=(1−λ)baseline_slack/2+λh，其中h是旧加强残差；不能仅以数量级预计成功代替全部准确行检查。任何负行不换λ、不改rhs、不放宽为仅≥0。

通过才复用原51574浮点缓存行作一次新参数筛查（K26604/A24970），不重建B/L或AS、不回算旧q；核旧manifest及缓存/准确点序hash。保留原目标和ε1/2旧基线slack目标，固定K1e−9/A1e−8；非finite或旧baseline负slack停止。一个缓存值同时作两种比较，不虚增点评价，也不把缓存筛查称严格生成元包围或独立完整系数回代。

将12个方向的准确细网格maps乘t′得到T，再准确a=s′*a1+T，检查canonical/N窗口、c′>0/s′≥0、每非零T索引和≤48（另记录实际支撑≤3/4）及逐项结构恒等式。所有其他T仍为0。新核等于(1−λ)K1+λKq，但Kq未全域认证，故不能由凸组合直接宣称新核全域可行。

最多一次已审整数40次根半积分成本包围：复用pin的round10/enclose_cost.py，仅调用structure与enclose，不执行其main、不调用模型/Arb。质量按α=β=47/40准确公式，主要比较旧已存B1准确下端和用户目标17/5；成本上端须严格低B1下端，且低于主侧本候选checker固定claim35909/10000，才导出plain candidate。3.5909只是执行claim配置，不替换用户3.4目标。没有新LP、修补反馈、积分或全域checker。

输出parameters.json（新旧参数/λ/hash）、structured_kernel.json（beta=alpha和structure={s,T,baseline_sha256}）、exact_cost.json和result.json。只有参数盒、全部准确halfslack行、51574缓存筛查、结构及准确成本均通过才写candidate.json，字段仅alpha/R/h/N/c/a；它与structured_kernel的实际核字段逐项相同，记录双hash。失败保留中断/实际计数，不输出候选，不重试。

唯一90秒hard含导入、全部输入/缓存hash、JSON读取、准确行运算、缓存乘法、精确展开、成本整数根及保存。内部时钟在循环/调用前后检查；0新模型、0LP、0新点函数。全域其他区域/轴/原点/截断界面/无穷远仍需主侧独立两套fullchecker及已审解析接口，候选成本低不等于新上界。当前已证B1保持不变，目标3.4尚未达到。

# 保存真实后验：方波/CDF与径向输运精确守卫

2026-10-07。独立新预登记，仅后处理 posterior_radial_age_probe_20261007 已存完整来源与连续窗口真实赢家后验；没有导入或重跑来源、arrival、winner oracle。原解析终稿与上游文件未修改。

登记规则：每轮按原 case/receiver/candidate_L 顺序，先为尚未选取的 n 取首个合法 pair，再以原顺序补足至多3对。资格固定为L>R且mL>mR，未按gap、测量值或检验结果筛选。round1的n1无合法pair，因此保留该未选择事实；三轮实际n为2/3/2、4/8/4、16/64/256。

|round|n|R|L|ΣD_i|全局方波平均|首个有效k|
|---|---:|---:|---:|---:|---:|---:|
|1|2|1|3/2|9/56|9/112|2|
|1|3|1|15/8|5/8|5/24|3|
|1|2|1|7/4|9/56|9/112|3|
|2|4|1|3/2|21/272|21/1088|6|
|2|8|1|9/8|4213/4898|4213/39184|3|
|2|4|1|7/4|21/272|21/1088|6|
|3|16|1|35/32|765/224|765/3584|2|
|3|64|1|17/16|3961/960|3961/61440|5|
|3|256|1|51/32|32757/704|32757/180224|9|

全部操作用Fraction：按实际坐标排序积分CDF差；phase按全部实际atom坐标mod2b断点完整积分，验证每坐标J_i=4D_i和全局均值2ΣD_i/(nb)。物理coincident来源的posterior合并，source signs只能依物理位置，完整原来源仍保存在每pair输出。

固定搜索s=R+(L−R)/2^k，k=1,…,32，首个1−c(s/L)^n>0处验证ΣD_i≥(s−R)(1−c(s/L)^n)/2，并从已存D直接核该半径的真实后验质量帽。全部attempt保留，没有修改输入或扩大搜索。直接距离下界ΣD_i≥ΣπL(D−R)_+/2也全部通过。

终态一次正常exit0，耗时0.172943秒。共9 selected pairs，359 coordinate identities，9 global identities，9 direct distance inequalities，9 positive local-radius inequalities；phase5155区间、CDF2067区间，无empty轮或未找到s。没有浮点/超越函数、Monte Carlo或oracle调用。

保存registration/code/results/receipt及9个pair JSON。每pair包含完整原source、保存receiver与winner/candidate posterior、物理合并posterior、全部坐标CDF区间、phase区间与signs、local search attempts；receipt记录全部上游及新输出hash。

这些检查只认证固定真实来源/receiver的精确算术身份与不等式。未积分receiver Lebesgue体积、未检验nearmax或确定K、未提供actual FIRST/history资格，也未证明一般维数阶数。

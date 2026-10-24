# 第44轮solve静态审读

结论：GO，仅在本轮builder完整成功且输入hash吻合后可唯一执行；未运行任何数学。

- `solve.py` SHA `e9074068c755b7635f7fcdfb2acf262e6a1d48e49143b4490e5d581978c753cd`。
- `solver.m` SHA `7c40dec75ac7eff1706ab8b96e8ab1a0a5d6bd0d0c06316bd597bbde95d646b8`。
- `run_stage.py` SHA `256d8d134fea5b76485489014b2fe979957db501e2803ac015b0ee9b6c90148c`。

- 26列为绝对c/s、固定12个旧column和12个z；原/加强K/A点行符号、连续midpoint/负penalty、24绝对值行及c−s*c1≥0均正确。总88319行，连续36720行；旧方向、矩阵/rhs/目标/bounds及MATLAB输入输出hash绑定。
- MATLAB独占batch一次linprog，20秒soft，180秒wrapper自有进程组hard；无fallback或basis warmstart声明。非exitflag1、非finite或状态异常拒绝，wrapper终态负责清理自有子进程。
- float_slacks复用首批51574点行；有限门槛已补，最终Q再查第二批51574，未加入第三批缓存点值。K1e−9/A1e−8只为数值筛查，不能解释为严格非负证书。
- 返回float先转准确Fraction，以预定1023/1024向baseline收缩，再唯一2^-50量化g/s/t；重建c=s*c1+g、z=|t|。准确原参数盒、结构及36720条half行逐项核验，不clamp、换lambda、修补或重试。
- c/s/T的真实半质量进入成本；导出准确s*a1+T不逐系数舍入，D4质量和T支撑≤3/4核验正确。唯一整数根成本调用以upper<B1lower为出口，17/5独立报告；成本本身不代替全域核验。
- 读取新的builder结果仅在complete后；细节schema及新旧方向一致性仍须与冻结builder联审。发现阶段候选即使通过局部准确行和缓存容差，也不是新C2上界。

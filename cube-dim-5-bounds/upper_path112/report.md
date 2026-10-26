# upper_path112：全树路径均值消费

固定五维 alpha=1、16条未单位化 body方向、共同Cauchy时间、Fourier (2pi)^-5，同一G。只消费已存在的H512质量和固定树，没有增加空间点或重新计算H/jet/log/exp/asinh。所有源目录只读。

完整97972节点已消费，94580物理叶、加权体积1。先256个真实加权self缺口最大的叶子pilot，严格全局收益后，另登记24个至多4096节点的有限批；最后一批小于4096。完整I上端由5.1807887336656288e-05降为4.7619359947475389e-05；可读严格有理cap为 `297621/6250000000`。同g105=21/200000的完整成本上端由36.4436413689652降为33.4972716580468。S=32(g/5)^(5/4)-I_U 的已证下端约-2.1285052332480712e-06，仍为负；本轮没有通过严格32门，也不改变正式C5。准确有理端点和根-free门余量见summary.json。

## 机制与固定接口

111的 J_a=integral W(rho)R(rho/a)d rho、W=E(max|U|-rho)+ 用真实四尾offset的精确CDF分段多项式构造。消费同节点凸弦的非负alpha/beta，初始每CDF片8panels，保留旧085/R右端分支取min。5种真实偏心几何、365个几何时间计划、19272正片的plan与alpha/beta hash逐项匹配原preflight。只补登记的5772个缺失R，192bit整数根外包；所需184旧R不重算，全部继承catalog保持来源记录。新mean factor精确有理结果再向上舍入至192bit，以保持全band求和高效且严格。

每节点从被冻结SHA的原packet读取完整73端点，逐band检查质量下端为正，末band至无穷用H末端完整支付。m_new=min(oldM,newmean,U)，强制L<=m_new<=U，原range L/U不变。x^(5/4)端点以224bit整数四次根向上包围，凸弦非负权组合乘真实2^(-4level)体积和orbit weight，再与old_selfU取min。所有ancestor保持min(self,childsum)，下侧保持原max及原值。仅97971节点mean严格改善，97971节点selfU严格改善；不能把局部因子改进当作同比例全局I改进。

完整范数预算使用076的原eps512短时E和g105，I_U=(UH^(4/5)+E)^(5/4)，384bit整数根外包，无tube。每成功节点prefix都支付其余节点旧界；UH单调不增。终态严格32多项式余量为负，故所有此前prefix均未过门；没有首正snapshot。

## 继承、运行与验证

旧owner源upper_path112的9文件逐字节复制，未重跑成功preflight；transfer_receipt记录每个SHA与字节数。旧scratch空且无存活计算。独立本地clone HEAD af57128518cf97aed4d5af8307b1ac7b1ae2cf6d，初次全checkout因源仓库无关项目缺失Git对象失败，随后仅cube-dim-5-bounds的sparse checkout成功。未修改源Git、其他项目、下界，未提交或使用MATLAB。

scalar实际PID16028完成5772新R；factor实际PID16113完成365计划，原接口hash冻结。pilot及24批的PID、有限key登记、source参数SHA和原始结果路径在各receipt中。production_v2相对consume仅把反复roots总和改为准确root delta累计，限定diff保留；数学公式不变。单worker足够，没有wall cap或总CPU预算。

作者核验完整97972结果、固定拓扑、上下ancestor账本、94580叶体积1和全部范围/均值次序，并做64个固定确定性实际73-band样本（首末和SHA优先），以不同384bit下根对照224bit上根凸弦。该抽样不冒充全量独立数值复核。独立audit_path121另行执行；此manifest为候选，最终理论及完整I验收仍由主代理负责。

## 复现与禁止重复

不要重跑prepare.py、consume.py scalars/factors或已完成生产对象。其成功输出是可直接消费的固定输入。read-only `python3 verify_author.py`可复核作者账本和有限样本，0H/0R。node_certificates.json.gz保存97972新M和selfU，numeric_hierarchy.json.gz保存完整树和最终预算；对照nodes.json.gz及其中冻结的source packet路径/SHA可做独立全量复核。原始输出日志在专用scratch upper_path112_astra，各批receipt记录SHA。外部source103/111/113和g105/E来源锁在continuation_registration.json。

本轮有限任务完成；未登记空间细化、新H或更密rho网格。

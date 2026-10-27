# 第020轮：新诊断及边界（接续记录）

1. 175生产确已完成195942新根、97971新总支付及97972节点完整树；15worker全exit0。暂为UH候选，等待独审与主Cn合成。producer identity.preflight_sha256因setup局部变量x复用，实际值为selfcheck SHA42f8a208...，不是preflight002的0b9ca61e...。真实preflight、注册、全部shard/topology仍被独立plan SHA及setup准确核对。独审与归档应分别typed绑定两个对象，保留历史误名原件，不改冻结、不重复生产；这是元数据别名缺陷，不是H或积分证书被替换。

2. 173原运行代码存在mass>=0与解析报告“非正停”的诊断不一致；v2将门改为mass>0。每个真实正时间band质量严格正，端点区间上包非正应报告来源不一致，而非裁剪。此次只修运行门，原173冻结件不变，无真实band/q/corner重跑。

3. 173原故障后无限join与fault时exit0已在v2加入故障专用有界清理、非零退出及成功前后全锁检查，正常数值批仍无wallcap。主读v2发现结果等待循环收到fault后仍等齐全部summary，使另一个失响应child阻挡finally；因此v2尚不批准真实运行。v3须在fault/stop后立即进入cleanup，并只测试此前未覆盖的“一个fault+一个忽略stop卡死”路径。故障专用宽限不等于正常程序时限。

4. 174实际来源绑定15worker均完成455对象后，小baseline join因独审器使用旧158错误枚举词evaluated（真实evaluated_single_source）而停。这是审计器schema错误，不是173候选反例。保留455成功prefix，仅新小join改枚举字段；不重跑455来源或旧数学。来源资格与173新运行版资格分开冻结。

主保留证明、来源消费及最终常数接纳责任。以上记录不代表新C5界或双侧突破。

截至用户暂停：173v3运行门已修复并完成原生产有限批，但其独立numeric020从未启动；只保存候选和1065pending。175历史alias经独立typed绑定闭合。正式界仍C_cert/U175，随后现有两侧审计和论文交付通过。所有以上故障与措辞修正均不提高下界或重启研究。

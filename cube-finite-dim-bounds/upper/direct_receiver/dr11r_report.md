# 11R完成报告：直接接收权对偶的有限准确诊断

2026-10-25。只upper/direct_receiver/，未改旧锁源/账本/其他项目，未提交推送或扩招。用户要求本轮结束后等待，故本批已经停止，不追加实验。

## 结果和真实范围

完整机制与证明在dr11r_mechanism.md。真实接收集E、可测共同半径r(x)、非负h满足Th≤1时，Lebesgue Tonelli准确给λ∫E h≤||f||1。若另有对任意见证数据的一致mass≥|E|/U，才得到一般L1、全部正尺度的Cn上界。本批没有该mass引理。

两disjoint receiver boxes、每box一个h常数时，Th为交叠长度乘积之和。每坐标断点l±r,u±r划分后multi-affine，所以全空间最大值确在断点笛卡尔积；不是仅查centers或任意grid。由此有限2变量LP可用准确可行h及相同目标的非负input原子prices证明最优，不用无限strong duality。

|每维配置|n3最优h|n3局部receiver比值|n4最优h|n4局部receiver比值|
|---|---|---|---|---|
|volume_critical|(3/4,1)|5/4|(3/4,1)|5/4|
|face_adjacent|(3/4,1)|8/7|(3/4,1)|8/7|
|vertex_adjacent|(63/64,1)|128/127|(255/256,1)|512/511|

这些比值是固定E/r、constant-per-box类的|E|/max∫h，不是Cn。每个h的连续global load准确≤1并取得1。critical例是此限制类准确损失，不排除更细h(x)。equal-face/vertex的非uniform weights确改善uniform two-box baseline，但不推出任意E。

中心检查有准确反例：h=(1,1)在两个receiver centers负载都1，全局却5/4或1+4^-n，断点处超载。不能将receiver centers当全域证明。

## 哪部分已覆盖、桥接缺口在哪里

h≡1/2的两尺度U≤2只是计数；固定r=1,2的2−2^-n是已有finite averaging-kernel envelope。对最多两disjoint axis boxes还可用分离轴证明uniform U≤5/4。本批不把这些有限类基准称新Cn能力。

实际严格bad set是open，各r见证集亦open，因此可以选countable内部dyadic receiver boxes、每盒选一个rational见证radius，finite prefixes在measure上逼近。故真正可检验的缺口可缩为：任意有限数量disjoint boxes、任意assigned positive radii的uniform dual mass下界，常数不随数量/尺度漂移。两盒成功不能成对相加，否则重复花input容量。

若坚持保留任意既给measurable E/selector，不可假设full boxes内逼近（fat Cantor例）；此例不否定上述open bad set自由选择的合法归约。机制文件另给容量误差Σ a_j epsilon_j/(2r_j)^n及mass误差Σ a_j epsilon_j的准确转移条件；只是预算，不是权重存在性证明。无限尺度的关键uniform mass引理仍未闭合。

## 实际运行与最小代码

固定六配置在计算前已锁dr11r_configs.json SHAe5c1c21eb66f61765016d8351a204d0a9f5cc88e88e845c1972e23f7476ae6d3。最初未运行的separated草案是平凡无重叠负载，按主要求缩小为critical-volume例；实际仍只有六配置，没有第七例/新grid/优化扩批。

Python Fraction以6进程并行独立配置，总8232断点（每n3例343、每n4例2401），完整Pareto、LP feasible h、primal=finite dual mass及center反例全PASS。总运行约0.11秒，未为CPU占用重算。最小源dr11r_check.py SHA18220eba65135b50775962ddf020bebb968d01d186d8a2febb18c1f14b0d938b；准确结果/finite atoms在dr11r_results.json。

主串行MATLAB实际exit0/session4782，脚本SHAfaa3c07f5faec736fc60020ce8c63f4c674ba20de63fd5575143ff88cccc7dc2。独立vectorized六配置8232断点全部maxTh=1、primal=dualmass、比值匹配；已实际读stdout，log/hash与最新主账本读取记录在dr11r_run_receipt.json。变量upDR11R，无clear/path。它是交叉实现，不冒充独审。

|类别|本批状态|
|---|---|
|准确证明/证书|Tonelli条件接口、有限连续最大值定位、六例可行/最优LP和center反例；open bad set有限盒归约及误差预算|
|实际工具|Python6进程Fraction全断点PASS；主MATLAB六例PASS|
|未闭合|一般接收集/任意多尺度的uniform mass下界；没有无限strong duality结论|
|未宣称|无新的Cn上界，无任意h或全尺度模型的普适障碍|

固定Lebesgue实际用于receiver交叠体积、(2r)^n正规化、输入/输出Tonelli、mass与|E|；不是一般背景覆盖论。

## 主账本同步及交付

初次输入锁定main为U3=8、U4严格14753/1000；交付前真实main账本SHA43d7d6e54f5e96ade5f26aff9b65c095932763a999c1f8989a305ed914178c67已由036/037更新为严格U3<1979/250及U4<14371/1000。精确代数comparison thresholds按该账本；本11R有限结果未参与或改变这些主上界，不能从历史8/16重新宣称突破。

白名单10文件：dr11r_inputs.json、dr11r_preregistration.md、dr11r_configs.json、dr11r_check.py、dr11r_results.json、dr11r_mechanism.md、dr11r_matlab.m、dr11r_run_receipt.json、本报告、dr11r_manifest.json。来源/配置、源码、准确数据和实际回执均锁hash；不归档额外缓存。当前工作完成，按用户指令等待主复审/论文及用户审核，不启动后续研究批。

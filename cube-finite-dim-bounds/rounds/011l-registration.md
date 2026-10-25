# 11L：四维纯body积分认证独立工作面

2026-10-25。新增原生continuous_integral_n4，请求GPT6.1-sol/high。独立clone /Users/zhengzhihao/Desktop/T/Math-Research-finite-integral4，HEAD d5f647caec30f0006d3dee7c95855b604dbdd05b，仅upper/continuous_integral_n4/写入。稀疏checkout提示其他项目缺blob，但本项目/AGENTS已实际可读，git status/diff干净；不修其他项目。

固定4D纯body，不搜tau、不重复g任务。目标I4=∫_[0,1]^3 G(1,y,z,w)^(4/3)≤19/100000。主数值10x11x12给.000183520694694，只用于候选发现，不能作证。所有点/系数归一使用u11_green.py SHAfdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d、proof63f759e850f3864390cf2af6bf39c57ead78130421389a2cfc5e853335063ea5；帽界读主proofs/round011_face_caps.md与audit030。完整表示独审pending，不自动由026继承。

输入hash另存round011_integral4_snapshot.json。原upper_theory负责g与按sign稳定的连续G评价，双方应共享确定的evaluator并披露，不各自重复同一化简。新代理职责是积分域分解、联合积分上包围、真实方向线帽及全部误差预算；可以研究piecewise rational G系数的polytope分区/高阶求积或严格正time-cutoff误差。不得用有限采样/阶差给上界，不删除node-zero层或把它们全当奇线。

先选定一条可检验机制并锁参数，做≤4096初始cells/区域的有限诊断，最多一次根据明确缺口的表示修订。禁止机械扩网格。若能认证目标，给全域覆盖/体积账、最小checker及证书送原audit；否则给真实上下界包围与未结算项/可复用诊断。不能仅交方案而未调用数学工具。

Python/Arb主动本机实际运行，MATLAB脚本由主串行。你不是唯一工作者，不能改其他项目/总账/其他代理源码，不提交推送或递归扩招。n3已有minimum与J代理，不重复；n4g由upper_theory，不重复；主本人负责真实模型/全L1/allr/最终成本整合与准确比较，U4保持16直到所有要求严格通过。

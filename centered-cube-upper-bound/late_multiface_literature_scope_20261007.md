# 晚参数多面路线：噪声半群文献范围及共享 mask 身份

2026-10-07。定向检索只为核对潜在半群端点捷径，未引入新的文献定理作为主证明。

读 Harrow–Kolla–Schulman, Dimension-Free L2 Maximal Inequality for Spherical Means in the Hypercube（Theory of Computing 10(3),2014）原PDF的 proof overview、Senate 定义、§4开头以及附录Theorem A.1证明：
https://www.theoryofcomputing.org/articles/v010a003/v010a003.pdf

实际读到的弱(1,1)定理作用于时间平均族 Sen(A)，不是未经平均的noise operators。本文的主L2比较不能直接迁移到本项目原ordered弱端点。没有据摘要中的“semigroup maximal”字样宣称已解决当前receiver选时。未全文重证其Krawtchouk部分，也未采纳它为本项目新增已付工具。

一个由本项目原核直接推出的身份：R=(I+S)^−1=∫_0∞e^−τ Hτdτ。原Gi连续，在Hτ内部实际活跃mask每坐标成功概率p=1−e^−τ。替换变量后e^−τdτ=dp，因此R的实际mask推前与∫_0¹Krdr一致。每个A的质量为

∫_0¹p^|A|(1−p)^(n−|A|)dp=1/[(n+1)binom(n,|A|)]。

因此活跃坐标数的推前均为0,...,n上的均匀分布。这只是标签边际身份；R在给定mask内仍保留Poisson重复跳与Gamma时间，∫Krdr则每个活跃坐标只有一次原Gi。不得据共享mask说两个空间核相等、点态可比，或让新退出测度与原FIRST一致。此身份未新增任何空间费用；无需为已知Beta积分重复数值。

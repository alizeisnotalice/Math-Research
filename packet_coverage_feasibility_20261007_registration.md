# Packet coverage feasibility：执行前登记

2026-10-07。新解析对象是任意预固定 fractional packets 的平均覆盖上界，不是 cube 弱上界的新阶数，也不是 actual FIRST/gates 数值样本。只执行本 prefix 的脚本一次，不重跑旧诊断。

三轮维数依次为 n=64,256,1024；η=1/4，θ=1/sqrt(n)，a=1,b=2，D=2n³，ε=1/n²，c=1−1/(4n²)，τ=1−1/n。两种完整来源为曲率近均匀盒，以及均值1、高度 n、周期1/2条纹乘同一曲率。只核解析参数：c>τ、nc>τsqrt(n)、源正性范围、两个 fixed endpoint phase、packing cap 的整数代数、两种密度 cap 的对数半径严格上包，以及 boundary 与 coverage 总上包。

对数半径 H=log(2 Hcap/(η² θ c)) 不用浮点决策：以正 Taylor 有理和给 exp(Hupper) 的严格下界，从而证明 H<Hupper。三轮 Taylor 项数为24,32,40；近均匀半径上包6,7,8，条纹半径上包10,12,14。所有主比较均使用 Fraction 和 factorial；记录 U_n<4/n²、两输入上包随三轮下降。浮点只展示近似量。

此 guard 不进行空间积分、packet 优化、Monte Carlo、连续 winner 搜索，也不认证原 softFIRST/fullfuture/CPGP/LCA。真正覆盖上界和所有 packets 的量词由解析证明给出。执行后保存本 prefix results 与 receipt，不修改主总账。

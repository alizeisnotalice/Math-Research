# 11D：主模型三维混合方向成本pilot

2026-10-25。cb11准确96锥代数及固定36参数点normalization已实际PASS；主已读源，尚未独审其完整表示。锁module SHA ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849，复制相同bytes到scratch不可变文件，不修改代理源。此为数值发现批，非全域上界认证。

固定tau=1/10,1/3,1,3，各face[0,1]^2上9x9格点记录有限min，随后各仅一个SLSQP最多256函数评价（超限记录，使用已见最小有限值）；不对奇性inf点启动优化。面积分分别8、16阶Gauss-Legendre。费用式3^(3/2)*I_est/g_sample^(3/2)，g_sample只是真minimum上界，I_est未认证，绝不写成合法C3上界。只一组固定四tau/两积分阶，不自适应扩网格。

使用本机Python，串行四tau避免并发过大fraction开销。首次异常保留点并停止相关case，不以报错当障碍。批末汇总成本、min位置、阶差和异常；如出现<8迹象优先独审/全域认证，如无收益调整生成元机制而非机械加密。MATLAB待代理固定校准脚本串行。主负责src/round011_mixed_cost.py/results/round011_mixed_cost.json，sourcehash锁在结果，旧证书不复验。

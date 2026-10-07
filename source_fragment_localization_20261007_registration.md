# 实际高势组件局部化：三轮空间组件守卫预登记

2026-10-07；登记先于执行。仅新 prefix，不重跑旧守卫、不改主账。

解析已推导：ξ=μ_hi|{原V>H} 的支持扩厚 S+(-b/2,b/2)^n 的可数连通组件给固定 ξ_j。任一原 query 的正 ξ 捕获只有一个标签，故 local-fullness 切线只收费 Σ_j I_Fj≤I_R 与 Σ_j w_j=w；替换旧 global-fullness 拆分，不再次抵扣 Hw。未证明长组件 fragment 预算。

本次数值是两个真正 L¹、有限尺度 {1,2} 的一维空间族的精确 Fraction 几何核验，n=1 固定，不是一般阶数实验。三轮组件数 J=12,24,48。全部来源 packet 宽 d=1/8、每 packet 质量1、完整密度8；τ=1/4、η=1/4、δ₀=1/8；μ_hi=μ。

远离族：centers=4j，原 E_R 为每中心两侧距离 (9/16,15/16) 的 intervals。原 R=2 是唯一 finite winner，R=1 响应0，R=2 完整捕获恰一个 packet；原 annular V 在来源上恒3/16，H=1/8，因此 B 覆盖全部来源。global-fullness 失败但 J 个局部组件各完整捕获，local-fullness 成立。

重叠链：centers=3j/2，原 E_R 为相邻中心 midpoint 两侧长度共1/8的 intervals（半宽1/16）。原 R=1 响应0，R=2 完整捕获恰相邻两个 packets；end-source V=1/64，interior-source V=1/32，H=1/128，因此 B 同样覆盖全部来源。扩厚组件只有1个，全部原 rows 保持 local-fragment。此例只否定“真实高势集合必有 local-fullness”的无阈值限制推断；H 很小，不反驳 H≈√n 的目标尾。

每轮核：receiver interval 上两尺度的全输入捕获公式、严格阈值/唯一 winner/annular 门、非相邻来源 halo、扩厚连通与远离分离、原 source profile 与真实 B、组件/fullness/fragment 质量门、总来源质量抵扣一次的代数。区间端点用有理数、等式用精确算术；不抽 MC，不增加背景，不把 E_R 子集当完整 E 弱比，不检测原 soft FIRST/CPGP/LCA/history。

脚本执行一次，保存 results/receipt 与 SHA256。任何失败保留，不自适应换参数。没有新普适估计阶，停止于这三轮空间组件证书。

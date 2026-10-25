# cb7_dual：连续B-spline对偶有限批预登记

2026-10-25。只新增本目录cb7_dual_*，不改cb7既有独审016证据、他人文件/主总账，不扩招或提交推送。已读主proofs/continuous_laplacian_dual.md与rounds/007-upper-priority.md，输入绝对路径和SHA见cb7_dual_sources.json。

独立解析核：C_c²测试可经共同有界支撑的正光滑卷积在C²一致逼近；ψ_epsilon(0)不必精确0，但趋ψ(0)=0，固定有限W及ν允许极限。因此∫KΔψ=∫ψdν≥0，q=1−Δψ≥0，cube障碍给∫K≥D=1−(1/8)∫cubeΔψ。连续证明无发现缺口；本批为该Laplacian正有限源核类的成本障碍测试，不是C3下界，也不依赖旧离散核LP。

有限模型（n3）：h=1/q，中心k∈[-K,K]^3∩Z³。中心三次cardinal B-spline B≥0、C²、support[-2,2]；ψ=Σ_k a_k Π B(x_j/h−k_j)，a_k≥0。仅允许max|k_j|≥2，故每个basis在0至少一个因子为0，ψ(0)=0准确成立。非负符号/置换轨道平均，对本有限系数模型保留约束/目标；不能把有限模型弱值当整个连续ψ类失败。

全部支撑胞元ℓ_j∈[-K−2,K+1]，空间区间h[ℓ,ℓ+1]。反射将ℓ→−ℓ−1，故检查代表0≤ℓ1≤ℓ2≤ℓ3≤K+1的64个tensor cubic Bernstein系数；每个代表包含所有64系数，反射/置换只重排系数。单cell矩阵必须把所有邻近中心逐个加入其轨道列，不把轨道重数漏掉或额外重复。目标按完整中心轨道求和。

B的degree3 Bernstein整数分子(分母6)，相对cell m=ℓ−k为−2,−1,0,1：
[0,0,0,1]; [1,2,4,4]; [4,4,2,1]; [1,0,0,0]。
B''线性升degree3后分子(分母3)：
[0,1,2,3]; [3,0,−3,−6]; [−6,−3,0,3]; [3,2,1,0]。
Laplacian系数为q²/108乘逐轴DD×B×B整数和，全空间逐点≤1由每系数≤1充分保证。C²拼接使面/棱/角无奇异项，support外Δψ=0。

目标准确用cube边通量：I_k=∫[-q,q]B(t−k)dt，J_k=B'(q−k)−B'(−q−k)，每中心∫cubeΔbasis=h Σ_j J_kj Π_(l≠j)I_kl；D=1−1/8乘轨道和。另用逐cell积分复算目标作为建模检查。

仅两次solve：q1,K3；q2,K6。各一次SciPy/HiGHS，不设置人为60秒求解上限；实际资源按本机负载协调。每解仅一次有理化修复：非负系数舍入至10^−10格，再准确计算全部Bernstein最大系数；若>1，整体缩至准确最大系数=1。不再求解/换分母。准确记录全部a、逐cell验证、ψ0、D及原/修复差。若严格D>8，立即证书送独审；否则只有限类诊断，不宣称连续障碍。

本机Python=/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python（numpy/scipy/sympy）。MATLAB只交cb7_dual前缀脚本给主串行，不直接调用，不clear/path。严格核验主要依赖Python整数/Fraction；MATLAB独立实现仍共享数学模型，明示交叉证据边界。停止条件两次solve及两次修复与准确checker后，无第三档/额外solve。

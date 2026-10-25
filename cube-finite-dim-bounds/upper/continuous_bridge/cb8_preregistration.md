# cb8：signed连续对偶三LP预登记

2026-10-25；复用原独立clone，只新增本目录cb8_*。不改旧证据/主账本/其他工作者文件，不扩招、提交推送；不是唯一工作者。已读主rounds/008-registration.md、round008_input_snapshot.json及已审signed_dual_interface/019；来源绝对路径与SHA见cb8_sources.json。

有限模型：n∈{3,4},h=1/q，全部中心[-K,K]^n，轨道为0≤k1≤...≤kn≤K，系数允许signed。每完整支撑胞元ψ的degree3 tensor Bernstein系数≥0；Deltaψ系数≤1。ψ0使用准确消元：若轨道仅含s个1、n−s个0，其乘3^n的原点系数为binom(n,s)2^(n−s)，其他0；消去全1轨道(系数1)。n3关系8/12/6/1，n4为16/32/24/8/1。origin等式不能留作浮点容差。

沿用已经独审的B分段定义，但本轮完整signed矩阵重新生成。ψ Bernstein整数矩阵分母6^n；Delta整数矩阵分母3*6^(n−1)，链式导数另乘q²。所有支撑cell ell=-K−2,...,K+1；符号反射ell→−ell−1与坐标置换约化cell，每代表保留全部4^n系数。轨道中心逐个累加，目标按完整中心通量积分并以全cube胞元积分逐列准确核。

仅三个HiGHS solve，按顺序n3 q1K3、n3 q2K6、n4 q1K3；变量free，准确origin消元后各一次，无重启/第三网格。旧019的n3/n4 q1严格例只作为新矩阵可行性回归，不重跑旧checker、不作为新成果。

每LP最多两次准确重建：第一次取数值ψ Bernstein≤1e−8的活跃非负facet，准确有理RREF/零空间投影，在自由参数10^−12网格重建；所有被选facet及origin必须准确0。若全ψ非负不通过，仅第二次可扩大活跃阈值到1e−6并加入第一次负facet，自由参数改10^−14，仍准确投影。不会以整体缩放修ψ负项或origin误差；只有全ψ Bernstein≥0且origin0通过，才按准确Delta最大值C>1整体缩1/C，并按1+(D−1)/C重算目标。失败如实记录，不以更多solve/重建兜底。

活跃行先以mod素数提取独立行，再用准确有理RREF；重建后逐个检查全部选中活跃行，如果漏秩则加入非零行继续同一次线性代数重建，直至活跃约束全部准确0（不调用优化器，不改变阈值）。所有最终ψ/Delta胞元系数都按Python任意精度整数/Fraction准确判号。新checker不导入浮点LP矩阵，核完整signed lattice。

严格D>8/16立即保存候选并送原独审，较低D仍只为合法连续ψ及有限模型诊断，不声称Cn界或方法障碍。Python实际=/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python。MATLAB仅交便携cb8前缀脚本主串行-batch，不自调、不clear/path。

停止：三次solve，每次至多两次准确重建，完整证据和准确比较；无额外网格/参数搜索。模型/算术共享披露，作者检查不代独审。

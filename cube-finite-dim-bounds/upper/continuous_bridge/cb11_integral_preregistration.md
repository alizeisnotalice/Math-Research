# 11H 初始有限诊断与一次有因表示修订预登记

2026-10-25，只cb11_integral_*，不改旧锁源/其他分支/总账，不扩招提交推送，非唯一工作者。已读主011h、minimum_snapshot、030准确预算；源hash见cb11_integral_sources.json。固定tau=1/3，I=∫_[0,1]^2 G(1,y,z)^(3/2)，目标57/100000，不攻g或tau优化，不重复minimum的4095 direct-log诊断。

表示选择：H_epsilon=∫epsilon∞p_t，以全部96锥（48反锥对）真实time表达，t>0各因子实部正，zero/repeated knots均安全。用Acb/Arb矩形时间区间包围，绝不在generic DD上除跨零的小差；不把中心值/Gauss阶差当上界。主cb11表示尚未完整独审，这一依赖如实保留。

空间覆盖：从4x4 dyadic方格起，对三个positive-quadrant奇性角(1,0),(0,1),(1,1)各把触角子格递归26次，去掉边长2^-28的三个最末角square。剩247个不重叠boxes；每个原始角square被030半径1e−8帽覆盖（sqrt2*2^-28<1e−8）。full-face 48帽预算除24换single-positive-quadrant，即用030 all-face-cap上端/24，不漏全face因子，也不把E下端当actual integral下界。

初始：epsilon=1/64，T=16。10个dyadic time段从2^-6到2^4，先5段各二分，得15有限time段；第16项为精确peak大尾上界。247*16=3952空间/时间项≤4096。每box用030 r_Q下界，短time差A log(1+400epsilon²/r²)，A=55200/pi³。直接形成真实G上界并求I上界，逐项记录H区间宽、短time、大time、空间cell和帽预算。

预见的失败理由：030统一短time常数非常保守，epsilon=.01量级即使远离奇点也可能付数十Green误差，不与minimum的lower-H用途混淆。因此仅允许一次表示修订：改为全域L^(3/2)误差的Minkowski分账，而非逐box加短time误差。对同一空间247 boxes，epsilon=2^-26，15个time段ratio4从2^-26到2^4，加大尾同为16项，不增加分区数。全域短time界可准确积分成K epsilon²，K=1600pi A^(3/2) Gamma(5/2)zeta(3/2)/2；quarter/full-face换算明示。若新H区间仍过宽，停止并给所有remaining boxes及预算，不追加网格/精度。

另可准确检查条件目标：若最终可证I_H≤9/16000，这一Minkowski分账能否推I≤57/100000；这是条件式，不把未证I_H当前值代入。无新空间采样/成本优化。两表示均Arb/Acb128bit，全部有理endpoints；不依赖float判真。实际Python本机venv；MATLAB只交cb11_integral前缀便携解析预算校准脚本主串行，不自调、不clear/path。

输出：有限计算真上界/未闭合box列表、cap和各tail换算、严格条件证书或明确未达。只在完整I上界严格过目标才优先送独审；否则不得称突破或方法障碍。单批结束不结束研究目标。

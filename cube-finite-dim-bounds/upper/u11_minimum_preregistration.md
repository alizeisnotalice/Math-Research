# 11K minimum预登记

2026-10-25，只upper/u11_minimum_*，固定4D纯body，G>=21/25000于完整cube，parent负责I/tau/cost。不提交推送、不改他人/总账、不扩招。输入hash见u11_minimum_inputs.json。

先自核DD3 sign的mixed1/3及2/2正和分母简式，提供240 node interval API供积分代理共享；不能用point函数作box。完整cube按齐次和群对称归到x=(1,y,z,w),1>=y>=z>=w>=0，边界由真实G连续/方向line∞处理，不能预设min位置。

初始有限诊断最多2048个longest-edge dyadic tetrahedra，直接从ordered simplex剖分；若积分代理的准确node-sign polytope分区先已交且<=2048，可只读复用其几何代替重复建arrangement。Arb对node的affine hull取真实顶点范围，known-sign混合式给每cone enclosure；sign不定用branch hull，正负cross-sum含0必须报unresolved，不能填0。

只允许一次有明确失败原因的表示修订，预期若出现多个zero-projection导致单cone无限相消，则合并同spectral vertex的exact residues/符号正和式；不机械扩box。修订仅代数/表示，不启动I/cost/tau或任意点搜索。若不能闭合，完整给未覆盖域/失败cone和下一准确接口，不称g通过。

Python Fraction/SymPy+Arb实际运行固定公式与有限cell诊断。MATLAB交主串行，只给相同DD/geometry诊断。初始cell count、Arb精度、修订次数在结果逐项保存。finite samples不作为全域证书，任何新g_low必须完整证明。

初诊实际后唯一表示修订登记：共享CI4的351准确sign polytopes，Arb128得到248单cone lower pole、103 finite lower不足、0cert，cap由独立one-big-jump lower>21覆盖。完整准确失败域存initial_results，非g反例。修订只处理两个真谱顶点代表B,D的各48cone simple-residue star：canonical linear因子合并、exact多项式消因子后评估。最多两star，common degree若>32或无法形成可核因子就停止并报缺口；不扩原351cells、不重复I或参数搜索。

同一次表示修订内的固定完成步骤：两个star common degree26均已消因子（一次JSON fmpz序列化实现修补重跑），但star自身仍有artificial projection pole。为完成全和取消，固定12个B line reps+8个D line reps，分别按物理排序的3个body sign区域合并为3个global rational式；不增加点/box/权重。global common degree上限96，三式各一次gcd/factor，超过或无法形成可核表示即停。非顶点pole消失以polytope vertex triangulation/layer-cake证明，不擅加其他star模型。

# 第11轮最小复核入口（本轮仍在进行）

以下命令从cube-finite-dim-bounds目录运行，python指安装python-flint和SymPy的解释器。本轮实际Python3.15；准确版本见results/round011_runtime_versions.json。NumPy/SciPy用于候选生成，四维证书的验证路径不调用LP，也不要求安装它们。MATLAB实际运行记录保存在各批receipt；它们不替代连续证明。

## 已完成的严格接口

三维全cube最小值的独立完整fresh复核：

    python audit/033_check.py

该命令不使用缓存，重算全部59976节点，默认不写结果；--save才覆盖审计输出，通常不必使用。采用独立点值与多项式构造，但共享已审032的ray/rate与Arb内核。完整下端保存在033_results.json，向外可取g3=2621323989/500000000000；它不是真实minimum精确值。

四维已保存积分证书的完整验证：

    python src/verify_round011_n4_integral.py

入口检查锁定checker/engine/certificate哈希，将输入位置设为本项目中的全和表示与接口。它重核2440四面体、多项式系数、体积和角帽；不生成候选、不调用LP。--inputs-only仅作便携入口/来源/布局检查，不能作为积分通过证据。原始ci4b_affine.py保留作者绝对路径作为不可改写的历史锁源，外部复现请使用此便携入口。

准确目标成本和三维g结果集合核对：

    python src/round011_assembly_check.py

它不重算Green节点，也不独自证明一般L1/全部尺度。四维C4<14.753已由035完整独审接受并进入有效账本；该算术脚本本身仍不是完整证明。

## 三维完整积分：全批与最终独审均已完成

固定8256区域、degree32、epsilon9/5000、384bit，所有参数与源hash锁定。现有12part已全部exit0，8256块作者selfcheck与aggregate通过，不应再启动重复任务。

外部从零复现可将以下三个文件复制到独立输出目录再运行，避免覆盖存档证据：upper/continuous_bridge/cb11_cones.json、cb11_integral_cheb_check.py、cb11_integral_cheb_selfcheck.py。四块缓存完全可选，不是证明前提。示例串行全批：

    python cb11_integral_cheb_check.py --stage part --part 0 --parts 1
    python cb11_integral_cheb_selfcheck.py --parts 1
    python cb11_integral_cheb_check.py --stage aggregate --parts 1

亦可用固定整数P拆分part0至P−1，但必须全部结束后核全覆盖并聚合，不能只累加已完成分片冒充全积分。结果写入输出目录，源文件不变。主已导入精简证据，可用 python upper/continuous_bridge/cb11_integral_cheb_final_verify.py 无节点重算地核全部区域/端点/chord与准确聚合；已实际通过。最终036已完成独立凸弦/覆盖/聚合审计，主按清单导入。

## 分析证明链与范围

比较和正部源：audit/026、031；三维真实表示032、全域g033；短time与插值积分接口034已独审通过；四维完整候选035已经完成独审并归档。主组装见proofs/round011_candidate_assembly.md。核主导、全部奇异源、Lebesgue背景及一般L1/共同尺度传递不能由有限代数检查代替。

失败路线原始大量记录与可再生节点缓存不作为最终归档入口；必要诊断、准确输入和最终代码按白名单保留。所有有效界以results/bounds_ledger.json为准。

## 四维完整独审与锐化成本

    python audit/035_check.py

该入口独立重建1152锥、全域g系数及2440单元的4,932,474个非零残差；继承已审026/031分析链。实际约63秒，不调用作者engine或LP；共享SymPy/Fraction/flint/Arb与锁定多项式输入，明确不是完全独立算术内核。主已核10个交付文件及43个来源哈希，逐节阅读完整证明，MATLAB有限交叉检查实际exit0。

    python src/round011_sharp_rational_cost.py

这一标准库Fraction脚本用Machin交错级数证π<355/113，准确给出三维7.916及四维14.371候选成本幂的严格差额。它只验证成本；三维036及四维037已经完成完整链/锐化验收，有效上界依据这些完整审计而非该算术脚本单独成立。

## 最终三维与四维锐化审计

    python audit/036_check.py
    python audit/037_check.py

036默认只读仓内compact证据和已审依赖，独立核全部8256凸弦、面积、准确聚合、短时预算和平方成本，不读取外部raw或重算Green节点。--raw-archive是可选增强；本次独审实际使用它核对原12份保留字段。037独立逐单项式展开三个四维case，不调用作者compose/LP；继承035已完整接受的表示/积分/一般输入证明链，核全域minimum锐化和准确成本。二者Matlab回执均实际exit0。默认复核不覆盖已锁结果。

当前严格简单结论C3<7.916、C4<14.371；准确比较应使用账本的exact_comparison_threshold。用户要求的第二次全文复审尚待本轮归档后启动；上述命令是本轮初审入口，不伪称已经完成第二次复审。

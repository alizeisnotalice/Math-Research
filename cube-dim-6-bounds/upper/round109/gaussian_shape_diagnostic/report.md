# 第109轮：12进程正Gaussian shape诊断

仅原body32等权Cauchy上界模型。106封存文件未改。正式C6上界仍64；本轮无全域g、g0=1/65000、完整I或正S证书。下面全部实验数值是固定有限输入的MC诊断，不是区间积分证书或有限方差定理。

## 实际运行与输入

点冻结预程序运行一次，wall0.774554秒，0MC/旧kernel调用。Sobol5固定seed20261027的32候选中，首16均合格：分母4096、六个绝对坐标严格递减且全正；与8个指定旧表中544个canonical点无重复。这些点均在真实S2外，证明与点表见method.md、points.json。

主MC由主代理全文静审后唯一运行。65536个iid Dirichlet(1/2,...,1/2) shapes，32个固定seed分块，每块2048；12个spawn worker均实际完成chunk，64个符号反射对同shape求平均。总wall含worker清理0.430003秒，child CPU2.705421秒；chunk内CPU和1.342225秒（其余包含spawn/import等开销）。NumPy2.5.4，macOS Accelerate；每个child在导入NumPy前将自身BLAS环境限制为1线程。PID54732与12个worker54741–54752已退出，无常驻session。未设程序时限、未扩大样本。

两个等价常数表达式的浮点相对差为5.55e−15，K=0.08507935783685838；uniform shape A=1024I的解析sanity在每worker通过。所有shape保留。2个condition>10^10的shape在同一原始Gamma输入上做80位mpmath回算；无clip、丢样或重抽。最大condition约7.82294e10；最大保存solve residual约1.62176e−7。失败目录为空，32块checkpoint完整，state=COMPLETE。

## 经验结果

| 诊断量 | 保存结果 |
|---|---:|
| 新16点中最小的反射平均经验均值 | 2.336064748881957e−5，point4 |
| 逐shape有限16点最小值的经验均值 | 1.6702334882300774e−5 |
| 上两者之比（有限表的取min顺序损失） | 1.3986456177198623 |
| 原逐shape全cube corner floor经验均值 | 3.511394793926783e−7 |
| 16点经验方差缩减比plain/sym范围 | 约2.38至35.77 |

point4为(4096,4002,3950,3212,3079,759)/4096。其identity-sign均值2.3687613711881633e−5，反射平均均值2.336064748881957e−5，经验方差缩减约20.30倍。单样本最大贡献份额从约1.123%降到0.2604%，top16份额从约8.775%降到1.465%。32块反射均值范围[2.1690697600350547e−5,2.527528261178241e−5]，identity-sign块均值范围[1.701940212547937e−5,3.374330197320234e−5]。这些是组间波动记录，不是置信区间。

condition的中位数约89.68、99%经验分位约3.63403e5、99.9%约3.91504e7；超过10^10的2个shape已高精回算。在point4，这一最高condition组的保存贡献为1.32205e−14，当前有限批并未由它们主导。该观察不能证明未采到的边界shape尾部可忽略。

反射平均确实改善了这张新点表的数值稳定性；提前逐shape取有限点min仍损失明显。原corner floor经验值也远低于106所证约4.28601e−6的严格天花板。没有因此排除或认证真实原模型。

## 证明接口与一个明确反例

设F_P是64反射平均后的正shape integrand。真正的新floor需要K E[inf_cube F_P]的严格下包。本轮m_grid(P)=min_j F_P(x_j)满足inf_cube F_P<=m_grid(P)，方向是上包。因此1.67023e−5既不是g下界，也不是B_orbit下界；16点经验min同样不证明全域g。其高于g0不能作突破门槛。

反例无需新程序：在实际cone中的uniform shape P_i=1/32，A=1024I，F_P(x)=1/(32|x|^5)。真正方体最小值在全1 corner；本轮16点每个自由坐标都严格小于1，所以每个点值都严格高于真正最小值。连续性还给出一个正概率shape邻域，有限点min在其中仍严格高估。这说明不能用有限点或采样精度代替连续x的min认证。

可审计的下一接口是明确空间cell上的矩阵下包：对每个反射的二次型，在cell上用其顶点最大值M_{cell,epsilon}(P)，则

L_cell(P)=det(A)^(-1/2) * average_epsilon M_{cell,epsilon}(P)^(-5/2)

严格下包该cell上所有F_P(x)。完整cell覆盖后可使用K min_cell E L_cell作为全域g下包，或在shape内取min_cell得到B_orbit的下包。若只用一整cube作为cell，反射后的各个二次型最大值全相同，恰好退回106已排的原corner-floor表达式；细cell的作用因此有明确数学位置。此时仍需严格的空间覆盖与正shape积分下估；为了势的下包可以把未认证shape事件贡献置零，但必须准确支付保留事件的概率/积分。当前经验方差改善和condition分层只帮助选择这些事件及cell，不支付这些证明条件。下一步应先量化这个连续空间误差，而非盲目增加MC样本。

## 保存与复现

run_shape_diagnostic.py SHA256：17c2b5218632fa6ebc20bec02b7dbe372f5a27c6f3cfde7642fff1220b1a2555。
mc_inputs.json SHA256：10eaa573518b5f0cff06f268ded1b54ac8746f55543cb5646eae2b651d9dc9f4。
points.json SHA256：a084d8cd13d93f942cd7bf18d51e1124850eb7a7f5e636445bd2a55d027b2327。
results.json SHA256：4203a3ea95d97b3b1d508f05494108ad5ae117685a6496d0d646f6024bd44c96。

primary_mc_release.json绑定运行前源码/输入。run_receipt.json保存工具chunk25d2e3、CPU/wall与清理事实；32组NPZ+JSON保存moments、top16、condition分层、seed、输入哈希及3个极端shape输入。它们约1.3MiB，不保存65536×16×64全量反射值。workers目录保存12个child依赖/线程环境。程序拒绝覆盖已存在run；不能把本报告的复现说明当作第二次运行授权。未调用旧kernel、091基线、下界、MATLAB、Astra或子代理，未stage/commit/push。独立clone仍透明保留106时partial-checkout状态。

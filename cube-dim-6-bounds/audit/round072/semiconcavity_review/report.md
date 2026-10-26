# B072：normal deficit、六列族半凹性与整盒接口独立审计

**结论：PASS，保留作者的严格范围。** 三维 normal deficit 的范数身份、有限交叠的负加权范数充分条件，以及任意可逆六列 Cauchy frame 在真 S2 外的局部半凹性，证明均闭合。5项新独立小配置一次通过。本审计另外明确线性变换和顶点整盒下证书的系数。

没有证明一般 body32 半凹性，没有发现 Cauchy-frame 反例，没有全域 g/I/Cn 结论。一般 frame/body32 推广仍由 F076 负责，本审计未接管或重复。作者8列例只作交叠几何；抽象 log 例只反驳 leading-cusp 推断，均不越级。

## 1. 冻结来源与本轮实际工作

只读作者目录 `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round072/green_semiconcavity_probe/`。实际核对主指定 manifest 和 theorem，全部13个最终产物 pin 与8个上游 source pin 一致；科学运行前后再次核对未变化。

|原件|SHA-256|
|---|---|
|manifest.json|f0cf4a06301a2bfe4e72f1c4beb81935693658ad481a45887edbd0f66418355d|
|theorem.md|47ef45f244c361409a2182d55a2ab6a9bd78e537f684d39a4fa82ecc08ff55c9|
|report.md|823c6b5f7c539e9e3bf98aec2709a17bdc038985791a761654272677f92e1170|
|results.json|4798ae3ca8b096870fc8057e0cc2d1a480025ee76b4485e3033de8ad8886526d|
|checked_results.json|e95230955721e101ad7b59aab9f9fde84452e2df0254a5d72aef09d37c1fb814|

完整 theorem、report、两脚本、登记及保存审计均已读取。作者9配置、43 pair残差和唯一 axis Green 路径只作已保存证据绑定；本轮不重算它们，也不运行069旧配置。

本轮新配置使用各向异性 radial 三维模型、新联合权重、带不同 rates 的余项参数、新非正交六列矩阵和一个三维参数线性盒。实际0 pair残差、0 Green路径积分、0 body32 kernel。科学 PID81095、tool chunk `3cc407`、exit0；CPU0.595275秒，程序墙钟0.744805秒，guard含启动/回收0.848637秒，子进程已reap。CPU30秒、总wall180秒限制均满足，无重试。

## 2. 三维 deficit 的准确常数与范数性质

对真实三维对称 Cauchy stable 密度，Fourier差使用非负 integrand
e^(−tψ(ξ))(1−cos(ξ·z))。因此先 Tonelli 积 t，再转极坐标是合法的；没有交换条件收敛的有符号积分。t³积分给6/ψ⁴，三维径向积分为

\[
\int_0^\infty\frac{1-\cos(ar)}{r^2}\,dr=\frac{\pi|a|}{2},
\]

其非振荡 Laplace 表示同样全程非负。归一化严格为

\[
\frac6{(2\pi)^3}\frac\pi2=\frac3{8\pi^2},\qquad
N_\psi(z)=\frac3{8\pi^2}\int_{S^2}
 \frac{|\theta\cdot z|}{\psi(\theta)^4}\,d\theta.
\]

球面积分 ∫|θ·z|dθ=2π|z|，故 c≤ψ(θ)≤C 给作者的两侧常数 3/(4πC⁴)、3/(4πc⁴)，没有丢失2π。正权 absolute linear functionals 的积分是偶、凸、一阶齐次且非退化的范数，具有有限 Lipschitz 常数；负号后的分布 Hessian 非正。

独立配置1从径向 Fourier/Laplace 公式导出三维 isotropic density t/[π²(t²+r²)²]，再作 T=diag(2,3,5) 变换。非振荡有理分解给

\[
\int_0^\infty\left(1-\frac{t^4}{(t^2+r^2)^2}\right)dt
=\frac{3\pi r}{4},\quad
N_T(z)=\frac{\|T^{-1}z\|}{40\pi}.
\]

这独立核对了 normal-deficit 系数与 Jacobian；不是重复作者 axis-product slice，也不是完整六维 Green 求值。

## 3. 负加权范数必须在联合变量上控制

对所有切向、法向变量合在一起的 x，平滑化后准确有

\[
D^2(-A N\circ P)[h,h]
=-A D^2N[Ph,Ph]-2DA[h]DN[Ph]-N D^2A[h,h].
\]

A≥0 只使第一项非正；后两项不能略去。对于正交 P 和单位 h，它们统一不超过 2L||∇A||+M||D²A||。用正 mollifier 保持凸性、非负性和 Lipschitz 常数，M仅增加趋零项；局部一致极限保留分布 Hessian 上界。故作者在完整联合变量中的充分条件正确。有限多个项逐项相加即可，重叠、重复 span、不同 normal 尺度并不要求 disjointness；但准确分解和统一余项必须先证。

若直接使用一般线性映射 T_j 而非正交 P_j，同一证明的交叉项应写成 2L_j||T_j||op||∇A_j||；M_j仍是实际 N_j(T_j(x−x_j)) 的上界。不能在一般线性像下默认为 ||T_j||=1。

独立配置2取 A=2+(y+z₁)²+z₂² 与 N_ε=√(z₁²+4z₂²+ε²)。在 |(y,z)|≤1、0<ε≤1 上，M≤3、L≤2、||∇A||≤4、||D²A||≤4，故统一半凹常数可取28。全联合 Hessian 乘积式精确通过，norm Hessian 的方向分子为

\[
\epsilon^2(h_1^2+4h_2^2)+4(z_1h_2-z_2h_1)^2\ge0.
\]

ε=0、z₁>0 的 (y,z₁) 截面，在新点 (1/2,1/3) 的 Hessian 行列式是 −25/9，存在正曲率方向。这说明半凹结论包含受控的联合交叉项，并不宣称负加权范数整体凹。

作者 f^q 的条件式也正确：在0<m≤f≤M、局部 Lipschitz常数L、D²f≤C_fI且C_f≥0时，内缩球上平滑化、链式法则及极限给

\[
C_{f^q}=qM^{q-1}C_f+q(q-1)m^{q-2}L^2,
\qquad1<q<2.
\]

未提供这些定量数据时不能把该公式当作 body32 的面积分证书。

## 4. 精确轴族分解与 C² 余项的两个时间端点

off-axis-S2 恰意味着至少三个坐标非零。选其中三个为 y，在闭邻域上 |y_i|≥δ>0，余三个 z 有界。作者精确分解

\[
G(y,z)=F(y)-A(y)N_{axis3}(z)+R(y,z)
\]

来自完整 product 密度恒等式，而非仅渐近 leading term。F(y)=G(y,0) 的积分及前两阶 y 导数在 t→0 有界、t→∞ 为 O(t⁻⁶)，故 F∈C²。

对 R，(1−H_y) 以及一、二阶 y 导数在0均为 O(t²)。乘 t³后，对 z 求0、1、2阶的 integrand 分别为 O(t²)、O(t)、O(1)，包括所有总阶≤2的 mixed 导数。大时间利用 q_t 偶性 ∇q_t(0)=0 和统一 ||D²q_t||=O(t⁻⁵)，Taylor 公式给 q_t(0)−q_t(z)、∇q_t(z)、D²q_t(z) 都是 O(t⁻⁵)；乘t³后统一为 O(t⁻²)。A及局部导数有界，(1−H_y)及其 y 导数有界，因此大端 mixed 导数也可积。

这些是局部统一支配，允许先对截断积分求导，再以支配收敛得到联合 C²；没有在 S3 上预设导数或使用069仅off-S3结论。R的连续 Hessian 在所选紧邻域有界，故第3节足以给轴族半凹性。

独立配置3取 tangential rates (1,2,4)、normal rates (2,3,5)，y∈[2,3]×[4,5]×[6,7]、|z|≤2。精确验证10种 y 导数的 (1−H)/t² 有限极限，保存小端联合导数的有理 dominator。normal density 的全局 Hessian 常数可取

\[
\|D^2q_t\|\le \frac{811}{27000\pi^3}t^{-5}.
\]

该常数来自完整正 product 的 score 和 operator-norm 控制；大端差/梯度/Hessian的常数相应为2倍、2倍、1倍。所有小端指数2、1、0和大端指数−2均准确可积。没有进行新的 Green 积分。

## 5. 可逆六列传递和顶点整盒下证书

对可逆六列矩阵 N，密度推前的 Jacobian 给

\[
G_N(x)=|\det N|^{-1}G_{axis}(N^{-1}x).
\]

真 S2 同样由 N 推前。若轴侧对应邻域有 D²G_axis≤C_axis I，则物理空间上可取

\[
C_N=|\det N|^{-1}C_{axis}\|N^{-1}\|_{op}^2,
\]

或以 Frobenius norm 替代得到明确保守值。Lipschitz 常数只乘一次 ||N⁻¹||和同一个 Jacobian；若直接传递 G^q，则 Jacobian 为 |detN|⁻q。不能把这些因子省略，也不能将非可逆矩阵当六维 density 变换。

独立配置4采用新的非正交 rational/integer 六列矩阵，|detN|=30030，N⁻¹x₀=(2,3,6,0,0,0)，x₀=(7,15,30,0,0,0)。唯一坐标表示证明其off真S2，无需pair枚举。二次函数的 Hessian pullback 精确通过；每单位 C_axis 的 Frobenius 保守 multiplier 是 `423631591/27081081027000`。这只是变换系数，不是某个 Green 已知的实际半凹常数。

一般线性盒令 x=c+Ls，|s_i|≤h_i，且**整个凸盒及所需邻域**均有同一 D²f≤CI。用各坐标的 product barycentric weights λ_v，使 ∑λ_vv=x，则

\[
f(x)\ge\sum_v\lambda_v f(v)
 -\frac C2\sum_i(h_i^2-s_i^2)\|L_i\|^2.
\]

证明是对凹函数 f−C|x|²/2 应用 Jensen，再用独立符号权重的准确协方差；交叉项消去。若顶点有严格下界 ℓ_v，整个盒因此有

\[
f(x)\ge\min_v\ell_v-\frac C2\sum_i h_i^2\|L_i\|^2.
\]

轴平行盒 L=I 对应系数 **C/2乘半宽平方和，或 C/8乘全边长平方和**。此处不需要 f 可微，也不要求 S3 避开盒子，只要求该函数的统一半凹性已经在整盒上成立。仅顶点off-S2、仅局部点wise常数或未证明的body32半凹性，都不能满足前提；不能让盒中穿过已移除S2的线段混入证书。

独立配置5在新的6×3非正交线性盒、8个顶点上验证完整协方差和二次函数等号。其最大 penalty 每单位C为 `7027/29645`。这个等号检验确认了1/2系数，未给任何实际 Green 顶点下界或全域g证书。

## 6. 必须保留的范围边界

作者8列 frame 的两个 rank3 spans 交维2且含一个off其真S2的点，严格说明有限frame可能有此交叠。该例没有计算此8列frame的Green/Hessian；几何本身不证明半凹性，也不证明body32具有同样交叠结构。

抽象函数 −|s|+s²log(1/|s|) 的负导数jump不能抵消邻近非零点的无界正经典二阶导数。因此它准确反驳“负leading cusp加off-stratum C²就足够”，并不构成Cauchy-frame反例。原040的正J₂发散也仍不决定有符号Hessian。

一般body32还缺实际多strata分解、非负且统一C²的振幅、所有交叠及相对normal尺度上的统一余项控制，以及可计算常数。独立六列 proof 的切向/法向 product 结构不能直接迁移到额外耦合列。F076继续处理这些问题。

本轮最终判定：**PASS 作者三项受限局部定理及上述条件整盒接口；无body32半凹性/globalg/I/Cn结论。** 所有新写入只在本审计目录；074、069及原作者来源冻结未改。0旧路径/43pair/069重跑、0body32kernel、0grid、0MATLAB、0递归、0stage/commit/push。完整准确结果、源绑定与真实收据已封存，主代理保留最终证明及集成责任。

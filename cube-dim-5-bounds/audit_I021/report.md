# I021：四子集完成项、真实变量 scope 与固定 type 权重独立审计

## 结论

**通过。** upper_I020 的分组 Cauchy–Binet、分数 Hölder/Finner 及
Mellin–Gamma 链条成立。主 round008 指定的唯一固定有理权重

    w = (1017/2000, 3703/10000, 303/2500)

也满足完整上界资格，并严格改善默认 I020 显式常数。可用严格有理界为

    I < 142142151/10000000000000 = 0.0000142142151.

该结论使用已审的同模型正 Green/轨道 Minkowski 接口，控制完整全时间、
全正 unit face 的辅助积分 I，没有截时、空间采样或删除退化形状。它没有
给出新的全域 g，因而不改变 C5<=32，也不是 C5 突破。

独立性的关键是：本审计没有重跑作者 1820 个四子集枚举，而从四阶
Hadamard minor 的整除性、法向零点集的协方差主子式及法向轨道入手，
完整重建所需的加权几何与每一个独立变量的 scope 占用。数值式也采用
未移位 Mellin Gamma 和 I 的七次幂/有理根包围，与作者计算路径不同。

## 锁定与依赖

等待生产锁后才读取候选哈希。upper_I020 锁时间为
`2026-10-25T18:12:03.746665+00:00`，manifest SHA-256：
`f6d6a51c5b5da809b92861a2b0dfcb955921ebf87ca68386eda23b91863dd32e`。
四项生产白名单均匹配 manifest；报告 SHA-256 为
`b186d9db778602e6e3c113080981657f050e9e580827eab3c361ebdad112ce63`。

固定权重源码和登记分别锁为
`e674656c1a4d21579e12d8e7b690fa017da38bbc90bdfb0d04e01879b5426049`、
`0df6f5f74278772fd3806a32d607545d440f96b04cc6871cf141586d1db48908`；
保存结果锁为
`cc129596e5ed3b037809689cb357d87c52baf437c397282b37aa32ae4091de98`。
14 项精确来源见 registration.json；来源不在运行中修改。

使用的上游接口是 round006/orbit_minkowski_proof.md 及已通过的
audit_upper019/report.md 中同一组归一化：16 条未单位化方向
`v=(1,+/-1,+/-1,+/-1,+/-1)`，独立标准正 stable 变量
`E exp(-lambda S_i)=exp(-lambda^(3/4))`，Lebesgue 测度，q=10/7。
本审计不重开旧轨道常数或旧 g 的计算，只核对新负矩如何接入该接口。

## 1. 按四子集分组的精确恒等式

令 V 的列为上述 16 个方向。符号求和直接给 VV^T=16 Id，故
W=det(VV^T)=16^5。对四子集 B 取其整数 cofactor normal n_B。
rank(B)<4 时所有 cofactor 都是零；其余情形定义

    w_Bi = (n_B dot v_i)^2,
    L_B = sum_i w_Bi S_i,
    C_B = sum_i w_Bi = 16 ||n_B||^2,
    theta0_B = C_B/(5W).

Laplace 展开给 `w_Bi=det(V_(B union {i}))^2`，选中方向本身权重为零。
每个五方向基恰由其五个四子集完成一次，因此

    det A = (1/5) sum_B (product_(j in B) S_j) L_B.

在 S_i=1 处，det A=W，右侧为 sum C_B/5，故 sum C_B=5W。
亦可不使用这次校准，直接由所有四阶 principal minors 的
Cauchy–Binet 求和得到

    sum_B ||n_B||^2 = sigma4(VV^T)=5*16^4,
    sum_B C_B = 5*16^5.

因此 theta0_B 是真实概率权重。置

    F_B = (product_(j in B) S_j) (L_B/C_B),

则 `det A = W sum_B theta0_B F_B`，加权 AM–GM 方向为

    det A >= W product_B F_B^(theta0_B).

没有遗漏 1/5，也没有将某些正完成项丢弃后冒用原归一化。

## 2. 不枚举四子集的完整 normal 分类

每个 4x4 +/-1 minor 都被 8 整除：减去一行再从其余三行各提出 2。
Hadamard 行范数界又给绝对值至多 16。因此非零 cofactor 分量只可能
是 +/-8 或 +/-16。除去最大公因子后，primitive normal 分量位于
{0,+/-1,+/-2}，且至少一个分量绝对值为 1。

只需检查 15 个排序绝对值模式。对每个候选 n，令

    H_n = {v in V : n dot v=0}.

一个 normal 可由四个独立方向产生，当且仅当 span(H_n) 的维数为 4。
所有非零候选都满足 rank(H_n)<=4，所以此判定既必要又充分。精确有理
行消元给出仅有三类：

| primitive normal 类型 | projective 法向数 | H_n 大小 | scope 大小 |
|---|---:|---:|---:|
| (1,1,0,0,0) | 20 | 8 | 12 |
| (1,1,1,1,0) | 40 | 6 | 14 |
| (2,1,1,1,1) | 80 | 4 | 16 |

法向数来自全部坐标置换及改号，最后识别整体正负。它们分别也等于
`choose(5,2)*2`、`5*8`、`5*16`。法向变换置换 projective 方向，
因此保留 dot-product 平方、rank 及所有加权包含关系。

## 3. 用 sigma4 重建每个 normal 的 C 总量及真实 scope

固定 primitive normal n，令 A_H=sum_(v in H_n) vv^T。所有 rank4
四子集 B⊂H_n 的 cofactor 都平行于 n；反之 normal 为 n 的四子集全
包含于 H_n。四阶 principal-minor Cauchy–Binet 因而给

    sum_(normal(B)=n) ||n_B||^2 = sigma4(A_H),
    sum_(normal(B)=n) C_B = 16 sigma4(A_H).

本审计用 Bareiss 整数消元计算每个 A_H 的五个四阶 principal minors，
得到以下精确表：

| 类型 | sigma4(A_H) | 单个法向的 sum C_B | 全类型 sum C_B | 总 theta0 |
|---|---:|---:|---:|---:|
| (1,1,0,0,0) | 8192 | 131072 | 2621440 | 1/2 |
| (1,1,1,1,0) | 3072 | 49152 | 1966080 | 3/8 |
| (2,1,1,1,1) | 512 | 8192 | 655360 | 1/8 |

总量为 5242880=5*16^5。这条重建不使用四子集数量或作者几何流哈希。
作者所报未加权数量 1160/480/80 未在本审计中重新计数；新界只需要上表
的加权总量，已被独立精确重建。

真实组函数只依赖

    A_B = B union S_n,   S_n = {i : n dot v_i != 0}.

B 与 S_n 不交。因此每组 scope 的大小恰为 4+|S_n|=20-|H_n|，
而不是只取 L_B 的非零支撑，也不是把不同 L_B 当作独立变量。

为独立检查 B 部分的加权包含量，对 i∈H_n 使用 rank-one 恒等式

    sum_(B subset H_n, i in B) ||n_B||^2
      = sigma4(A_H+v_i v_i^T)-sigma4(A_H).

右式对每个 principal determinant 都是 rank-one 更新的一次差。对三个
canonical H_n 中全部 18 个方向作精确核验，均等于
`(4/|H_n|)*sigma4(A_H)`；法向轨道变换将这些身份传给全部法向。
因此在一个 normal 组内，某变量若在 S_n 中就贡献全部该组 C，若在
H_n 中则贡献该组 C 的 4/|H_n|。

程序由此对 140 个 projective normal 逐一累计 16 个变量的加权 scope，
得到每个类型的条件占用分别恒为 12/16、14/16、16/16。故默认总占用
对每一个变量均为

    rho0 = (1/2)*(12/16)+(3/8)*(14/16)+(1/8)*1 = 53/64.

这是逐变量的精确恒等式，不是只检查一个平均占用。

## 4. 分数 Hölder 在原始独立 S 上合法

在独立产品概率中，若非负 f_B 只依赖 A_B 且 h_B>=0 满足
`sum_(B:i in A_B) h_B <= 1`，则

    E product_B f_B^(h_B) <= product_B (E f_B)^(h_B).

可逐坐标证明：固定其余坐标，对涉及当前坐标的函数用指数 1/h_B 的
广义 Hölder；不足 1 的份额用常数 1 补足。得到各函数对该坐标积分后
的同一乘积形式，再处理下一坐标。积分掉坐标不会扩大其余 scope。
对非负函数先逐个截断，再单调取极限即可处理无界情形。

取 h_B=theta0_B/rho0、f_B=F_B^(-b*rho0)，b=3/20，所有坐标指数和
恰为 1。不同 L_B 之间相关完全不妨碍该结论。本论证没有在倾斜形状
概率 R 下假设独立性；它只作用于原始 16 个独立 S_i。

## 5. Mellin–Gamma 与尺度抵消

对独立标准 positive beta-stable，beta=3/4，Tonelli 与 Mellin 积分给

    E S^(-s) = 1/Gamma(s) integral_0^infty lambda^(s-1)
                                      exp(-lambda^beta) d lambda
             = Gamma(s/beta)/(beta Gamma(s)),    s>0.

该公式同时证明所有所需负矩有限。若 X=sum_i p_i S_i，Laplace 变换为
`exp(-R lambda^beta)`，R=sum p_i^beta，则

    E X^(-s) = [Gamma(s/beta)/(beta Gamma(s))] R^(-s/beta).

在同一组内，B 的四个 S_j 与 L_B 独立，因为其完成项权重严格为零。
故 `E F_B^(-s)=m(s)^5 R_B^(-s/beta)`。三个 R 只依赖 normal 类型：

    R1=8^(1/4),
    R2=1+1/sqrt(2),
    R3=(7+8sqrt(2)+3sqrt(3))/2^(15/4).

独立 normal-dot 概率多重集分别为
`{(1/8):8}`、`{(1/16):8,(1/4):2}`、
`{(1/32):7,(1/8):4,(9/32):1}`，直接推出这些闭式。

由 T 的 Laplace 变换 exp(-16 lambda^beta)，

    Z=E T^(-beta)=Gamma(1)/(16 beta Gamma(beta))
                 =1/(16 Gamma(7/4)).

R 是以 T^(-beta)/Z 倾斜再推送的概率。d=det(A/T)=det A/T^5，
而 5b=beta，所以精确得到

    J=E_R d^(-b)=Z^(-1) E(det A)^(-b).

这里没有要求 T 与 S/T 独立。结合 AM–GM 和上节分数 Hölder，默认式是

    J <= 2 Gamma(7/4) m(b*rho0)^(5/rho0)
                         product_t R_t^(-theta0_t/5).

rho0=53/64 给 `b*rho0=159/1280`、`b*rho0/beta=53/320`，所以移位
Gamma 形式的参数恰为 373/320、1439/1280，外幂为 320/53。R 的指数
是 -theta0_t/5，因为 rho 在 `(-b*rho/beta)*(theta0/rho)` 中抵消。
作者上述指数均正确。

负矩有限性覆盖全部 stable 时间与任意退化形状邻域，同时闭合上游
Minkowski 可积性条件，不需要额外短时或长时预算。

## 6. 固定有理 type 权重扩展

设 theta0_t=(1/2,3/8,1/8)。主指定扩展是对每个四子集 B 取

    q_B = w_t * theta0_B/theta0_t,   B 属于类型 t.

这是 type 内仍按原 Cauchy–Binet 权重分配，绝不是按该类型的四子集
个数平均。Σq_B=1，且

    KL(q||theta0) = sum_t w_t log(w_t/theta0_t) =: KL(w||theta0).

对恒等式 `det A/W=sum theta0_B F_B` 换用 q_B 作 AM–GM，得到

    det A >= W exp(-KL(w||theta0)) product_B F_B^(q_B).

因此负矩前因子须乘 exp(b KL)，方向为增大；主公式符号正确。
type 内比例缩放保持对称性及条件占用，故每个原始变量的占用是

    rho(w)=(12w1+14w2+16w3)/16.

对指定 w，Σw=1、各项严格正，且每个变量的精确占用都是
66127/80000<1。对 h_B=q_B/rho(w) 重复相同 Finner 证明，得到

    J(w) <= 2 Gamma(7/4) exp((3/20) KL(w||theta0))
             [Gamma(1+rho/5)/Gamma(1+3rho/20)]^(5/rho)
             product_t R_t^(-w_t/5).

这个固定候选不需要全局最优性。本审计没有运行 SLSQP，也没有核算
保存结果中的另一个 floor 权重候选。

## 7. 完整 I 接口与实际独立数值证书

上游已审接口为

    I <= c^(10/7) (pi^2/60) J^(10/7),
    c=3/(2^(15/2) pi^(5/2)).

本审计先算未移位 Gamma，再用不同代数路径

    I^7 = 3^10 J^10/(2^75 * 60^7 * pi^11)

计算 512-bit Arb 球。最后将七次方球转换为分母 10^105 的有理根端点，
用整数七次幂比较验证下端在球之下、上端在球之上，没有使用浮点根。
默认与唯一固定权重结果分别包含于生产所保存的严格打印球。

| 显式上包围 | J 展示值 | I 展示值 |
|---|---:|---:|
| 默认 I020 | 1.51039608935354144429 | 1.42147007315294853246e-5 |
| 固定新权重 | 1.51035996772191164671 | 1.42142150920757661775e-5 |

精确有理端点比较确认固定新权重的 I 上端严格低于默认 I 下端；默认
I 上端又严格低于保存的旧 round006 球下端。旧 Gamma 基线没有重算。
另确认固定新权重的上端严格小于 142142151/10^13。完整端点、逐变量
占用、Gamma 参数和球均见 results.json。展示值不参与判定。

实际唯一审计执行 PID=28799，Python 为
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python`，
python-flint 0.9.0、Arb512 位，耗时约 0.0358 秒，退出码 0。
没有导入或执行作者程序，没有枚举任何四列子集、重启旧完整基线、
启动参数搜索、采样或占用 MATLAB。

## 保存、复现和剩余范围

本代理仅写独立 clone 的 audit_I021/。主和作者目录只读；不提交、
不推送，不修改其他审计结果。白名单为本报告、registration.json、
check_grouped_bound.py、results.json、manifest.json。

检查脚本只按显式 --project-root 解析登记的项目相对路径；
registered_origin 只作来源记录，不用于读取。默认只核 hash，不运行
数学计算；只有 --run 才执行本次登记的新审计。因此归档到具有已登记
来源目录的仓库后可直接复现，无需原来的外部 clone：

    python -B cube-dim-5-bounds/audit_I021/check_grouped_bound.py --project-root /checkout/cube-dim-5-bounds
    python -B cube-dim-5-bounds/audit_I021/check_grouped_bound.py --project-root /checkout/cube-dim-5-bounds --run

运行时临时 input view 只是指向已锁来源目录的 symlink 布局，不是归档
依赖；归档后使用完整项目根即可。

结论允许更新同模型完整 I 上包围，但没有证明 I<=11/10^6，也不能由
目前显式上包围大于该目标反推真实 I 大于目标。新的 g 或 C5 改善仍是
独立未完成任务，本批不触发 C5 双侧突破后的额外三轮。

# 152：150 时间分布积分机制的独立理论与实现准备审计

结论：在实际选定的 64 个未展开 root、同一 anchor 的完整 73 个正时间 band、已资格的 111/112 数据下，150 的机制方向正确。未发现公式反例。本目录目前只完成理论审阅、原始来源身份门和代码准备；没有执行新的时间 factor、mean、gain 或 UH 数值批。旧 145/141 及其他冻结目录未修改。

本审计直接读 029、111、112 的证明和实际 production_v2.py，并回读原 112 rows、certificate、numeric hierarchy、geometry 与 64 个 raw packet。完整 SHA 清单见 `input_locks.json` 和来源登记；不导入生产 consumer 或以作者通过标签替代公式审阅。旧 H/R 和完整 112 范围/积分资格由已冻结 audit121 继承，不重算。

## 1. 时间单调性、反向 Chebyshev 及无穷尾

029 的模型为五维 alpha=1、16 条未单位化 body 方向、共同时间、Lebesgue 和逆 Fourier `(2π)^-5`。其 §8 给同一个密度的正 Gaussian 混合

`p_t(c)=(4π)^(-5/2)t^(-5) E[(det A)^(-1/2) exp(-c^T A^(-1)c/(4t²))]`，

其中 A 的分布不依赖 t，且 A 几乎处处正定。对固定 c，每个非负被积函数随 t 非减，因此 `q_c(t)=t^5 p_t(c)` 非减。此步不交换导数与期望；极端接近退化的 A 也没有被删去。029 §2/§8 的正时间有限性及原点 Gaussian 期望给 `q_c(t)≤p_1(0)<∞`。结论包括 S2 上的 c；S2 的 Green 奇异性不影响正时间密度或 H512。

对 `0<a<b≤∞`，令 `dμ=t^-5 dt/Z`，`Z=(a^-4-b^-4)/4>0`，无穷端按零解释。设 `f(t)=F(1/t)`，F 非减，则 q 非减而 f 非增。对两个独立 μ 变量，

`2(Eμ[qf]-Eμ[q]Eμ[f]) = E[(q(s)-q(t))(f(s)-f(t))] ≤ 0`。

所以 `∫p_t(c)F(1/t)dt ≤ m Eμ[F(1/t)]`，其中 `m=∫p_t(c)dt` 是该 anchor 的真实 band 质量。这里的归一化是 t^-5，绝不能改成均匀 dt 或 p_t dt。q 有界，F 在 `[0,1/a]` 有界，且 Z 有限；故包含 b=∞ 的全部积分和双重积分均合法。没有额外 Fourier 或 Gaussian 常数，因为 q 的定义使常数在同一 m 中完全保留。

## 2. R/T 两条弦及端点方向

从 111 的逐固定 t Taylor 公式和完整平移比可分别得到空间均值上界 `p_t(c) F_R(1/t)` 与 `p_t(c) F_T(1/t)`，其中

`F_R(z)=R(hz)`，

`F_T(z)=1+μ∞ z+(25/24)∫W(ρ) z² R(ρz)dρ`。

R 是 111 §5 的正、递增、凸函数，W≥0，μ∞≥0。`z² R(ρz)` 的二阶导为 `2R+4ρzR'+ρ²z²R''≥0`，故两支均递增且凸；有限 h 允许在任何紧 z 区间用连续导数支配。R(0)=1，且 z²J(z)→0，所以无穷时间端两支准确等于 1。这里没有把 p_t 或固定 epsilon 的 H 当作凸函数或齐次函数。

换元 z=1/t 后，μ 的权为 z³dz。令 `za=1/a`、`zb=1/b`，则端点凸弦的左端系数为

`w_a = [∫zb^za z³(z-zb)dz] / [(za-zb)∫zb^za z³dz]`，

它等于 150 的 theta；另一权为 `1-w_a`，两权严格正且和为一。b=∞ 时代入 zb=0 准确得到 `w_a=4/5`，因而完整最后时间尾已支付。

必须先分别形成 `w_a R_a+(1-w_a)R_b` 和 `w_a T_a+(1-w_a)T_b`，再取二者及旧 factor 的 min。不能把端点 min 当成一个凸函数。每个端点的 R/J 上包都只作正权替换。右端上包还可与同支左端上包取 min：真实 F(zb) 同时不超过这两个上包；不要求保存的上包数值本身单调。192 位 ceil 后不增由旧 factor 处在同一 2^-192 网格保证，未来 checker 将逐项核该网格身份。

## 3. exact 弦与旧 224 位上包弦必须区分

令 `φ(x)=x^(5/4)`，旧真实范围为 `[L,U]`、`0≤L<U`。exact 端点弦 A 的斜率满足

`s=(U^(5/4)-L^(5/4))/(U-L) = U^(1/4) [1-(L/U)^(5/4)]/[1-L/U] ≥ U^(1/4)`。

因为 `r^(5/4)≤r` 对 0≤r≤1 成立。实际已核 `U_i≥2^144`、网格 N=2^160，所以 `U≥2^-16`，从而 `s≥1/16`。

**所需下界属于 exact 弦，不是分别向上舍入 224 位端点后所得弦的斜率。** 112 的实际代码是 `new_selfU=min(old_selfU, weight*ceil(N*rounded_chord/2^(4level)))`。因此只有实际 `new_selfU<old_selfU` 时，才可确定保存的新 selfU 来自 rounded chord，并得到 `saved_selfU/N ≥ orbit·volume·A(M_old)`。150 的全部 64 个所选项实际严格满足这一门。

若新的合法均值上端从 `M_old` 降为 `M_new`，对 exact 弦有

`E_B φ(H512) ≤ A(M_new) = A(M_old)-s(M_old-M_new)`。

所以可从保存的上包减去 `orbit·volume·(M_old-M_new)/16`。这个证明使用 old rounded chord 高于 exact chord 的事实，不需要计算或下界 rounded chord 的斜率，也不需要重算任何旧 power。

在 N 网格，`Δ_i=m_old_i-m_new_i≥0`，保守可减整数为

`G_i=floor[orbit_i Δ_i / 2^(4level_i+4)]`。

floor 方向使减量不大于已证安全减量。随后保留旧下界并核 `new_selfU≥old_selfL`；它是数据一致性门，不是本步减量证明的替代。

## 4. 根贡献、体积和完整 UH

直接全域减量需要本次项确为未展开根。如果是有 ancestor 的普通 child，ancestor 的 self arm 可能阻止全域下降，不能直接从全域减去其 self 变化。本次实际 64 项均为 112 根且无 child，`effU=selfU`。它们是旧合格完整 root 覆盖中的互异物理域，体积准确为 `2^(-4*5)`，orbit 由 indices 的重复度独立核为 `4!/∏mult!`。因此可从旧完整 UH 逐根减安全 G_i/N，其余域完整保留旧支付；没有新覆盖、空间细分或双重收费。

UH 指 `∫[0,1]^4 H_(1/512)(1,y)^(5/4)dy` 的完整上包，不是完整 G 的 I。此机制不改 029 模型、一般 L1 或所有坐标共同的正半径接口。未来即使 UH 严格改善，仍需另经主接纳与同模型短时合成；本审计不计算 M/JM、范数合成或 Cn。

## 5. 已实际完成的来源门

预登记 `source_registration.json` 后，实际 Python PID **66554**，exit **0**，约 **0.705 秒**。只读取原保存字段和元数据：

- 全 97972 个原 112 node/certificate/numeric hierarchy 的身份及 94580 物理 frontier；
- 实际 selected64 正是既定 eligibility 下按保存 gap 排序的前 64 项，没有换根；
- 64 项均为 level5/g002/无 parent/无 child 的 root；current effU=selfU=112 new_selfU<old_selfU；
- 保存 L≤M112≤U、L<U、U≥2^-16；旧 selfL/U 和 packet hash 对应原记录；
- orbit、由真实 indices 与 anchor 推出的全部 offset、h、μ∞、同一 73 时间端点、原 packet point/key/hash；
- 全部 H 输入是 dyadic 有理数，适合独立整数对齐求和；
- 97 个源 SHA 在检查前后相同，0 新 factor/mean/H/R/旧 power。

原完整 UH 同时锁到 112 hierarchy/preflight 与 audit121 whole-results。旧 H/R 数值资格继承明确的 audit121 成功来源，不重新运行其 7,151,956 bands 或旧 224 位/256 位 power。

## 6. 待主审的不同算术路径

`audit152.py` 不导入作者模块，`rational152.py` 使用显式整数/gcd 有理数，不使用 Fraction 或 root 函数。73 个新系数通过积分 z³ 与 z⁴ 后形成 barycentric 权，而非调用生产 theta。64 个均值把 raw H 的 dyadic 分母先统一成同一二进制网格，逐 band 用整数质量差乘 192 位 factor 整数，再一次整数除法确定 160 位向上量化；区别于生产 Fraction 全有理求和。gain 用整数位移核验，完整每根 prefix、最终 UH 和作者结果逐项比较。

目前 150 尚无本审计锁定的实际输出，`candidate_binding_template.json` 明确 execution_ready=false。未来必须用真实 producer PID、完整 factors/nodes/results/state/registration、最终 manifest 的精确 SHA 生成新的 ready binding，再由主授权准确 plan+binding；不会自动扫描 future、追加根、生成新 R 或重跑旧 power。准备只执行 AST 语法解析，不运行上述新算术。

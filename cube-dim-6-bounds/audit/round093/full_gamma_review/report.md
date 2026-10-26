# B093 独立审计：完整 Γ 账目与所选上界表达式的不可用性

**结论：A093 的完整局部上 Hessian 账目、固定盒/全部允许 T 的表达式障碍，以及条件正窗口比较均通过独立审计。建议主任务接受。** 未发现遗漏的076余项、空间mixed项、时间端点、归一化或适用域缺口。主任务保留最终接受权。

严格通过的是：在指定原始body32模型、固定rho/R/盒Q和0<T≤1/64下，一个完整且极粗的上界Γ_env有效；这套上界产生的指定顶点下证书不可能通过。**不是证明真实最小半凹常数很大，不是证明真实曲率很大，不是排除模型、候选g0或所有其他证书。** 实际G/F端点与梯度、全域g_L、完整I_U和最终S仍未取得。

本轮只在独占本地目录写入；未stage、commit、push或改动作者/主目录。按主任务转达的最新要求，仅完成当前已登记审计并封存，之后停止，不启动新轮次、子任务或代理。

## 1. 来源和有限审核登记

作者manifest锁为`e6bb3214a05c3e534a7e135c3a930f99cb4ff6fcc89c0b38c3a2ab6c8b41831b`，证明锁为`f0b558d683697b4133569918782e1ca7f5eee0f25b269bcf46da2274774456d7`。作者12文件（含manifest）及12项其他只读来源共24项完整绑定，具体路径/散列见`source_bindings.json`。

来源包括076、085、086、主unique_U附录、057，以及086最新独审report/manifest。086独审manifest为`03ad58a4cd0051164c0cd1f2e56f33926f36609357e22b9135b87cec2d6d6308`，其被审作者来源与093引用完全相同。它的独审通过可支持本次新账目审查，但不被伪称为当时已审过093，也不追溯修改作者原登记中的pending状态。主任务仍决定最终接受状态。

完整静态审查覆盖093 proof/report、账目checker、085提取、057原证明和适用域。先写`audit_plan.json`登记至多3个新纯算术配置，再冻结一个最小saved关系证书和guard。仅运行这一个新配置：验证保存数据之间的关系与精确不等式见证，不重建作者八配置，不执行作者checker，不求新T、不扫描、不调用核。

## 2. 精确分解及12项覆盖

076在本球上的精确分解是

    G=R4+R5−Rcut+Rquad−Rsigned−a_U N+Rnormal+H_T.

主unique_U附录保证只有一个incident rank3 span；其他rank3及rank≤2 pattern在整个球中为0。ρ=1/512、R=1/41472与085一致。球内dist(x,true S2)>1/4，可由sqrt(1/8)>1/3且1/3−rho>1/4直接认证。五维盒h=1/2048满足5h²≤rho²。全盒超出正单位面的部分不能计入I；作者明确使用全盒作更强局部证书，并只把正面交集体积计为2h⁵，未偷换积分域。

各预算覆盖如下：

| 账目项 | 实际被支付函数/项 |
|---|---|
| rank4_aggregate | 0到T的全部rank≥4 pattern |
| unexpanded_cutoff_removal | 对有界未展开a_ext去掉全部外cutoff的误差 |
| weighted_quadratic_remainder | 二次Taylor误差的完整物理Hessian |
| t5_inplane_correction | t⁵r_ext及其cutoff导数 |
| negative_cusp | −a_U(y)N(z)的上Hessian |
| signed_renormalized_F_curvature | 全时间减线性项后的angular profile的纯z Hessian |
| signed_upper_tail_curvature | 从T到无穷的signed重整化上尾 |
| signed_coefficient_cross_including_linear | Da_U·c_j的一阶系数导数与整个有限L_j的梯度 |
| signed_coefficient_second_including_linear | 系数二阶导数与整个有限L_j的值 |
| normal_layer_constant | 变化振幅a_U(y)乘Tq_1(0) |
| normal_deficit_upper_tail | a_U(y)乘虚拟三维deficit上尾的值/梯度/Hessian |
| original_full6D_t_ge_T | 原六维p_t从T到无穷的真实大时间尾 |

最后两种尾完全不同，不能互相替代；作者两者均支付。所有signed remainder都付绝对Hessian界，故原分解中的正负号不造成漏项。只有负cusp使用有效上半凹界。12项共20个非零系数，合并成六个严格正的T幂系数。

## 3. rank4、cutoff、quadratic和t5核查

**rank4统一body界通过。** ±1的6×6可逆矩阵det是非零32倍数；5×5cofactor是16倍数，Hadamard使其模小于64，因此≤48。故inverse每项≤3/2，四个tail行平方和≤54、两个completion行≤27。对所有至多35960个J及其086固定completion都成立，不需要样本矩阵或不完整records。代入086的完整sum定理，T、T²、T³系数分别为27(b2+b1²)、108a1b1、18(a2+3a1²)，共同因子35960/(32·729·R⁸)，与代码一致。

**cutoff-removal通过。** 每个telescoping项选含large标签j的固定商基，Λ≤9/4。外切向列模≤sqrt3，故C_J的Frobenius范数≤3；对应inverse的Frobenius范数≤sqrt2，切向transport≤3sqrt2<9/2。对单位完整物理方向(p,q)，振幅参数变化≤sqrt(1+(9/2)²)<5。作者V=5覆盖切向、法向及混合方向，不是只付法线。

未展开a_ext的两阶导数、一次振幅一次density和两次density三种Leibniz项，分别积分得到T³/3、T²和T。mixed项中的2与积分1/2相消，式(7)系数准确。large因子及其前两导数均保留tail支撑；把其不需要的t⁻¹补成更粗上界只会增大预算。此阶段没有weighted列与large列必须同时入基的错误要求。

**quadratic完整物理导数通过。** 在固定2I商基中，|p|≤1、|e|=|C_Jq/2|≤1。E的yy、yd、dd导数分别由E4|d|²/2、E3|d|和E2控制，完整二阶导数含两个mixed项。再用|d|²≤28|w|²、|d|≤6|w|，和纯selected Cauchy产品的相对导数界1/t、5/(4t²)，得到

    14E4m2 t² +(12E3m1+28E3m2)t
       +E2m0+12E2m1+(35/2)E2m2.

各14、12、28、35/2因子均已核到，式(9)时间积分无遗漏。固定微分basis与各正矩的含j估计basis不同是合法的：只是对相同正coarea积分换估计，不微分basis选择器。

**t5 correction通过。** 保留cutoff时使用绝对峰值而非可能在cutoff零点失效的相对score。固定基的l1方向范数≤sqrt3/2<1，所以作者s1、max(s2,s1²)足够；振幅参数导数≤2。三个项正好为m0[4Z2T³/3+2Z1s1T²+Z0h2T]。085的r_C2界在整个0<t≤T0统一成立，允许缩短T但不能未经新证明延长。

所有A/E/Z从085的π³-scaled数据除以27取上界，m0/m1/m2也分别正确除以27。各物理乘积的π因子没有重复或遗漏。冻结的完整28变量矩仍保留原标签，没有用7个商方向替换weighted law。

## 4. signed层、法向归一化与真实上尾

**signed层通过。** 12个axis标签和16个tetrahedral标签给D=11/2、D1=21/2、D3=753/16。对m_j的一阶导数，两种命中分别贡献D_jΛ_j/(2π³)和D_jΛ_j/π³，合为3D_jΛ_j/(2π³)。因此直接积分的L_j值和DL_j界正确，且包含重整化恢复的线性项。

纯z Hessian由076的 `L_j=F_j+T Dm_j(1,0)·z−Q_j` 支付：F_j angular界和Q_j的48D_jΛ_j³ρ/(π³T)均在账中。线性项纯z Hessian为0，但乘变化系数后并非全Hessian为0；作者用完整L_j、DL_j的直接有限时间界支付了这两个系数导数项。|c_j|≤2给2A1、2A2、2A3，全部mixed factor2也正确。

**normal层通过。** quotient coercivity取22<16sqrt2。三维球面积/径向积分给q_1(0)≤1/(π²c³)、||D²q_t||≤4t⁻⁵/(π²c⁵)，作者以π²≥9作有理化。evenness给零梯度，故D_T的值/梯度/Hessian分别≤Q2ρ²/(2T)、Q2ρ/T、Q2/T。乘a_U后的三个项和额外A2TQ0都已付；a_U Tq_1(0)不能因q_1(0)为常数而删除，作者没有删除它。

**原六维大时间尾通过。** AAᵀ=32I与max|v·ξ|≤sqrt6|ξ|给ψ6≥32|ξ|/sqrt6≥13|ξ|。S5面积π³、单位方向平方平均1/6、∫r⁷e⁻ᶜᵗʳdr=7!/(ct)⁸，结合原(2π)⁻⁶归一化给105/(8π³13⁸)t⁻⁸，积分后15/(8π³13⁸)T⁻⁷。作者的27有理化准确。它只积分t≥T>0，没有把040发散的正J2延伸到0。

## 5. 真G上界与全部T的表达式障碍

真G上界是独立于Hessian账目的正密度论证。阈值r_value=1/640满足r_value·Σ|a_i|<1/8，而球到true S2距离>1/4；故每个到达球的参数表示至少有三个独立large列。选三列的Gram行列式为正整数，Jacobian≥1，其余变量的三维quotient仍保留coercivity13。三个tail峰值与quotient峰值相乘，正union bound最多4960项，确得 `M=4960/(3⁵13³r_value⁶)`。先a.e.再由固定t连续性和邻域严格间距延拓到每点的步骤合法。

原六维peak是15/(8π³13⁶)t⁻⁶。固定T0=1/64分段积分得到

    G≤Gbar=11700693505150503878656/1172914587 <10¹³.

此为粗但真实的值上界，不是从曲率上界推出值。固定Γ_env产生的盒penalty约1.69249×10⁴¹，已超过Gbar。任何合法顶点lower endpoint不超过该顶点真G，故指定条件 `min l_v≥g0+penalty`无法满足。没有实际顶点也可证明这套测试不可用；这不否定真实G≥g0。

全部T结论的适用域是同一rho、R、h、冻结085界和0<T≤1/64。Γ_env至少为C1T+C_−7/T⁷。其全正轴最小值的八次方为 `8⁸ C1⁷ C_−7/7⁷`；乘盒系数的八次方后，精确见证给penalty下界≥10³⁵>Gbar。只用正数八次幂的单调性，不需定位极小点、第二个T或时间扫描。严格凸性公式也正确，但不是排除所需额外计算。

关键逻辑层级：这是 **chosen upper-bound expression** 的下界，不是实际最小Γ的下界。真实Γ可能远小于Γ_env。改变cutoff/盒几何、改良联合估计、使用不同证书或实际positive-window数据均不被此结论排除。不能把“粗界再大”改称“模型失败”。

## 6. 条件window比较与S状态

F=B_(1/250,1/2)是明确同点同模型的新比较对象，实际F/grad/vertices均未计算。057固定全32列right inverse Aᵀ/32满足l1方向范数≤1、平方l2范数1/32。标量log二阶上界1/(4t²)给 `D_uu p≤129p/(128t²)`，因此F上Hessian为129F/(128a²)，且|grad logF|≤1/a。这一上界推导已独立核对，不将057的lower-log界误当upper-Hessian界。

用同盒半径sqrt5 h≤9h/4、exp(s)≤1/(1−s)且s=1125/4096<1，可得window顶点penalty系数10078125/194707456乘**真实匹配中心upper endpoint**。中心/32顶点输入仍欠付，所以它不是已通过证书。

原mixture lower-log接口是Hess logF≥−I/(16a²)，Taylor损失带1/2，故q_box=5h²/(32a²)=78125/33554432。含gradient的保守interval条件以F_L作正分母下界合法；不能假设梯度零。零梯度例子明确只是尺度比较。该log损失与window的加性penalty单位不同，不能直接作为优劣的已认证margin。

条件gate `B<4/48750` 和 `I_U⁵<64⁵(g0/6)⁶`确实等价：B=I^(5/6)，64^(5/6)=32。g0仍为候选，实际g_L、I_U及S未付。局部Γ账目不改变当前全域端点，也没有使S>0。固定F不能使用full-G的空间齐次性；作者明确保留该限制。

## 7. 最小saved证书检查与执行收据

唯一新纯算术配置检查了：

1. 085三个完整导数数组、cusp系数和三个原28-label moment系数与093 extraction逐字对应。
2. 12个component名无重复，20个正系数减去保存的六项total后残差全0；保存的component值及比例与保存总额一致。没有重新生成12个预算公式或求新T。
3. 保存penalty的归一化、真G上界一致性、all-T八次幂见证恒等式及10³⁵至10³⁶的精确括号成立；没有decimal扫描。
4. 条件window与lower-log系数满足独立的无量纲恒等式 `P_F/q_box=(129/8)·F_box_multiplier`；B_gate⁶等于保存的I_gate⁵。全球端点和S状态仍为UNPAID。

所有运算使用本地Python3.15标准库Fraction。完整解析审计才建立每个公式的真实性；这个saved关系检查只防止提取、拼账和见证转写错误，不把自洽数据当作证明。

程序仅一次，1个新配置PASS；tool chunk `a96992`，guard PID20799、child PID20800，exit0。全guard墙钟0.028053458秒；guard CPU0.038576秒、child OS CPU0.021768秒，总0.060344秒。无活进程、失败或重发。guard允许165秒工作、3秒kill/reap，179秒硬alarm覆盖收据；child CPU27秒加guard CPU2秒，在CPU30及总墙钟180上限内。

计数：1新关系配置；0作者八项复跑、0新T候选、0时间扫描、0重新构造component formulas、0body矩阵/四子集枚举、0density/native/kernel/vertices/gradients、0MATLAB、0递归代理。所有代码、输入与执行来源已冻结。

## 8. 最终状态

093通过完整账目及表达式障碍审查；086新独审来源合法且相同。没有必须作者补证的缺口。该结论关闭的是“完整上Hessian预算有没有漏付”和“这套固定预算是否值得买顶点”的问题；没有关闭真正的全域g/I/S问题。

交付为本报告、有限审核登记、最小saved checker、单次guard、运行登记、完整来源绑定、结果、执行/工具收据及manifest。所有材料仅在本地独占目录；不单独归档远端，不改其他代理文件。按最新调度要求，完成本轮seal后停止，不启动下一任务。

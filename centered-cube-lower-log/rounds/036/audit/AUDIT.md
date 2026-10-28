# 第36轮mod来源独审

结论：round036_mod_theory/PROOF.md的宏观比较及加强到L≥c的Dirichlet论证通过，未发现明确数学错误；新增低频类别积分引理也正确，但结合补偿线时必须另证平滑基源已class-balanced，当前随机补偿样本有明确反例，见第7节。round036_mod_cell的来源、共同L循环DP、sampler、全pool边际RN与分半固定阈值CI主要接口通过。保存统计没有确认weak或log增益。本人只写本目录，未修改两条被审线或Git，未启动代理、重跑大MC。

本目录checks.py实际运行了252个Gaussian-rational roots小例、10个保存Fraction比较，以及保存holdout CI的独立标准库反演。checks.json记录全部结果及hash；这些有限检查不替代一般证明。

## 1. 强引理：抽象质量、L≥1与连续量词

假设每坐标/每表在每类别上质量恰1/R，保证每个产品律下类别独立均匀。取对应residue事件后每支质量1/R，外层b求和总质量1。没有将给定residue事件后的坐标误当独立。

共同L的(3)是有限正项展开。对L=q+ell≥1，一个periodic区间可分q整周期及余区间；余区间长<1，其总质量≤1。因此非零字符项/总质量≤E/(q+E)≤1/2，字符展开给epsilon=(R−1)2^-n。零频主项及最后外层b求和的1/R因子正确。cell质量假设下余区间最多两个端格，E≤min(1,ell+2/K)，所以I上的(7)也正确。

这些界逐点、对全部相应L成立，同一尺度下比较之后才取sup，没有坐标/表各自选L。非负性给max及weak系数的同一乘法比较；epsilon≥1时下方向取正部，正文已处理。非负周期L1密度的中心箱平均对每个固定x、L>0连续：在固定紧集以可积密度支配，箱边界Lebesgue零测，dominated convergence合法。有理稠密尺度恢复连续sup，closed端点亦可从内部逼近。此论证依赖密度，不应不加说明改成未平滑Dirac输入。

## 2. 加强到L≥c：常数逐式核对

N≥2，alpha=2pi min(h,R−h)/R∈[2pi/R,pi]。

- alpha≤pi/N时，2/N²乘距离和、1−cos u≥2u²/pi²及sum_(d=1)^(N−1)(N−d)d²=N²(N²−1)/12，确给1−Q²≥alpha²(N²−1)/(3pi²)。再用N²−1≥3N²/4与N≥cR/2，得到1−Q²≥N²/R²≥c²/4。
- alpha≥pi/N时，sin(alpha/2)≥sin(pi/(2N))。sin在[0,pi/4]上的弦下界给N sin(pi/(2N))≥sqrt2，故Q≤1/sqrt2≤sqrt(1−c²/4)，最后一步使用c≤1。

对c≤L<1，完整格数N≥KL−2≥Kc/2≥2。完整格沿圆周连续；R|K保证跨period时字符序列仍为同一等比序列，否则该步不能直接用。端格总质量e/K≤2/K，且N/(N+e)≥1/2。因此字符/总质量≤1−(1−sqrt(1−c²/4))/2≤1−c²/16≤exp(−c²/16)。完整格、端格及两种alpha情形的常数均闭合。L≥1的1/2帽更强，拼接为全部L≥c合法。

取c=4/K（K≥4）得到epsilon=(R−1)e^-n/K²。R≤K时n/K²−logK→infty足够使epsilon→0；这是对全部L≥4/K、相对本来源自己的prior的退化结论，不控制L<4/K，不比较任意非乘积来源，不是一般Cn上界，也不排除prior本身的收益。

实际252个R2/4、N2..64的有理平方几何和检查全部通过；主线1,389,888个浮点几何和只是附加有限探路。现有解析推导无需依赖这些扫描，也没有尺度数量费用。

## 3. 理论线脚本与证据范围

已读probe.py、rational_probe.py、short_probe.py及对应保存输出/注册。两份随机诊断的脚本/注册hash全部匹配。macro脚本按j mod R聚类，在每张表内先作坐标residue卷积，再取对应终residue；fixed表prod S与prior平均正确。Fourier_product_envelope和pointwise coarse envelope是相应字符余项的上方向。n63保存残差约1e−14而理论帽远小于此，脚本/证明正确解释为浮点底层，没有当成真实残差或反例。

short脚本R=K时沿cell轴FFT与residue轴相同；ifft在b位置取回概率，表似然softmax加权，并用R倍(prob−1/R)形成joint/prior−1，归一化顺序和相位符号正确。很小残差到1e−17层不构成区间认证。10个保存Fraction宏观比较的准确分式均满足其epsilon，checks.py重算通过；一般连续量词来自解析证明。

## 4. mod-cell真实来源、DP与提议

R整除K12，表b每cell一个正box，delta=1/(64KR)小于原子最小边距的两倍，故box留在原cell，每cell质量1/K。sum_i(j_i mod R)在均匀cell向量下均匀，实际prior=1/R，总质量1。

residue_weights将K个连续cell按K/R组、每组R类reshape后加总，恰取j mod R；table轴与residue轴没有混淆。cyclic[out,in]=(out−in)mod R实现正确循环卷积，最后对table b取终residue b的对角后求和。fixed表为prod_i sum_class z，prior为这些fixed的1/R平均。所有坐标及表调用同一个L，joint没有额外乘/除prior。

draw_proposal对joint先抽均匀cells，再算该cells的b；fixed固定b，prior才另抽均匀独立b。平滑噪声和中心箱位移逐坐标独立，同一行共享L。evaluate_queries的RN分母含全部R+2来源和全部J62尺度的平均，不是抽中latent项或max。因此W≤2；lambda-tail指标有帽min(2lambda,2(R+2)J)，代码的全阈值贡献assert与之相符。f_joint≤sum_b f_b=R prior的逐点检查也是正确的结构校验，不能反向当作gain定理。

小控制对12³个cell直接枚举，并对R2/3/4/6比较DP；Fraction在二维、三长度分别与直接有限和准确相等。这些小例及成本定标在MC前，N4096、12个seed及标签固定，未按收益扩N/R/标签。来源密度和cube积分是有限正和，但实际三角/box/DP浮点值未认证。

## 5. 分半验证、家族费用与稀疏候选

每case前2048点用ESS≥32的经验扫描选择阈值，nextafter向下处理经验值处的严格尾；后2048点只验证该冻结阈值。该调整只依赖discovery，未从holdout重新选择。每个选择阈值上的CI覆盖真实尾系数；它不是该来源所有阈值的weak上CI，也不能以上CI给整条曲线排序。

全部case模型数为3·(4+5+6+8)=69，双方向138项，alpha_each=.05/138。checks.py从保存均值/帽独立重算138个Bernoulli-KL端点，最大差0，alpha总和.05；脚本、PLAN、calibration、registration hash均匹配。跨模型共享样本不妨碍union，也不免除家族数量。全样本unrestricted/ESS64扫描仍是同数据描述，不附带另一层推断。

典型稀疏例：n15/R6的unrestricted joint系数1.42794，阈值约7.494e8，尾ESS仅13.93；ESS64筛选退回1.01595。独立half-holdout选择的是约1.0616阈值，系数1.00852、CI[.91571,1.10187]，不是对那条7.494e8候选的确认。n7/R6的unrestricted ESS约10.18、n7/R4约16.44、n3/R3约25.57，亦须保留稀疏限制。当前12个joint holdout下端全部<1，而L1本已给weak≥1；没有新增已确认weak下界。

宏观Fourier诊断使用另一组64个固定查询中心，对非零字符乘积按各表似然加权，与DP余项相符；它是描述核验，没有把细小残差当收益。全部weak探路只含有限L，且没有完整weak上围栏；目前既不能证joint最好，也不能证其真正weak更差。强引理限制的连续大尺度机制与有限短尺度探索分开，边界准确。

## 6. 残留数值层级

未做区间认证的内容包括box积分/DP、真实sampling浮点概率与RN对应、阈值分类、极小Fourier余项及CI舍入。保存摘要复算与hash支付的是一致性，不是这些输入误差。结论是解析强引理通过、实际有限模型实现方向通过，未出现可宣布的新weak/log/Cn阶。

## 7. 新增低频类别引理与跨格平滑的必要条件

已只读PROOF.md新增第6节及class_mode_probe.py/results，脚本hash与结果一致。对R|K的类别并集C_r，逐格积分的有限根单位和正确：

    ∫_(C_r)exp(2pi i m y)dy
      =[(exp(2pi i m/K)−1)/(2pi i m)] exp(2pi i m r/K)
        ·sum_(s=0)^(K/R−1)exp(2pi i mRs/K).

当非零m不为K/R的倍数时后一和为0，故1≤|m|≤M<K/R的任意零均值三角修正对每个类别积分均为0。若未修正的平滑基源在每类别已有质量1/R，修正后又另证非负，则抽象L≥1比较仍可用。它不保持逐格质量，故不能免费迁移L≥4/K强结论。11组浮点残差只是有限探路，根单位恒等式才是一般依据。

**应用范围修复：平滑前每格一个原子，不保证平滑后每格/每类别质量平衡。** 本人的compensation来源允许原子靠近格边，对称箱会跨格。低频修正的类别积分为0，只保留已有类别质量，不能把这种失衡自动修复。因而新增引理不能仅以M<K/R为由直接套到本轮全部随机补偿样本；必须同时检查平滑基源class balance，或先构造/修复它。

已实际运行自有class_balance_checks.py，以保存U的精确dyadic数和精确delta作Fraction边界转移；完整结果在class_balance_checks.json，不抽新随机数、不改原来源。对R3（整除全部三个K），9组中8组的mu_delta类别质量非1/3，只有K255/v128没有箱跨格。例：

- K63/v8有两箱跨右边界，mu类别质量最大偏差准确为33177325500527/4433230883192832，约.00748378019884；乘p=.5后，projected/B校准来源的相应类别质量仍有约.00374189009942偏差。
- K1023/v128仍有一箱跨右边界，mu类别最大偏差准确为1481767654713/5998935441145856，约.00024700510103。

两例都满足M<K/R；所以这是应用条件缺失，不能归咎于低频门槛。H重标的常数背景和峰质量同时变动，但类别偏差仅乘正的ps，亦不会自动归零。来源正性和单位总质量证书不受影响；被否定的是将该来源直接移植进class-balanced模分类定理的步骤。

MC mod_cell的箱中心离格边至少1/(2KR)，delta/2=1/(128KR)，严格留格，故每格及类别质量确实准确；此问题不影响第4/5节审查。收尾已重新读取父线第6节修订：它明确要求未加三角修正的平滑基源已类别平衡，并写明本轮随机补偿样本跨格反例，不能直接套用。修复后的相容性表述通过；本审计保留问题及修复边界，不修改任何被冻结脚本或数据。

## 8. 失衡来源的正确归一化与实际先验版本

已只读父线新增PROOF第7节、imbalance_probe.py/results。一般解析扩展通过，尤其两种prior使用相同乘法因子是正确的，不需额外再乘一次Z误差。

每个产品f_b质量1，设其整类别非零字符幅度统一≤d<1。整周期字符展开给|R m_b−1|≤zeta=(R−1)d^n，故|Z−1|≤zeta，没有多余R因子。正文明确要求zeta<1，保证m_b>0、Z>0；F*=F0/Z才是单位质量来源，不能继续将F0当质量1。

任意L=q+ell≥1的非零字符/总坐标质量≤(qd+E)/(q+E)≤(1+d)/2，使用q≥1、0≤E≤1及d<1。于是F0的cube平均与uniform prior的误差≤epsilon=(R−1)((1+d)/2)^n。除以Z得到正文(12)的[(1−epsilon)_+/(1+zeta),(1+epsilon)/(1−zeta)]。

actual-prior比较中，w_b=m_b/Z，令D_actual=sum_b m_b A_L f_b。该未归一化分母在(1±zeta)A_L bar f之间，而A_LF*/A_L bar f_actual=(A_LF0)/D_actual，Z完全消去。因此同一乘法因子准确成立。实际权重的比较混合使用unconditional f_b；它不是声称F*给定table后仍以unconditional f_b为分量，后者实际还有限制E_b。连续sup和weak系数继承该逐点比较，与第1节已审量词一致。

等测度类别及rho_i,b≥kappa时，减去常数kappa后的非负剩余质量为1−kappa，均匀项非零字符消失，所以d≤1−kappa。kappa独立n、logR=o(n)时两误差趋0；kappa_n减小时可用明确充分条件n kappa_n/2−logR→infty（同时控制zeta、epsilon），不能只用有限数值下端替代整个匹配参数族的资格。有限36个已认证来源及任意平移的单位质量/正下端确实保留，可取.01，但不证明未认证的任意K_n、v_n来源都有同一下端。

自有imbalance_checks.py只读保存的K63/v8精确类别分式，实际以Fraction重算R3、n1/4/7/15/31/63的全部m_b、Z及actual weights；全部界通过，源U hash相同，保存Z浮点显示与准确分式一致。n1时Z−1约−.00979881225，必须归一化；n15/63时这个特定循环表的Z恰为1，不应把它当一般恒等式。该脚本没有重新抽样或计算高维尾。

本扩展修复了「跨格仍假装类别平衡」的问题，并只恢复L≥1宏观比较。没有恢复逐格质量、4/K范围或全尺度weak结论；仍不提供新的Cn/log阶。未发现该扩展的明确数学错误。

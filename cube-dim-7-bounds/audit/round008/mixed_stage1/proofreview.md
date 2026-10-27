# R008 mixed StageI 独立 Sol 审计：接受失败门槛，停止StageII

2026-10-27。新独立clone起点 `22208695024e443ac5307e242d5a11a29593bd60`，仅写 `audit/round008/mixed_stage1/`。理论源 `upper/round007/theory/mixed_moment_design.md` SHA256=`ad6d9275b0facc4db882c98142e7214fab2139dd132022a59488f9ce90a68af1`。全文读该设计，以及Terra StageI registration、run、assemble、check、certificate、report、receipt、sources、manifest和预运行hash。本审仅接受StageI及其所需理论；未计算或授权StageII。

**结论：mixed Gamma/Mellin、cross-kernel误差和正ℓ7 Minkowski下界推导成立；StageI的正增益gate严格失败。** 独立准确signed账本和768-bit标量重组确认

`L_H < 0.901274565714358557 < 1`。

这是固定下界surrogate本身的严格上包，因而该surrogate不能证明μ>0；不只是数值区间尚未分离。真实H≥1仍成立，LH<1暴露的是下界路线松弛。接受“StageI完成、停止StageII”，没有新I预算、g或Cn结果。

certificate里仍保存 `optimistic_surrogate_M_lower` 形式字段。由于LHupper<1，该算式代入负候选μ，不属于Radau非负均值可行域，**不得作有效乐观lower，更不是真实M lower**。Terra最终report/checker已明确此限制。本审只核它的算术值，拒绝其字段名的lower语义；positivity gate本身足以决定停止。

## 1. 原模型、pair matrix与Cauchy–Binet系数

A、B为两个不交8方向Hadamard块的raw covariance，TA=detA^(1/7)、TB=detB^(1/7)，L=TA+TB，R=det(A+B)^(1/7)，H=R/L≥1。s1/8倾斜dν=L^-s dP/q，q=EL^-s来自完整R004证书。于是 `E det(A+B)^(-1/56)=q Eν H^-1/8`，指数及时间尺度正确。

四个独立matrix pairs经determinant Minkowski和四项AM–GM给 `M≤(1/2)Π J_pairmat,type`。若所有pair type有同一上界J*，才可写J*^4/2；作者没有假设不同相对frame的pair matrix law必相同。每个8-column frame的Gram8I及独立同权raw variables保持原body64共同Cauchy时间，不是新生成元。

对7-subsets J，Cauchy–Binet系数cJ=det(VJ)^2≥0。按A中选j列分组，有

`Cj=Σ_(|J∩A|=j)cJ=8^7 binom(7,j)`，

因为这些正好是det(t VA VA^T+VB VB^T)=8^7(t+1)^7的系数。它适用于任意两tight blocks，无需枚举11440 minors；零minor无贡献。此处是binom7j，不是binom8j，也不能将power coefficient sum与sum的power混淆。

## 2. F_j正weighted Mellin公式和解析域

对同一块选j个S_i，`F_j(z;e)=E[(Πselected S_i^e)T^-z]`。块determinant identity对八个S完全对称，所以只依赖选取数j。写Y_i=Rp_i，R~Gamma4、p~Dirichlet(1/2)^8独立；T=16^(1/7)R^-1(Πp)^-1/7。权ΠS_i^e带4^-je R^-je Πselectedp^-e。直接积分得到

`F_j=2^(-4z/7)4^-je Γ(4+z−je)/Γ(4+8z/7−je) · Γ(1/2+z/7−e)^j Γ(1/2+z/7)^(8−j)/Γ(1/2)^8`。

radial power是z−je，simplex选中参数减e，所有4和7因子均正确。对e≥0、j0..7，Re z>7e−7/2保证选中及未选中simplex参数正；radial参数>1/2+(7−j)e>0，denominator参数也是正simplex参数总和。故这是一条足够common domain。j0时可放宽，但作者使用较强统一域安全。

在该域正weighted积分绝对有限；局部用附近两个实指数的有限矩控制差商/对数因子，给holomorphy。由positive representation有 `|F_j(v+it;e)|≤F_j(v;e)`；Hölder给log F_j(v;e)凸，即使未加T幂的原权可能不单独可积，也不妨碍在有限矩域内应用Hölder。Gamma分母按reciprocal Gamma使用，无隐含除零分支。

StageI e1/7、b9/8、a9/16；rho27/64使a−rho=9/64>max(0,−5/2)，a+rho=63/64。全部j0..7和三个实envelope节点的Gamma参数均正，本审准确有理检查通过。没有误用StageII较窄mixed-pole strip。

## 3. Signed Mellin/Fourier公式，不是modulus-square

对b>0、a=b/2，beta/Fourier identity把(A_scalar+B_scalar)^-b分为两侧T^-a和Γ(a±it)因子。对分别带j和7−j个S^e权的独立块期望，得到

`I_j(e,b)=[πΓb]^-1∫_0∞ Re[Γ(a+it)Γ(a−it)F_j(a+it;e)F_(7−j)(a−it;e)]dt`。

全线原系数是1/(2πΓb)，因为K_j(−t)=conjugate K_j(t)，半线为real part并乘2。F_j和F_(7−j)一般不同，因此不能写|ΓF_j|²；虚部不应被强制为0，real Taylor terms也不可取绝对值。

绝对Fubini成立：先对随机权/T幂取绝对期望得F_j(a)F_(7−j)(a)，Gamma product在compact t上有限，尾由recurrence可积。原I_j是正weighted expectation，故最终完整积分下端才可与0相交，不能clip signed panel。交换两块给I_j=I_(7−j)，这只把八个j压到四个，并未改变cross-kernel。

Terra run正确构造G_j(z)=Γ(z)F_j(z)，乘G_(7−j)的coefficientwise conjugate series，而非自身conjugate。在实panel变量上，这正是holomorphic partner Γ(a−it)F_(7−j)(a−it)的Taylor coefficients。pilot额外核direct constant、负频率conjugacy和点Gamma recurrence；没有imaginary-zero断言。所有65个integrated even real terms及signed panel result保留。

本审准确读取发现3115个term intervals严格负、256个constant imaginary intervals均排除0。所有256 panel积分区间本次恰正，但这不证明一般kernel非负；signed算法和±tail仍是必要条件。每个panel区间都包含其65项endpoint intervals准确和，本审未重算任何series coefficient。

## 4. Crossed Cstar、Cauchy remainder及双侧无限尾

H_j(v)=Γ(v)F_j(v;e)正且log-convex。对complex frequencyw=t+iy，Gamma/moment bound给 `|K_j(w)|≤H_j(a−y)H_(7−j)(a+y)`。该右端的log是y的凸函数，故在|y|≤rho的最大值在两个端点之一：

`Cstar=max{H_j(a−rho)H_(7−j)(a+rho), H_j(a+rho)H_(7−j)(a−rho)}`。

两个crossed products必须都取，不能沿用R004对称单乘积；无需Gamma单调性。Terra先取两product的向外upper再取max，正确。

64 panels覆盖[0,32]，半宽1/4、centers(2k+1)/4，eta=h/rho=16/27。degree128遗漏129，raw aggregate remainder是

`32 Cstar eta^129/(1−eta)`。

factor32是总区间长度，不是面板数。Cauchy圆严格位于Gamma/mixed moment可积域。给定recurrence N16，完整normalized absolute tail是

`F_j(a)F_(7−j)(a)Γ(a+16)^2/[πΓ(b)31·32^31]`。

它支付t≥32全部未积分部分，符号未知。正确lower为 `(poly_lower−raw_error)/(πΓb)−tail_abs`，upper为 `(poly_upper+raw_error)/(πΓb)+tail_abs`。Terra两侧均已支付，不会把旧positive-tail规则误用到这里。四个完整I_j区间下端均正，实际未触发zero clipping。

本审768bit用log-Gamma独立重评F_j及两cross，证明declared Cstar dyadic≥两个独立cross上端；再按同一declared安全Cstar重组raw error，独立重组signed tail和四I_j，均被certificate端点包含。

## 5. Positiveℓ7 Minkowski下界与pair-type适用性

对每个7-subset J，xJ=cJ^(1/7)Πselected S_i^(1/7)≥0，R=||x||ℓ7。以正权L^-9/8的integral triangle inequality给

`E[R L^-9/8] ≥ ||E[x L^-9/8]||ℓ7 = [Σ_J cJ I_j(1/7,9/8)^7]^(1/7)`。

这是正确下界方向。有限各I_j来自上面的绝对weighted Mellin；有限多个coordinates的ℓ1 bound也给所需norm integrability，未循环假设下界对象有限。固定系数cJ及完整exchangeable expectation保证没有随机data-dependent weights被移入期望。

代Cj和j↔7−j后，

`EνH ≥ LH =8[2(I0^7+7I1^7+21I2^7+35I3^7)]^(1/7)/q`。

四个I_j取complete positive lower、q取upper得安全LH lower；反向端点给的是这个**明确surrogate表达式LH**的upper，不是真实EνH的upper。这个区别决定下节gate解释。由于F_j只依赖j、Cj只依赖tightness，该surrogate统一适用于不同pair orientation，不需要矩阵pair law相同。

## 6. Gate严格失败与31/500的正确量词

本审从declared exact I_j端点及完整q端点，分开向外算LH下/上端，被certificate的两个准确dyadic值包含；并准确验 `LH_upper<1`、`mu_lower=LH_lower−1<0`。因此固定LH真值本身也<1，重新提高同一积分精度不会把它变成>1。这条L7 lower surrogate无法通过μ_L>0 gate，StageII不应开始。

真实H≥1和真实μ=Eν(H−1)≥0仍成立；它们只给本证书的零修正baseline，不能从LH<1推出实际matrix gain不存在。

设计中的optimistic极限推导需μ_candidate≥0：在ν最小可行值μ²下，Radau给M_budget=M_pair/sqrt(1+μ_candidate)。若LHupper≥1，用LHupper限制此特定surrogate可认证的正μ大小，才可给certificate可达性lower。当前LHupper<1，直接代M_pair/sqrtLHupper属于负μ候选，约.08536；它没有非负Radau可行域内lower含义。certificate字段名仍沿旧名称保存，Terra最终report和checker已明确“formal only”；本审拒绝其有效bound语义。

按已有H≥1把可认证μ lower归0，只能获得δlower0、J_pairmat≤q和原pair预算，不能给该路线的新严格提升。原pair确定预算表达式的已认证lower仍大于目标31/500=.062；这描述本证书在该gate失败后的能力限制，不是actual M lower，更不是Cn反例、model exclusion或证明目标不可由其他方法达到。positivity失败已足够停止，不需要用无效形式gate支持结论。

## 7. 实际本审运行与文件责任

minchecker.py为stdout-only可携带入口，支持--numeric-dir、--design、--r004-dir；默认合并后项目相对目录。没有写结果逻辑，shell仅将stdout保存为本审results.json。一次执行exit0，0 assertion failure、0 repair，0.421287792007206秒；Python3.15.0、python-flint0.9.0、Arb768bit。共享基础库，非第二算术内核。

准确核256唯一终态、65terms/panel、term signs和imaginary data，仅读ledger并汇总，不执行run.py/check.py，不新增节点或MC。R004 certificate逐字节与原件相同，复制/原manifest的certificate项一致，reused_q对象一致。所有理论和Terra输入在本次前后hash相同，prehash中的run/assemble/理论/旧证书也匹配final bytes；Terra manifest逐项通过。全源hash保存在results及sources。

数值面板系数仍由Terra generator及certified acb_series算术负责；本审检查其数学规格、signed实现和标量接口，不称独立第二mesh。Terra最终报告已清楚区分StageI失败与物理/算子结论，主集成人保留最后采用和证明责任。

## 8. 结论范围与停止

接受F_j与StageI signed MB、cross envelope、ℓ7 surrogate数学，以及“positivity gate严格失败、StageII停止”的证书。没有StageII K7或ν/δ计算，没有新完整I上预算，更没有g/Cn结果。R001/R002的complete boundary和一般L1/all-common-r接口没有被本轮改写或重新认证。

此独立审计到此停止，不自动进入II、不修改surrogate、不追加panel、指数、优化或新策略。旧封存及其他工作者文件未动；无MATLAB、子代理、commit或push。

# Positive gated polynomial 的McDiarmid upper接口

本稿先完成本机小例计算，再给有限n自含证明。没有Gaussian、saddlepoint或外部未核验渐近。当前固定选择 `r0=10^−6`、`C=81/20=4.05`；C来自确定数学预算，不按尾收益调参。原40轮source/pilot字节及CI均未修改。

## 1. 条件query law及非空正项

先固定ordered query-bin标签，各隐藏坐标X_i独立；本来源对坐标置换对称，所以给定counts可用规范ordered标签代表，不声称histogram条件下原坐标本身独立。

固定表b、共同cube长度L，记两类正核为z_{i,b,0}(X_i)、z_{i,b,1}(X_i)。在macro L≥1时，完整周期保证每类核至少其source class质量除L，因此0<η<1时两类均严格正。box宽δ>0又保证有限上界。

对固定的非空count集合E_b⊂{0,…,n}，表内gated项是

`F_b(X)=Σ_{s∈E_b}[t^s] ∏_i(z_{i,b,0}(X_i)+t z_{i,b,1}(X_i))`。

它是有限正multilinear sum，且在全部隐藏位置严格正。固定其余坐标后，可写成 `a_0 z_{i,b,0}+a_1 z_{i,b,1}`，a_s≥0且a_0+a_1>0。因此对任意x_i,x'_i，

`F_b(x_i)/F_b(x'_i) ≤ max_s z_{i,b,s}(x_i)/z_{i,b,s}(x'_i)`。

空E_b给F_b≡0，保留标签但从log/concentration项中去掉；不对0取log。fixed产品是E_b全支持的情形。prior项为 `π_b∏_i(z_0+z_1)`，π_b>0这个常数不改变log振幅；π_b=0则同样为零项。没有把gated项额外乘π。

## 2. 旧quarter-cap给坐标log振幅

设每类核在第i个bin内的振幅不超过v_{i,b,s}，并有正均值下包l_{i,b,s}。使用40轮旧quarter-cap定义

`r_{i,b}=max_s v_{i,b,s}²/(4 l_{i,b,s}²)`，`q_{i,b}=2sqrt(r_{i,b})<1`。

令m为该核真实均值。z_min≥m−v≥l−v，且v/l≤q_i，所以

`z_max/z_min ≤1+v/(l−v)≤1/(1−q_i)`。

第一节于是给H_b=log F_b的坐标bounded-difference常数

`c_{i,b}=−log(1−q_{i,b})`。

假设所有r_i≤r0<1/4，q0=2sqrt(r0)。积分1/(1−q)给

`c_i≤q_i/(1−q0)`。

再用 `log(1+r)≥r/(1+r)≥r/(1+r0)`，有

`S_b^osc:=Σ_i c_i² ≤ C0 Σ_i log(1+r_i)`，

`C0=4(1+r0)/(1−2sqrt(r0))²`。

令 `Rold,b=∏_i(1+r_i)−1`，则可用确定上帽

`S_b^osc ≤ S_b:=C log(1+Rold,b)`，只要C≥C0。         (1)

count表示中log(1+Rold,b)=Σ_q N_q log(1+r_qb)，正是旧source stats的第4列之和，不需要n² PB或逐坐标查询。

这里r必须是旧range²/4 quarter-cap。若以后把variance帽换成range²/12，不能仍把q写成2sqrt(r_new)；应保留旧quarter r用于(1)，或重新调整系数。方差改善与log振幅预算是两个接口。

## 3. 当前来源r0的确定预算

40轮source全部R、表、class的振幅≤B/L，其中

`B=29879/44945280`。

共同source类质量≥.49，故真实class均值m≥.49/L。假设real-arithmetic kernel近似误差预算e≤10^−8，代码均值下端是tilde m−e。注意

`tilde m−e ≥ m−2e`，

因此L≤2时，其L倍至少 `.49−4e`，不能遗漏第二个e。于是旧quarter-cap满足

`r_i ≤ B²/[4(.49−4·10^−8)²] = 4.6016316918…×10^−7 <10^−6`。 (2)

Fraction核验保存完整有理比较。source解析中点误差即使只用L≥1的粗界也≤11/1283637196800，另列的经验guard10^−13相加远小于10^−8。但这不认证实际全部FP误差；(2)对理想核计算/误差预算成立，数值执行仍携带原FP资格。

r0=10^−6时q0=1/500，

`C0=1000001/249001 <81/20`，

最后的整数比较为20000020<20169081。因此固定C=81/20可直接用于当前旧source stats，覆盖全部R=4/8/16/32，不依数据/收益选择。

## 4. 自含Doob–Hoeffding推导

先证明有界中心变量D的Hoeffding引理。若D落在区间[a,b]且ED=0，令ψ(t)=log E exp(tD)。指数倾斜后的D仍在[a,b]，所以

`ψ''(t)=Var_t(D)≤(b−a)²/4`。

ψ(0)=ψ'(0)=0，积分两次给 `E exp(tD)≤exp[t²(b−a)²/8]`，任意实t成立。variance帽可由与区间中点之差平方≤(b−a)²/4直接得到，无分布近似。

对H=log F_b，取Doob martingale `M_i=E(H|X_1,…,X_i)`。独立性意味着固定prefix时，作为X_i函数的未来积分仍有振幅≤c_i。故D_i=M_i−M_{i−1}的conditional均值0、conditional range宽≤c_i。逐步应用上述conditional mgf帽并迭代：

`E exp[t(H−EH)]≤exp[t² Σ_i c_i²/8]≤exp(t² S_b/8)`。

对t>0作Markov并优化，任意a>0得到

`P(H−EH≥a)≤exp(−2a²/S_b)`，S_b>0。             (3)

若S_b=0，则所有c_i=0，H在product law下恒定，不使用(3)的零分母。非空正项的log在紧bin product上有有限上下界，mgf和EH均合法；空项已单独处理。

## 5. 真实均值上包与单项upper

若F_b非负、真实均值≤u_b，则u_b=0意味着F_b=0 a.s.，可以去掉。对正项、u_b>0，Jensen给

`E log F_b ≤log E F_b≤log u_b`。

所以λ>u_b时，(3)给

`P(F_b>λ)≤exp[−2 log(λ/u_b)²/S_b]`。          (4)

边界规则：

- 已知零项：λ≥0的strict tail为0。
- S_b=0：F_b恒定且≤u_b，因此λ≥u_b时strict tail为0；λ<u_b时一般只能用1。
- S_b>0、λ≤u_b：用1；不把负log gap平方后误作小upper。

对λ>u_b可把(4)再与Markov min(1,u_b/λ)取小。strict λ=u_b、S_b>0时1保守有效，不声称概率0。

joint表项可用 `u_b=A_b^+ PB_b^+(E_b)`，空段的PB严格设0；prior真实项用实际固定π_b乘产品，均值上包可用 `u_b=π_b^+ A_b^+`。π区间上端不需要归一化，Σπ_b^+可以大于1，仍是逐项合法均值上包。actual π_b正时其log只是常数，坐标振幅不变；actual π_b0为零项。fixed产品直接u=A^+。这说明posterior概率与prior概率的误差包只进入u，不进入新的随机系数。

## 6. Sum与共同尺度

对总和F=Σ_b F_b，令U=Σ_b u_b。U=0时F=0 a.s.。若λ>U>0，对所有u_b>0分配

`t_b=λ u_b/U >u_b`。

只要F>λ，至少一个F_b>t_b，否则各项求和≤Σt_b=λ。无需表间独立。由union和(4)，令g=log(λ/U)>0，

`P(F>λ) ≤min(1,Σ_{b:S_b>0} exp(−2g²/S_b))`。  (5)

S_b0的正项在该阈值下事件为空，零项不参与。令Smax=max_b S_b；Smax>0时更简洁地

`P(F>λ)≤min(1,R exp(−2g²/Smax))`。             (6)

Smax0且λ≥U时tail为0；λ<U时一般用1。若实现统一对λ≤U返回1也保守，只在Smax0、λ=U处舍弃这个更紧的零界。空段/零prior权重可继续计入R使(6)更松，但不可计算log0。

每个L先用同一个L形成全部表和F，再用(5)/(6)。最后对有限共同尺度集取各L upper之和、截到1，是合法union bound；不允许每坐标/每表单独最大化，也不假设尺度独立。连续窗口max不由此自动受控。

## 7. 流式接口及旧贡献帽

当前source stats中每b,L有old log product variance cost，故S_b=C·stats[3,b,L]，只作O(NRJ)代数。joint/prior各自的 `log U=logsumexp(log u_b)` 不同；Smax可以共用同一组表内旧r代价。fixed0用Sfixed而不是R表union；coarse形式中的R取1。数值应在log域计算g及exp尾，单独处理−inf均值和S0。

McD upper本身不保证≤U/λ，因此**不能独立替换旧upper后继续套原Markov贡献帽**。为保留原 `min(2λ,2RJ C_Z C_E)`，可以逐L取

`min{旧Cantelli/Markov upper，variance/3 Cantelli upper，McD upper}`，

再作L union，或以原aggregate upper作global min保底。该新upper≤旧upper，原importance权重/cap随之保留。没有对新upper单独使用lower active条件。

本接口只证明上包及数据需求。任何同seed恢复source充分统计的执行需独立核对原count/stream hash，并作为旧数据重放后处理登记；本稿没有启动高维MC、没有生成新CI或宣布新的weak/C_n结论。

# 第44轮独审：表内镜像标签、精确质量及消失的先验谱价

本证明在本机执行真实周期窗的Fraction小例与数值倾斜积分后整理；实际结果见 `CHECKS.json`。先执行了340个镜像核例、5个维数的质量例和100个均匀求积倾斜例；随后提出下述一般证明及二类别增强帽。之后另补100个Fraction成对查询律的精确倾斜核验。没有新增C_n数值估计。

## 1. 来源与table-specific labels

θ=(0,1/8,1/4,3/8)，d=(1/8,3/16,1/4,5/16)，δ=1/20，eps=2/25。ρ_b为正周期密度

\[
\rho_b(t)=\varepsilon+\frac{1-\varepsilon}{2\delta}
\sum_{\sigma=\pm1}1_{[\theta_b+\sigma d_b-\delta/2,\theta_b+\sigma d_b+\delta/2]}(t\bmod1).
\]

每个峰按周期延拓解释。令C_{b,+}=[θ_b,θ_b+1/2) mod1，C_{b,−}为补集。反射R_b(t)=2θ_b−t mod1保持ρ_b、交换这两个class（边界忽略），故每个class的ρ_b质量恰为1/2。

这不是把同一个class划分在各table间共享。实际Fraction窗检查中235/340个kernel例的own-label与common-label积分不同，n=1,2,3的整体输入平均也有精确有理差异。当前特定四个table在固定common class [0,1/2)上恰也各有质量1/2，但这不使其后验或输入相同，也不是任意移动峰的普遍性质。

## 2. 质量1不需要逐点winner

设w_b(k)≥0，对每个k=0,...,n满足Σ_bw_b(k)=1，且全部查询/尺度使用同一组w。令N_b(y)=Σ_i1_{C_{b,+}}(y_i)，定义

\[
f_w(y)=\sum_b\rho_b^{\otimes n}(y)w_b(N_b(y)).
\]

在table b自己的乘积来源下，N_b∼Bin(n,1/2)。因此

\[
\int_{[0,1)^n}f_w
=\sum_b\sum_{k=0}^n\binom nk2^{-n}w_b(k)
=1.
\]

定义π_b=Σ_k Bin(n,1/2){k} w_b(k)，则π_b为每个table项的真实质量，Σπ_b=1；matched product prior Σ_bπ_bρ_b^{⊗n}也有周期质量1。

注意N_b(y)随b改变，故即使w是每列one-hot，也不保证Σ_bw_b(N_b(y))=1。精确例：n=2，w_b(k)=1_{b=k}（b=0,1,2,3；k=0,1,2），y=(11/20,4/5)，有N_0=0、N_1=1，两个table同时开启，权重之和为2。这是合法的非负多项叠加输入；不能继续称它为逐点唯一winner。质量仍由上式精确等于1。

本归一化依赖每个table自己的class质量均为1/2。任意更改labels、source或权重后，应重算Σ_bE_{Bin(n,p_b)}w_b；不能因为每列权重和1就忽略不同p_b。这里的结构性反射使这个重算闭合。

## 3. 真中心窗的反射恒等式

对任意全长度L>0，令

\[
z_{b,\pm,L}(x)=L^{-1}\int_{x-L/2}^{x+L/2}\rho_b(y)1_{C_{b,\pm}}(y)dy,
\quad a_{b,L}=z_++z_-,\quad p_{b,\pm,L}=z_\pm/a.
\]

eps>0使a≥eps，因此这些量处处定义。换元y↦R_b(y)保持周期Lebesgue测度，反射中心窗并交换class，故精确有

\[
z_\pm(R_bx)=z_\mp(x),\qquad
a(R_bx)=a(x),\qquad p_\pm(R_bx)=p_\mp(x).
\]

任意实t的score倾斜律

\[
Q_{b,L,t}(dx)=a_{b,L}(x)^{1+t}dx/m_{b,L}(t),
\quad m_{b,L}(t)=\int_0^1a^{1+t}dx
\]

保持R_b不变。ρ_b有正上下界，所以m有限且为正，包括负t和短窗。于是

\[
E_{Q_{b,L,t}}p_+=E_{Q_{b,L,t}}p_-=1/2.
\]

这里不声称任意查询律ν都保持先验；需要ν对R_b不变。score倾斜满足这一条件，而只保持class质量、破坏source反射对称性的control一般不满足。

## 4. 原先验γ价消失；长窗二类别帽可加强为9/8

对固定table、L、t，其倾斜class先验μ_t=(1/2,1/2)恰等于来源η。故第43轮的χ²(μ_t||η)=0、原先验Gram leading eigenvalue γ_t=1。使用原Bin(n,1/2)分母与正确新先验分母一致，没有原先验改变引起的γ_t^n费用。

若L=q+ell且整数q≥1，0≤ell<1，canonical残余窗表示给

\[
p_+=\frac{q/2+r_+}{q+r_++r_-},\qquad 0\le r_+,r_-\le1/2.
\]

该函数对r_+递增、对r_-递减，因此

\[
\frac q{2q+1}\le p_+\le\frac{q+1}{2q+1}.
\]

反射给Ep_+=1/2。二类别Gram非平凡特征值为

\[
\tau=4\operatorname{Var}(p_+)\le\frac1{(2q+1)^2}=:c_q.
\]

这也直接来自p_+离1/2至多1/[2(2q+1)]。二类别histogram trace的精确对称张量公式为

\[
\mathcal T_n=\sum_{j=0}^n\tau^j
\le\frac1{1-c_q}.
\]

所以q1可用严格统一帽9/8，比直接接第43轮的2更强。没有声称该具体正背景box族能取等。整数L时class积分恰L/2，p_+=1/2，实际τ=0，脚本对整数L设定kernel=.5正确。

q=0短窗不享受上述小于1的帽：仅有0≤p_+≤1、τ≤1，故通用hist trace≤n+1。源虽有正背景，但窗口完全落在某个class内部时另一class积分仍为0。随L→0，τ可能趋近1；不能把q≥1的9/8或第43轮2移用于所有短L。主数值mirror短L出现τ≈.998，与本结论相容。

## 5. 表内反射不能自动转成整个多表输入的无增益结论

对单项F_{b,w}(y)=ρ_b^{⊗n}(y)w_b(N_b(y))，逐坐标同时使用R_b可得

    F_{b,w}(R_b^n y)=F_{b,w^rev}(y),  w_b^rev(k)=w_b(n−k).

其中心窗平均也有对应query反射恒等式。这证明单项的query分布在权重反序后相同。各table的R_b不同，不能把不同坐标变换分别套到一个共同query下的Σ_bF_b及其共同尺度max；这样会改变不同table项的联合关系。

本组特殊参数确有一个有限的全局闭合关系，不能笼统说不存在任何whole-f反射。全局S(y)=−y使ρ_0、ρ_2保持，交换ρ_1与ρ_3；labels中C0交换正负，C2保持，C1与C3互换。因此f_w(S^n y)=f_{w'}(y)，其中

    w'_0(k)=w_0(n−k), w'_1(k)=w_3(k),
    w'_2(k)=w_2(k),    w'_3(k)=w_1(k).

但Σ_bw'_b(k)=1+w_0(n−k)−w_0(k)，一般不再满足原逐列simplex约束。反射输入本身仍质量1；该等式只说明原受限参数化不对这次全局反射封闭（除非额外满足w_0反序对称）。即使有封闭子族，也只是生成等价输入，不推出最优门控等于product prior。

因此本轮移除了每个table在score倾斜下的原先验漂移谱价，仍没有支付m_{b,L}(t)^n、不同table倾斜查询律、权重依赖、共同尺度max或真实weak积分。固定channel的trace帽不推出C_n没有log，也不推出新的C_n下界。

## 6. 数值与归档范围

实际执行340个Fraction真周期窗反射例，100个Fraction对称有限query law倾斜例，5个维数质量核验，100个16384节点的均匀求积倾斜例，及336个主脚本kernel的独立Fraction交叉检查；最大主核差约3.34e−15。主脚本扫描与高维训练的资格见 `REVIEW.md`，没有把数值误差小或guard当向外误差认证。

## 7. 一个共同长尺度下，保留实际score的多表联合界

这项结论比前面的mirror接口更一般：不要求source反射，也不使用t≠0倾斜。每个table b只需两个class子密度ρ_{b,+},ρ_{b,−}，各周期质量1/2。对于几何labels，子密度就是ρ_b乘各自class指示；也可用一般非负子密度及潜在class标签的等价展开。不同table的class划分可以不同。

固定同一个全长度L=q+ell，q≥1。令a_b(x)为该table的真实中心窗平均，Q_b(dx)=a_b(x)dx，p_{b,+}(x)为该窗的class posterior。周期Fubini给∫a_b=1及E_{Q_b}p_{b,+}=1/2；这个均值恒等式不要求mirror。上一节的范围估计仍成立：

    q/(2q+1) ≤ p_{b,+} ≤ (q+1)/(2q+1),

所以单坐标原先验Gram的谱为1、λ_b，且0≤λ_b≤r²=(2q+1)^−2。

令β_n(k)=binom(n,k)2^−n，w_b(k)≥0且对每个抽象count k有Σ_bw_b(k)=1，π_b=E_{β_n}w_b。将π_b=0的table删除：因β_n(k)>0，这时w_b(k)对全部k均为零，不存在除零问题。定义

\[
A_b(x_1,\ldots,x_n)=\prod_i a_b(x_i),\quad
T_b=E[w_b(N_b)\mid X^n],
\]
\[
F=A_Lf_w=\sum_bA_bT_b,\qquad
G=A_Lg=\sum_b\pi_bA_b,\qquad
g=\sum_b\pi_b\rho_b^{\otimes n}.
\]

F、G均非负周期密度、积分为1。q≥1还使每个a_b≥q/L>0，故G>0。逐点加权Cauchy给

\[
\frac{F^2}{G}
=\frac{(\sum_b\sqrt{\pi_bA_b}\,\sqrt{A_b}T_b/\sqrt{\pi_b})^2}{\sum_b\pi_bA_b}
\le\sum_b\frac{A_bT_b^2}{\pi_b}.
\]

这里已经保留真实source score A_b；没有把它删除或换成代理。对共同Lebesgue查询积分后，每项分别是

\[
\int A_bT_b^2=E_{Q_b^n}T_b^2.
\]

不同Q_b只在Cauchy分离后逐项使用，不需要在原多表和内对每项另作不同的query换元。

### 共同源端count基

令Z_i为独立公平±1标签，N为plus count。对d=0,...,n定义

\[
\psi_d(N)=\binom nd^{-1/2}\sum_{|S|=d}\prod_{i\in S}Z_i.
\]

这些量只依赖N。不同子集的Rademacher乘积正交，故E_{β_n}ψ_dψ_e=1_{d=e}。n+1个这样的函数构成count空间的正交基；这就是这里使用的normalized Krawtchouk基。特别ψ_0=1。它在各table间共用，因为各source端count先验都是同一个β_n；不需要各table的query律相同。

展开w_b=Σ_d ξ_{bd}ψ_d，其中ξ_{bd}=E_{β_n}w_bψ_d。binary单坐标Gram在共同源端基(1,Z)下为diag(1,λ_b)，独立n坐标的tensor Gram使次数d的对称分量对应λ_b^d。因此精确有

\[
E_{Q_b^n}T_b^2=\sum_{d=0}^n\lambda_b^d\xi_{bd}^2.
\]

对每个d再用一次加权Cauchy：

\[
\xi_{bd}^2=(E_{\beta_n}w_b\psi_d)^2
\le\pi_b E_{\beta_n}(w_b\psi_d^2),
\quad
\sum_b\frac{\xi_{bd}^2}{\pi_b}
\le E_{\beta_n}\left[(\sum_bw_b)\psi_d^2\right]=1.
\]

degree0实际取等，因为ξ_{b0}=π_b且Σπ_b=1。将这些非负项组合，得到完整的有限n联合界

\[
\boxed{\int_{[0,1)^n}\frac{(A_Lf_w)^2}{A_Lg}
\le\sum_{d=0}^nr^{2d}
\le\frac1{1-r^2}.}
\]

由于∫F=∫G=1，进一步

\[
\boxed{\int_{[0,1)^n}\frac{(F-G)^2}{G}
=\int\frac{F^2}{G}-1
\le\sum_{d=1}^nr^{2d}
\le\frac{r^2}{1-r^2}
=\frac1{4q(q+1)}.}
\]

q1即9/8和1/8，常数独立n、table数量及来源形状。本结论真正覆盖一个共同长尺度下的actual source score与多表联合；它不要求mirror，镜像条件的额外作用仍在非零score倾斜时消除先验漂移。

### 证据与边界

`JOINT_CHECK.py/JSON`实际执行284个Fraction Krawtchouk正交检查、44个跨表权重Cauchy检查、48个精确条件期望谱展开、12个具有不同score/Q_b的有限模型完整联合积分及28个真实box来源的逐点Cauchy检查，全部通过。有限模型仅验证恒等式，不冒称它们全部来自窗口；一般窗口结论由上述证明给出。root先前真实n1/2/3 grid扫描的F²/G及分离积分均小于9/8，但grid积分误差没有向外认证，不能替代此证明。

这只控制固定同一个L。q0短窗只能在该推导中用r²=1，费用至多n+1；共同连续尺度sup可能选取依赖query的尺度，不能直接把逐L积分帽移过sup。尚未推出C_n无log、C_n新阶或多尺度weak排序。

## 8. 条件性共同正算子接口及有限K的显式费用

第7节的跨表付款可以抽象为共同源端Hilbert空间上的frame不等式。设各table的source type先验都是同一个π，Hilbert空间为L²(π)，count posterior条件期望映射的Gram为G_b≥0。w_b为该空间中的非负函数，Σ_bw_b=1，π_b=E_πw_b。π_b=0的项照旧删除。对任意h∈L²(π)，加权Cauchy给

\[
\sum_b\frac{|\langle w_b,h\rangle|^2}{\pi_b}
\le\sum_bE_\pi(w_b|h|^2)=\|h\|_2^2.
\]

即正rank-one算子S=Σ_b(w_b⊗w_b)/π_b满足0≤S≤I。若已经有共同正算子0≤G_b≤G_0，则在逐点F²/G的同一次Cauchy分离后，合法逐表积分得到

\[
\int\frac{F^2}{G}
\le\sum_b\frac{\langle w_b,G_bw_b\rangle}{\pi_b}
\le\operatorname{tr}(G_0S)\le\operatorname{tr}G_0.
\]

这里不要求不同G_b交换或具有共同特征基；需要的是实际共同L²(π)坐标上的Löwner支配。

若单坐标class先验共同为η，u=√η，且已证明

\[
0\le H_b\le H_0:=uu^T+D,\qquad D\ge0,\quad Du=0,
\]

则正张量的Löwner单调性给H_b^{⊗n}≤H_0^{⊗n}；压缩到共同的对称source-type子空间得到G_b≤G_0。设D在u^⊥上的特征值为d_1,...,d_{K−1}，若每个d_j<1，则

\[
\int\frac{F^2}{G}\le h_n(1,d_1,\ldots,d_{K-1})
\le\prod_{j=1}^{K-1}(1-d_j)^{-1}.
\]

不必要求trD<1；若trD<1则另有较粗的1/(1−trD)帽。该接口的假设不能从逐表trace帽自动变成一个维数无关的低trace D。

不过，第42轮原Q_b下的严格帽tr(H_b−uu^T)≤β_q=1/[4q(q+1)]已经无条件给出

\[
H_b-uu^T\le\beta_qP_{u^\perp}.
\]

因此可选择明确的共同D=β_qP_{u^⊥}，得到任意有限K的严格联合界

\[
\boxed{\int\frac{F^2}{G}\le(1-\beta_q)^{-(K-1)}.}
\]

它独立于ambient n和table数R，但随K指数增长；若K=K_n随n增长，不能把该费用描述为维数无关。K2的第7节直接谱帽给9/8，比此处通用帽8/7更强，故二类别结论仍单独保留。更好的各向异性D若来自数值有限扫描，仍需要其连续来源族的严格认证。

root另运行了真实K4低维grid探路，数值F²/G最高约1.01019、分离积分约1.09918。其部分粗网格结果因质量求积偏差小于1，不能解释成负χ²或严格积分下界。本节依据共同正算子及第42轮一般定理，而非这些浮点网格数。共同尺度sup及q0短窗限制继续保留。

## 9. 稀疏共享长尺度的weak比较推论

这里保持第7节的同一固定二类别门控输入f及实际prior g，周期质量都为1。令𝓛⊂[1,∞)为有限或可数的共享长度集合，F_L=A_Lf、G_L=A_Lg，M_𝓛f=sup_{L∈𝓛}F_L。对0<c<1、λ>0定义

    E={M_𝓛f>λ, M_𝓛g≤cλ}.

每个x∈E至少有一个见证L满足F_L(x)>λ而0<G_L(x)≤cλ。设v=G_L/λ≤c，利用(1−v)²/v在0<v<1上递减，得到

\[
\frac{(F_L-G_L)^2}{G_L}
>\lambda\frac{(1-c)^2}{c}.
\]

因此对有限𝓛，用覆盖和积分界直接得

\[
\boxed{\lambda|E|\le\frac c{(1-c)^2}
\sum_{L\in\mathcal L}\frac1{4\lfloor L\rfloor(\lfloor L\rfloor+1)}.}
\]

对可数𝓛，可取E_j={max_{i≤j}F_{L_i}>λ, M_𝓛g≤cλ}。这些集合单调增加到E，且每个E_j使用同一前j项覆盖预算；测度单调极限给相同不等式。这里不能把两个有限max同时放进E_j后声称单调，因为有限prior max也随j增加；另一种合法方法是直接对单尺度见证集合取可数并集。

若每个[q,q+1)至多含一个L，则预算至多

\[
\sum_{q=1}^\infty\frac1{4q(q+1)}=\frac14.
\]

这个预算独立于稀疏尺度的数量。分解{M_𝓛f>λ}⊂{M_𝓛g>cλ}∪E，定义周期weak指标W_𝓛(h)=sup_{λ>0}λ|{M_𝓛h>λ}|，便得

\[
\boxed{W_\mathcal L(f)\le c^{-1}W_\mathcal L(g)
+\frac{c}{4(1-c)^2}.}
\]

若右边prior weak为无穷，该式按扩展实数理解；本推论不另行证明它有限。任意有限𝓛还可使用实际预算代替1/4。root实际Fraction预算探针得到前Q个单位区间总预算Q/[4(Q+1)]，以及(1,2)内J个不同长度的逐项预算J/8，代数与本推论吻合。

对周期L1输入，f∈L1_loc(R^n)，任意固定x的L↦A_Lf(x)在L>0连续：变化窗口包含在共同有界盒内，窗口对称差测度趋零，局部积分的绝对连续性及L^−n连续性给结论。因此连续长尺度sup可用可数有理长度表示，确保可测；但稠密单位区间的逐项预算发散，所以本推导不能控制连续sup。不能把∫sup替换成sup∫。

这是特殊归一化门控族与其实际prior在稀疏长尺度上的比较，既不是所有L1输入的统一强L1误差界，也不是一般HL上界。周期模型本身不是欧氏L1；既有有限尺度截断归约仍可用于合法下界转移，但此推论不自动支付完整欧氏连续最大算子的全部尺度。

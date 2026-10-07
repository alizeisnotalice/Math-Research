# 完整 nearmax 的平滑与输运：一阶障碍边界及双平移有限差分

2026-10-07。父任务指定 6.1-sol high。只新增本稿，不改他人文件或主账；本轮先解析，无数值执行、阶数拟合或大 MC。

**严格所得：** 任意完整保质量来源变换有原赢家冻结的有限 log 差分合同；其一阶平滑条件在精确极大点只是既有 S≤K、来源接触的推论。另有一条新的双平移 pooling 合同：两个完整响应的空间幅度差给正 logcosh 项，必须保留真实 cross-winner regret 与阈值跨越项。它没有 local-test→global 的量词漏洞，但仍未证明 regret 小，也未由径向窄峰推出Φ正增益或 n⁻¹⸍² 峰宽。一般浅壳/source-once/geom 未闭合。

使用已读 [M03](</Users/zhengzhihao/.codex/skills/math-m03-near-extremizer-stability/SKILL.md>)、[E02](</Users/zhengzhihao/.codex/skills/math-e02-nonlocal-hjb-supersolution/SKILL.md>) 的亏损、对称性、生成元与全空间条件审计；已读 method/cube-interface/provenance。这里自证有限差分，不引用外部稳定性、HJB比较或未核 subharmonicity。

查重并读 [log-max 乘性变分](log_max_variation_20261007.md)、[加性 obstacle](log_max_additive_obstacle_20261007.md)、[nearflat 几何](nearflat_cube_geometry_20261007.md)、[fixed jump obstacle](fixed_jump_obstacle_interface_20261007.md)、[resolvent repacking](resolvent_obstacle_repacking_20261007.md)。source S 的一阶 KKT、弱Lebesgue异常和旧固定生成元 odometer均不是新成果。目录检索未找到本文§4的完整双平移 pooling/threshold-crossing 合同。旧 nearflat 双向乘性 hinge 使用一个全局 v；本稿扰动则是两份整来源的同一空间平移，不能相互冒充。

## 1. 同一完整问题、原赢家和源势

完整有限正 μ，质量W>0；允许尺度为共同有限族或全连续[a,b]，a>0。闭cube全边长约定，M_μ=sup_R h_R*μ，原可测真赢家R，E={M_μ>τ}，I=τ|E|，

\[
 \Phi_\tau(\mu)=\tau\int\log_+(M_\mu/\tau)dx,
 \quad S(y)=\tau\int_E\frac{h_{R(x)}(x-y)}{M_\mu(x)}dx.
\]

记Z=1+nlog(b/a)，则S≤Z、I≤ZW、∫S dμ=I。K为同一完整来源/阈值/尺度族的Φ/W上确界。假定当前Φ≥(K−ε)W；不假定极值存在、紧性、来源形状或 nearmax 后验产品化。

加性稿已证：精确全局极大点有S≤K逐点、S=K在supp μ上；nearmax一般只给输入依赖的点态误差或无容量参数的弱Lebesgue异常体积。不能把后者提升成统一点态obstacle。本文所有平移/卷积/transport均作用于整份μ，原R、E仅用于下界，不把新赢家强行等同原赢家。

## 2. 任意保质量 smoothing 的真实有限 log 合同

令T为一份预先固定Markov来源变换，ν=Tμ≥0、ν总质量W。T可为卷积κ*μ、热半群P_hμ、固定跳核或完整几何pushforward；须保证ν属于定义K的合法来源类。完整有限正Borel类满足这个闭合；若坚持结构子类则另核变换资格。

对0<t<1，来源μ_t=(1−t)μ+tν，质量仍W。在原E冻结R，定义

\[
 r_T(x)=\frac{h_{R(x)}*\nu(x)}{M_\mu(x)}\ge0,
 \quad p_T=r_T-1\ge-1.
\]

新最大值至少原query对新完整来源的响应M_μ(1+t p_T)。E外新log项非负、E内log_+≥log，准确得到

\[
 \boxed{\tau\int_E\log(1+t p_T(x))dx
       \le\Phi_\tau(\mu_t)-\Phi_\tau(\mu)
       \le\epsilon W.} \tag{1}
\]

不要求新winner不变、不求winner导数，亦不删除阈值外的可能新收益。r_T≤W/(τaⁿ)，1+t p_T≥1−t，所有积分有限。

若需一阶形式，令

\[
 D_T=\tau\int_Ep_Tdx=\int S\,d(T\mu)-I,
 \quad Q_T=\tau\int_Ep_T^2dx<\infty.
\]

log在[1−t,∞)的Taylor下界给

\[
 \boxed{D_T\le\epsilon W/t+
                 tQ_T/[2(1-t)^2].} \tag{2}
\]

粗容量仅Q_T≤(1+d)ZW，d=W/(τaⁿ)；其中r_T²≤d r_T，负交叉−2r_T可保留或删除。没有统一的nearmax/Q_T小量。不能把ε固定而让d、热核峰值或输入norm隐入常数。

对卷积κ，准确D_T=∫[T* S−S]dμ，其中T*S(y)=∫S(y+z)κ(dz)。固定源µ只领一次W；没有每个z重新领W的平均weak合同。对可积概率族κ也可直接用其平均ν一次代入(1)。

## 3. 一阶非局部/扩散条件没有额外次调和性

若精确全局极大点真的存在，则由(2)令t↓0得D_T≤0。但既有S≤K和来源接触已经逐点给

\[
 T^*S(y)\le K=S(y)\quad(y\in\operatorname{supp}\mu),
 \qquad\int(T^*S-S)d\mu\le0. \tag{3}
\]

这对任意Markov核成立，只是上障碍与接触的推论；不说明(T*S−S)在整个空间非正，也不提供新的HJB/Green odometer。对nearmax，Lebesgue异常不能自动控制Tμ给异常的质量；卷积密度或测试容量需显式支付，正是加性稿的限制。

如要取热生成元极限，必须先给领域。例：f≥0、f∈W²,¹，P_h=exp(hΔ)，则(P_hf−f)/h→Δf于L¹；若需要直接对Φ差分取一阶，下界中相对响应需由该L¹误差和a⁻ⁿ/τ控制。精确max只能推出来源加权∫SΔf≤0，而非ΔS≤0或≥0全空间。若再假f光滑紧支撑，旧接触S=K在supp f上、supp Δf⊂supp f直接使∫SΔf=K∫Δf=0；这不构成独立的扩散刚性。

几何transport同样如此。对光滑紧支撑向量场v和合法小步diffeomorphism，推来源的一阶为−div(vf)；只有在这份L¹变分已证时才可写−∫S div(vf)。未经S的弱导数资格，不能继续写f∇S，更不能用各receiver自选方向。精确接触使这类支撑内一阶条件仍属旧KKT。

### 3.1 真正exactmax的平滑可严格降低Φ

用合法退化窗口a=b、n=1说明不能对所有nearmax免费宣称smoothing正增益。固定µ=Wδ₀、τ=W/(ea)。只有h_a，完整响应W/a在Q_a上、外为0。标量τlog_+(u/τ)≤u/e及∫h_aµ=W给K=1/e，当前µ准确达值Φ=W/e。

取0<z<a/2，来源ν_z=W(δ_{−z}+δ_z)/2，两份整来源平移混合。响应在重叠interval长a−2z上为W/a，在两边总长4z上为W/(2a)，故

\[
 \frac{\Phi_\tau(\nu_z)}W
 =\frac1e\left[1+\frac{2z}{a}(1-2\log2)\right]<1/e.
 \tag{4}
\]

原score精确为S(y)=(1/e)(1−|y|/a)_+，满足旧pointwise obstacle/contact，但分布二阶导数
S''=(ea)⁻¹(δ_{−a}−2δ₀+δ_a)同时有正负部分。不能把它当全空间subharmonic或superharmonic。这个完整source exactmax例只审计普适平滑/生成元推断；它不是b>a大K浅壳的反例。正宽化L¹ nearmax会逼近该固定尺度现象，不宣称δ模型已具有原soft/history资格。

## 4. 新双平移 pooling：幅度增益、真实regret与阈值跨越

固定一个全局z∈Rⁿ，整来源µ_+=T_zµ、µ_−=T_−zµ，ν=(µ_++µ_−)/2。它们各质量W；平移协变精确给Φ(µ_±)=Φ(µ)、M_±(x)=M_µ(x∓z)，R_±(x)=R(x∓z)是各自原truewinner。

写A=M_+(x)、B=M_−(x)、T₀=(A+B)/2，并定义两个原winner对另一完整来源的regret

\[
 d_1=B-h_{R_+}*\mu_-,\quad
 d_2=A-h_{R_-}*\mu_+,\quad
 d=\min(d_1,d_2),\quad\beta=d/(A+B).
\]

所有d_i≥0；0≤d≤min(A,B)，所以A+B>0时β∈[0,1/2]，A+B=0置β=0。将这两个合法原尺度用于ν，得到

\[
 \boxed{M_\nu\ge\tfrac12(A+B-d)=T_0(1-\beta).} \tag{5}
\]

该选择只取两个**完整输入平移后的真赢家**，没有为每receiver选择不同来源扰动，也没有把sup和卷积交换。新完整M_ν可更大；只用(5)作下界。

记F₀={T₀>τ}，c_β=−log(1−β)≤2β。正部对向下log平移的1-Lipschitz界给

\[
 \Phi(\nu)\ge\tau\int\log_+(T_0/\tau)dx
                  -\tau\int_{F_0}c_\beta dx.
\]

F₀有限体积，β≤1/2，regret成本可积。定义scalar pooling差
Γ(a,b)=log_+((a+b)/2)−[log_+a+log_+b]/2。

当a,b>1，Γ=logcosh((log a−log b)/2)≥0；当a,b≤1，Γ=0。单边跨阈值时，假a>1≥b≥0，最坏b=0：1<a≤2给Γ≥−(log a)/2，a≥2给Γ≥(log a)/2−log2，所以统一Γ≥−(log2)/2。只在两份严格E的对称差上可能出现负项，阈值plateau并未删掉。

令E_±={M_±>τ}=E±z。由两份平移的Φ都等原Φ、ν的全局nearmax上界≤KW，严格得到

\[
 \boxed{\tau\int_{E_+\cap E_-}
  \log\cosh\left(\frac{\log A-\log B}{2}\right)dx
 \le\epsilon W+\frac{\tau\log2}{2}|E_+\triangle E_-|
                  +\tau\int_{F_0}[-\log(1-\beta)]dx.} \tag{6}
\]

这是独立有限差分约束，不要求source平坦来推导。nearmax的同一输入K与质量W只用一次；来源µ_±质量各W不表示另付两份原input。两个原Φ的算术平均仍是Φ(µ)，而扰动ν质量W。

若两个平移winner尺度在一个区域相同，则该处β=0；幅度变化产生的logcosh增益只能由阈值跨越与其它区域regret抵销。一般尺度不相同，β不可删除。阈值项也不可删除：固定尺度例(4)恰有真实边缘损失，虽winner完全相同。

因为|E|有限，z→0时|E_+△E_-|→0（L¹平移连续）。这没有统一速率；除以|z|²前还需真实perimeter/threshold层控制。也没有证明β积分比幅度增益小。故(6)不能被口头升级成平滑正gain、Dirichlet能量下界或highK峰宽。

## 5. 径向窄峰的来源一次信息，与平滑缺项的准确连接

全连续[a,b]且b>a时，固定概率核混合

\[
 \overline U(x)=\frac1{\log(b/a)}\int_a^bU_s(x)\frac{ds}s,
 \qquad\int\overline U=W.
\]

原完整weak joint P_X=uniform E，m=I/W，有

\[
 \mathbb E_P[\overline U/M]\le1/m. \tag{7}
\]

这只用M>τ和固定核的一次来源质量，属于已知aggregate fixed-reference机制。high m确实迫使原receiver平均相对响应的log-scale宽度小，不能再用nearuniform宽峰例子否认这份资格。

原posterior age A_age=nlog[R/max(a,D)]（避免与§4的响应A混用）对任意有限正Borel来源有准确Stieltjes/layer-cake身份

\[
 G_R(x):=\mathbb E_{\pi_R}e^{A_{\rm age}}
     =1+n\int_a^{R(x)}\frac{U_s(x)}{M(x)}\frac{ds}s.
\]

closed/open中间cube的差只影响径向参数的可数质量原子，ds积分同值，不假captureCDF可微。由固定核Fubini还得一份真实joint矩

\[
 \boxed{\mathbb E_Pe^{A_{\rm age}}\le1+B_0/m,
 \quad\int e^{A_{\rm age}}dJ\le I+B_0W,
 \quad B_0=n\log(b/a).} \tag{8}
\]

这没有角独立或receiver采样替换Lebesgue，也不是critical moment的uniform常数；high m时它比逐row的1+B₀更强，却仍含未知m。熵身份的方向不能从该上界推出m小。

要把径向窄峰变成(6)的严格正Φ增益，仍需一条实际空间桥：为一份固定全局平移/卷积尺度，high m输入必须产生足够幅度差，同时两平移winner的cross-regret与阈值跨越总费更小。相对径向积分(7)不自动控制β：µ(x∓z)在相同receiver的完整捕获改变了来源位置和query中心，各自最佳尺度可能不同；这些正是(5)保留的量。

若使用来源dilation y↦cy，须另记准确换窗/阈值身份

\[
 \Phi_\tau^{[a,b]}((D_c)_\#\mu)
      =\Phi_{c^n\tau}^{[a/c,b/c]}(\mu).
\]

质量W保持，但同τ、同窗口下不是原平移对称。不能把各receiver径向最优微调L当作一次source dilation；新尺度越出原窗口或阈值改变均需结算。本文未找到从(7)一般地强迫width≥n⁻¹⸍²的合法变换。

## 6. 正则化、旧固定障碍及actual迁移边界

对完整L¹来源，源TV扰动给|M_ν−M_μ|≤M_|ν−μ|及

\[
 |\Phi_\tau(\nu)-\Phi_\tau(\mu)|
                    \le Z\|\nu-\mu\|_{\rm TV}.
\]

因为u↦τlog_+(u/τ)是1-Lipschitz。所以截断/正compact mollifier可构造任意精度的光滑紧支撑nearmax来源（允许归一化损失显式随截断小量调整）。这只消除正则性门槛，没有输入统一的近似速率、平滑正gain、密度上界或source几何稳定性。

旧fixed-jump/resolvent障碍从固定自伴生成元构造u≥0及ν−Lu=μ_cap≤κ，需真实源、域外符号、有限占用与共同算子。当前S是原全赢家列，既有source KKT不提供这些odometer数据；不能把(3)称成同一个Lu≥ν−κ，亦不能将旧fixed W_b与当前完整W混同。E02的特定stable HJB比较也不赋予S新的全空间超/次解性质。

完整Φ扰动改变original FIRST、fullfuture、CP/GP、birth、LCA与history。这些门的导数/输运残量未出现于(1)或(6)，不能将它们直接认作当前R_angle或heavy paid费。当前严格新接口为双平移(6)，配套未付字段是其真实regret和threshold跨越；(1)是一份通用有限扰动合同，(3)明确一阶平滑没有超出旧obstacle。一般highK浅壳重复占用仍需新的空间证明。本轮没有运行数值，也没有新的普适阶数结论。

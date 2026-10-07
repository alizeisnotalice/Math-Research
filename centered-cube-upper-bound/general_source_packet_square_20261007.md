# 完整原来源的 packet 平方预算：fullness、dominance 与真实覆盖

2026-10-07。只新增本稿，不改其它文件或主账。按父任务 6.1-sol high 要求工作；工具没有公开实际运行模型 ID。独立核验一般完整 μ 的 source-square 引理，不添加背景来源，不执行数值或 toy 守卫。

**引理通过：预固定不重计来源分配的 fullness/dominance incidences，具有 4/(ηθ) 倍原 weak receiver 质量的 source-square 预算。** 它不需要 packet 小直径、输入产品结构或 winner 驻点。θ=n^{-1/2}、η为常数时，只支付这个被明确截出的平方子支；剩余 incidences和它们的来源交叉仍须实际覆盖或支付。

使用 [A01](/Users/zhengzhihao/.codex/skills/math-a01-fractional-normalization/SKILL.md) 的原分母/齐次量核对及 [G01](/Users/zhengzhihao/.codex/skills/math-g01-carleson-measure-box-and-tent-test/SKILL.md) 的来源重计审计，读取 SKILL、method/cube-interface及 provenance入口。下文是自含正重排/Tonelli/Jensen证明，不引用未验树或 Carleson 嵌入。

## 1. 查重与一般化范围

已全文读取 [source_partition_dominance_extension](source_partition_dominance_extension_20261007.md)、[shift_averaged_source_dominance](shift_averaged_source_dominance_20261007.md) 和 [single_capture_source_rearrangement](single_capture_source_rearrangement_20261007.md)，并检索旧 source-operator/dominance稿。

前两稿以固定 small-box packet 和输出 dominance支付完整 hardband交通与 Γ 的径向能量；它们不自动支付本稿 S(y) 的来源平方。single_capture稿 §4.2 已给大份额原子与 fullness packet 的同 1/(ηθ) 阶来源平方，那里原完整模型为 spike/packet来源加 β1_D，分母 βR^n+捕获 spike质量是 additive。

本稿明确移到**任意原完整 μ、无新背景**：两条实际下界只先给 max，而非 additive分母，因此用 half-sum比较得到常数4。来源合同、profile量和未付覆盖必须保持；这不是另一个一般 √n定理，不给同一旧实际事件重复领费。

## 2. 原输入、分母及预固定 packet

固定完整有限正 Borel来源 μ，质量 W>0，原可测 finite truewinner side length R(x)>0，

\[
 Q_x=x+[-R(x)/2,R(x)/2]^n,\quad
 m(x)=μ(Q_x),\quad M(x)=m(x)/R(x)^n,
 \quad E=\{M>τ\},\quad I_1=τ|E|.
 \tag{1}
\]

τ>0。同一完整 μ、τ、R、E 全程冻结。finite尺度集使 |E|≤N W/τ<∞（仅为有限性，不采用 N作为费用）；bounded continuous window的已知正包络也给有限性。证明实际上仅用 |E|<∞、捕获及 m>τR^n，不用最大性；可以套原已经保存的 selected rows，不重判 winner。

预先选至多可数的 measurable来源分配

\[
 μ_i=a_iμ,\quad a_i\ge0,\quad\sum_i a_i\le1\quad μ\text{-a.e.},
 \quad M_i=μ_i(\mathbb R^n)>0,\quad\sum_iM_i\le W.
 \tag{2}
\]

零质量项忽略。这可以依完整 μ选择，但不能在积分的 receiver x处重组、归一化或重新分配。盒可以重叠；必须以 (2) 控质量而非仅数支撑盒。未分配来源 μ_0=(1−Σ_i a_i)μ继续属于原 m、M、E和未付 profile。

给 η,θ∈(0,1]，定义真实 capture quantities及合格 rows

\[
 m_i(x)=μ_i(Q_x),\qquad
 E_i=\{x\in E:m_i(x)\ge θm(x),\quad m_i(x)\ge ηM_i\}.
 \tag{3}
\]

两条等号均给该 paid子支，补集使用严格失败；若主账选相反等号规范须整体一致。分母都是完整原 m。来源分配固定、(x,y)↦1_{Q_x}(y) measurable，故 m_i、E_i measurable。

仅保留该 packet 的来源-receiver incidence，定义

\[
 S_i(y)=τ\int_{E_i}\frac{\mathbf1_{Q_x}(y)}{m(x)}dx,
 \qquad S_{dom}(y)=\sum_i a_i(y)S_i(y).
 \tag{4}
\]

S_i不是以 μ_i重启 maximal/FIRST 后的 kernel，也没有把 Mi换成当前 m_i。

## 3. 两条真实下界及半和比较

在 x∈E_i、y∈Q_x 时，R(x)≥2||x−y||∞。由 μ_i≤μ、fullness及原超阈门，

\[
 m(x)\ge m_i(x)\ge ηM_i,\qquad
 m(x)>τR(x)^n\ge τ(2\|x-y\|_\infty)^n.
 \tag{5}
\]

因此（令 V_y(x)=(2||x−y||∞)^n）

\[
 \frac τ{m(x)}\mathbf1_{Q_x}(y)
 \le\frac{2τ}{ηM_i+τV_y(x)}\quad(x\in E_i).
 \tag{6}
\]

这是 max(A,B)≥(A+B)/2；不能无背景地写 m≥ηMi+τR^n。捕获以外左边为0，故 (6) 可在全部 E_i上积分。

## 4. 逐原来源的 cube径向重排

|{x:V_y(x)≤u}|=u，对每个固定 y准确成立，n不产生额外系数。对递减 φ(u)=τ/(ηMi+τu) 和任意 Lebesgue集合 A、|A|=v，层饼给

\[
 \int_Aφ(V_y(x))dx
 \le\int_0^vφ(u)du
 =\log\left(1+\frac{τv}{ηM_i}\right).
 \tag{7}
\]

具体是每个 level set的交集体积≤min(v,level-set volume)，再正 Tonelli；中心 cube {V_y≤v}取等号。A无需为 cube，y不必是 atom。由 (6)，每个原 y均有

\[
 S_i(y)\le2\log\left(1+\frac{τ|E_i|}{ηM_i}\right).
 \tag{8}
\]

|Ei|=0时两边为0，Mi>0排除零分母。以 2√u≤1+u 给

\[
 \log(1+t)=\int_0^t\frac{du}{1+u}\le\sqrt t
 \quad(t\ge0).
 \tag{9}
\]

故独立于 packet 位置/形状/原子数/维数的 labelled square为

\[
 \boxed{\int S_i(y)^2μ_i(dy)
 \le4M_i\log^2\left(1+\frac{τ|E_i|}{ηM_i}\right)
 \le\frac{4τ}{η}|E_i|.}
 \tag{10}
\]

这里没有 β1_D、背景质量或来源坐标独立假设。partial capture只通过真实 m_i≥ηMi使用，未升级为 full capture。

## 5. 每行 packet数与 fractional Jensen

从 (2)，Σ_i m_i(x)≤m(x)。每个 E_i命中至少 θm(x)>0，所以

\[
 \sum_i\mathbf1_{E_i}(x)\le1/θ\quad(x\in E).
 \tag{11}
\]

Tonelli与 (10)得

\[
 \boxed{\sum_i\int S_i^2dμ_i\le\frac4{ηθ}I_1.}
 \tag{12}
\]

fractional profiles不能先无权相加后平方。正确逐来源 Jensen/Cauchy为

\[
 S_{dom}(y)^2\le\left(\sum_i a_i(y)\right)
                    \sum_i a_i(y)S_i(y)^2
 \le\sum_i a_i(y)S_i(y)^2,
 \quad
 \boxed{\int S_{dom}^2dμ\le\frac4{ηθ}I_1.}
 \tag{13}
\]

对 countable indices先取 finite partial sums，再正单调极限。Σa<1时把其余质量作为零接受标签也直接给 Jensen。

η为常数、θ=n^{-1/2}时，(12)/(13)是该部分 source-square≤(4/η)√n I1。右边为原 receiver质量 I1，不是自动≤√nW，也不是已经支付 current R_angle。

## 6. 实际覆盖量、未付 incidence及一个可核覆盖证书

原全profile S(y)=τ∫_E1_{Q_x}(y)/m(x)dx，∫S dμ=I1。精确写

\[
 S=S_{dom}+S_{res},\quad
 S_{res}(y)=τ\int_E\frac{1_{Q_x}(y)}{m(x)}
 \left[a_0(y)+\sum_i a_i(y)\mathbf1_{x\notin E_i}\right]dx,
 \quad a_0=1-\sum_i a_i.
 \tag{14}
\]

残余包含未分配来源、dominance失败以及fullness失败。receiver Ei重叠是允许且已由 (11)收费；source上dom/res相加产生交叉，只有

\[
 \int S^2dμ\le2\int S_{dom}^2dμ+2\int S_{res}^2dμ
 \tag{15}
\]

等合法 Young/Minkowski合同，不能按 incidence不交就删除该 square cross。

可直接保存/核验的 row覆盖量为

\[
 q_{dom}(x)=\sum_i\mathbf1_{E_i}(x)\frac{m_i(x)}{m(x)}\in[0,1],
 \quad I_{dom}=\int S_{dom}dμ=τ\int_Eq_{dom}(x)dx.
 \tag{16}
\]

若独立证明 receiver平均覆盖 τ∫E q_dom≥κI1（κ>0），则 Cauchy与 (13)严格推出

\[
 \kappa^2 I_1^2\le I_{dom}^2
 \le W\int S_{dom}^2dμ\le\frac{4W}{ηθ}I_1,
 \qquad
 \boxed{I_1\le\frac4{ηθ\kappa^2}W.}
 \tag{17}
\]

|E|=0无需除I1。pointwise q_dom≥κ是充分但更强的证书；只需原真实 receiver平均捕获覆盖。**κ一般保证尚未证明。** (17)是明确的 conditional覆盖回代，不是由 packet分解存在就获得的上界，也不据此登记新阶或数值。

每row存在一个合格packet只保证 q_dom≥θ，不能据此把κ当维数无关常数。θ=n^{-1/2}时，所需常数覆盖必须由合格packets的实际总捕获份额另证，或者用旧空间small-box输出覆盖支付相应row。

## 7. 固定合法空间 packets 如何真的覆盖

小直径 packets、固定空间 cells、输入自适应 clusters都可以合法选择，只要其 a_i及Mi在 receiver积分前冻结且满足 (2)。若来源 μ_i集中于固定 B_i，Q_x完整包含 B_i，则 m_i=Mi，fullness自动成立（η≤1）；仍须计算 Mi/m≥θ。若只部分切过 B_i，盒很小也不保证 m_i≥ηMi，必须在原 Q_x实算。

旧 small-box dominance覆盖是另一合法工具：side≤a/n时只要 m_i≥θm，就能按 Mi一次覆盖支配输出并支付完整 hardband交通；不要求本稿 fullness。它与本稿能互补使用，但量和资格不同，不能把同一已经删除输出交通又计一次source-square实际费用。

现有 multires_packet_network保存输入可用原八点packet标签以及生成公式的小盒支撑；source_partition稿已经核其 Σμ_i=μ。但必须重新**读取/核算原行数据**的 m_i/m、m_i/Mi及q_dom，才能知道本稿截出多少交通；packet中有八个atoms、平均捕获数或固定微格协同标签均不能代替这些比值。本稿未跑该核算，未宣称现有样本已覆盖。

small份额不等于几何diffuse：把同一空间packet人为细分或分成多个同位置 fractional标签，会让每个 m_i/m低于θ而没有改变任何原 cube几何。Mi阈值和partial fullness一起是来源分配的真实统计量，失败补集仍不能仅凭标签大小叫nonconcentration。真正空间分散必须由固定支撑/距离/packing字段另证。

多尺度全部树节点的 μ|P通常重复来源，不能直接放入 (2)。可选来源反链，或预固定 a_i=w_i1_{P_i}且每原来源路径Σw_i≤1；后者的 dominance测试必须使用实际缩小后的 m_i，不保留未缩放的旧份额。多个独立分解scheme若逐个重领W、最后取任何scheme成功的并集，亦不满足本稿合同；需共享来源权或另付scheme重数。

因此当前有一个可执行的一般certificate审计：固定原合法packets后，保存每个原 row的完整m、全部正m_i、Mi、η/θ等号规范，计算q_dom及原source weighted residual。下一解析/probe应针对 (14)残余的真实重用，或证明 (16)receiver覆盖，而非按x拼出刚好被捕获的临时packet。

## 8. 完成状态与实际迁移边界

无背景的 (8)–(13)严格成立，原任意正输入及 finite winner保留；只用 capture/fullness/dominance和原M>τ。查重确认旧dominance交通/Γ与background packet square相关但不等同，当前常数4扩展已独审。

没有一般合法packets覆盖原所有source incidence的证明，也未支付 S_res及current soft FIRST/fullfuture/early/CPGP/LCA/allhistory核。若用于后者，必须给原kernel与本稿 restricted hard incidence的真实支配/映射，不能凭source-square名称迁移。未新增主账paid费，不添加任何background质量；没有新猜测阶数、toy检查或旧守卫重跑。

# E01：继承规格引用的冻结来源定位

继承规格所指的旧审计文件未在输入中找到。这里逐卡列出冻结批次已记录的来源、假设和定位；本索引没有新增独立数学审核。未读外引、私有接口、被隔离原式和两轮缺失保持原状态。

## P-0297391ab7ac5c7b

[冻结阅读卡](papers/TEAM_EF-P-0297391ab7ac5c7b.json)，输入SHA-256 `d8cdcae3be480c34293629cb21594680450afd9b71cdb8cb2fd9cd622faf755b`。

Theorem 2: full weak LP and modal weak LP equivalence

完整occupation LP(8)与m个modal measure LP(10)最优值相等。该等价发生在两个弱LP之间；不默认原始非凸switching OCP无gap或存在最优switch。

假设：["T和x_0,x_T固定；X紧；f_j,l_j连续，f_j对x Lipschitz。", "U为有限单位向量集合，动态与cost对mode coordinate u_j线性。"]

定位：{"start_line": 409, "end_line": 511, "pdf_pages": [9, 10]}

证明骨架：["正向取μ_j=(π_{t,x})#(u_jμ)，弱流与cost逐项保持。", "逆向μbar=Σ_jμ_j；ū_j=dμ_j/dμbar≥0且Σū_j=1 μbar-a.e.，明确构造 μ(dt,dx,du)=Σ_j ū_j(t,x)δ_{e_j}(du)μbar(dt,dx)，保持每个mode矩。", "这是TEAM_EF核验的有限mode投影/合并构件；正文Theorem2证明首两段的(8)/(10)标签颠倒且δ1(du_j)写法含糊，按上述有类型构造读取。"]

缺口：["不从此推原switching controls recoverability；p_S=p_W的Theorem1仅略述并引用Vinter。"]

Primal attained; HJB duality (Lemma 2, Theorem 3)

可行modal LP有有限总质量T，连续cost下弱星紧给primal minimum；dual为sup_v[v(T,x_T)−v(0,x_0)] s.t. l_j−(∂_tv+∇v·f_j)≥0 on K。作者陈述primal/dual无gap，未断言C¹ dual supremum attained。

假设：["上述compact/continuous数据；modal LP可行（作者Theorem3）。"]

定位：{"start_line": 518, "end_line": 664, "pdf_pages": [11, 12, 13]}

证明骨架：["v=t固定总质量；Alaoglu在紧K上给measure subsequence，弱约束和连续cost保极限。", "dual slack积分给对每个可行trajectory/measure的lower certificate。", "作者借Anderson–Nash cone closure criterion证no gap；本次保原文依赖，不独立重证其拓扑conic dual theorem。"]

缺口：["Corollary2以τ_j=spt π_tμ_j而非加权mode marginal积分；同页把τ_j再写为一个时间积分存在集合/标量类型冲突，暂不将该support-unweighted criterion列为可执行认证。"]

Theorems 4–5: polynomial hierarchy converges to weak modal LP

无限moment formulation等价weak modal LP；有限阶p_d单调非降并渐近趋p*，每阶为最小成本lower bound。不自动证明p*等于原始switching infimum。

假设：["f_j,l_j及K defining functions均为多项式；K basic semialgebraic，描述中含对全部(time,state)变量的ball constraint。", "所有moment/localizing PSD约束和允许test degrees随阶增加；有限SDP只能作为松弛。"]

定位：{"start_line": 667, "end_line": 812, "pdf_pages": [14, 15, 16]}

证明骨架：["将弱流对单项式测试写成线性moments，cost亦线性。", "Archimedean Putinar criterion给全阶PSD序列的representing measures（外引Lasserre）。", "截断约束扩大可行域，阶增加缩小，因此p_d≤p_{d+1}≤p*；compactness/Putinar给作者asymptotic result。"]

缺口：["核心representing-measure与conic finite-SDP no-gap外引未独立重证；numerical positivity/residual仍需认证。"]

Joint modal reconstruction (Section 4)

按共同时间格点约束Σ_jΣ_k μ̃_j{(t_i,x_k)}=Δt_i，再同时拟合modal moments，估state barycenter与mode duty cycles；是candidate recovery，需真实dynamics integration与feasible cost来验收。

假设：["已知有限modal moments；作者解释假设optimal pair unique且scalar state以便说明。"]

定位：{"start_line": 942, "end_line": 1036, "pdf_pages": [19, 20, 21]}

证明骨架：["分开重构各modal moments未保Σmodal=dt的时间边缘，会破坏duty cycle。", "共同时间质量约束保证每个时间cell conditional state/mode mass规范化。", "候选经local solver/refinement后比较feasible upper cost与certified lower bound，不能把有限moment matching直接称真实轨迹。"]

缺口：["网格非零atoms不自动形成unique trajectory；作者experiments未由TEAM_EF重跑。"]

常数：无dimension-free norm定理；总质量T精确，mode数m与state维数n影响moment block大小。§3.4是rough real-arithmetic complexity estimate；不作为实际CPU保证。

端点与限制：["Young/convexified p_S可能严格低于原switching OCP，原文p6明确保state/endpoint constraint gap和inward-pointing外引条件。", "在切换小例p24，打印dynamic为ẋ=a_σ x(t)，但后文x(t)=1/2−t及cost1/24对应constant drift a_σ；本次PDF已确认该矛盾，不能把该例直接用于实算。", "p24和p25数值只是作者实验；不代表有限阶SOS已普遍精确。", "积分约束、closed-loop、impulse extensions多为sketch；每一种重新验证可实现性与对应no-gap假设。"]

## P-078f65ff156d77ea

[冻结阅读卡](papers/EG-P-078f65ff156d77ea.json)，输入SHA-256 `3364b59e511f0d70cb415e13a4a716d3a906ec04e1d581843cb0105f6192d725`。

定理1：有限 N 代理目标与经验占用测度的等价表达

N 代理目标精确等于经验测度对上的线性运行项、线性终端项和相同时间下的二次交互项；同一经验测度还满足初始经验分布 ρᴺ₀ 对应的弱 Liouville 等式。交互项中 δ(t−t′) 只是相同时间积分的形式写法。

假设：["各代理轨迹满足共享的连续动力学，f 关于状态局部 Lipschitz 且对控制一致；µᴺ、νᴺ 是单轨迹占用测度的经验平均。", "运行成本 ℓ₀ 与终端成本 Ψ 有界连续；交互核 W 可测且对应积分有限。"]

定位：{"start_line": 206, "end_line": 250, "pdf_pages": [4]}

证明骨架：["逐代理使用 µ[ωᵢ] 定义换写运行成本，使用 ν[ωᵢ] 换写终端成本。", "将经验占用测度定义代入 W 项，双重求和恰还原 ∑ᵢ(W*ρᴺ_t)(xᵢ(t))/N；再把 t′ Dirac 表示为同一时间交互。", "沿每条轨迹对 v(t,xᵢ(t)) 用链式法则积分并对 i 平均，得到弱 Liouville 关系。", "我视觉核对了 PDF 第4页公式 (7)–(9)，检查双重求和系数与严格时间积分对应；未对有限代理近似任何连续控制策略作推导。"]

缺口：["这是精确的有限代理表示，不是 N→∞ 的 mean-field 收敛定理。", "相互作用目标是 quadratic measure functional；整体优化仅在核 PSD/交互权重符号合适等条件下凸。"]

定理2与推论1：PSD 核下目标凸、可行集弱*紧且存在最优解

若 W 为正半定，则交互二次项凸，J 在相应正测度空间上凸；Δ 由线性弱 Liouville 约束定义，是弱*紧凸集。J 在 Δ 上弱*下半连续，故存在最优 (µ*,ν*)。

假设：["X 与 U 紧，Σ=[0,T]×X×U 紧；f 连续；初始 ρ₀ 固定。", "J 中运行成本 ℓ₀ 与终端成本 Ψ 有界连续。", "W 是对称正半定核；另外，正文的交互权重 λ 应非负或并入 PSD 核，尽管定理段落未明写符号条件。"]

定位：{"start_line": 364, "end_line": 372, "pdf_pages": [6]}

证明骨架：["线性运行/终端项自动凸；把交互项写为 signed-measure 二次型，插值余项为 −α(1−α)∫∫W d(µ₁−µ₂)d(µ₁−µ₂)，PSD 给非负曲率。", "Δ 为正测度对空间与连续线性 Liouville 等式的交集，因底空间紧、质量固定而弱*紧凸。", "对最小化序列取弱*收敛子列，线性项弱*连续，PSD 连续核二次泛函弱*下半连续，极限点达到最优值。", "我视觉核对了 PDF p6 定理2/推论1和凸性恒等式；我没有独立重证弱*紧性或二次泛函的下半连续通用结论。"]

缺口：["convexity 的交互项外权重符号需要 λ≥0；在所读正文定理的显式假设列表中未见该符号，应用时需把 λ 的约定写清楚。", "弱*紧性使用固定总质量与 compact 状态/控制/时间域；非紧域不能直接使用此结论。"]

定理3–4：方向导数与凸问题的一阶最优性变分不等式

Gateaux 导数是线性运行项加终端项及相互作用梯度：δJ=⟨g_µ,δµ⟩+⟨Ψ,δν⟩，其中 g_µ=ℓ₀+2λ(W*ρ_t)(x) 为同一时刻的状态卷积。若 J 在 Δ 上凸，则可行点最优当且仅当对所有可行 (µ,ν)，⟨g_{µ*},µ−µ*⟩+⟨Ψ,ν−ν*⟩≥0。

假设：["(µ,ν)∈Δ；方向为有限 signed measures (δµ,δν)，保持总质量并满足线性化 Liouville 约束。", "定理4另假设 J 在 Δ 上凸。"]

定位：{"start_line": 374, "end_line": 444, "pdf_pages": [6, 7]}

证明骨架：["对二次相互作用在 µ+εδµ 上展开，二次余项除以 ε 后趋零，得到双线性主项的两份相等贡献，形成系数 2。", "凸可微函数的一阶支撑不等式给必要充分 VI；必要性由方向导数，充分性由凸性。", "我直接由对称二次型展开检查因子 2 与 VI 方向符号；未将可行方向说成 Δ 的全部代数方向，非负性边界仍由凸可行集最优条件处理。"]

缺口：["interaction derivative uses same-time interpretation of δ(t−t′);严格实现应依赖式(9)的 dt 形式。", "不显式处理非光滑成本或不可微核的次梯度版本。"]

定理5：精确 FW 线性预言机由初态参数化经典 OCP 实现

把每一初始状态的最优经典轨迹占用测度按 ρ₀ 聚合，所得 (µ̄ᵏ,ν̄ᵏ)∈Δ 是整体 FW LMO 的解。由此每个 FW 线性方向可以分解为互相独立的、按初态索引的 classical optimal-control 子问题。

假设：["固定迭代 µᵏ 与其线性化成本 g_{µᵏ}。", "对 ρ₀-a.e. 初态 ξ，经典轨迹问题 min∫₀ᵀg_{µᵏ}(t,x,u)dt+Ψ(x(T)) 具有最小解；且存在 ρ₀ 可测的最优轨迹选择 ξ↦ωᵏ(ξ)。", "若最小值不达，正文说明以 ε-optimal 轨迹理解，需相应采用 inexact-oracle 收敛结果。"]

定位：{"start_line": 456, "end_line": 501, "pdf_pages": [7, 8]}

证明骨架：["附录 Lemma1 将 µ disintegrate 为 dtρ_tλ_{t,x}；对可行 Liouville 式应用 continuity-equation 表示与 superposition，得到路径空间概率 η。", "对 η-a.e. relaxed trajectory，把平均状态速度与线性化成本速度增广成一个 convexified differential inclusion；附录 Lemma3 引用 chattering 定理构造均匀逼近的经典轨迹，端点和代价一并逼近。", "由此每条 relaxed path 的线性化成本不小于其初态 ξ 的 classical OCP 值 V(ξ)；任一测度解总体成本下界为 ∫V(ξ)dρ₀。", "聚合 measurable selection 中每个 ξ 的经典 minimizer 达到该积分下界，因此聚合测度是 LMO 最小解。", "我阅读了 Appendix A–C 的完整论证链，并核对下界与达到两侧如何衔接；continuity-equation superposition/chattering 与 measurability 定理引用外部来源，未独立重证。"]

缺口：["定理要求逐初态存在解和可测最优选择；作者承认非紧设置可能只有 ε-optimal 解，此时精确 oracle 结论应改用 inexact 版本。", "区间轨迹松弛和点态控制核可由 chattering 近似，不意味着任一 measure-optimal solution 是单一 classical feedback trajectory。"]

定理6：L-smooth convex 目标下的 FW 速率与 fully-corrective 变体

对 k≥1，J(µᵏ,νᵏ)−J*≤8L(T+1)²/(k+2)，为 O(1/k)。若采用 FCFW 并精确解累计凸组合权重的修正问题，其目标值不高于同一字典上的普通 FW 更新，所以相同下降率成立。W 非 PSD 时目标可能非凸；文章仅援引标准 FW 在适当条件下趋向 VI 驻点，不给全局 gap 率。

假设：["J 在 Δ 上凸并对积 total variation 范数 L-smooth。", "算法1 每步使用精确 LMO，步长 αₖ=2/(k+2)；有限 horizon T，µ 总质量 T，ν 总质量 1。"]

定位：{"start_line": 548, "end_line": 609, "pdf_pages": [8, 9]}

证明骨架：["定义 L-smooth 为 J(z′)≤J(z)+δJ(z;z′−z)+(L/2)‖z′−z‖²。", "每个线性化子问题最优性与最优可行点 z* 给 δJ(zₖ;sₖ−zₖ)≤J*−J(zₖ)。", "单轨迹测度与 ν 的总变差差分别至多 2T 与 2；按积 TV 求和 ‖sₖ−zₖ‖≤2(T+1)，得 Δₖ₊₁≤(1−αₖ)Δₖ+2L(T+1)²αₖ²。", "用 αₖ=2/(k+2) 的归纳得 8L(T+1)²/(k+2)。p9 视觉复核了 Theorem6 的步长与系数，附录 D 中我复算了递推常数。", "FCFW 对历史聚合占用测度权重做 simplex QP；其中包含以当前点进行普通 FW 步的可行候选，故精确校正目标不增。"]

缺口：["L 未在本文主定理内进一步显式界为 λ、‖W‖∞ 等模型数据的函数；速率需以已知 smoothness 常数为前提。", "数值 LMO 由网格控制序列的 Adam 梯度优化近似，FCFW QP 也用投影梯度；实验没有 oracle/QP 误差认证，不能直接把 exact-theorem gap 数值当作实际运行的认证界。", "速率仅适用于凸和 L-smooth 的情形；PSD 失败时不适用。"]

松弛速度反例与轨迹型算法的限制

OM-MFC 可行连续状态松弛允许 λ=½δ₋₁+½δ₊₁，得到平均速度 F=0；该“ghost relaxation”不是原始单体动力学可一步实现的速度。附录 chattering 只证明可由快速切换经典轨迹一致逼近；FW LMO 若限于经典 Γ 内则所生成迭代是经典轨迹占用测度聚合，不会在有限一步直接选取此平均速度。

假设：["一维受控系统仅允许 u∈{−1,+1}；两个控制向量场在当前位置给出速度 −1 或 +1。"]

定位：{"start_line": 303, "end_line": 332, "pdf_pages": [5]}

证明骨架：["由 relaxed control 的平均速度公式计算 (−1+1)/2=0，不属于 U 的两个可取速度。", "正文以 differential inclusion 的 convex hull 解释该松弛轨迹；通过 chattering 近似的论证延后到附录 Lemma3。", "我核对了例子的输入输出；未将极限的快速切换轨迹称作同一条经典轨迹。"]

缺口：["该例解释 feasible measure relaxation 会扩大轨迹类，并非宣称整体 mean-field 弱式问题无解或必有正 relaxation gap。", "对群体概率边缘而言，多条经典路径的混合可精确产生分布变化；single-trajectory ghost 不应等同于 population-level empirical measure 的可实现性问题。"]

数值案例：二维 PDE 对照、三维障碍与卫星轨道（经验性）

作者报告 2D 单初点及分布初点情形中 FCFW 与 PDE 的群体演化定性接近；3D 十障碍实例约 20 分钟运行；卫星群由共同低轨释放点飞向高轨并沿目标轨道扩散。所有具体轨迹/时间是有限次数值实验；软障碍并非硬可行约束。

假设：["UAV 例使用单积分器、软避障势与 Gaussian 排斥核；2D 部分将轨迹法结果与栅格 PDE 密度作定性比较。", "3D 情形使用 10 个障碍且不作 PDE 对照；卫星例采用三维 Kepler 引力加有限推力控制。"]

定位：{"start_line": 611, "end_line": 895, "pdf_pages": [9, 10, 11, 12, 13, 14]}

证明骨架：["逐段读完 6.1 UAV 各初始分布/维数设置、迭代参数与 PDE 对照说明。", "逐段读完 6.2 卫星场景的 Kepler 动力学、推力界、权重和算法迭代。", "我只记录文章报告，不重跑控制序列、PDE 基准或运行时间。"]

缺口：["没有提供 finite-N mean-field approximation theorem、模型误差或数值不确定度。", "障碍通过罚函数处理，不能从图像宣称形式化 collision-free 保证。", "连续初始分布聚合是对 ρ₀ 的积分；有限粒子/轨迹样本只对离散支持精确，连续情形依赖 Monte Carlo 近似。"]

附录 A–C：测度可行性、连续方程表示与 classical LMO 证明

Lemma1 以 v(t,x)=ψ(t) 从 Liouville 式识别 µ 的时间边缘为 dt，再 disintegrate 成 dtρ_tλ_{t,x}；Lemma2 得到 continuity equation 并通过 superposition 得 path-space η；Lemma3 对状态与附加成本做 chattering 逼近；附录 C 用 path measure 推出任意 LMO 目标下界与聚合可测最优经典轨迹的匹配，证明定理5。

假设：["紧致 X,U、连续且满足文中 Lipschitz/可测要求的动力学；固定初始分布。", "松弛测度与线性化成本满足附录假设。"]

定位：{"start_line": 963, "end_line": 1055, "pdf_pages": [15, 16]}

证明骨架：["Lemma1：对只依赖 t 的测试函数 h，弱关系强制 µ 的时间投影 dt；空间控制 disintegration 得 dtρ_tλ。", "Lemma2：把 f 对 λ 平均成 F，弱关系即连续性方程；引用概率表示定理将解表示为 F-积分曲线分布。", "Lemma3：把 (x, cost) 增广；松弛轨迹速度在经典速度集合的 convex hull 内；套用 Filippov/chattering 引理得到统一逼近。", "附录 C：由任一可行 µ 的路径表示将目标写成路径成本平均；逐初态下界 V(ξ) 给整体下界；按 ρ₀ 聚合逐初态 minimizer 达到下界。", "我全读附录证明并核算逻辑次序；所引 Bogachev/Vinter/Filippov 定理未独立证明。"]

缺口：["superposition 与 chattering 结论依赖紧性、boundedness/continuity 和值函数选择假设；不可直接迁移到任意状态约束或不可测动力学。", "ρ₀ 连续支持时聚合定理是对概率分布积分，不自动产生有限条离散轨迹样本。"]

附录 D：定理6 的下降递推

从 L-smoothness、线性 oracle 最优性与范数界推出 Δₖ₊₁≤(1−αₖ)Δₖ+2L(T+1)²αₖ²，再在 αₖ=2/(k+2) 下归纳给定误差率。

假设：["定理6 的 convex/L-smooth、exact LMO 与步长假设。"]

定位：{"start_line": 1057, "end_line": 1098, "pdf_pages": [16]}

证明骨架：["将 zₖ₊₁=zₖ+αₖ(sₖ−zₖ) 代入光滑性不等式。", "用 oracle 最优性把 δJ 替换成 J*−J(zₖ)，用 ‖sₖ−zₖ‖≤2(T+1) 控制二次项。", "得到递推并按 k 归纳；我视觉核对 p16 递推系数，且手算 TV 界与 2L 常数一致。"]

缺口：["这是标准 exact-oracle FW 证明，不包含算法实现中的数值离散/控制优化误差。"]

常数：Theorem6 明示 gap≤8L(T+1)²/(k+2)，显式依赖时间跨度 T、积 total-variation norm 下的 L-smooth 常数 L 和迭代 k。总变差几何因 µ 质量 T、ν 质量 1 而给 2(T+1)。本文未给 L 关于 λ、W 或具体系统参数的显式统一表达。

端点与限制："有限时间 [0,T]；x,u 空间紧致是弱*紧性/存在性的基础。一般 positive semidefinite W 的定义针对对称状态核及所有有限 signed measures；目标外乘 λ 的符号需要非负或吸入核中的约定。线性 Liouville 可行集是 convex relaxation，可含 ghost averaged velocity；定理5 仅说 linear oracle 的值可由经典轨迹序列实现，不能推出 finite-N/共同反馈策略收敛。连续初始分布需要对初态积分，Monte Carlo 仅为近似。FW 全局 O(1/k) 只对凸且 L-smooth、精确 LMO；软障碍数值图不提供硬碰撞保障。"

## P-22e8303b48a27276

[冻结阅读卡](papers/EG-P-22e8303b48a27276.json)，输入SHA-256 `a6937bb233451e3d1840b2e7011a064a2c93d0e69068d6cd481884946e90d3e3`。

定理1：矩阵值相互作用半正定蕴含可行目标凸；公共标量核下的权重判据

可行集 Δ=Δ₁×Δ₂ 上的目标 J 凸。证明将交互项写为关于状态边缘的对称双线性二次型 Q；对两个可行点的插值，Q(ρθ)=θQ(ρ¹)+(1−θ)Q(ρ⁰)−θ(1−θ)Q(η)，其中 η=ρ¹−ρ⁰；半正定性使 Q(η)≥0。若为两人口且公共标量核 φ 正定，则足够的显式权重条件是 κ₁₁κ₂₂≥((κ₁₂+κ₂₁)/2)²。

假设：["X 紧；各 Kₚq 在 X−X 上有界且 Borel 可测。", "矩阵核 K 对每一组有限有符号 Borel 测度满足文中式 (12) 的二次型非负。", "要使用显式权重判据时，另要求所有 Wₚq 均等于同一个标量正定核 φ，且 κₚq≥0；此时仅是充分条件。"]

定位：{"start_line": 153, "end_line": 210, "pdf_pages": [3]}

证明骨架：["将原有序交互项按人口指标交换与变量替换精确改写为 K 的二次型；Kₚq(z)=K_qp(−z) 给出对称双线性型 B。", "状态边缘映射 µₐ↦(ρₐ,t) 在可行集上为仿射，因此交互项的凸组合恒等式只需检验 η。", "对公共 φ 情形，把 K 写成 φ 乘以 ½(κ+κᵀ)；后者为 2×2 半正定矩阵的行列式条件即上述充分判据。", "我视觉核对了 PDF 第3页的式 (10)–(14)，并独立检查了 2×2 行列式条件；定理证明依赖弱测度分解等文内事实，但没有独立重证所有测度论前提。"]

缺口：["一般有向核时不能只检查权重矩阵；必须验证对所有有符号测度的核二次型半正定。", "给出的公共核权重判据是充分条件，并未宣称对异质核或一般有向核必要。", "本文不提供有限代理人数逼近、松弛间隙消失或 directional 非凸情形的全局收敛结论。"]

定理2：Frank–Wolfe 线性极小化预言机按人口分离

冻结由当前点确定的 gₐ 后，整体 LMO 的目标是两个人口的线性目标之和、可行集又是直积，故整体极小化等价于分别对每个人口解 min_{(µₐ,νₐ)∈Δₐ}{⟨gₐᵏ,µₐ⟩+⟨Ψₐ,νₐ⟩}。这只保证结构分离，不等于已给出具体子问题的解析解。

假设：["当前点 zᵏ∈Δ 固定；可行集是直积 Δ₁×Δ₂。", "定理仅涉及一阶线性化问题；不要求原目标凸。"]

定位：{"start_line": 217, "end_line": 240, "pdf_pages": [4]}

证明骨架：["对交互项求一阶变分，分别得到待优化测度作为第一变量与第二变量时的两组卷积项；其和构成 gₐ。", "固定 zᵏ 后 gₐᵏ 为已知函数，线性化目标可拆成按 a 求和。", "利用 Δ=Δ₁×Δ₂，逐人口取线性极小即为全局 LMO 的极小化。", "我直接按目标与可行集直积检查了该分离推理；式 (15) 的两组有向项亦与原双线性交互目标的一阶导数一致。"]

缺口：["每个人口的连续时间控制最优问题仍可能难解；定理没有给数值 oracle 误差或离散化误差控制。", "该定理不推出非凸情形的下降或全局最优保证。"]

引理2：附加内点与可测选择假设下的经典轨迹实现

将每一初态的经典最优轨迹测度按初始分布 ρ₀ᵃ 聚合，得到解人口 a 线性 oracle 的 (µ̄ₐ,ν̄ₐ)；分别聚合两个人口即实现整体分离 oracle。作者给出 proof sketch：superposition 表示松弛可行测度为松弛轨迹混合，增广状态–成本系统的 chattering 逼近经典轨迹，内点条件确保逼近仍在 X 内。

假设：["Uₐ 紧、Ψₐ 连续，fₐ 与当前 gₐᵏ 连续，且在 X 的邻域内关于状态局部 Lipschitz，并对 (t,u) 一致。", "从 ρ₀ᵃ-a.e. 初态出发的所有可行松弛轨迹留在某紧集 Dₐ⊂int X。", "每个初态的经典控制问题存在可测选择的经典最小化轨迹。"]

定位：{"start_line": 241, "end_line": 268, "pdf_pages": [4]}

证明骨架：["先对几乎处处初态的经典控制问题选取可测最小化轨迹，再对初始分布积分生成占用与终端测度。", "由 superposition 与 chattering 说明经典轨迹 oracle 值与松弛 oracle 值相同；内点假设提供可行性余量。", "作者明确将论证结构归于参考文献 [4, App. B–C]，并调用文献 [12]、[13] 的 superposition/chattering 构造。", "我已读完作者写出的 proof sketch 与其列出的条件；未独立证明所引用 superposition 与 chattering 定理，也未验证一般系统必满足这些强假设。"]

缺口：["所需紧性、内点条件、存在性与可测选择均是显著限制，不能从占用测度约束本身推出。", "作者给的是证明提纲并依赖前作附录及教科书定理，本文未逐项重证完整测度控制实现结果。"]

定理3：精确 oracle 下的曲率界与 O(1/k) Frank–Wolfe 率

曲率常数满足 C_J≤32T maxₚ,q‖Kₚq‖∞；并且对 k≥1，J(zₖ)−J*≤2C_J/(k+2)。因此在有限时域和有界核下得到精确预言机的标准 O(1/k) 上界。

假设：["定理1与引理2条件在各次迭代均成立。", "使用 γₖ=2/(k+2) 且每步线性极小化精确；可行目标有有限曲率常数。"]

定位：{"start_line": 233, "end_line": 268, "pdf_pages": [4]}

证明骨架：["交互目标为关于状态边缘的二次型；凸组合线段上的二阶余项恰为 γ²∫₀ᵀQ(ρˢₜ−ρᶻₜ)dt。", "两可行状态边缘之差是两个概率测度之差，故每个人口 TV 范数至多 2；每个 (p,q) 积分项绝对值至多 4‖Kₚq‖∞，共四项得到 16 max‖K‖∞。", "曲率定义含因子 2，积分时域长度为 T，于是 C_J≤32T max‖K‖∞。", "剩余迭代率由标准精确 Frank–Wolfe 曲率定理推出。本人独立复核了 16 与 32 的 TV/计数常数；标准下降率本身引用 Jaggi 定理，本文未重证该通用优化结论。"]

缺口：["迭代率依赖精确 LMO；文中数值实现用时间网格、梯度下降/Adam 及有限凸组合校正，并未为这些近似误差给出理论界。", "常数对 T 与核的 sup 范数依赖为显式线性关系；若耦合核随维数或参数变大，其 sup 范数仍可能增长。", "非凸 directional 例仅报告数值现象，定理3不适用。"]

数值例：软避障与 directional 人口排序（非安全/最优性定理）

作者报告了两个满足凸性充分条件的二维交叉避障例的 fully-corrective FW 目标下降；第三个 directional 例报告 sampled ordering margin 非负且均值 0.373、障碍净距 +0.20 至 +0.64。避障仅由软罚项实现，作者明确说明数值轨迹不构成硬安全保证。

假设：["UAV 例采用单积分器、有限离散时间网格及 Gaussian 公共核；文中给出的权重满足充分凸性判据。", "3D 搜索救援例使用有向核 W₁₂(z)=φ(z)(1+εtanh(β_d dᵀz)) 与 W₂₁(z)=φ(z)(1−εtanh(β_d dᵀz))，该例不在文中认证凸性的公共偶核充分条件内。"]

定位：{"start_line": 282, "end_line": 380, "pdf_pages": [5, 6]}

证明骨架：["第一个例的权重为 [[1,.5],[.5,1]]，Gaussian 核正定且 1≥.25；第二个权重 [[1,.8],[.1,.3]] 满足 .30≥.2025，均通过式 (14)。", "directional 例取 ε=.6、β_d=1.5、κ₁₁=1、κ₂₂=.8、κ₁₂=κ₂₁=.9；作者不声称 Theorem 1 的 convexity certificate 成立。", "我视觉核对了 PDF 第5页权重矩阵与第6页方向核、报告数值及软障碍限制；未重跑数值实验，也未将图示观察视为可认证保证。"]

缺口：["报告依赖离散化与数值优化，不提供数值误差或最优性 gap 的认证。", "软障碍罚项不能排除碰撞；3D 方向排序仅对被采样时刻/数值解报告。"]

常数：显式曲率上界为 32T maxₚ,q‖Kₚq‖∞，仅显示线性依赖有限时域 T 与矩阵核分量的统一上界；FW gap 再乘以 2/(k+2)。该理论界要求每轮精确 LMO。未给有限代理数逼近误差、离散化/求解器误差或非凸情形的常数。

端点与限制："定理限于有限时域 [0,T]、紧致状态/控制域及所列连续性条件；未给 T→∞ 极限、有限代理人数的 mean-field 收敛率、一般松弛间隙结论。凸性需全测度空间的矩阵核半正定；公共标量核时的式 (14) 只是易检验充分条件。有向 directional 数值例未获凸性/收敛证明。数值避障是软罚且只记录计算轨迹观察，不是硬约束。"

## P-3b44b82c06eccd85

[冻结阅读卡](papers/EG-P-3b44b82c06eccd85.json)，输入SHA-256 `71ea62a266b6c73790b91d912fd4ca2b0f794d8496c6748e39dc7df083d2c293`。

定理1：紧半代数状态/控制域上 MSOS 界收敛

第 d 阶 SOS 次 HJB 及对应矩松弛值 J_d 单调趋于原弱 OCP 值 J*；因此极限上弱对偶无 gap。定理不称每个有限阶已精确。

假设：["扩散 OCP 数据满足 Assumption 1 的联合多项式条件；X,U 可分别表示为多项式不等式集并含足够大的球约束 R_X−‖x‖²≥0、R_U−‖u‖²≥0。", "有限多项式矩/SOS 层级按作者给定的 MSOS 松弛/限制构造。"]

定位：{"start_line": 318, "end_line": 333, "pdf_pages": [6]}

证明骨架：["作者先由球约束得 [0,T]×X×U 紧；以测试函数稠密性将弱约束从所有光滑函数缩至多项式函数。", "对测试函数 w≡1 与 w(t,x)=t，弱等式分别固定 ν 的质量为 1、ξ 的质量为 T。", "最后引用 Tacchi 2022, Corollary 8 得层级收敛；我核对了质量约束推导与作者引用位置，但没有独立证明外部一般收敛定理。"]

缺口：["作者特别指出无界 X/U 时界仍有效，但本收敛无 gap 结论需要紧致假设；引用 LPZ06 的无界状态反例说明可能保留严格 gap。", "基础多项式与 SOS 证书层级的完备性来自外部结果，并未在本文重证。"]

推论1：全局 HJB 次解由 Dynkin 公式给值函数下界

对任意 (t,z)，w(t,z)≤V(t,z)。

假设：["w∈C^{1,2}([0,T]×X)，Aw+ℓ≥0 且 w(T,·)≤ϕ；任意 admissible control 下随机过程满足文中有限矩/积分可积条件。"]

定位：{"start_line": 280, "end_line": 290, "pdf_pages": [5]}

证明骨架：["沿任意 admissible policy 用生成元与 Dynkin 公式，将 w(t,z) 写成终端期望减去 ∫Aw。", "用 Aw+ℓ≥0 及终端 w≤ϕ 得 w 不超过该 policy 的期望总成本。", "对 admissible policies 取下确界。本文给出此短证；对随机积分的可积性依赖假设/通常 Ito 理论。"]

缺口：["只给可行次解的下界；不给达到 value function 的次解存在性或最优控制构造。"]

推论2：跨时间/状态分区的分片 HJB 函数仍为价值下界

按所在时间区间和状态片拼成 w(t,x)=w_{i,k}(t,x) 后，对任意 (t,x) 它不超过价值函数 V(t,x)。随机扩散可能双向穿越空间界面，因此要求连续；当 g≡0 时，作者引述已有工作指出可用沿漂移法向的单向界面条件放宽。

假设：["分区满足 Assumption 3：覆盖且互斥；闭包和空间公共边界为基本闭半代数集。", "每片 w_{i,k} 足够光滑，片内满足 Awi,k+ℓ≥0；时间界面条件 (7)、随机扩散的空间界面连续性 (8)、终端下界条件 (9) 均成立。"]

定位：{"start_line": 391, "end_line": 418, "pdf_pages": [7]}

证明骨架：["把受控路径分解成停时区间，在每一时段/状态片内应用推论1的 Ito/Dynkin 论证。", "时间条件保证倒推经过时间节点时次解不向上跳破坏下界；空间连续保证随机路径跨片时拼接值相同。", "Appendix A 对空间切片的进入/退出停时逐段求和并令 N→∞；我阅读了附录逐步论证，未独立核验所有停时可积性。"]

缺口：["确定性单向界面放宽只在 g≡0 且动力方向确有单向跨界时可用，不能移植给含扩散的双向穿越。"]

局部 primal transport：分片弱问题由通量项平衡，连续层级比较需严格表示条件

按分片弱 Dynkin 等式和 π 界面通量构造的任何可行点可聚合为全局弱 OCP 可行点且目标相同；将全局光滑 HJB 次解限制到各片也保可行。进一步的有限阶 MSOS bound 在上述严格表示条件下至少与传统分片前方法一样紧，常可更紧。无该表示条件时，不能直接援引文中的有限阶 tightness 比较。

假设：["时间/状态分区满足前述适配条件；分片 ξ_{i,k} 为非负占用测度，ν_{i,k} 为节点测度，π_{i,j,k} 为反对称界面有符号测度。", "比较 SOS 松弛紧度时，每个闭片 X̄_k 的表示使用的多项式不等式需严格多于全局 X 的表示。"]

定位：{"start_line": 441, "end_line": 503, "pdf_pages": [8, 9]}

证明骨架：["作者用停止时刻之间的 Dynkin 公式解释局部弱等式；离开/进入状态子域产生带符号的 π 通量。", "把全部 ξ_{i,k} 求和、终端 ν_{n_T,k} 求和后，内部片界面通量反对称抵消，得全局弱可行关系且成本保持。", "全局 HJB 次解限制后片内可行；作者据更丰富的每片约束推出 SOS 证书包含关系。Appendix B 具体推导了二片且不分时间的通量公式；我已读全文及视觉核验第22页。", "按原文条件记录其 tightness 断言；没有独立证明一般 SOS 模块的表示层级嵌套。"]

缺口：["作者只在适配分割以及闭片严格多一组表示多项式的条件下作有限阶松弛比较。", "Appendix B 详细算例只展示二分片/不分时间，复杂分区推广被称 straightforward，本文未逐种拓扑展开。"]

跳过程推广：分片次 HJB 可构造，但无限状态的多项式 SDP 有额外过近似

使用跳生成元 Aw=∂_t w+Σ_i a_i(w∘h_i−w) 以及按可跳入邻域定义的接口条件，可形式化对应分片下界。作者指出闭基本半代数集若要恰好等于跳过程离散 X，只在 X 有限时成立；无限或很大的 X 为得到 tractable MSOS 松弛需过近似，而局部化能减轻、不能消除由过近似带来的保守性。

假设：["Assumption 4：状态为离散可数集，跳映射/强度、成本关于 (x,u) 联合多项式，控制集合基本闭半代数。", "跳生成元多项式封闭，以便对多项式测试函数应用 MSOS。"]

定位：{"start_line": 802, "end_line": 870, "pdf_pages": [14, 15]}

证明骨架：["从跳过程 generator 直接写出成本次解不等式和局部入邻域匹配。", "将有限集合与无限离散集合的半代数表示约束分开，说明无限状态时需对每个片区闭包以半代数集过近似。", "我逐段读了第7节；作者给出的是框架推广和可表示性限制，不是一条新的普适跳控制收敛定理。"]

缺口：["无限跳状态空间中的 tractability 依赖近似几何及多项式证书，未给统一误差界。", "数值基因调控案例是有限实验，不代表无限状态控制的全局认证。"]

数值验证：扩散与跳系统的分割收益是实验观察

图表显示分区可在相同近似阶下改善作者估计的 bound quality / computation trade-off，并缓解高阶矩病态问题；这是有限组数值例和图示报告，不是误差率或每个模型保证。

假设：["多项式随机 Lotka–Volterra 扩散模型与受控基因调控跳模型按论文设置；作者比较层级、时间/状态分区与计算成本。"]

定位：{"start_line": 543, "end_line": 653, "pdf_pages": [10, 11]}

证明骨架：["读完两节模型、代价、分割与结果描述；并逐页视觉检查相关排版图表所在页。", "记录报告的对比结论，不将作者估计的 relative optimality gap 当作具备独立误差证书的统一定理。"]

缺口：["实验求解器及耗时受硬件/实现影响；未重跑。", "状态集合、初始分布和格点/控制近似的图示不能外推到一般系统。"]

附录 A–B：停时拼接证明与局部 transport 公式

Appendix A 将一条过程路径按穿越子域的进入/退出停时分段，对每段使用 Ito 公式，求和后利用边界连续及非负漂移项得到 value-underestimator；Appendix B 对二分割构造进入/退出 measure，分解为边界通量及初末端 measure，并由停时 Dynkin 关系得到片内 transport 等式。

假设：["扩散过程满足文中 finite moments/非爆炸性假设；片内函数足够光滑；停时定义适用于被考察的分割。"]

定位：{"start_line": 1254, "end_line": 1381, "pdf_pages": [21, 22]}

证明骨架：["Appendix A 对停时之间的 Itô 公式取期望，鞅积分项用平方可积性取零，再以空间连续、时间单调界面条件拼接；对 N→∞ 并到 T。", "Appendix B 将跨边界 entry/exit measure 分解为 π 与 ν，利用反对称通量抵消并以非爆炸性识别终端局部 measure。", "视觉复核 PDF p21–22；Appendix B 原文只详细推导二片情况。"]

缺口：["停时鞅项的平方可积性和非爆炸条件依赖 Assumption 2 与过程正则性；通用过程理论未独立重证。", "转换文本有公式排版错位，关键停时关系已通过视觉查看原 PDF p22 核对。"]

常数：本文没有给新的显式 sharp 常数。MSOS 层级收敛是渐近结论而非有限阶误差公式；论文的结构表报告变量和 LMI 数关于分区数 n_T、n_X 的 O(n_T)、O(n_X) 线性尺度且 LMI 维度 O(1)，阶数 d 对 moment/SOS 规模带来组合增长。统计/实验耗时不是统一复杂度界。

端点与限制："扩散 value-underestimator 需要片内 C^{1,2} 与可积 Dynkin/Ito 条件；传统 MSOS 渐近无 gap 需紧状态/控制集合的球约束表示，弱测度 LP 本身可能严格松弛原轨迹控制问题。随机扩散跨空间界面要求连续性；仅确定性 g≡0 可讨论单向法向条件。有限阶局部/全局 MSOS 紧度比较附带“闭片表示多于全局表示”的条件。跳过程无限离散态需半代数过近似，额外保守性无法由文章普遍消除。"

## P-610d8ca4e35b8a23

[冻结阅读卡](papers/EG-P-610d8ca4e35b8a23.json)，输入SHA-256 `c8b5421e72f5efb9e7ad3b4bc332bb1d96cdf5c16b50d63a21fd939f0aa4dfcd`。

命题3.1/定理4.5：VaR tail-bound chance-peak 的凸 measure SOCP 上界

逐时刻 VaR 由 Cantelli/VP moment bound 支配后，chance-peak objective 可上界为 sup over stopping distributions of ⟨p,µτ⟩+r√(⟨p²,µτ⟩−⟨p,µτ⟩²)。替换 SOC 变量后得凸 measure SOCP (24)，其有限维 SOC 是统计矩之间的方差锥约束。

假设：["以固定初始分布 µ₀ 和 Markov generator L 表示从 0 到停止/指定时间的过程；X 与 [0,T] 紧、边界首次接触时停、p 连续。", "Cantelli 可任取其有效 r；VP 额外要求输出分布单峰及 ε≤1/6。"]

定位：{"start_line": 380, "end_line": 514, "pdf_pages": [7, 8, 9]}

证明骨架：["由 tail inequality 把 VaR 目标替换成均值加标准差 r 倍，随后用 stopping-measure Dynkin 关系扩展到 measure program。", "引入 c≤√(b−a²)，其最大化目标促使最优 c 取根号；把 c²+a²≤b 用标准 Lorentz cone 恒等式改写。", "PDF p8 视觉核对了 (21)–(24) 与 SOC 坐标，手工展开确认 (1+b)^2−(1−b)^2=4b。", "Problem 3.1/4.1 明确将 t* 作为 [0,T] 中的确定性标量终止时刻；固定时刻当然也属于停止时。Theorem 4.2 取 τ*=t*∧τ_X，并声称这一标量时刻轨迹家族与 measure feasible set 一一对应；上界构造只需正向映射，该 one-to-one/任意可行流的反向覆盖未由所读证明建立（见 gaps）。"]

缺口：["VP 只对单峰分布及 ε≤1/6 有效；一般状态分布应退回 Cantelli。", "Problem 3.1/4.1 的 sup_{t*∈[0,T]} 明确是确定性标量终止时间，不是该记号含糊。Theorem 4.2（PDF p8，lines439–448）从固定时刻诱导流给出上界，并声称 one-to-one；Lemma4.3（p8，lines450–456）把 every feasible (µτ,µ) 说成 source trajectory 并外引[20]，未证明其仍属于同一个固定标量时刻轨迹家族。缺口限此反向覆盖/对应证明；不把作者用 stopping-time 字样当作已证允许一般随机化或适应停止，也不把后一解释当唯一修复。"]

定理4.7及附录B：tail SOCP 的函数对偶与强对偶

对偶在连续 v 与有限维 u 上最小化，满足 Lv≤0、v+u₁p²−2u₂p≥p、([u₁+u₃,−r/2,u₂],u₃)∈L₃；作者声明该无限维 measure SOCP 与函数/SOC dual 在 A1–A4 下强对偶。

假设：["弱对偶要求 A1–A3；作者的强对偶声明另用 A4（初分布概率质量）、紧性与连续 p。", "Appendix B 以 SOC 的 LMI 表示、可行轨迹、µ/µτ 有界和 SOC 矩阵 Φ² 的 trace 有界验证抽象测度 SDP 强对偶定理条件。"]

定位：{"start_line": 528, "end_line": 541, "pdf_pages": [9]}

证明骨架：["正文由 (25) 的 measure equalities 配对 v,u 并写拉格朗日函数得到 (26)。", "Appendix A 建立含 measure 与有限 PSD/SOC block 的抽象 cone strong-duality 定理；Appendix B 将三维 SOC 嵌入 PSD 矩阵并验证可行、measure 有界与平方迹有界。", "我独立展开了作者式 (60) 的 trace：6+6⟨p⟩²+6⟨p²⟩+6⟨p²⟩²。文中 (61) 的界用 Π₁=max p 替代 |⟨p⟩|；若 p 的最大值小而终值测度集中在负值区，这不能控制均值平方。有效但较粗的替代是 ⟨p⟩²≤⟨p²⟩≤Π₂，因此 B²可取 6(1+2Π₂+Π₂²)。故结论存在可修复证明缺口，不能照录 (61) 作为已核通过。"]

缺口：["附录B 的 (61) 有上述符号/绝对值界问题；强对偶结论可用改正后的有限界补救，但本文给出的具体证明界不成立于所有 A1–A4 数据。", "函数空间和抽象 duality 结果仍依赖 [49] 一般凸对偶定理；我核对推演结构，未重证外部定理。"]

定理5.3–5.6：Expected Shortfall domination LP 与停止问题

通过 ε ν+ν̂=p#µτ，把标准fractional ES写成对ν的线性目标；measure LP (29)上界标量时刻chance-peak ES（Theorem5.4）。作者Theorem5.5声明A1–A7下(29)与其Problem5.1/(27)无gap，但该字面固定标量t*定义下声明为假：下面checked-local双分支ODE给固定时刻最优2/3、measure LP最优1。强对偶是该measure LP与其函数dual的不同主张，保外引依赖，不能修此模型gap。

假设：["0<ε≤1；µτ 的状态边缘为概率测度且 p 连续；A1–A7 确保作者援引的 occupation measure stopping 表示。", "注意原子分位点处必须采用允许分割 quantile atom 的标准 CVaR/ES 约定，不能照用事件指示式 (3) 的整粒子尾质量。"]

定位：{"start_line": 582, "end_line": 642, "pdf_pages": [10, 11]}

证明骨架：["从 ν≪ψ 且 dν/dψ≤1/ε 构造 slack measure ψ−εν，得到 domination equality。", "把 stopped process 构造为 µ、µτ，并取 pushforward ψ=p#µτ，得到上界方向。", "无 gap 证明将每个 Dynkin 可行测度表示成过程轨迹，再引用 ES measure representation；证明依赖文献 [20] 的表示定理。", "独立端点反例：V=1 (概率1/2), V=0 (概率1/2), ε=1/4。按正文式(3)，VaR=1 且 ε⁻¹E[1_{V≥VaR}V]=2；按式(4)最小化 λ+4E[(V−λ)_+] 得 1，按式(5)/(28) domination LP 最优也是1（取 ν=δ₁）。因此式(3)含整颗分位原子而与 LP/标准 CVaR 不等价；需分割原子后才可保持后两式。", "2026-10-06 checked-local source counterexample（由audit_ef_16与audit_fg_16独立核约束）：X=[0,2]∪[3,5],T=1，平滑cutoff drift在两初态轨迹上为1，初测½δ_.5+½δ_3.5；连续p为.75与4.25附近两个.1半径bump。固定t时正reward两窗不交，标准ES_.75≤2/3且可达；初态可测停时.25/.75诱导正flow和terminal lawδ1，ν=δ1、νhat=.25δ1使(29)目标1。完整模型、每个A1–A7及全C¹ chain-rule flow在source_counterexamples单列，不写成作者定理。"]

缺口：["打印的 Definition 2.1 式(3) 与 Lemma 2.1 / measure LP (5) 在 quantile atoms 上不相容。作者在 Remark 3 承认 atom splitting，但未修正式(3)本身；后续 ES theorem 适用于 standard fractional-tail ES。", "Problem 5.1/式(27)（PDF p10，lines566–574）的 t* 明确为确定性标量。Theorem5.4（p10，lines609–618）固定时刻也可称 stopping time，正向构造给 upper bound；Theorem5.5（p11，lines623–626）从 every feasible flow 引 Lemma4.3 便声明 no gap，没有证明任意 feasible flow 的 terminal law 由一个共同固定 t* 产生。缺口是此 measure-flow 与固定终端轨迹家族的反向覆盖，不是 t* 记号歧义；[20]未独立全文核。", "强对偶依赖 A1–A4 和 [50, Thm 2.6]；附录C核对了三个条件链，未独立证明引用定理。", "Theorem5.5 no-gap在显示的Problem5.1确定性标量t*语义下已被上面明确可行反例推翻；(29)无共同t* hyperplane或terminal-time marginal constraint。Theorem5.4固定时刻→measure LP的upper-bound方向仍可用。把原问题改为一般适应停止是一种更改问题的迁移，不是原(27)记号的自动解释。"]

定理6.3/6.7：多项式支撑上的 moment-SOS 上界收敛

Cantelli/VP moment-LMI 上界随 d 趋近 tail measure SOCP 最优值；ES LMI 上界单调递减趋于 ES measure LP 的 chance-peak 值。故在前述 assumptions 全满足时，数值 SDP 序列渐近收敛到无限维凸问题值。

假设：["A1–A7 的 generator/stop-process regularity；A8 要求 X₀、X 基本半代数并含球约束；A9 要求 generator 将时空多项式映为多项式；A10 要求 p 多项式。", "相应 measure 解有界；ES 的 ν、ν̂ 支撑在 p(X) 的紧区间。"]

定位：{"start_line": 700, "end_line": 848, "pdf_pages": [12, 13, 14]}

证明骨架：["把生成元作用于单项式测试函数，Dynkin 约束成为仿矩线性方程；对状态/时间支持加 localizing PSD 矩阵，对方差项保留 SOC。", "Theorem 6.3 先用 µ、µτ 紧支撑和测试函数 1,t 控质量，继而引用 Tacchi [43] 强对偶/收敛结果。", "Theorem 6.5 对 ν,ν̂ 用 ⟨1,ν⟩=1、⟨1,ν̂⟩=1−ε 及 p(X) 紧支撑证明有界；Theorem 6.7 再引用 Archimedean Lasserre convergence 推出 ES LMI 单调收敛。", "作者 Theorem6.6 用最坏矩分布下的 Cantelli 界推出 tail/Cantelli 上界不低于标准 ES；我只核对表达式及假设角色，不独证一般 Lasserre 结果。"]

缺口：["收敛为 d→∞ 渐近结论；文章不给任意有限 d 的精度率。", "紧致性或 ball/Archimedean 条件缺失时，上界仍有效但层级未必收敛至 measure LP（Remark 6）。", "Theorem6.3 对非凸的 VP objective? 其变分 SOCP对给定 r仍是凸形式，须保持 VP 适用的分布单峰假设。"]

扩展：切换系统的模式占用测度、概率安全距离

ES/SOC stopping模型可将 generator balance 改成 µτ=δ₀⊗µ₀+ΣL_ℓ†µ_ℓ；针对 unsafe set 的距离目标可用扩展状态和 pushforward ES measure 改写，但变量/PSD 尺寸增加且 mixed powers 会阻碍 correlative sparsity。

假设：["有限个 generator L_ℓ；总占用 measure 可分成非负 µ_ℓ，不加 dwell-time 约束。", "距离扩展的点集函数通常不 polynomial，矩-SOS实现需增加辅助状态或采用联合 measure/边缘相等约束。"]

定位：{"start_line": 975, "end_line": 1085, "pdf_pages": [16, 17, 18]}

证明骨架：["逐式读完切换模式分解 (45)–(46) 与距离联合分布/边缘匹配 (47)。", "作者以 generator decomposition 得切换关系；距离 LP 的边缘约束强制 η 和峰值测度具有相同状态边缘。", "该节属于模型扩展表达，没有另给统一 gap 或距离精度定理。"]

缺口：["距离点集函数一般不多项式；所述维数和计算复杂度增加，具体效果取决于 unsafe set 几何。", "无 dwell-time 的 switched occupancy relaxation允许任意快切换，未给有限切换速率/离散实现误差界。"]

附录 A–C：一般 measure-SDP 强对偶与应用

Appendix A 将 measure 与有限 PSD/SOC 块纳入同一锥对偶框架，借 compact base 论证强对偶及最优可达；Appendix B 将其应用到 tail SOC，Appendix C 验证 ES LP 的有界性、可行性和连续性条件。

假设：["可行 measure/matrix 变量有界，至少一可行点；抽象 cone 有合适紧截面。ES 证明另用紧支撑、有界 measure mass 和连续数据。"]

定位：{"start_line": 1450, "end_line": 1752, "pdf_pages": [25, 26, 27, 28, 29]}

证明骨架：["读完 Theorem A.2：由可行面上的二次矩/质量统一界验证 A2′，取 ψ=(I,1,…,1) 建立正锥紧截面，再引用 Barvinok 凸对偶定理。", "读完 Appendix B SOC-LMI embedding (55)–(61) 并视觉核对 pp.28：得到 trace 展开及 Π₁ 界错误；用概率均值二阶矩不等式给出上文修正有限界。", "读完 Appendix C pp.28–29：根据 Theorem6.5 的质量界、轨迹诱导可行 measure 及连续 generator/p 检验强对偶条件。", "上述一般拓扑/锥对偶结果依赖 [49] 与 [50]，未独立重证。"]

缺口：["附录B的定量 trace 上界如前所述需要修正；不能按原式(61)复用。", "Theorem A.2 的闭锥截面/弱*紧性论证依赖 compact support 与质量/矩约束；离开论文设定不可直接套用。"]

数值案例：连续、离散、切换系统与安全距离

五类演示包括二维 SDE、三维 stochastic twist、离散参数随机映射、switched SDE、unsafe set 距离。作者报告 Cantelli/VP/ES SDP bound 与 MC estimate；对于概率距离实验，VP 下界很保守，ES 有改进。数表仅是所列系统与近似阶的实验。

假设：["结果来自作者 MATLAB/YALMIP/MOSEK 实现，Monte Carlo 用 50,000 条路径及 Δt=10⁻³ 估 VaR/ES。"]

定位：{"start_line": 1085, "end_line": 1398, "pdf_pages": [19, 20, 21, 22, 23, 24]}

证明骨架：["读完所有系统动力学、支持域、p、时间范围、表1–17与图说明，并视觉核对主要定义和优化公式页。", "实验数值按作者报告，未重新执行求解器/50,000条 Monte Carlo。"]

缺口：["MC 样本和步长不是概率置信证书；表中 SDP 为上/下界但具体有限阶数值精度不可由对比采样推出。", "概率距离章节的“负下界截为0”是报告后处理，不是原始 LMI bound 的性质。"]

常数：Tail statistic 的 r_C=√(1/ε−1)，VP 的 r_VP=√(4/(9ε)−1)，故显式随尾概率 ε 变差；无新尖锐常数。SOCP strong-duality Appendix B 以 Π₁=max_X p、Π₂=max_X p² 给出 trace 上界 (61)，但 Π₁ 有符号且不能普遍控制 |E[p]|，该显式界有缺口；可用 |E[p]|²≤E[p²]≤Π₂ 修为 6(1+2Π₂+Π₂²)。moment-SOS 结果只有随 d→∞ 的收敛陈述，复杂度依赖变量数/多项式度数/矩阵阶数。

端点与限制："VaR 取 ε∈[0,1)；ES 公式应限定标准分数尾 CVaR，且 ε>0，VaR 原子处需要按质量拆分。VP 要求 unimodality、ε≤1/6。stop-measure 表示假设 compact [0,T]×X、过程在首次触边时停、generator domain 符合 A5–A7；Problem3.1/5.1 的 t* 是确定性标量；Theorem4.2 的 one-to-one 与 Lemma4.3/Theorem5.5 的任意 feasible flow反向覆盖到该固定终止轨迹家族不能由(29)保证；本卡checked-local模型给Theorem5.5 literal no-gap反例，upper-bound不受损。moment-SOS 渐近界另需 ball-representable semialgebraic sets 与 polynomial generator/objective。数值例不提供复杂度/精度统一界。"

## P-b8164dc1ce4dcf5d

[冻结阅读卡](papers/EG-P-b8164dc1ce4dcf5d.json)，输入SHA-256 `295c7964b0cc57b5590434f4acff2062246bcec0a22e634760fd33707beb502c`。

定理3.1与推论3.2：凸数据下仿射与完整占用测度松弛均无间隙

仿射测试函数测度松弛值 M_aff=M。由于 M_r 位于 M_aff 与原值 M 之间，推论3.2进而给出 M_r=M_aff=M；即在这些联合凸性假设下，无论使用 affine 还是文中完整的非线性弱 Liouville 测试，均无松弛间隙。

假设：["Ω 为有界连通开集，边界分片 C¹；p∈[1,+∞] 且占用测度满足相应矩有限/紧支撑条件。", "Y、Z 闭且凸。对每个 x，体内 L、Bᵢ、Cᵢ 关于 (y,z) 凸，Aᵢ 关于 (y,z) 仿射；边界 A∂,ᵢ 关于 y 仿射、B∂,ᵢ 与 C∂,ᵢ 关于 y 凸。允许 I 为不可数索引集，含点态与积分约束。"]

定位：{"start_line": 261, "end_line": 382, "pdf_pages": [5, 6, 7]}

证明骨架：["任一原可行 Sobolev 函数诱导 (µ,µ∂)，故 M_aff≤M。", "由连通域常值引理识别 µ 的 x 边际为 dx，再 disintegrate 并对每个 x 取条件均值 y(x)、z(x)。", "把仿射测试函数代入 Liouville 式，推出 y∈W¹,p 且 Dy=z；边界条件显示 µ∂ 的条件均值是 Sobolev trace y∂。", "仿射等式约束在条件均值下保留，凸不等式与积分约束由 Jensen 保留，闭凸 Y 也包含均值；目标成本由 Jensen 不增。因此由任一 affine 可行测度构造原问题的可行函数，得到 M≤M_aff。", "我按原文检查了两边值比较方向、条件均值与 Jensen 各自承担的步骤；Sobolev trace/Stokes、Radon–Nikodym/disintegration 及可测性定理均是文中引用或标准结果，未独立重证。"]

缺口：["无间隙依赖体内和边界的联合凸性/仿射条件；仅有 L 在梯度变量上凸一般不满足此定理。", "完整测度类较 affine 测度类更小，但从 affine 等于原值和 M_aff≤M_r≤M 推出同值；非凸时不能沿用此夹逼结论。", "该定理给值相等，不提供原可行函数的唯一性、测度解由单条轨迹生成，或数值 SDP 的有限阶精确性。"]

定理4.2：无积分约束时 affine 松弛等于逐点凸包络问题

定义 Mhat 为在 conv K(x) 及边界 conv K∂(x) 上、以 L 在 K(x) 上的逐点凸包络 Lhat 为成本的凸化变分问题，则 Mhat=M_aff。故无积分约束时有“凸包络值 = affine occupation-measure 值 ≤ 完整 occupation-measure 值 ≤ 原始值”。

假设：["本节取 p=∞、Y 与 Z 紧；点态可行纤维 K(x) 闭，且 x↦K(x) 关于 Hausdorff 度量连续。", "L 对每个 x 关于 (y,z) 连续；除特别展示反例的例子外，本节主定理排除积分约束。", "边界约束通过 K∂ 的相应闭凸包络处理。"]

定位：{"start_line": 435, "end_line": 529, "pdf_pages": [8, 9]}

证明骨架：["从凸包络可行函数出发，逐点用 Choquet 表示其 (y,Dy) 为 K(x) 上概率测度的重心，并使成本等于 Lhat；通过紧弱拓扑、Prokhorov 和 Kuratowski–Ryll-Nardzewski 选取可测核，生成 affine 可行测度，得 M_aff≤Mhat。", "反向从 affine 可行测度取条件重心，沿定理3.1论证得 y∈W¹,∞、Dy=z 且边界均值为 trace；重心位于 convK，凸包络的 Jensen 不等式给 Mhat≤M_aff。", "我核对了两方向构造及值不等式方向；Choquet 与可测选择定理本身依赖被引的外部结果，未独立重证。"]

缺口：["Hausdorff 连续纤维与紧性是构造可测 Choquet 核的关键；不能仅由逐点凸包定义自动省略。", "对一般积分约束，该等价可失败；作者的例4.3/4.4分别给非凸 Y、非凸 Z 的例子，其中 Mhat=−1 而 M_aff=0。", "Mhat=M_aff 不意味着完整非线性测试松弛必等于 convex envelope；M_r 可严格更强。"]

例4.3–4.4：积分约束下凸包络松弛可严格弱于 affine OM

两个例子都算出凸包络变分值 Mhat=−1，而 affine occupation-measure 值 M_aff=0，故有严格差距。例4.3的 M_aff 下界由 L=−C+‖z‖² 且 ∫C≤0 推得非负，零由 z=0 的可行均匀测度达到；例4.4同理由 L=−C，零值由 y=0 的可行均匀测度达到。

假设：["例4.3取非凸 Y={−1,1}²、凸速度集 Z 为单位圆盘；L=y₁y₂+‖z‖²，C=−y₁y₂≤0 为积分约束。", "例4.4取 Y 为单位圆盘、非凸 Z={−1,1}²；L=z₁z₂，C=−z₁z₂≤0 为积分约束。"]

定位：{"start_line": 530, "end_line": 601, "pdf_pages": [9, 10]}

证明骨架：["例4.3中 Y 四角点上 y₁y₂=±1，conv envelope 形成 Lhat=|y₁+y₂|−1+‖z‖² 与 Chat=|y₂−y₁|−1；凸包络问题沿连接两反对角角点的常值曲线达到 −1。", "原 affine 测度问题因 ∫C≤0 而有 ∫L=−∫C+∫‖z‖²≥0，取 z=0 且平衡正负 y₁y₂ 得到 0。", "例4.4交换非凸性到速度集 Z，得到 Lhat=|z₁+z₂|−1、Chat=|z₂−z₁|−1；凸化问题用导数沿反对角线段的可行曲线取 −1。", "affine measure 问题再次由 L=−C 与积分约束得到非负值，取 y=0 的均匀测度达到 0。PDF 第9–10页视觉核对了两个构造与值方向。"]

缺口：["这两例仅证明在给定数据下凸包络值严格低于 affine OM；不推出所有带积分约束的问题均有严格差距。", "例子中的凸包络值、可行曲线和测度显式构造是专用反例，不能抽象替代定理4.2的 Hausdorff 连续性假设。"]

推论4.1：指定凸下估计量形成占用测度下界链

M_under≤M_aff≤M_r≤M；若某凸下估计问题恰好 M_under=M，则 affine 与完整占用测度松弛均与原值相等。

假设：["体内与边界等式约束仿射，使它们无需进一步凸化。", "选取的成本/不等式/积分约束下估计量在原可行点上不超过原数据；凸化问题按式(14)定义。"]

定位：{"start_line": 390, "end_line": 420, "pdf_pages": [7]}

证明骨架：["对凸化下估计问题应用推论3.2，得到其占用测度松弛值等于 M_under。", "原 affine 可行测度对凸化问题仍可行，且下估计成本不高于原成本，因此凸化值不超过 M_aff。", "结合始终成立的 M_aff≤M_r≤M 得结论。", "我复核了松弛集合包含与成本方向；该推论不会自动构造达到 M 的某个下估计函数。"]

缺口：["需验证候选下估计量对所有原可行点成立；任意凸函数替换不是自动合法的下界。", "若等式约束非仿射，作者指出需进一步构造其凸化，而文中推论未直接覆盖。"]

定理5.2：带控制的 variational occupation-measure formulation 无间隙

在 OC1–OC3 下，控制问题 M_oc=M_r^oc。作者通过平均掉占用测度中的状态、导数、控制后，再按均值状态选择逐点最优控制；可测选择与凸性确保得到原问题可行解，Jensen 保证成本不增。

假设：["采用第5节的 p 阶矩/紧支撑设定，原始状态 y∈W¹,p、控制可测，且 µ 的 Liouville 等式只对 y 仿射测试函数。", "OC1：点态控制可行纤维非空凸，积分约束密度凸；边界对应条件也成立。", "OC2：给定 (x,y,z) 后控制成本最小值可达且 argmin 非空闭；边界相同。", "OC3：约化成本 L̄、L̄∂ 对 (y,z)/y 凸且局部有界。"]

定位：{"start_line": 610, "end_line": 828, "pdf_pages": [11, 12, 13, 14]}

证明骨架：["原函数与可测控制诱导 occupation measures，故 relaxation 值不高于原问题值。", "反向对 relaxation 测度按 x disintegrate，取 y,z 和边界 y 的条件均值；引用定理3.1相同弱式论证给出 Sobolev 函数与其 trace。", "OC1 与 Jensen 保留等式、不等式、积分及边界约束；OC2 加 measurable graph 假设，使逐点 argmin 可用 measurable selection theorem 选出控制。", "OC3 使约化成本凸，Jensen 比较平均控制状态成本与原测度成本，得到 M_oc≤M_r^oc。", "我逐项读了可测选择、约束保留和成本链；selection 与 trace 定理是外部引用，未独立证明所有一般 Polish 控制空间下技术条件。"]

缺口：["OC1–OC3 是真实的结构假设；一般非凸可行控制纤维、控制极小值不达或非凸约化成本均不由该定理覆盖。", "该控制推广仍是变分域 Ω 上的弱占用测度结果，不是任意动态 mean-field/game 系统的通用无间隙定理。"]

命题6.1：二维 micromagnetics 的 affine 占用测度形式复现原问题下确界

式(31) 在 O_aff² 上的线性 occupation-measure 问题，其最小值等于原 micromagnetics 问题(24)的下确界。作者利用凸包络问题与 Young-measure 凸化值等于原 infimum 的既有结论，再套用定理4.2。

假设：["二维刚性铁磁模型：磁化 m∈L²(D;R²)、|m|=1 a.e.；磁静场满足 div(−∇u_m+mχ_D)=0 的弱式。", "作者把 y=(m₁,m₂,u_m)、z=(∇m₁,∇m₂,∇u_m) 作为变分状态/梯度，单位圆条件作为支持约束，磁静场弱式用一族线性积分约束表达。", "目标中对圆周的各向异性能密度不必凸；使用 affine occupation relaxation 和定理4.2，另引用先前 Young-measure/convex-envelope 结果。"]

定位：{"start_line": 936, "end_line": 980, "pdf_pages": [16]}

证明骨架：["将 m 与磁势并成 y，梯度放在 z 中；用支持集编码 |m|=1，用 C_φ 表示 div 约束的测试积分。", "该圆周约束不凸，所以一般不能直接套定理3.1；作者改用 affine 松弛，由定理4.2与问题(27)的凸 envelope 相连。", "式(26)的 Young measure/凸包络等价来自引用的 [6, Thm.3.4, Remark 3.7]，本文没有重证；据此把式(31)值等同于原始 infimum。"]

缺口：["结论是 infimum 的等价，不保证单位圆约束下原始问题取得 minimizer；文中指出序列可振荡而弱极限落在单位球内。", "式(26)中凸包络极限等价依赖外部 micromagnetics 文献；本文未独立证明。", "所谓可计算依赖后续 moment-SOS 阶层近似；本文没有给具体数值实验、截断阶误差或显式复杂度界。"]

附录 Lemma A.1：连通域上满足常向量通量弱式的质量边际是 Lebesgue

则 µ 等于 Ω 上的 n 维 Lebesgue 测度；此引理用于由 occupation-measure Liouville 式的质量约束识别 x 边际。

假设：["Ω 为有界连通开集；µ 是正且紧支撑的 Radon 测度，并且 µ(Ω)=|Ω|。", "对所有紧支撑光滑测试函数，∫Ω∇φ dµ=0。"]

定位：{"start_line": 988, "end_line": 1033, "pdf_pages": [17]}

证明骨架：["先取 Ω 内小平行多面体 R 及可行平移 τ(R)，沿坐标轴分解平移。", "构造支撑在 R 与其平移凸包内的积分测试函数并用光滑逼近，将弱式变成 µ(R)=µ(τ(R))。", "所有允许小平移下测度一致，调用文献 [24, Thm.2.20] 得 µ=c dx；质量条件给 c=1。", "我读了构造和外部结论用途，未独立证明“所有小矩形平移同质量推出常数倍 Lebesgue”的测度唯一性定理。"]

缺口：["连通性确保可由局部平移连接域内位置；非连通情形正文 Remark 2.6 要求直接规定 x 边际。", "证明依赖将指标函数测试以光滑函数逼近及 Rudin 的测度识别结果。"]

附录 Lemma B.1：增长受控的弱 Liouville 测试函数积分可积

边界积分 φ(x,y)n(x) 对 µ∂ 可积，体内项 ∂ₓφ+zᵀ∂ᵧφ 对 µ 可积，因此 Definition 2.7 的完整弱 Liouville 等式有意义。p=∞ 时由紧支撑及连续性直接得到。

假设：["p∈[1,∞)；µ 与 µ∂ 满足 Def.2.4 中 y,z 的 p 阶矩条件。", "测试函数 φ∈F_p，使 ∂ₓφ、∂ᵧφ、φ 满足文中对 |y|^p 与 |y|^{p−1} 的增长界。"]

定位：{"start_line": 1035, "end_line": 1090, "pdf_pages": [17, 18]}

证明骨架：["p=∞ 由紧支撑与连续 integrand。", "p<∞ 时对边界项用 y 的 p 阶矩控制。", "体内 z·∂ᵧφ 采用 Hölder，配合 z 的 p 阶矩与 ∂ᵧφ 的 |y|^{p−1} 增长，得到有限积分。", "我阅读并检查了所用矩条件与指数互补关系；没有逐行独立证明每个泛函分析可积性细节。"]

缺口：["引理只保证特定增长类测试函数的积分定义良好，不提供测试类稠密性或松弛等价本身。"]

常数：无显式收敛率或误差常数；核心是精确值等式 M_aff=M 或 M_hat=M_aff，分别在凸性或 compact/Hausdorff 连续纤维条件下成立。文章介绍 moment-SOS 层级作为有限维近似背景，但未在本文给可计算截断阶数对应的误差常数。

端点与限制："区分 p∈[1,∞) 的矩有限测试/测度类与 p=∞ 的紧支撑版本。非凸主定理4.2限于紧 Y,Z、连续紧值点态纤维、Hausdorff 连续且无积分约束；含积分约束的例4.3/4.4显示 convex-envelope 问题可严格低于 affine OM。凸定理3.1需关于 (y,z) 联合凸，不足以只假设对梯度 z 凸。带控制定理另需 OC1–OC3。没有有限 SOS 阶误差或动态 mean-field 模型定理。"


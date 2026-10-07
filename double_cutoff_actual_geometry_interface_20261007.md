# 双计数截断后的实际来源几何容量接口

2026-10-07。本稿给出任意输入上的一个有偿来源几何分层：来源固定的高质量核心，由最粗标记祖先清算；核心外来源合成一份固定小质量测度，用原移动硬尺度列支付。它支付明确的当前实际子交通，未证明所有实际余项都属于该子支。所有原 FIRST、共同完整 future cap、两赢家、出生、严格 CP/GP、唯一 LCA、重捕获、原容量币、共享种子和历史门保持；诊断 Bernstein 掩码不等于 early continuation 活跃坐标。

## 1. 原合同与查重

沿用 `actual_parent_occupancy_route_20261007.md` 的原来源及父组记号。固定原带、共同 q>0 和森林种子 θ，令

\[
 M_P=\mu_{\rm hi}(P),\qquad h_P=\mu_{\rm hi}(P\cap B),\quad
 W_{\rm hi}=\|\mu_{\rm hi}\|\le W.
\]

B 是原来源固定出生资格集合。原完整软行、硬行分别满足 S(x)≤q、H(x)≤C_hq，C_h≥1。当前实际正交通先合并唯一输出带、唯一 LCA 和所有互斥标签，仍由 q⁻¹ S(x)H(x) 支配，密度≤C_hq。逐原来源核删门时只作正支配；不反推被放大的对象继承原门。

当前 R** 还保留 lowS、lowcoin、辅助 K<min(κ_C,k_ε) 和相对于原 L_s 的 outside 计数 O_s<k_ε。K=0 子核就是 (1−σ)ⁿh_{L_s}，而 K>0 不具有该硬支撑。严格 CP/GP 为

\[
 v_jh_P>q\min(L_s,R_h)^n\sqrt n,\quad 0<v_j\le1,
 \qquad h_P>q(L_sR_h)^{n/2}. \tag{1}
\]

原 tex 2960–3015、4766–4840 给这些门；2608–3015 的 CA/CP、`high_marked_receipts/large_parent_capture_route.md` 已检查过捕获份额乘回后只恢复行界。旧共享覆盖 tex 4100–4135 与旧最粗祖先反链可付 q|U|；父盒加固定 b 的密度判据已在 `actual_parent_occupancy_route_20261007.md` §5 记为旧工具应用。本稿新分层用来源固定 CP 半径和可丢少量来源的核心，不重新领取历史整账。

## 2. 原 CP 的来源固定半径及可覆盖两序

只标记 h_P>qaⁿ√n 的父组；实际 (1) 的父组均相关。定义

\[
 \widetilde\rho_P=(h_P/(q\sqrt n))^{1/n},\qquad
 \rho_P=\min(b,\widetilde\rho_P). \tag{2}
\]

于是 a≤ρ_P≤b（a<b 时 a<ρ_P），qρ_Pⁿ≤h_P/√n，且 P⊂A 蕴含 ρ_P≤ρ_A。由 (1) 只有

\[
 \min(L_s,R_h)<\widetilde\rho_P,\qquad
 \min(L_s,R_h)\le\rho_P. \tag{3}
\]

ρ_P<b 时才有第二个严格不等式；当根半径>b且实际 min=b 时等号可发生。下文一律用闭方体 Q(0,ρ)，保留原闭面规范。

令互斥的同一输出选择域为

\[
 D_h=\{R_h\le L_s\},\qquad D_s=\{L_s<R_h\}. \tag{4}
\]

本节可覆盖的当前实际子交通为 D_h 上全部原 K 标签，加 D_s 上 K=0 标签。在 D_h 上取原 hard 来源端点 z；在 D_s,K=0 上取原 soft 来源端点 y。分别有 x−z∈Q(0,ρ_P) 和 x−y∈Q(0,ρ_P)。等尺度只归 D_h。不用低 K 断言其它掩码也有硬支撑。

## 3. 可执行的任意来源核心容量

固定 0<ζ<1、T>0。令 𝒦 是所有全空间有限有理闭轴盒并集的可数枚举；不要求 K⊂半开来源 P。标记相关 P，当某 K∈𝒦 满足

\[
 \mu_{\rm hi}(P\setminus K)\le\zeta M_P,
 \qquad q|K+Q(0,\rho_P)|\le T h_P. \tag{5}
\]

从枚举取第一个满足者 K_P。来源质量零者不标记。也可用固定细分格闭盒的有限并；必须有任意精细的有理尺度。

这是完整输入上的来源固定判据，不使用当前捕获 μ(Q(x,R_h))，不按 x 重分配来源。有限盒并的 Minkowski 并体积是有限矩形并的可测函数；质量及 ρ_P 随原种子可测，故每个候选判据和第一个满足者可测。原森林有限深度、可数节点保证祖先选择可测。不是把所有候选核心的体积或所有种子的覆盖免费并在一起。

对于任意紧核心 K，若 μ(P\K)≤ζM_P 且有严格体积余量 q|K+Qρ_P|<Th_P，则闭有理外覆盖 K_j↓K 可取有限盒并。质量损失不增加，K_j+Qρ_P↓K+Qρ_P，有限体积连续从上使某 j 满足 (5)。薄原子面保留在闭外覆盖中。纯等号只给 T+η 的近似判据，不能声称必有同 T 的有限枚举证书。

## 4. 最粗祖先清算：一个真正的空间费用

取所有最粗标记祖先 A。它们在原来源森林中两两不交，Σ_AM_A≤W_hi，Σ_Ah_A≤W_hi。若原 LCA P 的祖先链含标记组，则它有唯一所选 A⊃P。核心 K_A 不必包含子组核心 K_P；下文重新检查实际端点是否属于 **K_A**，不继承 child-core 好状态。

记

\[
 U_\theta=\bigcup_A(K_A+Q(0,\rho_A)),\qquad
 \mu_{\rm bad,\theta}=
 \sum_A\mu_{\rm hi}|_{A\setminus K_A}. \tag{6}
\]

这是原输入的固定子测度，不是重启 FIRST 的新输入。由于来源 A 不交，

\[
 \|\mu_{\rm bad,\theta}\|\le\zeta W_{\rm hi},\qquad
 q|U_\theta|\le T\sum_Ah_A\le TW_{\rm hi}. \tag{7}
\]

在 §2 的实际两序子交通中，原选定端点位于 K_A 的 good 部分由 (3) 和 ρ_P≤ρ_A 支撑在 U_θ 内。先合并原唯一标签再付行密度，得到

\[
 \mathcal R_{\rm good,\theta}\le C_hq|U_\theta|
 \le C_hT W_{\rm hi}. \tag{8}
\]

对 bad 部分：D_h 上 hardborn 坏来源满足 μ_hi|_{B∩A\K_A}≤μ_bad，完整软行≤q，故点态被 h_{R_h}*μ_bad 支配。D_s,K=0 上，原 P^∅_{σ,L_s}≤h_{L_s}，完整硬行≤C_hq，故被 C_h h_{L_s}*μ_bad 支配。其它门均保留在原子交通中，仅在证明该正上界时删去。

关键是 (4) 是原同一 x 的互斥选择，两项对 **同一** μ_bad 成为

\[
 \mathcal R_{\rm bad,\theta}(x)
 \le C_h\sup_{a\le r\le b}(h_r*\mu_{\rm bad,\theta})(x).
 \tag{9}
\]

如果另一个接口的尺度顺序依 pair 或额外标签，(9) 不能直接套用，需保守加两份列费用。本稿的 L_s(x),R_h(x) 是原共同赢家，因此无此变化。

原硬列包络常数

\[
 N_h=\int\sup_{a\le r\le b}h_r(u)\,du
 =1+n\log(b/a) \tag{10}
\]

给 ∫sup_r(h_r*ν)≤N_h‖ν‖。有限原 J 只减小 sup，连续列 (10) 是合法正上界。故

\[
 \boxed{\quad
 \mathcal R_{\rm marked,\theta}^{D_h\cup(D_s,K=0)}
 \le C_h(T+\zeta N_h)W_{\rm hi}.
 \quad} \tag{11}
\]

种子平均保持同一上界。取 ζ=n⁻¹ᐟ²、T=√n L_n，b/a≤2 时，(11)≤C_h√n[L_n+log2+1/n]W。没有 J、atom count、父组数量或选核数。它只删除当前未付交通的标记祖先子支；不是把整个空间 U_θ 里所有历史账再付一次。

## 5. 可选全 K 扩展：用严格 GP 控制 hard 半径

不经 mark 均值，(1) 和 L_s≥a 还直接给

\[
 R_h<(h_P/q)^{2/n}/a,
 \qquad R_h\le r_P:=\min\{b,(h_P/q)^{2/n}/a\}. \tag{12}
\]

r_P 沿来源祖先单调；相关 strictCP 父组有 r_P≥ρ_P≥a。将 (5) 中 ρ_P 换 r_P，重新作一个来源固定最粗祖先反链，所有原 K/两序都以 hard z 为选定端点。good 费用 C_hT W；bad 由 soft 行≤q 给 h_{R_h}*μ_bad，无需 C_h。这给独立可选版本

\[
 \mathcal R_{\rm GPcore}\le(C_hT+\zeta N_h)W. \tag{13}
\]

该条件更严：r_P≥ρ_P，不能宣称 (11) 的全部标记组自动满足 (13)。若同时采用两个版本，要在当前实际交通取依次互补交集，并显式付两份费用；不能把两个不同核心反链暗合为一份。旧 GP-COLUMN 给的是已付大几何平均尺度列，本节 (12) 是其严格补集的来源固定上半径应用。

但 (13) 的未饱和部分并未推进困难空间交叉。ζ<1、M_P>0 使合格核心非空，故 |K+Qr_P|≥r_Pⁿ。若 r_P<b，标记必有

\[
 h_P/(qa^n)=q r_P^n/h_P\le T,
 \qquad R_h<aT^{2/n},\quad
 \beta_h=n\log(R_h/a)\le2\log T. \tag{13a}
\]

于是任意该输出的原 hardsource-once 表示，直接使用旧低 β 硬列工具即可付
\[
 [1+\min\{n\log(b/a),2\log T\}]W.
\]
这不需 core 或 GP 核费。T=√n k_n⁴ 时只是 β_h≤log n+8log k_n 的易短窗口；不将它作为新困难空间分支。已饱和 r_P=b 时 (13) 就是下面的 trimmed-b 支持覆盖，与 GP 无关。

### 5a. 主账推荐：任意来源的 trimmed-b 核心分层

为避免给易短窗口换名字，主账只推荐如下更直接的来源固定判据。对每个 M_P>0 的原父组，标记当存在规范 K 满足

\[
 \mu_{\rm hi}(P\setminus K)\le\zeta M_P,
 \qquad q|K+Q(0,b)|\le T M_P. \tag{13b}
\]

这里可用完整 M 而非 h，覆盖更强，费用不变。保留 h 版也正确但更严格。(13b) 不需要 CP/GP 推半径；原 Rh≤b 已够。取它自己的最粗标记祖先反链 A 及其 K_A，不使用 §4/5 的其它核心继承。对所有实际 K 和两序，重新检查 hard z∈K_A 的 good 部分；这时 U=∪_A(K_A+Qb)，q|U|≤TΣ_AM_A≤TW。bad hardborn 来源仍被 μ_bad=Σ_Aμ_hi|_{A\K_A} 支配，‖μ_bad‖≤ζW。原软行≤q 给单份硬列。因此

\[
 \boxed{\mathcal R_{\rm trimmed\text{-}b}\le(C_hT+\zeta N_h)W.} \tag{13c}
\]

这是旧 b-density source-cover 工具的鲁棒核心升级：新的付款内容是允许丢少量来源，并以一个固定 badsource 列清算，不是新的 GP 交叉核定理。ζ=n⁻¹ᐟ²、T=√n k_n⁴、k_n=⌈log₂(n+2)⌉，给 [C_h√n k_n⁴+N_h/√n]W。只在当前 R** 删除该实际来源子支；§4/5 的其它版本保留工具，不另加费用。当前主账移项系数 8192/49 乘新费用一次。留下 R***，其原 LCA 的每个祖先均无 (13b) 核心证书，原 R** 所有门仍在。没有把旧 128/63 改写为当前系数。

## 6. 为什么 bbox/随机 mark 不给免费闭合

若 ρ_P<b 且把核心取整父盒，边长 ℓ_P，那么 qρ_Pⁿ=h_P/√n，(5) 必须满足

\[
 (1+\ell_P/\rho_P)^n\le T\sqrt n.
 \tag{14}
\]

T=√nL_n 时 ℓ_P≤b[(nL_n)^{1/n}−1]。原实际 far 要求 ‖y−z‖₁>H_mark，而 y,z∈P 给 nℓ_P>H_mark。于是若

\[
 nb[(nL_n)^{1/n}−1]\le H_{\rm mark}, \tag{15}
\]

整父盒标记在该实际 far 域为空。旧 tex 2354–2385 的 H_mark 包含 (1664/3)a log(n+2) 及正项；b≤2a、L_n=n^{o(1)} 时 (15) 渐近成立。ρ_P=b 的整父盒判据是旧 b 覆盖判据的子类。GP 未饱和 r_P<b 时 q r_Pⁿ/h_P=h_P/(qaⁿ)>√n，整父盒要求 (1+ℓ_P/r_P)^n<L_n，更严。因此本稿不能拿整父盒或有限网络特别族冒充一般新覆盖率。

K>0 的另一尝试是冻结诊断掩码 A 和 dimensionless mark v，将 soft source 管扫成 K_P+Qρ_P+{tv:a≤t≤ρ_P}。但原 cap 在 x 仅约束平均核 P_{σ,L_s}*μ(x)≤q；冻结 mark 后的 h_{L_s}*μ(x−L_s v) 没有同样 q 上界。L_s=L_s(x) 时也没有 x↦x−L_s(x)v 的已给 Jacobian。后验 mark 分配不能换成先验概率。不管扫管体积是否已证，该条件行上限是独立缺口。

此外，即使先补上条件行控制，σ(x) 可随输出选择，全 mask 正包络要付已证

\[
 D_n=\sum_{k=0}^n\binom nk(k/n)^k(1-k/n)^{n-k}
 \le1+\pi\sqrt n/2,
 \quad\sum_A|A|c_{|A|}=nD_n/2. \tag{16}
\]

其中 c_k=(k/n)^k(1−k/n)^{n−k}，端点按 0⁰=1。因此不能用 E K≤n 免费平均一个输出依赖的 mask/mark selector；先 supσ 会多 √n。这是旧 `signed_band_geometry/allsoft_angular_review.md` 的已核工具，不重跑其系数。本稿 (13) 用原 hard 端点避免这次交换，仍需实际核心容量足够小。

## 7. 支付后的真实补集

若登记 (11)，在 D_h 或 D_s,K=0 的剩余实际源对中，其原 LCA 的每个相关祖先 A 均无 (5) 的可执行证书，即对每个规范核心 K，

\[
 \mu_{\rm hi}(A\setminus K)\le\zeta M_A
 \ \Longrightarrow\ q|K+Q\rho_A|>Th_A. \tag{17}
\]

任意满足该质量条件的紧核心至少有 q|K+Qρ_A|≥Th_A；等号只是外近似的边界。若 Borel 核心覆盖≥(1−ζ/2)M_A，有限测度的内正则性可取紧子核心覆盖≥(1−ζ)M_A，同样给非严格体积下界。不能无余量把 ζ/2 去掉。若登记 (13)，(17) 换 r_A，约束所有 K/两序。

主账若采用 (13b)，其准确补集改为对原 LCA 的每个祖先 A、每个规范 K：
\[
 \mu_{\rm hi}(A\setminus K)\le\zeta M_A
 \ \Longrightarrow\ q|K+Qb|>T M_A. \tag{17b}
\]
紧核心和 Borel 核心的边界范围依上文相同，只需把右侧换 TM_A。这里每个实际父组的 M_A>0，避免零质量判据；祖先不要求另有严格 CP 资格才作 (13b) 检查。

这是真正来源几何的鲁棒容量限制，不因给支持添极小质量的远点立即失效。它与 (1)、原捕获小响应、lowS、lowcoin、双 cutoff、短 hard 壳和全部原历史门同时保留。CP/GP 给的是 h_A 与物理尺度的关系，双 cutoff 给原 soft 标签响应；(17b) 给高质量来源核心的 Minkowski 体积下界。目前没有从三者联立推得剩余列小或所有祖先必标记。将未标记的来源换成每个 x 当前捕获子来源仍非法。一般 √n 目标尚未闭合。

更具体地，h_A≥h_P 使原输出的两条数值严格不等式沿祖先仍成立（不宣称祖先继承其它历史资格）。令 m_x=min(L_s,R_h)、V_K=|K+Qb|，(17b) 对每个合格规范 K 给
\[
 \frac{V_K}{m_x^n}
 >\frac{T M_A}{q m_x^n}
 >\frac{T\sqrt n}{v_j}\frac{M_A}{h_A},
 \qquad
 \frac{V_K}{(L_sR_h)^{n/2}}
 >T\frac{M_A}{h_A}. \tag{18}
\]
推荐 T 下第一项大于 n k_n⁴ M_A/(v_jh_A)。这是严格的来源鲁棒空间扩散合同，保留 v_j 和 h_A/M_A；它不是当前捕获体积，也不是输出集合测度下界或上界。原 lowS/双 cutoff 约束同一 source y 的加权核响应及相对于 L_s 的坐标位置，没有已证步骤把 (18) 的整来源核心体积与不同 x 的来源固定接受列配成可积尾。source-only entropy 很大时，这个交换依然是未付缺口；不能用捕获子测度替代 μ_hi|A 以制造小核心。

## 8. 三轮守卫的预登记范围

新脚本前缀 `double_cutoff_actual_geometry_`，三轮 n=4,16,64（√n=2,4,8），只核本稿新一般接口的有限有理证书：任意正来源权及出生子质量；最粗反链源质量；非嵌套 ancestor/child cores 的端点重新分类；闭面及未饱和/饱和半径幂比较；有限盒并精确体积；同一 x 两序互斥的原边际支配与同一个 μ_bad 清算。逐 x 的 soft/hard 值是预注册正行代数，不模拟原 φ、actual FIRST 或完整历史。不会据其报告实际覆盖率、实际余项样本或新阶数。

执行已终态：`double_cutoff_actual_geometry_exact_guard_20261007.py` 与同前缀 `_results.json`。三轮各 510 项，共 **1530 项 Fraction 精确 PASS**，无随机种子、无 live handle。每轮用 8 份不等正有理来源权与不同出生子质量、三种有限标记配置；有实际 child-core good 但 ancestor-core bad 的重新分类检查。高维盒并体积由共同后 n−2 坐标截面因子乘精确二维并面积，二维面积同时用坐标扫与独立 inclusion-exclusion 比对，绝不用总 bbox 代替并体积。

半径幂检查另覆盖未饱和、恰达 b、根半径>b 后 min=b 的合法闭端点、单尺度 a=b。边际组件预注册了正软/硬行数组与一个精确归一的有限来源列常数 1，验证 GP 单 hardbad 列，以及 CP 同 x 互斥两序的单个 C_h 包络。此有限列常数是组件证书，不替代连续原 N_h；数组不是原 φ，也未认证 fullfuture/FIRST/far/CP 联合实际样本。几何组件实际比较的是较强 qV≤Th，并未另运行 M 版标记比较；h≤M 使这些保存证书解析蕴含 qV≤TM。保存的最粗反链、ΣM≤W、坏来源质量检查直接对应 M 预算。主版本 (13c) 的证明不依赖该有限输入或正数组，升级 M 判据未重跑旧守卫。

当前主账仅推荐 (13c)，(11)、(13) 是范围工具；不为 GP 未饱和易短窗或旧全 mask 系数另跑守卫、不另领费用。一般空间补集仍是 (17b) 联立全部原实际门，尚未闭合。

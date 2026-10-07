# 原 hard 半径的来源固定核心：可和预算与同价扩用

2026-10-07。一般平方根空间余项仍未闭合。本稿给出一个可证明的来源预算升级：覆盖证书可冻结在不同物理半径，每个来源节点只取一份证书；沿每条来源路径的收费权重总和不超过 1。覆盖并集和一个固定坏来源集合因此同时有总预算，不按每个 β 带重复领取 W。保留原 FIRST、whole μ、fullfuture、原赢家/出生/容量币/LCA/history/continuation 和当前所有低门。

本轮实际阅读 `math-n04-summable-error-entropy-descent/SKILL.md` 及 references 的 provenance、method、cube-interface；只使用 method 中基础望远镜与正误差可和步骤，不引用其中待验文献为 cube 定理。下文势能是来源收费余额，不是未经证明的实际交通熵下降。

## 1. 查重及一个同价的根端增强

旧主笔记“任意输入的大质量核心与一次坏来源清算”、`double_cutoff_actual_geometry_interface_20261007.md`、`actual_parent_occupancy_route_20261007.md`、`gated_occupancy_bridge.md` 的来源一次/固定窄窗工具已经核对。旧固定 b 核心判据为

\[
 M_P=\mu_{\rm hi}(P)>0,\quad
 \mu_{\rm hi}(P\setminus K_P)\le\zeta M_P,
 \quad q|K_P+Qb|\le T M_P. \tag{1}
\]

ζ=n⁻¹ᐟ²，T=√n k_n⁴，k_n=⌈log₂(n+2)⌉，a≤b≤2a。q 是原共同 FIRST 高度。核心取全空间有限有理闭轴盒并的规范可数枚举；允许退化闭盒，严格体积余量也可用非退化闭外逼。原森林是有限深度、可数节点，原半开来源格仍唯一分配来源面上的质量。

固定 b 支的证明仅需 **原 hard z 属于所选来源节点**，不需要 soft y 属于该节点，也不需要该节点是实际 pair LCA 的祖先。取 (1) 的最粗标记节点反链 𝒜_b，它们在来源上不交；S_b=∪𝒜_b。对所有实际 hard z∈S_b 的交通，以其唯一所属 A∈𝒜_b 分类。

若 z∈K_A，原 Rh≤b 给 x∈K_A+Qb；合并原唯一标签后，整个子交通密度≤C_hq。若 z∉K_A，出生子来源被固定 μ_bad=Σ_A μ_hi|_{A\K_A} 支配；原完整软行≤q 给 h_Rh*μ_bad。故同一费用 (C_hT+ζN_h)W 支付所有 hard membership 子交通，包括 y在A外、原LCA未marked 的情况。原两端历史/LCA没有重定义，只是在证明正上界时删除它们。不能逐节点重新付一份总行。

因此旧 LCA-ancestor 覆盖是这个根端支的子集。同价替换后的剩余 hard z 在它自己的来源路径上没有任何固定 b 核心证书。原保存1530项守卫不认证该增强，它们确实先要求 y,z 同在所选祖先内；本轮另注册 hard-child-only 组件。

## 2. 来源路径的预算合同

令 𝒫 为允许认证的原来源森林节点族。明确记

\[
 D=\sup_z\#\{P\in\mathcal P:z\in P\}<\infty. \tag{2}
\]

若旧 J 是边深度，则 D≤J+1；不可把 root/leaf 的一层漏掉。零来源质量节点不用。取来源固定、输出无关的非负 w_P，满足

\[
 \sum_{P\ni z}w_P\le1\quad(\mu_{\rm hi}\text{-a.e. }z). \tag{3}
\]

例如所有节点取 1/D；也可在一个来源反链上取 1、其它节点取 0。不假设输入原子数、来源标签独立性或每个父组孩子数均匀。

w_P 是本稿的收费预算，完全不同于原容量失败权 (1−p_Cy)(1−r_Cz)。它只进入核心证书与收费余额，不乘回、归一化或替换原实际交通；原失败权及其全部历史保留。

按可数节点枚举定义余额

\[
 H_m=\int\Big[1-\sum_{j\le m}w_{P_j}1_{P_j}(z)\Big]d\mu_{\rm hi}(z),
 \qquad d_m=w_{P_m}M_{P_m}. \tag{4}
\]

H_m≥0，H_0=W_hi，逐步 H_m=H_{m−1}−d_m。所以 N04 基础求和严格给 Σ_P d_P≤W_hi。无迭代次数费用。d_P 是预分配来源收费，未把原输入削减为一个新 FIRST 输入。

## 3. 每节点仅冻结一份最大半径证书

预先固定有限物理网格 a=r_0<⋯<r_L=b，设 β_ℓ=nlog(r_ℓ/a)。可取 r_ℓ=min{b,a(1+1/n)^ℓ}，正步长 β≤1，含 b 闭端点。理论允许原有限尺度集或连续窗口中的 Rh；网格只用于 **源证书**，不更换原 winner。判断实际覆盖必须 Rh≤证书半径，不能拿 floor(Rh) 代替 Rh。

对 w_P>0，某半径 r_ℓ 可认证，当存在规范 K 满足

\[
 \mu_{\rm hi}(P\setminus K)\le\zeta w_PM_P,
 \qquad q|K+Qr_\ell|\le T w_PM_P. \tag{5}
\]

在可认证半径中选最大的 r_P，再从该半径的可数枚举取第一份合格 K_P。一个节点最终只保存一份 (r_P,K_P)。没有证书则略去，不令 r_P=a 冒充合格。因为网格有限，最大半径和首核心均可测。ζw_P 和 Tw_P 是来源固定预算；不随 actual x、Rh、σ、history 或捕获质量改变。

一份在 r_P 认证的核心覆盖它所需的所有更小 Rh，且其体积一直按 K_P+Qr_P 付一次。可认证半径存在性关于 r 单调，但核心不必跨节点、跨候选半径嵌套；最大半径选择后只使用保存的那份 K_P。不得每个 β 带再取一份核心而保留同一 (5) 预算。

## 4. 固定并集与固定坏来源的总清算

令 𝒞 为最终持有证书的节点。定义

\[
 U_\theta=\bigcup_{P\in\mathcal C}(K_P+Qr_P),\quad
 F_\theta=\bigcup_{P\in\mathcal C}(P\setminus K_P),\quad
 \nu_{\rm bad,\theta}=\mu_{\rm hi}|_{F_\theta}. \tag{6}
\]

ν_bad 是 **集合并的限制**，不是把重复节点的坏测度直接相加后冒充 μhi 子测度。来源节点可以嵌套，核心可以不嵌套。由 (4)–(5)、Tonelli 和并集的次可加性，

\[
 q|U_\theta|\le T\sum_{\mathcal C}w_PM_P\le TW_{\rm hi},
 \qquad\|\nu_{\rm bad,\theta}\|\le\zeta\sum_{\mathcal C}w_PM_P\le\zeta W_{\rm hi}. \tag{7}
\]

这里两类正误差分别由 Td_P、ζd_P 支付，总和已经核对，不能从“每步小”直接跳到总量小。这正是 N04 在本问题的实际用途。

本次实际可覆盖判据为：原 hard z 属于某认证 P，且原 Rh(x)≤r_P。可在 z 的原节点路径上选最粗的半径合格节点；不要求它包含 y 或是 LCA 祖先。若 z 属于该节点的保存核心，x∈U_θ；否则 z∈F_θ。必须检查本次保存的 K_P，不能继承其它节点的 good 状态。所有原门仍在当前实际子交通内。

good 总行≤C_hq。bad 保持原 hard source-once 接受权≤1，出生坏来源≤ν_bad，原完整软行≤q。因此

\[
 \boxed{\mathcal R_{\rm adapt,\theta}
 \le C_hq|U_\theta|+N_h\|\nu_{\rm bad,\theta}\|
 \le(C_hT+\zeta N_h)W_{\rm hi}.} \tag{8}
\]

N_h=1+nlog(b/a)；原任意 Rh∈[a,b] 的一源硬列被它支配。种子平均无新倍数。此处没有 β 带数、D、J、atom count 或候选核心数的费用；D 只使具体 (5) 判据更严格。所有 K/两序共用原 hard 端点，完全不需要把移动 σ 的掩码后验换成先验平均。

## 5. 与旧固定 b 支的同价拼接

为了严格保留旧所有固定 b 资格，而不新增一份费用，先保留 §1 的 𝒜_b、原 K_A 和 r_A=b。令

\[
 w_P=\begin{cases}
 1,&P\in\mathcal A_b,\\
 1/D,&P\text{ 与每个 }A\in\mathcal A_b\text{ 来源格不交},\\
 0,&\text{其它节点}.
 \end{cases} \tag{9}
\]

对旧反链不重选核心；它们满足 (5)。对新的不交节点作 §3 最大半径认证。z∈S_b 时其路径只在唯一旧 A 处有预算 1；z∉S_b 时至多 D 个新节点各有 1/D，故 (3) 成立。尤其

\[
 \sum_{A\in\mathcal A_b}M_A+
 \sum_{P\cap S_b=\varnothing}\frac{M_P}{D}\le W_{\rm hi}. \tag{10}
\]

没有把 outside S_b 来源重新作为 FIRST 输入；新 P 的 M_P、原 μ、原 soft 选择和历史都保持。只是把尚未消费的来源收费余额分给这些原节点。

用该证书族时，(8) **同价包含**旧所有 fixed-b hard membership 支，并增加 outside S_b 的 adaptive 半径资格。它替换现有核心费用，不额外增加 (C_hT+ζN_h)W。现有条件主账系数 8192/49 对这份费用仍只乘一次。其它既有费用保持原互补登记。

若没有采用 (9) 的拼接而另选一组 w，(8) 仍正确，但不能免费宣称覆盖旧支；如果先付旧支再付一个无拼接的新系统，则必须另外登记新费。两种情况不能混写。

## 6. 与固定 b 的严格范围比较

对规范有限盒并 K 和 a≤r≤b，有

\[
 |K+Qb|\le(b/r)^n|K+Qr|. \tag{11}
\]

自含证明：在一个坐标切片，K+Qr 是有限并的闭区间，每个非空连通分量长≥r。加长 b−r 至多每个分量增加 b−r，重叠只减小长度，故该坐标方向体积倍率≤b/r。已扩大其它坐标后剩余坐标切片仍是 r 长区间并。依次 n 个坐标用 Fubini 得 (11)，不把高维并体积换成 bbox。

在 fixed-b 补集里，若 (5) 认证 r_P，则同一个核心满足旧质量容限（ζw_P≤ζ），但旧体积判据失败。因此

\[
 T M_P<q|K_P+Qb|\le(b/r_P)^nT w_PM_P,
 \quad n\log(b/r_P)>\log(1/w_P). \tag{12}
\]

即 adaptive 新证书必须离 b 的 β 上端至少 log(1/w_P)。w=1/D 时须 β_P<B−log D，B=nlog(b/a)。这个是靠近 b 的已有覆盖重叠检验，**不是** β≤O(log n) 的 GP 易短窗。半径不是由 h/q 的 GP 幂定义，不能将旧 GP 未饱和的 β≤2logT 套到它。

在纯 source-contract 容量层，扩用可以严格发生：取 V_b/V_r=2^d、qV_b/(TM)=2、w=1/D，则旧 b 不合格，而新 qV_r/(TwM)=2D/2^d<1（2^d>2D）。这些是容量比例证书；没有声明实现了真实 winner/fullfuture/FIRST/far/history，所以不能把它叫 actual 余项反例或实际覆盖率。

旧 source-cone 对 source-fixed 宽窗给 (1+宽度)W，覆盖整个窗回到 O(n)W。单靠逐 β 窗求和不能给 (8)。本稿区别在于同时证明了实际 good 输出并集的 q|U|预算和一个 **全窗固定** ν_bad 的小质量，依靠 (5) 与来源路径余额，而非只重新命名窄窗列。仍未证明剩余所有来源都能取得这些几何证书。

## 7. 准确剩余与下界表检验

同价拼接后的未付交通保留原所有门。在其原 hard z 路径上：不存在 fixed-b 证书；对每个持有新 (5) 证书的原节点 P，Rh(x)>r_P。未持证书的节点不能用默认半径虚构新约束。

等价可执行表述：对每个 w_P>0、P∋z 和网格 r_ℓ≥Rh(x)，所有满足 μhi(P\K)≤ζw_PM_P 的规范 K 都有 q|K+Qr_ℓ|>Tw_PM_P。只需检查第一个 grid ceiling(Rh)，因为体积随 r 单调。紧核心的体积只有非严格下界，Borel 核心需质量余量；与前稿的闭外逼边界范围相同。

这仅是本轮注册有限网格的补集。不能由网格未认证断言连续最优半径或连续半径证书不存在；未证明免费插值 transfer。主账可记更强补集 R_dagger：每个原 hard-source 路径上的 **已保存** 证书 P 都满足原 Rh>r_P，未保存的节点忽略。旧 R***保留为历史，不把这个新判据当新 FIRST 定义。

该补集比旧 b 覆盖具有 source-only 半径信息，但是没有推出实际 Lebesgue 接收空间小。κ_C/低诊断 K 不使 (5) 自动成立，也不让不同 x 的捕获子源可以重启 FIRST。全 future cap 与 lowS 仍作用在原完整核响应，没有从它们证明来源路径预算能大范围认证。

已复读下界构造总表 A/B：非均匀联合峰质量、Gibbs/XOR/隐变量、增长 Cantor 峰型都允许高度复杂的正来源标签。这里 (3) 只按原有限森林路径分配预算，独立于标签数与是否产品；不把 Bernstein K 当下界源占用 k。不援引固定峰数或最细短窗上限排除一般输入。新数值不重跑这些旧构造。

## 8. 预登记：三轮新接口守卫

新前缀 `adaptive_core_source_budget_exact_guard_20261007`。三轮有限组件 n=4,16,64，并分别附压缩超大维 n=2⁴⁰、2⁶⁴、2⁸⁰ 的精确指数/有理预算检查。无随机种子。有限组件用任意不等正来源权、原层级节点及非嵌套核心，核 (3)/(7)/(9)/(10)、每节点只留最大认证半径、grid 上闭端点与严格剩余、bad **集合并**质量、原 hard membership 边际。另加 old b 的 hard-child-only marked、y 在节点外且 LCA 未marked 的资格收据，不称原1530已验此增强。

几何组件只用有限原子/闭点盒的精确并体积；未计算原 φ/FIRST/future/history，可能已被其它实际 paid 门排除。压缩大维只核容量比与预算；不构造超大维 cube 输入。必须报告 T 相对 n 的有限范围：n=2⁴⁰ 时 √n k⁴仍可能比 N_h 大，不能用小 n PASS 外推有用渐近改进。大维比例检验还核必要 CP/GP 标量兼容性，但不冒称完整联合合同。

执行已终态：`adaptive_core_source_budget_exact_guard_20261007.py` 与同前缀 `_results.json`，无 live handle。三轮分别 92、298、991 项，共 **1381 项 Fraction 精确 PASS**。未导入或重跑任何旧脚本/保存模型。

有限来源几何组件用 11 份不等正有理质量、来源固定出生子质量及一棵可实现为轴盒层级的有限源树，D=5。核心为有限闭点盒，r≤b 时不同点的 padding 方体互不相交；同一点多个核心的 padding 方体同心嵌套，所以真正全并体积可精确按最大半径算，不用 bbox。旧 b 判据的点数最小覆盖检查认证只有 hard-child 节点2合格。新每节点最大 radius 选择没有逐 β 重建 core。

三轮各有严格新增的 hard-child-only 收据：y 在节点2外、z在节点2内、原LCA=1未marked，但新 hard membership 合格，并在独立正行数组中有正交通。此证明所需 row/column 是未归一化原来源权的代数组件，不是原 φ/FIRST/fullhistory。n=4 没有新增 adaptive 节点；n=16、64 各有9个，ancestor-core bad / child-core good 重新分类分别3、37次，坏来源 **集合并**质量均1/1000。没有把这些比例当 actual coverage。

压缩参数的 T/n 如下；只保存 cancelled 指数和有理比，不产生超大维几何数据：

| n | T/n | T<n/2 |
|---|---|---|
| 2⁴⁰ | 2825761/1048576 | 否 |
| 2⁶⁴ | 17850625/4294967296 | 是 |
| 2⁸⁰ | 43046721/1099511627776 | 是 |

压缩证书旧 b 比值2被拒，新 adaptive 比值为2D/2^d<1；另核 h=M 下 CP/GP 必要的标量尺度关系，且其 β 不属于旧 GP 的 β≤2logT 易短窗。此处没有注册原 σ 带的实际 v_j 或 φ/wholefuture/FIRST，因而不认证原严格 CP、其它实际 paid 门或完整 actual 残余。n=2⁴⁰ 的 T 本身仍大于 n，小 n组件PASS当然不能证明渐近策略已优于原 N_h；只报告已证渐近预算 O(√n log⁴ n)，不由数据拟合阶数。

新脚本压缩源分账用不交来源区的2/5与3/5收费，合计1；没有把两个单位路径都算在同一 source 上。保存结果以此澄清后的最终版本为准。

重要解析结论为同价包含式 (8)–(10)。它给一个可和的实际空间子支付款，并使旧 fixed-b 的 branch 资格严格可扩用；严格扩用数值只认证 source-contract，可能早已被其它实际门排除。完整一般空间余项仍未闭合。

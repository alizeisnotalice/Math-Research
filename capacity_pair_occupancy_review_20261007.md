# Paid-coin amplification 的独立核查

2026-10-07。数学结论接受：对原 HR 的 sameTop union/failure 分数，可以直接放大旧更广 receipt 费用，支付当前带全部实际门的高 union-fraction failure 子支。它是旧工具的可选 corollary，不是新的空间列定理；本审稿不登记主账、不新增模型、不运行数值。

读取依据：原 `概率接口阶段证明.tex` 2690–2760（HR-COST/HR-PARTIAL）；`high_marked_cover/next_spatial_capacity.md` 的 SC1–SC5；主提示词 112–185 的最新条件账。初审时作者稿尚未生成，先直接审原合同；最终又全文只读 `capacity_pair_occupancy_20261007.md` §§1–5 与 root `current_joint_budget_20261007.md` §§1–4。以下采用最终阈值，最终对照范围见第6节。

## 1. 原分数、广义交通与严格门

令 D_θ 为 sameTop 且两个来源在不同最细格的合法 distinct-child/LCA 域；原 actual far 子项在此域内。以下 sameTop 指 D_θ。在此域的原源对有唯一 LCA P，原来源固定分数为

\[
 p=p_{C_y}=\min\{1,m(C_y)/(\sqrt n\,h(P\setminus C_y))\},
\]
\[
 r=r_{C_z}=\min\{1,h(C_z)/(n\,m(P\setminus C_z))\},
 \quad \ell=p+r-pr,
 \quad w=(1-p)(1-r)=1-\ell.
\]

零分母按原约定取 1。域外令 G=ℓ=w=0；w=1−ℓ 只在 D_θ 使用，不能在无合法 LCA 的域外套它。m、h 始终是原完整 source/出生合格 sibling 质量，不能换成实际捕获质量、low-S 热源质量、首落点或 terminal 来源质量。两个端点的一份 receipt 共用一个 Bernoulli；union 的期望是 ℓ，不能平方 p 或 r。

用 dΠ 表示原更广 normalized HMP 交通，含原 FIRST、两个赢家、来源资格等。HR 已证明

\[
 R_{\rm union}:=\mathbb E_\theta\int
 1_{\rm sameTop}\ell\,d\Pi
 \le J A_n W_{\rm hi},
 \quad A_n=N_h/\sqrt n+C_hJ_s/n. \tag{1}
\]

当前 actual failure 子项须保留原 w：若所有后加 CP/GP、历史、future、low-S 等门的接受权为 0≤G≤1，则被

\[
 \mathbb E_\theta\int1_{\rm sameTop}wG\,d\Pi \tag{2}
\]

支配。接受权只能保留或丢去以求正上界，不能重新定义来源或 FIRST。若某个 reference 放大丢失了原 w 或改变了交通归一化，不能仅凭它仍叫 failure 就套此审查。最新 conditional 账的回代系数也不能从旧名称推定。

## 2. 放大不等式及端点

预先固定 0<δ≤1。因为 (1−t)/t 在 t>0 递减，在 ℓ≥δ 上逐点有

\[
 w=1-\ell\le\frac{1-\delta}{\delta}\ell.
\]

由 (1)–(2)，

\[
 \boxed{R_{{\rm current},\ell\ge\delta}
 \le C_\delta J A_n W_{\rm hi},
 \quad C_\delta=(1-\delta)/\delta.} \tag{3}
\]

甚至不需要广义 union 与当前所有门的精确相容分割；先将当前接受权正放大到旧基测度即可。因此新 branch 由旧更广已付分支的 certificate 直接推出，没有证明新的 kernel column。

等号 ℓ=δ 计入 paid，补集严格为 ℓ<δ。δ=1 时高支仅 ℓ=1、failure 为零，(3) 的系数为零；δ=0 不适用。p=1 或 r=1（包括原零分母约定）也使 failure 为零。最终预设 k_n=ceil log₂(n+2)、δ_n=k_n⁻⁴，n≥512 时在 (0,1)，确定的是输入维数阈值，不事后选择有利 pair 阈值。它替换初审候选的自然对数 coin 阈值，仅使端点有理可复现；原 low-S 的 ε_n=log⁻⁴(n+2) 不变，δ 与 ε 不能混用。

既有 J=O(log n)、A_n=O(√n) 给 (3) 为 O(√nlog⁵(n+2))W，符合已有 log⁶ 目标余量。这不是把 δ⁻¹ 隐藏进常数。

## 3. 旧 HR/MM 的账：只适用于旧账

旧 HR-PARTIAL 与 MM 回代为

\[
 R_{\rm HMP}\le2J A_nW+2R_{\rm bal},
 \quad\lambda|E|\le P_{\rm MM}+(64/63)R_{\rm HMP}.
\]

所以在旧 R_bal 的当前补支新增 (3)，其交通转换系数是 **128/63 一次**，不得再乘 HR 的 2 或 direct-output 的 4/3。原 HR union 费 J A_n 已在旧已付系数中；加 CδJ A_n 后总计恰 δ⁻¹J A_n。等价地可以把这一个既有 HR 项的系数从 1 替换为 δ⁻¹，并保留 low-ℓ failure。不能“替换”后再额外叠加同一 Cδ 项。

只在仍未付的 current 子交通应用 (3)，先与所有旧 paid 支取准确互补。引用旧 union certificate 作为较广上界是合法的，但不得把该 union 本身再当一份新的实际交通收费。SC 的非层级币与 HR 的层级币定义不同，本稿不能混用两份 ℓ 或把 SC/HR 两个替代整账相加。

## 4. 最新主提示词的条件账：保守使用新系数

原提示词明确给

\[
 \lambda|E|\le2\mathcal B_{\rm paid}+(4096/49)R_{\rm heavy}. \tag{4}
\]

若将 (3) 在 R_heavy 的实际 failure 子项中使用，未吸收前新增费为 **(4096/49)CδJ A_nW 一次**。不能因旧 receipt 用了 HR 就把 (4) 的系数改回 128/63；新旧余项的名称不证明归一化相同。

若当前另外已有

\[
 R_{\rm heavy}\le F_nW+(49/8192)\lambda|E|+R_{\rm lowS},
\]

那么先代回 (4) 再精确移项得到

\[
 \lambda|E|\le4\mathcal B_{\rm paid}
 +(8192/49)F_nW+(8192/49)R_{\rm lowS}. \tag{5}
\]

在 R_lowS 上加 (3)，最终新费是 **(8192/49)CδJ A_nW 一次**，剩余 R_lowS,ℓ<δ 的系数也是 8192/49。这个最终倍 2 来自已公开的吸收，不是另一份 HR 顶格覆盖倍数。若 F_n 已含某费，它也随这次移项统一倍 2，不能只倍新费或保留旧主账系数混拼。

本审查只核上述条件代数；上游最新三个尺度原件尚未核得，故不宣称已证明旧 R_geom 与新 R_heavy 完全等同，也不擅自登记历史费用替换进最新版主账。

## 5. 剩余获得的约束及实际限制

在严格补集 ℓ<δ<1，必有 p<δ、r<δ，所有相关 min 截断不活跃，于是

\[
 m(C_y)<\delta\sqrt n\,h(P\setminus C_y),
 \qquad h(C_z)<\delta n\,m(P\setminus C_z). \tag{6}
\]

它们是来源固定容量比的较强必要条件；两个分母是不同 sibling 质量，不能推成两个 endpoint 质量的简单互比。更精确的联合条件仍是 p+r−pr<δ，(6) 只是其必要弱式。

(6) 不限制 moving winner/capture 的空间占用，也不证明 source-once 的低成本共同覆盖。不把固定 parent h 改成 h_x，不将已有 source-cone/CA 表示重报为新费用。新筛选可以合法地让未付 pair 更难被原 receipt 覆盖；没有因此支付全部低-ℓ actual 分支或证明它为空。

最终审查：原 coin 放大、严格端点及旧工具继承均接受；最新账保留 (4)–(5) 的独立系数口径。仅作为可选已证删支，不新增主账登记，不声称新的空间几何进展。未运行实验。

## 6. 最终作者稿与当前合并账的只读对照

作者最终稿的合法 LCA 域、域外零定义、双币失败条件下接受率 G∈[0,1] 的表示均与上文一致。即使后续 actual 门依赖原辅助币，也可在“双币失败”条件下取接受率；此处并未假定门与币独立。必须仍保留原 failure 先验质量，不能把 prompt 的任意 A≤1 正放大对象自动视为 G w。

root `current_joint_budget_20261007.md` §1 的

\[
 F_n=4e^6\eta_0^{-1}\sqrt n+(1+v_*)+24C_h\sqrt n
       +2C_hK_{n,\alpha}/\epsilon_n
\]

在其引用的旧优先互补接口下口径通过。最后 **两份** C_hK/ε 分别支付完整输出高平均与逐原 soft 来源高强包络：先删高平均输出，再在其补集删高来源；平均高不等于每份来源高，所以不能只留一份费。强包络来源费替换旧 selected-Z 来源费，不能另加第三份。W 是完整输入质量，不重新归一化低-S 或热来源，也没有给每个历史／parent 副本重新领取 W。其余三项仅沿 root 引用的既有微盒、低径向层及 tiny-firstmark 等已证合同继承，本审没有重做这些分支。

root 当前 §3 的 (4) 与本审 (5) 完全一致：吸收乘积 (4096/49)(49/8192)=1/2，最终为

\[
 \lambda|E|\le4\mathcal B_{\rm paid}
 +\frac{8192}{49}\,[F_n+(k_n^4-1)J A_n]W
 +\frac{8192}{49}R_{\rm lowS,\ell<k_n^{-4}}.
\]

没有额外乘 C_h、HR 顶格的 2、旧 128/63 或 direct-output 4/3。若采用这份最新版条件账，不能再叠历史 HR/MM 整账；其缺失最新三个尺度原件的上游依赖仍保留。

作者 §4 的 512-child chain 与乘积上界只说明 source-tree 代数不能给额外 ancestor 衰减，未认证 actual FIRST/CPGP/far/history 反例；这一范围标注正确。作者三轮 55,674 项 Fraction 守卫与 root 独立保存字段重构是已经记录的有限代数证据。本审仅只读核作者终态范围及 root 通知的 receipt，不重跑、不重新重构，也不把它们提升为实际空间认证。最终接受新增可付子支及合并账的条件代数；完整 lowcoin∩lowS 的 actual 空间费仍未证。

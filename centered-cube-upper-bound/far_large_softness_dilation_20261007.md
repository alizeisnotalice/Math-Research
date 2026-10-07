# 大软度远源的 dilation：uniform Fisher、共同参数费与逐来源退出桥

2026-10-07。只新增本稿；不重新核算已有支付。本文给出直接解析的新工具推导，待 root 全文独审；不把它们视为一般实际远源已闭合。完整 FIRST、future cap、实际门、原来源一次归一化始终保留，删除门只在明确的正上界或障碍中注明。

## 1. 当前 residual 与原已付近全软带

沿用 `nonconcentrated_actual_far_interface_20261007.md` 的实际 `R_far,new`、`μ,W,q`、原 `σ∈[0,1]`、`L_s,R∈[a,b]`、`a≤b≤2a`、`n≥512`。当前早首跳时间 `0≤t<T_n=(12+log(n+2))/n<1/2`，`c=1−t`；当前 hard source 子测度质量 `≤W`，所有 sourcepair 门 `A_f≤1` 已合并唯一 LCA 与分数，不重复来源或父层数。

大软度分支指 `σ>S0=1/(16√n)`。原资格另有 `σ>1/n`，非集中微盒、固定 hard 短壳、首跳大 mark 与 selected-vertex-far 都保持。下述核工具不通过删除这些条件扩大实际输入类；付正支时可以用放松核上界。

原 tex 3225–3253 **已经支付** `E_T={σ>1−d_n}`，其中 `d_n=1000√(log(n+2)/n)`，在适用大维数域用完整 `(σ,L)` FTC 得 `O(√nlog(n+2))W`，其余维数有已登记回退。3255 后原账明确删除该输出带。故仅证明这份 near-full 带的 `m+√n` dilation 不是新主账支，不能重加。

最新上游三个尺度原件仍缺；如它们已经覆盖下面扩大的参数带，所得费用只能替换当前未付交集，不能再加历史整账。

## 2. 原 positive spectrum 与 soft continuation Fisher 常数

原 tex 4879–4896、5096–5103 给 `G_c` 正概率密度及其谱表示（只需 `1/2≤c≤1`）：
\[
G_c(x)=\int_1^\infty\rho_c(u^2)e^{-u|x|}du,
\quad
\rho_c(v)=\frac{cv}{[(1-c)\log(v-1)-cv]^2+(1-c)^2\pi^2}.
\tag{1}
\]
用 `c=1` 的端点读作 `ρ_1(v)=1/v`。下文只使用原真实谱，不用 finite quadrature 来定义 `G_c`。

### 2.1 解析率矩

对 `1<v<2`，`log(v−1)<0`，分母至少 `c²v²`，所以 `ρ_c(v)≤2/v`。对 `v≥2`，有 `0≤log(v−1)≤v/2`，且 `1−c≤c`；因此括号 `≤−cv/2`，分母至少 `c²v²/4`，得
\[
\boxed{\rho_c(v)\le8/v\quad(v>1,\ 1/2\le c\le1).}
\tag{2}
\]
`log(v−1)≤v/2` 可对差函数求导：其最小值在 `v=3`，为 `3/2−log2>0`。

令 `ell_u(x)=(u/2)e^{-u|x|}`，则 (1) 是 Laplace 概率混合
\[
G_c=\int ell_u\,d\Lambda_c(u),\quad
d\Lambda_c(u)=2\rho_c(u^2)u^{-1}du,\quad
\Lambda_c([1,\infty))=1,
\]
\[
\boxed{\mathbb E_{\Lambda_c}u
=2\int_1^\infty\rho_c(u^2)du\le16.}
\tag{3}
\]
概率质量一来自原 `G_c` 的归一化及 Tonelli。这项率一阶矩界在谱公式下已由 (2) 解析证明，不是诊断尾界的无证迁移。

### 2.2 单个均匀加 Laplace 核

记 `h=1_[−1/2,1/2]`、`b_u=h*ell_u`。对 `x≥0`，
\[
b_u(x)=\begin{cases}
1-e^{-u/2}\cosh(ux),&0\le x\le1/2,\\
e^{-ux}\sinh(u/2),&x\ge1/2.
\end{cases}
\tag{4}
\]
其 log-scale score 是 `a_u(x)=−1−x b_u'(x)/b_u(x)`。正核及指数尾保证分部积分 `∫a_u b_u=0`，并令
\[
J_u=\int x^2(b_u')^2/b_u,\quad
I(b_u)=\int a_u^2b_u=J_u-1.
\tag{5}
\]
在内部，`b_u≥(1−e^{-u})/2≥1/4`（`u≥1`、`e>2`），且 `|b_u'|≤(u/2)e^{-u(1/2−|x|)}`。所以
\[
J_{u,\rm in}\le u/4.
\]
外部直接积分得
\[
J_{u,\rm out}=(1-e^{-u})(u/4+1+2/u)\le u/4+3.
\]
于是
\[
\boxed{I(b_u)\le u/2+2.}
\tag{6}
\]

### 2.3 真实 `b_t=h*G_(1−t)` 的混合得分

固定本轮 `t`，置 `b_t=h*G_c=∫b_u dΛ_c`。dilation 不改变 `Λ_c`；混合 score 是 latent `u` 条件下 score 的 conditional mean。对平方用 Jensen，再 (3)、(6)，严格得到
\[
\boxed{\int b_t=1,\quad\int a_tb_t=0,\quad
I(b_t):=\int\left(-1-xb_t'/b_t\right)^2b_t\le10.}
\tag{7}
\]
这里只是核本身的概率混合，绝非任意来源 posterior 的坐标独立性。积分微分可先在有限谱截断、有限 `x` 上操作，再由上述正矩界及 weak derivative 极限；`b_t` 的一阶导数也由 `b_t'(x)=G_c(x+1/2)−G_c(x−1/2)` 直接给出。

## 3. 固定一个 mask 的 `m+√n` dilation 身份

真实 continuation 令
\[
r=\frac{1-\sigma}{1-t},\quad
F_{r,t,L}=[rh_L+(1-r)(b_t)_L]^{\otimes n}.
\tag{8}
\]
固定 active mask `A`，其硬／inactive 坐标数 `m=n−|A|`，完整正形状为 `M_A,t=h^⊗m⊗b_t^⊗(n−m)`（坐标排列保留原标签）。硬核 log-scale 导数是 `−h+(δ_−1/2+δ_1/2)/2`，不能删掉面原子。

因此 `D_logL M_A,t,L` 的 regular 部分 score 是 `Σ_active a_t−m`，均值 `−m`、方差 `≤10(n−m)`；正硬面总质量 `m`。由独立性**仅在这份固定 product 核**成立，
\[
\boxed{\|D_{\log L}M_{A,t,L}\|_{TV}
\le m+\sqrt{m^2+10(n-m)}
\le2m+\sqrt{10(n-m)}.}
\tag{9}
\]
导数总质量零，取正变差并积分物理 `logL` 可得
\[
\boxed{\int\sup_{a\le L\le b}M_{A,t,L}
\le1+\log(b/a)[m+\tfrac12\sqrt{10(n-m)}].}
\tag{10}
\]
这是固定 mask 的正列费；不能直接用它宣布实际 `σ(x)` 或 mask posterior 也只付该常数。

## 4. 共同 `r` 的选择费：一个明确峰值障碍

若先把每个 mask 系数都放宽至 `0≤r≤r0` 的最大值，再逐 mask 求和，正质量是
\[
D_n(r_0)=\sum_{m=0}^n\binom nm
\sup_{0\le r\le r_0}r^m(1-r)^{n-m}
\le1+\sqrt n\arcsin\sqrt{r_0}.
\tag{11}
\]
单个共同参数被每个 mask 各自最大化，这笔费用不能写成一。

还可严格看见其增长：对整数 `1≤m≤floor(nr0)`，取 `r=m/n`。Binomial(n,m/n) 的 `m` 是众数，方差 `≤m`；Chebyshev 给距 `m` 小于 `2√m` 的至少 `3/4` 质量，所含整数数至多 `5√m`。因此对应众数概率 `≥3/(20√m)`，得
\[
\boxed{D_n(r_0)\ge\frac3{20}\sqrt{\lfloor nr_0\rfloor}
\quad(nr_0\ge1).}
\tag{12}
\]
例如 `r0~n^{-1/2}` 仍有 `D_n(r0)≳n^{1/4}`。把 (10) 的 `√n` 项逐 mask 相加会收到 `√nD_n` 这类额外预算；这是该松弛／收费方法的障碍，**不是**原共同参数最大核的 `n^{3/4}` 下界，更不是 actual FIRST 反例。

## 5. 保留共同参数的两参数 FTC 列工具

这一节避免逐 mask 独立峰值。原 `T_n<1/16`：在 `n=512` 可用 `log514<10`，之后 `T_n` 递减。这里 `φ(1/2)>3/8` 也可自含核验：原 `w` 的正混合给
\[
\phi(1/2)=\int_0^1s(1-e^{-1/s})ds
=\tfrac12-\int_0^1s e^{-1/s}ds.
\]
对 `0<s≤1`，`e^{-1/s}≤s/e`（令 `y=1/s≥1`，`ye^{-y}≤e^{-1}`），故末积分 `≤1/(3e)<1/8`；最后可用 `e>65/24>8/3`。由原 resolvent `φ=c b_t+t b_t*w`、`b_t*w≤1`、`φ(1/2)>3/8`，在硬区间有
\[
5/16\le b_t(x)\le1\quad(|x|\le1/2),\qquad |b_t'|\le2.
\tag{13}
\]
上界和导数界用 `G_c≤w/c≤2`。只考虑共同 `r∈[0,r0]`、`r0≤1/2`，记一维 `j_r=rh+(1−r)b_t`、softness score `B_r=(h−b_t)/j_r` 与 regular scale score `a_r=−1−x(1−r)b_t'/j_r`。直接得到
\[
\mathbb E_{j_r}B_r=0,\quad\mathbb E B_r^2\le64,
\quad\mathbb E a_r=-r,\quad\operatorname{Var}(a_r)\le11,
\quad\int j_r|\partial_ra_r|\le7.
\tag{14}
\]
证明细节：内部 `j_r≥5/32`、`|B_r|≤32/5<7`，外部 `B_r=−1/(1−r)` 绝对值 `≤2`；`a_r` 是 latent 硬分量的 regular score `−1` 与 soft 分量 score `a_t` 的 conditional mean，latent 方差 `r(1−r)+(1−r)I(b_t)≤11`。最后 `∂_ra_r=x b_t'h/j_r²`，内部 `|x b_t'|≤1`、长为一，积分至多 `32/5<7`。

对 `F_r=j_r^⊗n`，保留完整面原子并展开 product 导数，得
\[
\|\partial_rF_r\|_1\le8\sqrt n,
\quad\|D_{\log L}F_r\|_{TV}\le nr+\sqrt{n^2r^2+11n},
\]
\[
\boxed{\|\partial_rD_{\log L}F_r\|_{TV}
\le8\sqrt n\sqrt{n^2r^2+11n}+8n+8nr\sqrt n.}
\tag{15}
\]
regular 混合导数为 `(Σa_r)(ΣB_r)F_r+(Σ∂_ra_r)F_r`，前项用 Cauchy–Schwarz。硬面导数每面贡献系数导数质量一及切向 product 导数 `≤8r√n`，全部 `n` 面质量保留。

从 `(r,L)=(0,a)` 应用两参数 FTC／BV，定义
\[
\begin{aligned}
J_n(r_0)={}&1+8r_0\sqrt n+\log(b/a)\bigl[\sqrt{10n}\\
&+r_0\{8\sqrt n\sqrt{n^2r_0^2+11n}+8n+8nr_0\sqrt n\}\bigr].
\end{aligned}
\tag{16}
\]
严格列工具是
\[
\boxed{\int\sup_{0\le r\le r_0,\ a\le L\le b}F_{r,t,L}\le J_n(r_0).}
\tag{17}
\]
闭核面采用 BV 的相应上跳规范；`r` derivative 是普通密度，`L` 与混合 derivative 的面原子通过 TV 控制。这里已经支付共同 `r,L` 的选择，不另乘 (11)。常数对每个固定 `t<T_n` 相同；没有取 `sup_t`。

## 6. 当前 residual 上可用的参数带费用和旧支重叠

实际映射仍是 `r=(1−σ)/(1−t)`。若
\[
\sigma\ge1-\delta_n,\qquad
\delta_n=(1-T_n)r_0,\quad r_0\le1/2,
\tag{18}
\]
则对全部原早首跳时间 `t<T_n` 都有 `r≤r0`。固定首跳参考 (上一稿 (14)) 的总质量 `≤2W`，原 hard／soft 行因子 `u/q≤C_h`，遂由 (17) source-once 支付
\[
\boxed{R_{\rm near-full\ continuation}\le2C_hJ_n(r_0)W.}
\tag{19}
\]
实际门、非集中、短壳、far 等保持在起点，之后只为正支配删掉。不是按每个 mask 或每个输出重领一份 `W`。

取 `r0=L_n/√n≤1/2`，则 (16) 给
\[
J_n\le1+8L_n+\log2\,\sqrt n
[\sqrt{10}+8L_n\sqrt{L_n^2+11}+8L_n+8L_n^2].
\tag{20}
\]
所以 `L_n=n^{o(1)}` 时费用是 `√n n^{o(1)}W`。可用 `r0=min(1/2,L_n/√n)`，公式对所有维数照常成立。

这只有当 `δ_n>d_n` 且当前未付交通尚含该扩展带时才新增输出支付。旧 `L_n≈1000√log(n+2)` 的近全软区已付，不能重新加。若选择 `L_n=log²(n+2)`，它**并非在所有维数**超过旧常数 1000 带；仅在 `log^{3/2}(n+2)>1000/(1−T_n)` 且未触发 `r0=1/2` cap 时，扩展才真超过旧带。选择其它 subpower `L_n` 也必须逐范围比较，不凭渐近口号宣称当前维数已有新增行。

还有一个范围事实必须明确：原 tex 3225–3253 的三项 `(σ,logR)` 导数界本来就对全部 `σ≥1/2` 成立，其 FTC 中的 `d_n` 可直接换成任意 `0≤δ≤1/2`。因此原同源 future cap 已给整个输出带 `E_δ={σ≥1−δ}`
\[
\begin{aligned}
q|E_\delta|\le \mathcal T_n(\delta)W,\quad
\mathcal T_n(\delta)={}&1+10\sqrt n\,\delta+\log2\,[2n\delta+\sqrt{28n}\\
&+\delta\{10\sqrt n\sqrt{28n+n^2\delta^2}+6n+10n\delta\sqrt n\}].
\end{aligned}
\]
代 `δ=L_n/√n≤1/2` 同样得到 `√n n^{o(1)}W`。这是原公式的通用参数推论，完整输出可直接用它支付，无需先改成 continuation。因此“渐近扩大 near-full 带”**不是本轮新主账进展**。本轮新工具是统一 `t` 下的 (7)、固定 mask 的 (9) 及共同 `r` continuation 的 (17)；它们可用于以后需要 continuation 结构的比较，但不能重复登记旧整行工具的后果。

因此 (19) 是一般输入的明确可用 continuation 参数带工具；它不是 `σ>S0` 整支的支付。其余 `S0<σ<1−δ_n` 还没有来源一次费。

## 7. 全部剩余 `σ>1/n` 的逐来源退出证书

读取 `far_source_geometry_bridge_20261007.md` 的辅助噪声定义：对 hard 来源 `z∈Q(x,R)`，用**辅助**原噪声 `N_{σ,R}=[(1−σ)δ0+σw_R]^⊗n` 定义退出概率 `e_x(z)`。它不是原 `L_s` firstjump 的同一历史。

原 tex 2095–2097 已有 `φ(0)<3/4`；故每个内部坐标 `ψσ(v)=1−σ+σφ(v)≤1−σ/4`。包括闭面、任意 hard 来源和完整混合输入，
\[
\boxed{e_x(z)=1-\prod_{j=1}^n\psi_\sigma((x_j-z_j)/R)
\ge1-(1-\sigma/4)^n
\ge\frac{n\sigma}{4+n\sigma}>\frac15.}
\tag{21}
\]
最后只用原资格 `σ>1/n`。对大软度 `σ>S0`、`n≥512`，有 `nσ>4/3`，所以还可用 `e>1/4`。无需任意来源 posterior 独立。

因此任意实际 `0≤g(x,z)≤1` 及硬子源 `ν_h≤μ` 都满足
\[
\int g e\,dr_x\ge\tfrac15\int g\,dr_x,
\quad
R_{\rm far}\le5\int h_{R(x)}(x-z)g(x,z)e_x(z)dx\,d\nu_h(z).
\tag{22}
\]
这解除“完整行退出均值不能绑定低质量 gate”的问题：这里是逐 `z` 下界，绝非以 `∫e dr>1/2` 推选中门也有同样平均。原 complete future cap 仍给合法辅助同源行证书 `P_{σ,R}*μ≤q`，但 (21) 的 uniform 下界更强于仅作门质量分层。

## 8. 退出不等于 source-once 空间费：最小缺口

仍不能把 (22) 右侧直接支付。内部 auxiliary exit 核是
\[
h_R(x-z)e_x(z)=\mathbf1_{Q(x,R)}(z)
[h_R(x-z)-P_{\sigma,R}(x-z)].
\tag{23}
\]
它包含完整 hard 面的上跳，不是只含 `b_t` soft scale score 的平滑核。对删除 FIRST/futurecap 与实际门后的自由 fixed-source 列，(21) 给
\[
\int\sup_{a\le R\le b}h_R(x-z)e_{\sigma,R}(x,z)dx
\ge\tfrac15[1+n\log(b/a)]
\tag{24}
\]
（可固定一个 `σ>1/n`）。所以不能将 (17) 的 near-full smooth continuation 费用挪给这个 hard exit 最大列。这只否定删掉共同实际资格后的强 Schur 方案，不是实际输入／hardband／FIRST 反例。

真正需证的是在同一原 `μ`、actual `g`、连续或合法原物理 winner 和完整 future cap 下，对 (22) 的**门控**退出交通建立 `√n n^{o(1)}W` 列费，或者有偿比较它与 (19) 已覆盖的共同参数带及其补集。uniform exit 可补权，却没有让辅助噪声退出测度变成原 first-exit 参考来源；也不能把其条件化来源重新当作原 `μ` 来测试 future cap。

当前严格结果是 (7)、(9)–(10)、保留共同参数的 (17)–(20)，及适用于所有剩余软度的 (21)–(22)。已有近全软支付不重复，扩展参数带的新增性必须核实际范围，余下大软度远源列仍开放。没有把未找到费当反证。

## 9. 三轮常数守卫与数值诊断范围

独占脚本 `far_large_softness_fisher_guard_20261007.py` 和同名 `_results.json` 已实际运行，三轮 `n=512,1024,4096` 全部 PASS，无随机种子、无随机采样。纯有理部分逐轮核 `T_n<1/16` 的安全上界、自含 `φ(1/2)>3/8`、(13) 的 `5/16` 下界、score 常数、五个固定 mask 的 TV 三角预算、三个带截断 `r0≤1/2` 的实际 `r=(1−σ)/(1−t)` 映射和共同参数 FTC 安全有理上界。另用精确 Binomial 众数质量核 (12) 的平方形式，未以浮点概率认证峰值下界。

独立标为浮点诊断的部分共九个 `u`：三轮依次 `1,2,8`、`1.5,4,16`、`3,32,64`。它以 Laplace CDF 差核对 (4)，以内部边界层变量的 `8000/16000` 双 Simpson 网格分别计算得分平方与 `J`，核质量一、得分均值零、`I=J−1`、`J_in≤u/4` 及 (6)。最高 `u=64` 的初始 `4000/8000` 双网格差超过预设 `10^{-8}`，故统一加密后才登记 PASS；未通过放宽阈值掩盖网格误差。

这些只核解析常数、离散系数和标量核实现。谱率矩与真实 `I(b_t)≤10` 的证书是 §2 的解析推导；脚本没有将有限谱 quadrature 当 `G_c`，没有测试原实际 FIRST 或共同输入，更没有证明 whole-large-σ 空间列费。

# 一般非集中输出后的实际远源接口：固定短壳、微小首跳与有限顶点收费

2026-10-07。root 已全文核对第 2、4–6 节及原 tex 4870–4902、5515–5585 的相关核合同，通过本文限定范围的解析独审。主账登记采用“仅在当前未付交通取依次互补交集”，不重复加历史整账。最新上游三尺度原件尚未取得，具体接口一致性仍单列依赖。本文另有三轮精确离散系数守卫，不模拟原核或实际 FIRST 样本；一般 geom 主目标仍未闭合。

## 1. 一般输入已经取得的非集中合同与实际交通

原完整正输入 `μ`、`W=μ(ℝ^n)>0`、窗口 `a≤b≤2a`，`n≥512`。原可测硬赢家 `R(x)`、软度 `σ(x)∈[0,1]`、物理尺度 `L_s(x)∈[a,b]`、共同 FIRST 阈值 `q∈[3λ/4,λ]` 均冻结。若原接口包括 `σ=1`，使用原核的 `σ↑1` 端点版本；(11) 的系数此时按 `p=1` 理解，后文形状包络也覆盖该端点。不会把有限时域停止合同未经核验扩大到 `σ=1`。
\[
u(x)=R(x)^{-n}\mu(Q(x,R(x))),\quad
E=\{2\lambda<u(x)\le4\lambda\},\qquad C_h=16/3.
\tag{1}
\]
原未截断 FIRST 给 `P_{σ,L_s}*μ(x)=q`。全部 future cap、GOOD、score、owntrace、出生、森林、唯一 LCA、共同重捕获、原分数、严格 CP/GP 与删除门保留。原早期第一跳即退出子核完整保留后续传播。

令 `η0=49/65536`、`η_n=η0/√n`、`h=2a/n`。root 与 tensor 已核验的 greedy witness cover 在任意输入中直接支付输出 `E_paid`，
\[
\int_{E_{\rm paid}}u\le\frac{4(1+3h/a)^n}{\eta_n}W
\le\frac{4e^6}{\eta_0}\sqrt n\,W.
\tag{2}
\]
互补 `E_res` 对**每个**闭轴盒 `B`、全边长 `h` 满足
\[
\boxed{\mu(B\cap Q(x,R(x)))\le\eta_n\,m(x),
\qquad m(x)=\mu(Q(x,R(x))),\quad x\in E_{\rm res}.}
\tag{3}
\]
输出存在任意中心盒的事件按已核 greedy 接口的 Lebesgue 可测规范处理。这个合同不是特殊来源分区假设，而是 (2) 支付后的确定性全中心约束。所有完整来源、未分配来源均继续留在 `m,u,E,R` 中。

完整 soft 与 hard 行概率为
\[
p_x(dy)=q^{-1}P_{\sigma,L_s}(x-y)d\mu(y),\quad
r_x(dz)=m(x)^{-1}\mathbf1_{Q(x,R(x))}(z)d\mu(z).
\tag{4}
\]
实际 high／birth 端是这两者的子概率，未重新归一。原唯一 LCA 与全部分数已合并为 `0≤A≤1`；早期核比值可以留在 `A`，或在下文显式展开早期核，但不能两次计入。

基本近源子支是 `||y−z||∞≤h/2`。由 (3)、soft 行质量一而得
\[
R_{\rm near}\le\eta_n\int_{E_{\rm res}}u
\le4\eta_n\lambda|E_{\rm res}|.
\tag{5}
\]
该证明不假定两个实际后验独立；它是对每个固定 soft 来源 `y` 的 hard 条件质量上界，再用 `A≤1`。

## 2. 优先一般约化：连续赢家给固定短壳

记 `β(x)=nlog(R(x)/a)`、`s(x,z)=nlog(2||x−z||∞/a)`，硬捕获意味着迟延 `d=β−s≥0`。原 prompt 的 `G0=sup_{a≤R≤b}h_R*μ` 为连续窗口，实际硬赢家最大于这个窗口。取任意固定 `v>0`。

### 2.1 低 `β` 输出 source-once 支付

在 `β<v` 的输出，半径只在 `[a,min(b,a e^{v/n})]`。完整 actual source-once 合并给 `0≤g(x,z)≤1`、来源 `ν_h≤μ`；固定列包络严格积分为
\[
\int\sup_{a\le R\le c}h_R(x-z)dx=1+n\log(c/a).
\]
故包括全部实际门的这批输出交通有
\[
\boxed{R_{\beta<v}\le[1+\min(B,v)]\,\nu_h(\mathbb R^n)
\le(1+v)W,\quad B=n\log(b/a).}
\tag{6}
\]
不要求连续赢家最大，只有输出半径范围；不是每个输出再取得一份 `W`。

### 2.2 `β≥v` 输出的长迟延吸收

此时内缩半径 `R_-=R e^{-v/n}≥a` 仍在连续允许窗口。真实最大条件给
\[
\mu(Q_{\rm open}(x,R_-))
\le\mu(Q(x,R_-))\le u(x)R_-^n=e^{-v}m(x).
\tag{7}
\]
事件 `d>v` 精确对应开内 cube，故
\[
\boxed{R_{\beta\ge v,\ d>v}
\le e^{-v}\int_{E_{\rm res}}u
\le4e^{-v}\lambda|E_{\rm res}|.}
\tag{8}
\]
用的是完整 `μ`；high／birth 子源及全部门只缩小交通。`d=v` 留短壳，`β=v` 的内半径恰为 `a`，两种等号均合法。

取
\[
v_*=\log(65536/49),\qquad4e^{-v_*}=49/16384.
\tag{9}
\]
于是剩余真正硬来源满足 `β≥v_*` 和 `0≤d≤v_*`，固定厚度与 `n` 无关。无需再令 `v~log n` 来支付空间全窗。

**有限尺度限制。** (7) 在原 `𝒥` 含 `R_-` 的行也成立；任意有限 `𝒥` 不含这个内缩尺度时不成立，不能由有限集合的赢家冒充连续赢家。该情形应保留先前 `(1+B)e^{-v}W` 长迟延费用，或另证传递。主 prompt 的连续窗口资格允许本节 (8)。

## 3. 完整首跳与逐时间参考形状的明确接口

原初始格分配 `Σ_Cν_C≤μ_hi≤μ`。记 `T_n=(12+log(n+2))/n<1/2`、`c_t=1−t`。对首跳坐标 `i`、mark `r`、首落点 `w=y+r e_i`，原早期 first-exit 测度为
\[
c_t^{n-1}\mathbf1_{\{t<\min(\sigma,T_n),\ w\notin C\}}
G_{c_t,L_s}(r)dt\,dr\,d\nu_C(y).
\tag{10}
\]
首跳前无事件率、每坐标跳率与后续核正是原非齐次过程，没有用条件后验 `1/n`。

完整 continuation 满足
\[
S_{\sigma,t,L}=h_L^{\otimes n}*K_{\sigma,t,L}
=\sum_{A\subseteq[n]}p^{|A|}(1-p)^{n-|A|}M_{A,t,L},
\quad p=\frac{\sigma-t}{1-t},
\]
\[
M_{A,t,L}(x)=\prod_{j\in A}(h*G_{c_t})_L(x_j)
\prod_{j\notin A}h_L(x_j),\qquad\int M_{A,t,L}=1.
\tag{11}
\]
这是原真实 transition 的正掩码分解；`A` 是首跳之后至少发生一次跳的坐标集合。它在**固定参数、尚未端点条件化**时为 Bernoulli 权，端点 posterior 不再独立。下文使用正系数包络而非把后验当 Bernoulli。

所用固定尺度形状工具的出处是原 `概率接口阶段证明.tex` 5522–5582：
\[
N_h=1+n\log(b/a)\le n,\quad
D_n=\sum_{k=0}^n\binom nk\sup_{0\le p\le1}p^k(1-p)^{n-k}
\le1+(\pi/2)\sqrt n\le2\sqrt n.
\tag{12}
\]
对固定 `t`，定义完整 Bernoulli 系数峰值形状 `H_t=Σ_A c_|A| M_A,t,1`，则 `∫H_t=D_n`。每份因子偶、坐标单调，混合的 level set 为 star-shaped。已核 dilation layercake 因而给
\[
\int\sup_{\sigma\in[t,1],\ a\le L\le b}
S_{\sigma,t,L}(x)dx\le N_hD_n.
\tag{13}
\]
**逐时间性。** `G_{c_t}` 和包络可以随原 `t` 改变；(13) 对每个固定 `t` 有同一质量常数。随后在已固定的参考首跳测度上积分 `t`。没有额外取 `sup_t H_t`，也没有把 `N_hD_n` 移用到包含不同 `t` 的联合 sup。

固定参考 first-exit 测度取原已核比较 `G_{c_t,L}≤2G_{c_t,b}`：
\[
d\bar\kappa_{C,i}=2c_t^{n-1}\mathbf1_{\{0<t<T_n,\ w\notin C\}}
G_{c_t,b}(r)dt\,dr\,d\nu_C(y),\quad w=y+r e_i.
\tag{14}
\]
对所有 `L_s(x)∈[a,b]` 同时正支配，但只比较一个首跳坐标，没有 `2^n`。固定 `(C,i,t,r,y)` 后 `w` 与空间输出 `x` 无关；这是 (13) 平移并积分输出的合法位置。它不是为每个移动赢家重复制一份来源。

## 4. 新可支付子支：微小外部首跳 mark

原软核 `w(r)=∫_0^1 exp(−|r|/s)ds` 满足 `w(0)=1`。已核 resolvent 恒等式
\[
w=cG_c+(1-c)G_c*w
\]
给 `G_c≤w/c≤2w≤2` 于 `c≥1/2`；出处另见原 tex 4879–4896 与 5744。因此 `t<T_n` 时
\[
\int_{|r|\le a/n}G_{c_t,b}(r)dr
\le\frac{4a}{nb}\le\frac4n.
\tag{15}
\]
保留 first-exit 条件并限制这份 mark，(14) 的来源总质量有
\[
\sum_{C,i}\|\bar\kappa_{C,i}^{\rm tiny}\|
\le\frac8n\,n\int_0^{T_n}c_t^{n-1}dt\,W
\le\frac8nW.
\tag{16}
\]
实际 tiny mark 交通先保留全部门，再删除门作正支配；完整 hard 行积分至多 `u`，用唯一 `u/q≤C_h`，逐固定参考项由 (13) 支付 continuation 的 `σ,L` 选择。由 (16) 得
\[
\boxed{R_{\rm tiny}\le\frac{8C_hN_hD_n}{n}W
\le16C_h\sqrt n\,W.}
\tag{17}
\]
这批首跳仍可真实退出原格：源可靠近原格面；它不是已经被删除的“首次跳仍在格内”交通。也不要求 `σ≤n^{-3/2}`，不会复活旧空支。它只删除当前 residual 中 `|r|≤a/n` 的部分；是否已被更新上游三尺度覆盖须另核，不能把两套费用叠加。

## 5. 非空小软度区的真实 active-coordinate 尾

令
\[
S_0=\frac1{16\sqrt n},\quad K=\lfloor\sqrt n\rfloor-2.
\tag{18}
\]
原 hardwinner 与同 `σ` 全物理 future cap 已严格迫使 `σ>1/n`（原 tex 5788–5817）。`1/n<S0` 当且仅当 `n>256`，故在原 `n≥512` 范围内，区间 `1/n<σ≤S0` **未被该必要资格排空**；这里仍不声称存在完整 FIRST 合法输入实现每个参数。

在此区域，(11) 给 `0≤p=(σ−t)/(1−t)≤σ≤S0=p0`。当 `k≥K+1`，有 `k>np0=√n/16`，故 `p^k(1−p)^(n−k)` 在 `[0,p0]` 单调增加。于是多数 active 坐标的正 continuation 子核被
\[
H_{t,\rm many}=\sum_{|A|\ge K+1}
p_0^{|A|}(1-p_0)^{n-|A|}M_{A,t,1}
\]
支配，其质量严格为
\[
\tau_n=\Pr\{\operatorname{Bin}(n,p_0)\ge K+1\}
\le2^{-(K+1)}(1+p_0)^n
\le4\exp[-(\log2-1/16)\sqrt n].
\tag{19}
\]
这不是对条件后验使用 Chernoff；它是固定正形状包络的系数质量恒等式与纯代数估计。相同 star-shaped dilation 只扩大物理 `L`，对固定 `t` 给 `∫sup_L H_t,many,L≤N_hτ_n`。

(14) 求和首跳参考质量至多 `2W`。因此完整 many-active 子交通为
\[
\boxed{R_{\sigma\le S_0,\ |A|>K}
\le2C_hN_h\tau_nW
\le8C_h n e^{-(\log2-1/16)\sqrt n}W
\le8C_h\sqrt n\,W.}
\tag{20}
\]
最后用 `c=log2−1/16>1/2`、`z e^{-cz}≤1`。来源、首跳时间、坐标与 mark 各求和一次；没有掩码数 `2^n` 的额外费用，因为其系数和已完整计入 `τ_n`。

本节不支付 `σ>S0`；原 `σ>1/n` 与完整 future cap 始终保留。式 (20) 只对这份明确 `|A|>K` 正子核收费，不把“少坐标”后验概率认证成独立。

## 6. 新恒等与吸收：完整软 posterior 上的一条辅助折线

将 (10)–(11) 进一步展开为正历史参数 `ξ_j∼G_{c_t,L}`（`j∈A`）及初始 cube 偏移。具体未乘实际 pair 门前的 soft 子后验为
\[
\begin{aligned}
d\Pi_x={}&q^{-1}c_t^{n-1}\mathbf1_{\{t<\min(\sigma(x),T_n),\ w\notin C\}}
G_{c_t,L_s(x)}(r)dt\,dr\,d\nu_C(y)\\
&\times p^{|A|}(1-p)^{n-|A|}\prod_{j\in A}G_{c_t,L_s(x)}(\xi_j)d\xi_j
\;h_{L_s(x)}^{\otimes n}\!\left(x-w-\sum_{j\in A}\xi_j e_j\right).
\end{aligned}
\tag{21}
\]
求和原 `C,i,A` 并积分参数，恢复真实早期 `Q_x^early/q`；其总质量 `≤P_{σ,L_s}*μ_hi(x)/q≤1`。所有 soft 参数在给定 `x` 后按原目录固定，
`U=x−w−Σξ_j e_j` 是同一次初始 cube 偏移，未重新随机平均。原 sourcepair 门 `A_f(x,y,z,θ)` 原样乘入 (21) 与 hard 子概率；若实现存在额外路径门，保留其原 backward projection 的 `[0,1]` 密度，不能凭独立性删除它。

对一次 (21) 的 realization，按**一个预先固定坐标次序**排列 `A={j1<…<jk}`，只定义一条辅助 coordinate-aggregate 折线的顶点
\[
V_0=y,\quad V_1=w,\quad
V_{\ell+1}=w+\sum_{r=1}^{\ell}\xi_{j_r}e_{j_r}\quad(1\le\ell\le k).
\tag{22}
\]
这是原真实 transition 的 endpoint 参数化与 bookkeeping；它不是实际跳时顺序，也没有把重新排序的顶点当停时。**只选这 `k+2` 个点，绝不遍历 `2^k` 个坐标子集顶点。** 其它顺序或其它子集产生的点没有被本支支付。

更一般地，选择任何不依赖硬源 `z` 的可测顶点清单 `𝒱`，长度至多 `M`，可依赖原 `x` 和全部 soft realization。由 (3) 对每个固定 realization 的所有中心适用，
\[
r_x\left(\bigcup_{V\in\mathcal V}Q(V,h)\right)\le M\eta_n.
\tag{23}
\]
再用 `Π_x` 总质量 `≤1`，原 actual pair 门 `≤1`，硬 high／birth 子源 `≤r_x`，得到
\[
\boxed{R_{\rm vertex-near}\le M\eta_n\int_{E_{\rm res}}u
\le4M\eta_n\lambda|E_{\rm res}|.}
\tag{24}
\]
没有 soft／hard 后验独立假设；逐 realization 的确定性 hard 条件上界足够。该精确子支比原 `y` 邻近支允许首落点及经过完整后续传播后的 aggregate 顶点；仍不能支付整条折线、任意子集顶点或连续路径 tube。

在 `σ≤S0,|A|≤K` 中取 (22) 全部顶点，则 `M=k+2≤floor√n`。在其它仍未付区域，可只取 `y,w,w+Σξ_j e_j` 三点（`3≤floor√n`）。统一有
\[
\boxed{R_{\rm selected-vertex-near}\le4\eta_0\lambda|E_{\rm res}|
=49\lambda|E_{\rm res}|/16384.}
\tag{25}
\]
它已包括原 near-source `V0=y` 交通，故**替换** (5)，不能与 (5) 再相加。保留其它 pointwise path 门只会缩小该交通。若某实现没有提供 (21) 的额外 path 门 projection，则该额外历史门的确切参数化须另补；当前原 prompt 的 sourcepair 门位于 `Q` 外，(21) 正好适用。

## 7. 组合账与最小剩余接口

采用连续原 prompt 窗口，可以在当前 residual 内按顺序切分：greedy paid 输出 (2)；低 `β` 输出 (6)；长硬迟延 (8)；tiny 首跳 mark (17)；小软度多数 active 坐标 (20)；所选顶点邻近 (25)。每步只取前步补集。两项吸收 (8)、(25) 总系数正好
\[
49/16384+49/16384=49/8192.
\tag{26}
\]
条件登记的新交通总支配因此为
\[
R_{\rm current}\le\left[\frac{4e^6}{\eta_0}\sqrt n+(1+v_*)
+24C_h\sqrt n\right]W
+\frac{49}{8192}\lambda|E|+R_{\rm far,new}.
\tag{27}
\]
root 已独审第 4–6 节的解析推导。式 (27) 登记时只在当前未付交通取上述依次互补交集，不重复加历史整账；最新上游三个尺度原件尚缺，其与当前实际接口的一致性仍作为独立依赖，不凭本稿补造。该账保留 `R_far,new`，因此不是目标闭合。

`R_far,new` 保留全部 actual 门和 source-once 来源，满足非集中 (3)、`σ>1/n`、`β≥v_*`、固定硬短迟延 `0≤d≤v_*`、首跳真实退出且 `|r|>a/n`。其 `σ≤S0` 部分只有 `|A|≤K`，且硬源远离 (22) **这一条**辅助折线的全部所选顶点；其 `σ>S0` 部分只保证远离三个所选点。它可能接近其它子集顶点、接近折线中点，或具有高度条件相关的 first coordinate 与硬面，均仍未支付。

真正缺口是对这份保留 complete future cap 的短壳远源交通给空间列费；微盒非集中、少 active 坐标或端点 probability 单独没有推出该费。尤其不能从 `σ>1/n` 猜出 posterior `p_i≤C/(nσ)`，不能把最初共同输入的 future cap 当作条件化停止来源的上限。

## 8. 有限原子搜索不能认证这份残余

若完整来源只有 `N` 个正质量原子，任意非空赢家捕获至多 `N` 个原子，至少一个原子捕获份额 `≥1/N`。每个原子都可置于边长 `h` 的闭盒中，所以只要 `η_n<1/N`，该输出违反 (3)，已被 greedy paid 输出删除。原保存搜索 `N≤128`，在 `η0=49/65536` 就有 `η0<1/128`，而这里 `η_n≤η0` 更小；因此那些搜索不提供 `E_res` 或 `R_far,new` 的非空数值证据。

一般非原子输入或足够多、小份额原子仍可能留下残余。本稿没有执行新实例、没有以未找到或未证当反证。所有新增一般子支的费用都明确写出；协同实际远源仍是待付项。

## 9. 三轮精确离散系数守卫

独占 `actual_far_branch_exact_guard_20261007.py` 已执行并保存同名前缀 `_results.json`，状态 `PASS_EXACT_ACTUAL_FAR_BRANCH_COMPONENTS`。三轮 `n=1024,4096,16384`，`√n=32,64,128`；没有随机种子、没有原 `G_c` 数值模拟或实际 FIRST 样本。

- `p0=1/(16√n)=1/D`、`K=√n−2`，以公共整数分母 `D^n` 精确计算 `τ_n`：从 `D^n` 减去 `k=0,…,K` 的 binomial lower-complement 分子和。每步整数 recurrence 的除法余数均为零。
- 精确比较 `τ_n≤2^{−(K+1)}(1+p0)^n` 及 `nτ_n≤4√n`，全部通过；从而核验第 5 节选用的离散 tail 常数。
- 每轮 `k=K+1,2K,n`、`p=p0/4,p0/2,p0`，共 9 个系数链，以同一整数分母 `(4D)^n` 检查单调；另核 `kD−n>0`，明确系数峰值在受限参数区间右端。这里不是后验独立性检查。
- 核 `1/n<S0` 与 `n>256` 的必要软度区间条件；核 `a/b=1,3/4,1/2` 下 reference tiny-mass 代数上界 `8/n`；核 `M=3,√n−1,K+2` 时 `Mη_n≤η0`。完整来源归一化保持原 `W`，`M=K+2=√n` 的吸收系数精确为 `49/16384`。

很大有理数仅输出完整整数 numerator／denominator 的 unsigned big-endian SHA256、bit length 与诊断用 `log10` 近似；所有 PASS 比较实际使用完整整数，近似值不参与判定。shape envelope、resolvent 与实际门接口仍由解析证明及 root 原文核对认证，不能由这份离散结果代替。新 guard 不认证最小剩余 `R_far,new` 的预算或非空性。

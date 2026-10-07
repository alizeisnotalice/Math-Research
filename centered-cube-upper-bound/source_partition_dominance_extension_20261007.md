# 来源自适应小直径分配的支配输出：固定微格覆盖的范围强化

2026-10-07。本稿仅独审覆盖 lemma 的一般范围及现有 packet 字段，不新增实验、不改已有稿。结论是一般来源分配的 `O(W/η)` 输出支，与 `Γ_1` 的相应能量费；它不证明所有赢家都属于这支，也不关闭协同残余。

## 1. 合法来源合同：可重叠盒、不可重复质量

冻结原有限正 Borel 输入 `μ`、质量 `W>0`、窗口 `R(x)∈[a,b]`，以及原完整响应
\[
m(x)=\mu(Q(x,R(x))),\qquad u(x)=m(x)/R(x)^n,
\qquad E=\{2\lambda<u(x)\le4\lambda\}.
\tag{1}
\]
不要求 `R` 最大；若原接口给真实赢家，则保持原赢家不改。

在输出积分之前，选定至多可数的非负来源分配
\[
\mu_k=\alpha_k\mu,\quad \alpha_k\ge0,\quad
\sum_k\alpha_k\le1\quad\mu\text{-a.e.},
\]
\[
\mu_{\rm cov}=\sum_k\mu_k\le\mu,\qquad
W_{\rm cov}=\sum_kM_k\le W,\quad M_k=\mu_k(\mathbb R^n).
\tag{2}
\]
这允许分数分配；测度合同 `Σμ_k≤μ` 与 (2) 等价。每个正质量 `μ_k` 均集中于固定轴盒
\[
B_k=\prod_{j=1}^n[c_{k,j}-\ell_{k,j}/2,c_{k,j}+\ell_{k,j}/2],
\quad 0\le\ell_{k,j}\le a/n.
\tag{3}
\]
盒可以重叠、不来自同一网格；可以根据完整输入 `μ` 预先选择。关键是 `μ_k,B_k` 不随正在积分的接收点 `x` 重新分配来源。单纯每个输出重新挑一份小包、再各给它原质量，不满足 (2)。零质量项可忽略。`ℓ∞` 直径 `≤a/n` 的来源集合可取各坐标上下确界的轴盒，满足 (3)。

未分配来源为 `μ_0=μ−μ_cov≥0`；始终保留在原完整 `m,u,E,R` 与 profile 的来源积分中。不是把 `μ_cov` 归一化成一个新输入。

## 2. 支配判据与严格费用

固定 `η∈(0,1]`，定义
\[
E_{\rm dom}=\{x\in E:\exists k,
\mu_k(Q(x,R(x)))\ge\eta\,\mu(Q(x,R(x)))\}.
\tag{4}
\]
**分母是完整原捕获质量 `m(x)`。** 不能在部分覆盖时直接改成 `μ_cov(Q)`，否则原 hardband 下界不足以控制覆盖半径。若另已证明 `μ_cov(Q)≥ζm`，以 covered 捕获质量为分母的 `η` 支配只得到本稿参数 `ηζ`，费用为 `O(W_cov/(ηζ))`。

若 `x` 由 `k` 支配，则
\[
M_k\ge\mu_k(Q(x,R))\ge\eta m(x)>2\eta\lambda R^n.
\tag{5}
\]
捕获正质量保证 `Q(x,R)` 与集中盒 `B_k` 有交点，因而
\[
\|x-c_k\|_\infty\le R/2+a/(2n)
\le(1+1/n)R/2.
\tag{6}
\]
于是与半径严格上界合并得到
\[
E_{\rm dom}\subseteq\bigcup_{k:M_k>0}
Q\left(c_k,(1+1/n)(M_k/(2\eta\lambda))^{1/n}\right).
\tag{7}
\]
这是一份对输出的可数覆盖；盒相交不影响上界。体积求和用每个 `M_k` **一次**，由 `ΣM_k=W_cov` 得
\[
\boxed{\lambda|E_{\rm dom}|
\le\frac{(1+1/n)^n}{2\eta}W_{\rm cov}
\le\frac e{2\eta}W_{\rm cov},
\qquad\int_{E_{\rm dom}}u(x)dx
\le\frac{2e}{\eta}W_{\rm cov}.}
\tag{8}
\]
因此完整原 hard 交通，或任何门控且来源 `ν≤μ` 的实际子交通，在这批输出上都由 (8) 支付；没有要求这份子交通自身只使用 `μ_cov`。在原 residual 中只取输出交集，避免把已付交通再相加。

`≥η` 定义与严格 `>η` 都可用同一费用。前者的协同补集满足所有份额 `<η`，后者的补集满足 `≤η`。数值稿采用后者、将等号留在协同支；回代时必须沿用所选规范，不能重叠。

## 3. `Γ_1` 能量与保留未覆盖来源的残余

取 `β(x)=nlog(R(x)/a)`、`B=nlog(b/a)`，令 `Γ_dom` 是原 `Γ_1` 只将输出事件换成 `E_dom` 后的 profile，仍以完整 `μ/W` 积分来源。记 `s(x,z)=nlog(2||x−z||∞/a)`。已核 cone Jacobian 精确给
\[
\int_0^B\Gamma_{\rm dom}(t)dt
=\frac1W\iint\mathbf1_{E_{\rm dom}}(x)
\mathbf1_{\{0\le s\le B,\ 0\le\beta(x)-s\le1\}}
e^{\beta(x)-s}h_{R(x)}(x-z)dx\,d\mu(z).
\tag{9}
\]
逐来源 `dx` 的对角例外为零，Tonelli 适用于有限原子或任意有限 Borel 输入。由于 `Γ_dom≤1`，(8) 得
\[
\boxed{\int_0^B\Gamma_{\rm dom}(t)^2dt
\le\frac eW\int_{E_{\rm dom}}u(x)dx
\le\frac{2e^2}{\eta}\frac{W_{\rm cov}}W
\le\frac{2e^2}{\eta}.}
\tag{10}
\]
没有连续赢家假设、没有包数或原子数费用。

对固定这份分配，原输出精确分成 `E_dom` 与 `E_coop=E\setminus E_dom`。原 `Γ=Γ_dom+Γ_coop`，后者仍用完整 `μ/W`，含未分配 `μ_0`。一般 polylog 能量目标由 L² 三角缩减为 `Γ_coop` 的预算，不能由存在一份小盒覆盖推定 `E_coop` 为空。对每个协同行只知
\[
\mu_k(Q(x,R(x)))<\eta m(x)\quad\text{for every }k
\tag{11}
\]
（或严格支配规范下的 `≤`）。如果分配仅覆盖部分来源，未覆盖捕获可以占据全部或大部分 `m`；该输出当然可能仍在协同残余。若多个包共同捕获且每包份额都小于 `η`，同样没有被支付。

固定 `η` 给常数能量费，`η^{-1}=polylog(n)` 给 polylog 能量费。`η=n^{-1/2}` 可付实际平方根交通，但仅给平方根能量费，不能自动达到本轮固定 `v=1` 的 polylog 能量目标。

## 4. 已存 packet 搜索的适用字段核验

只读取 `adversarial_shell_pressure_20261007/cell_residual_shell_search.py`、保存 validation、scope receipt 与 root saved review；没有调用生成函数、重建点云或重新运行搜索。

五个 validation 保存行的 family 均为 `multires_packet_network`：candidate `6,20,34,48,55`，维数分别 `8,16,32,64,64`，原子数 `32,32,32,32,128`。脚本 74–85 行给确定的来源标签分组 `group=floor(atom_index/8)`，所以每份 packet 包含八个原子，不按输出改变标签。

其最终来源位置为
\[
z_i=A c_{g(i)}+\delta U_i,\qquad
\delta=\frac{\max(.08,\min(.42,2.5\,\mathrm{jitter}))}{n},
\quad U_i\in[-1,1]^n.
\tag{12}
\]
所有包均集中于以 `A c_g` 为中心、边长 `2δ≤.84/n` 的轴盒。原窗口是 `[1,2]`、`a=1`，所以这些盒满足 `side≤a/n`。这个适用性是由生成公式直接推出，**不依赖** root 记录的实际直径约 `.15–.62` 除以 `n` 的浮点测量。

脚本 112–121 行对所有原子给正权并归一化原总质量 `W=1`。令 `μ_g` 为原八点标签的权重子测度，则 `Σ_g μ_g=μ`，不丢点、不复制、不改变权重；packet 盒重叠也允许。因而这五个输入全部具备本稿 (2)–(3) 的来源分配资格。若未来保存行使用不同 packet 生成或后续修改几何，必须重新核验这些字段，不能只按 family 名套用。

既有搜索只测试 origin0、side `1/n` 的半开网格捕获。一个 packet 横跨网格边界可以使每格份额小，而整个固定 packet 的份额仍大；这个现象不违反搜索原判据，却会被本稿更强的 packet 支配判据删除。因此“固定格协同”不等于“对所有合法来源自适应小盒分配协同”。

**本稿没有检验每条输出的 packet 捕获份额。** 八点是 packet 的来源标签数，不表示每个赢家恰好捕获同一整包，也不表示各赢家都满足 (4)。要登记具体被删交通，需在原完整赢家 cube 上计算 `μ_g(Q)/m`，保留同一 `λ,R,E`，按约定等号规范切分；不能只用 winner 平均 captured atom count，不能把某八点包的全部 `M_g` 当作其当前捕获量。

**状态。** 一般 bounded-diameter、来源自适应、分数分配支配 lemma 及费用 (8)、(10) 通过直接解析独审。保存 packet 输入满足来源合同；具体支配输出比例未在本稿计算。余项须是对选定分配重新定义的完整输入协同交通，不据这项覆盖工具宣称一般赢家都已支付。

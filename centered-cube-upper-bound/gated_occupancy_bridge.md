# 实际门控硬来源的一次径向占用：精确身份与能量范围

2026-10-07。解析审计；只新增本文件，没有启动数值实验或修改主账。

**结论。** 原唯一 LCA、源对分数分配与已经积分完第一跳和完整 continuation 的子核，可以先合并再积分软后验，得到同一硬来源上的 `0≤g≤1`。全部重交通及其深补支因此有不含父层数 `J` 的精确 cube-cone 径向身份。这个身份给全窗 `O(n)W` 的既有粗费用，尚未给平方根费用。径向 `L²` polylog 是一个充分但更强的候选；抽象任意 `R,g` 长壳例只否定过松外类中的能量界，**该例连原硬幅度带也不满足**，不能用来排除保留真实 hard winner、hardband 与 FIRST 的候选。本文另给硬带下单原子、有限原子数与固定窄尺度窗的严格已付子支。

## 1. 读取与冻结范围

已读本目录 `fragmentation_budget_review.md`、`deep_gram_bridge.md`；上一层 `anisotropic_cube_flow/actual_profile_occupation_20261006.md` 与 `heavy_mass_fee_review_20261006.md`；原 `概率接口阶段证明.tex` 的 2696–2753（唯一 LCA）、5522–5582（first-jump 测度、晚时支付）、6095–6200（实际积分子核、固定尺度规范接触）段。也读取 2026-10-06《中心立方体_geom余项最终估计_Pro提示词》，采用其明确给出的完整行归一化与最新重捕获接口。提示词引用而本次未获得的三份更新上游文档未被独立认证。

已应用 G03/E01 的 SKILL、provenance、method、cube-interface 范围守卫：实际来源映射和 occupation 边缘必须给出，普通非负核不能冒称 Gram，源抽样不能替代 Lebesgue 输出。本稿具体身份由原核与 Tonelli 直接证明，不调用其外部定理、未定义 TC-A14 或占用 LP 无间隙结论。

冻结一份完整实际输入 `μ=|f|dx`，质量 `W`，原事件 `E`、软赢家 `L_s(x),σ(x)`、共同阈值 `q`、硬赢家 `R(x)∈[a,b]` 及硬幅度
\[
u(x)=h_{R(x)}*\mu(x),\qquad 2\lambda<u(x)\le4\lambda\quad(x\in E).
\tag{1}
\]
保留原 FIRST、全部未来上限、GOOD、score、owntrace、出生、时钟、严格 CP/GP、抽样删除、森林种子和早期第一跳资格。硬来源固定为
\[
\nu_h=\mu_{\rm hi}|_{B_0},\quad B_0=\{t_C\le b\},\quad
H=\nu_h(\mathbb R^n)\le W.
\tag{2}
\]
以下 `B=n log(b/a)` 只表示窗的对数体积宽度，与出生集合 `B_0` 区分。`H=0` 时全部目标交通为零，后面取 `H>0`。核以全边长定义：`h_R(v)=R^{-n}1_{||v||∞≤R/2}`。

## 2. 先合并所有父组，再积分软端：没有 J

完整原软行的高来源子概率是
\[
d\bar p_x(y)=q^{-1}P_{\sigma(x),L_s(x)}(x-y)d\mu_{\rm hi}(y),
\qquad \bar p_x(\mathbb R^n)\le1.
\tag{3}
\]
固定种子 `θ`，对每个实际来源对 `(y,z)`，同顶格时仅有一个 LCA `P(y,z,θ)`；不同顶格的余项门置零。将该唯一父组的重门、原分数权、所有上述资格以及实际早期子核比值合并为
\[
A_h(x,y,z,\theta)\in[0,1].
\tag{4}
\]
准确说，如果以父标签写系数 `a_P`，则 `A_h=Σ_P a_P`，且对每一固定 `(x,y,z,θ)` 至多一个父标签非零，故 **总系数**仍不超过一。不能只检查每个 `a_P≤1` 然后逐父积分一次软概率。

原第一跳子核是
\[
Q_x^{\rm early}(x,y)=\int_0^{\min(\sigma,T_n)}(1-v)^{n-1}
\sum_i\int_{y+re_i\notin C(y)}G_{1-v,L_s}(r)
[h_{L_s}*K_{\sigma,v,L_s}](x-y-re_i)\,dr\,dv.
\tag{5}
\]
它保留第一跳即退出与**全部后续传播**，并满足 `0≤Q_x^early≤P_{σ,L_s}`。因此 (4) 中可使用这个已经积分完历史的比值；`P=0` 处 `Q=0`，比值取零。若 `C(y)` 有原分数标签，须先按原 `Σ_Cη_C≤1` 的分配合并；不能为每个标签独立再取得一个一。固定 `(v,r)` 的历史密度比不具有 (4) 的概率界。

定义实际门控硬列份额
\[
g_\theta(x,z)=\int A_h(x,y,z,\theta)d\bar p_x(y),
\qquad0\le g_\theta(x,z)\le1.
\tag{6}
\]
把 `g` 在 `E` 外置零，并在 `E` 外显式将 `R(x)` 延拓为 `a`，从而 `β(x)=0`。这使后面全部 cone 点上的 `β` 都有定义，且不改变任何原交通。原硬后验是 `u^{-1}h_R(x-z)dν_h(z)`，所以因子 `u` 准确抵消，正 Tonelli 给
\[
\boxed{R_{\rm heavy}=\mathbb E_\theta\int_E\int
h_{R(x)}(x-z)g_\theta(x,z)d\nu_h(z)dx.}
\tag{7}
\]
如果起点只是原交通的正支配接口，则 (7) 左右改为 `≤`；不能把正扩大所得核称为原实际精确核。提示词明确给出的实际概率交通接口下 (7) 为等式。

深补支只需把 `1_{d_{P(y,z,θ)}>D}` 放进 (4)，得到 `g_θ^{>D}≤g_θ≤1`，并有同一个 (7)。浅深是精确交通分解，未复制硬来源。这里不需要节点门向祖先保持，也不需要先删除所有资格做 telescoping。`J` 只在此前逐父固定质量包络中出现；它不属于 (7) 的来源质量。

## 3. cube-cone 极坐标与真实尺度迟延

令 `ς_n` 为 `S∞={ω:||ω||∞=1}` 上归一化 cone measure：在 `[-1,1]^n` 取均匀点并径向投影至边界即得到它。对任何非负 Borel `F`，严格有
\[
\int F(x)dx=\int d\varsigma_n(\omega)\int_0^\infty
F\bigl(z+\tfrac12v^{1/n}\omega\bigr)dv.
\tag{8}
\]
这是因 `|{x:(2||x-z||∞)^n≤v}|=v`。它使用 cone measure，不是未经归一化的面面积测度。

取来源中心对数体积坐标
\[
t=n\log\frac{2\|x-z\|_\infty}{a},\quad
x_{z,\omega}(t)=z+\frac a2e^{t/n}\omega,\quad
\beta(x)=n\log\frac{R(x)}a\in[0,B].
\tag{9}
\]
于是 `v=a^n e^t`、`dx=a^n e^t dt dς_n`，精确核/Jacobian 乘积是
\[
h_{R(x)}(x-z)dx
=e^{t-\beta(x)}\mathbf1_{\{t\le\beta(x)\}}dt\,d\varsigma_n(\omega).
\tag{10}
\]
`t` 是来源中心径向坐标，不是 first-jump 时间 `v`、停止时刻或质量深度。`β(x)-t` 是真实硬赢家相对该来源距离的迟延；不能把它免费置零。

对每个来源、种子定义
\[
\Phi_{\theta,z}(t)=\int g_\theta(x_{z,\omega}(t),z)
e^{t-\beta(x_{z,\omega}(t))}
\mathbf1_{\{t\le\beta(x_{z,\omega}(t))\}}d\varsigma_n(\omega),
\]
\[
\phi(t)=H^{-1}\mathbb E_\theta\int\Phi_{\theta,z}(t)d\nu_h(z).
\tag{11}
\]
所有实际门通过 (6) 原样保留；软端可能有非紧尾也仍在 `g` 内。由 (7)–(10) 得
\[
\boxed{R_{\rm heavy}=H\int_{-\infty}^B\phi(t)dt.}
\tag{12}
\]
深补支对应 `φ^{>D}`，同样有精确式。原子对角线 `x=z` 对每个固定来源的 `dx` 为零，再由 Tonelli 处理，不删除来源原子；闭硬面规范在 (10) 保留。身份可先对有限空间盒截断再取非负极限，因此不预设目标事件体积有限。

由 `β≥0`、`g≤1`，
\[
0\le\phi(t)\le
\begin{cases}e^t,&t<0,\\1,&0\le t\le B,\\0,&t>B.\end{cases}
\tag{13}
\]
故核心径向交通 `≤H`，且
\[
R_{\rm heavy}\le(1+B)H\le[1+n\log(b/a)]W.
\tag{14}
\]
这是 source-once 的已有 `O(n)W` 粗天花板，不是平方根改进。它优于无必要逐父复制的 `J(1+B)W`，不声称新增已付分支覆盖此前主账。

## 4. L² 能量是充分桥，不能从恒等式免费得到

先注意目标按完整来源 `W` 收费，最自然的归一化是
`Θ(t)=(H/W)φ(t)`，此时 `R_heavy=W∫Θ`。若 `H≪W`，要求下面 (15) 的 `H` 归一化能量界，会额外强于目标所需的
`∫_0^B Θ²≤C log^A(n+2)`。后者已经足以得到 `O(√n log^{A/2}(n+2))W`。本文不假定 `H` 与 `W` 可比；以下先展示较强 (15) 的充分方向，第 8 节回到正确的 `W` 归一化。

若对**所有实际输入**能证明
\[
\int_0^B\phi(t)^2dt\le C\log^A(n+2),
\tag{15}
\]
则 Cauchy–Schwarz 给
\[
R_{\rm heavy}\le H\{1+\sqrt{BC}\log^{A/2}(n+2)\}
=O(\sqrt n\log^{A/2}(n+2))W.
\tag{16}
\]
这里 `B≤nlog2`。它不需要给每个来源、方向或父组单独证明 (15)。源质量/种子平均后的 `φ²` 较弱；Jensen 只能给
\[
\int\phi^2\le H^{-1}\mathbb E_\theta\int d\nu_h(z)
\int\Phi_{\theta,z}^2,
\tag{17}
\]
反方向不成立。要求每个来源的能量都 polylog，可能进一步过强。

由于 `0≤φ≤1`，原目标的 `L¹` 平方根预算只直接蕴含 `∫φ²≤O(√n n^{o(1)})`，并不蕴含 polylog。取 (15) 是额外平方校准假设；不能以“把交通写成一个平方”证明它。`φ²` 的数值非负也不证明原 LCA 源对核是 PSD，`deep_gram_bridge.md` 的实际非零 cross-child 核障碍仍有效。

更弱且正好足够的实际接口是
\[
H\int_0^B\phi^{>D}(t)dt
\le C\sqrt n\,n^{o(1)}W+\varepsilon\lambda|E|,
\quad\varepsilon\le49/8192.
\tag{18}
\]
它允许能量集中、来源间补偿与可吸收行费用，不强求 (15)。

## 5. 抽象 R,g 的长壳守卫：缺少硬带，不是 actual 反例

仅在外类 `R(x)∈[a,b]`、`0≤g≤1`、固定来源中，取 `ν_h=δ_0`，`0<L≤B`，
\[
E=\{0<t(x)<L\},\quad R(x)=2\|x\|_\infty\ (x\in E),\quad g=1_E.
\tag{19}
\]
在 `E` 外将 `R` 延拓为 `a`。于是 `β=t`，精确有
\[
\phi(t)=\mathbf1_{(0,L)}(t),\quad
R/H=L,\quad\int_0^B\phi^2=L.
\tag{20}
\]
取 `L=√n`，一次费用本身正好平方根，能量仍为 `√n`，说明 `L¹` 目标在这个外类并不迫使 polylog `L²`。取 `L=B=Θ(n)` 则说明只凭 `R,g` 范围根本不能证明目标。

**关键限制。** 这里单原子的硬响应是 `u(x)=R(x)^{-n}=a^{-n}e^{-t}`；对同一个 `λ`，`2λ<u≤4λ` 只能支持宽 `log2` 的 `t` 窗。因此 (19) 的长壳没有原 hardband，更没有原实际赢家/FIRST、共同重捕获、出生/CP/GP/早第一跳资格。它只排除“由 (13) 直接推得 (15)”这种外类论证，不排除保留真实 winner 与 hardband 后的 (15)，也不是一般中心立方体弱型反例。

## 6. 保留真实 winner 与同一硬带后的初步界

对一般完整输入 (1)，捕获质量不超过 `W`，故
\[
R(x)^n<W/(2\lambda),\qquad
\beta(x)<\log\frac{W}{2\lambda a^n}.
\tag{21}
\]
若右侧非正，事件为空；否则设
`B_*=min(B,log[W/(2λa^n)])`。真实 profile 支撑于 `t≤B_*`，所以
\[
\int_0^B\phi^2dt\le B_*,\qquad R_{\rm heavy}\le(1+B_*)H.
\tag{22}
\]
这保留了完整 hard winner/hardband，但仍有输入质量深度依赖，不能登记为统一 polylog 能量。当前没有从 hardband 的**上界**部分获得一般 diffuse 来源的维数统一能量改善。由已知 `φ≤1` 所得 (22) 只是明确初步界，不是硬带候选的失败证明。

### 6.1 每源固定窄尺度窗：严格 source-once 已付分支

若某笔子交通对固定 `(θ,z)` 只在
\[
c(\theta,z)\le\beta(x)\le c(\theta,z)+L
\tag{23}
\]
上发生，其中 `c` 在输出选择前固定，`L≥0`，则该来源所有 cone 方向的积分都至多
\[
\int_{-\infty}^{c}e^{t-c}dt+
\int_c^{c+L}1dt=1+L.
\tag{24}
\]
即使 `c` 位于原窗之外，(24) 仍为合法上界。积分来源/种子给子交通 `≤(1+L)H`，没有 `J`。特别地，任何**一个固定**对数体积宽 `L=√n` 的实际硬尺度带可 source-once 支付。覆盖整个窗的多带不能让每带再免费获得一份 `H`。若 `c` 随 `x` 移动，(24) 不再适用。

### 6.2 单原子完整来源：硬带已消除长壳

若完整 `μ=mδ_z`，捕获的所有实际硬行都有 `u=mR^{-n}`。由 (1)，
\[
\log\frac{m}{4\lambda a^n}\le\beta(x)
<\log\frac{m}{2\lambda a^n}.
\tag{25}
\]
这是每源固定宽 `log2` 窗。故任意实际 `g≤1` 的交通 `≤(1+log2)m`；而由 (13)，profile 能量也 `≤1+log2`。这里只在实际接口允许原子时使用；原 `L¹` 输入的单窄峰不自动等于单原子。该检查说明 (19) 不能用于保留硬带的候选。

## 7. 有限硬原子数：捕获份额分层给出明确条件费用

设硬来源 `ν_h=Σ_{j=1}^N m_jδ_{z_j}`，`m_j>0`，完整来源 `μ` 还可有其他质量。取 `0<η≤1`。对捕获的原子定义完整硬行份额
\[
\rho_j(x)=\frac{m_jR(x)^{-n}}{u(x)}
\quad(z_j\in Q(x,R(x))),
\tag{26}
\]
未捕获原子的份额为零。其和不超过一；不把硬高子源重新归一化。

把原交通按硬来源分别切成 `ρ_j≥η` 与 `ρ_j<η`，实际 `g` 和全部门不变。高份额交通满足
\[
2\eta\lambda<m_jR^{-n}\le4\lambda,
\]
所以每个原子都有固定尺度窗
\[
\log\frac{m_j}{4\lambda a^n}\le\beta(x)
<\log\frac{m_j}{2\eta\lambda a^n},
\tag{27}
\]
宽 `log(2/η)`。由 (24) 相加，高份额费 `≤[1+log(2/η)]H`，每个来源只用一次。低份额逐行用 `g_j≤1` 得
\[
\sum_{j:\rho_j<\eta}m_jh_R(x-z_j)g_j(x)
\le N\eta u(x)\le4N\eta\lambda.
\tag{28}
\]
从而严格有
\[
R_{\rm heavy}\le[1+\log(2/\eta)]H+4N\eta\lambda|E|.
\tag{29}
\]
取 `0<ε≤1`、`η=ε/(4N)`，得到
\[
\boxed{R_{\rm heavy}\le[1+\log(8N/\varepsilon)]H
+\varepsilon\lambda|E|.}
\tag{30}
\]
所有实际门可保留，高低两支精确覆盖原交通。给深补支同样切分也成立。若 `ε=49/8192`，可按提示词主账半吸收，费用有显式 `logN`。

这是**有限硬原子数的范围检验与既有窄窗工具推论**，不作为一般路线的下一目标。一般 `L¹` 输入无原子，任意离散逼近的 `N` 可无限增长；不能由 (30) 推出一般 `O(√n n^{o(1)})`。把原子换成固定质量片时，捕获份额只控制片的**当前捕获质量**，不保证片全质量都被捕获，(27) 的固定质量下端因而失效。一般路线必须继续保持 `N` 无关。

## 8. 尚欠的具体几何交叉估计：保留 hardband 与实际最大赢家

此节给候选精确对象，不声称已证明。令 `Θ=(H/W)φ`，并对实际重/深门 `g` 定义
\[
\psi_v(t)=W^{-1}\mathbb E_\theta\int d\nu_h(z)\int d\varsigma_n(\omega)
g_\theta(x_{z,\omega}(t),z)
\mathbf1_{\{0\le\beta(x_{z,\omega}(t))-t\le v\}},\qquad v\ge0.
\tag{31}
\]
`g` 含原 `E`。由 `e^{-d}=∫_d^∞e^{-v}dv`，非负 Tonelli 严格给
\[
\Theta(t)=\int_0^\infty e^{-v}\psi_v(t)dv,
\quad
\int_0^B\Theta(t)^2dt
\le\int_0^\infty e^{-v}\int_0^B\psi_v(t)^2dt\,dv.
\tag{32}
\]
第二式是相对于概率 `e^{-v}dv` 的 Jensen。这把能量预算转成**尺度迟延不超过 `v` 的共同径向壳占用**，保留真实 `R`；不是把 `R` 改成来源到输出的距离。

若先做仅删除实际门的 hard 放松，用完整 `μ` 代替 `ν_h` 且用 `1_E` 代替 `g`，得到
\[
\Gamma_v(t)=W^{-1}\int d\mu(z)\int d\varsigma_n(\omega)
\mathbf1_E(x_{z,\omega}(t))
\mathbf1_{\{t\le\beta(x_{z,\omega}(t))\le t+v\}}.
\tag{33}
\]
逐点 `ψ_v≤Γ_v≤1`。这里 `E` 仍由同一完整输入的 hardband 定义，`R(x)` 仍是原实际尺度集 `\mathcal J⊂[a,b]` 上的最大赢家，不在每个来源或方向上重选。`\mathcal J` 是原允许的物理尺度集，与森林层数 `J` 区分；原接口若只给有限 `\mathcal J`，此处不将其扩大为连续窗口。

这个放松候选要证明的是，对所有完整正输入、同一实际 `R`，有某些维数无关 `C,p,A` 使
\[
\boxed{\int_0^B\Gamma_v(t)^2dt
\le C(1+v)^p\log^A(n+2)\quad(v\ge0).}
\tag{34}
\]
式 (34) 必须区分两个量词版本：**连续窗口试探**使用 `\mathcal J=[a,b]` 的 hardband 与连续最大赢家；**原尺度集版本**使用原任意（包括有限）`\mathcal J` 的 hardband 与该集上的实际赢家，并要求常数对 `\mathcal J` 一致。只有后一版本与上述原接口直接相接。对连续版本的证明或压力检验，不能自动视为有限 `\mathcal J` 的上界：扩大尺度集可能改变赢家，并使原 `E` 上的连续硬最大值超过 `4λ`，破坏同一 hardband。

若原尺度集版本成立，(32) 中 `∫e^{-v}(1+v)^p dv` 为固定常数，就得到目标能量。它是充分的 **N 无关、保留原 hardwinner+hardband 外模型桥**，仍强于仅在 actual geom 上证明 (18)。从连续窗口试探传回原有限尺度集，须另证转移及全部资格费用；或从稠密连续极限重建原 FIRST、赢家、hardband 与门控接口，目前尚未完成。若某个 hard 放松版本 (34) 失败，应把 `g` 的 FIRST/早期退出等门放回 (31) 收紧，而不能据此否定原目标。

准确交叉项不是一个同输出的来源对 Schur 核。记
`χ_v(t,z,ω)=1_E(x_{z,ω}(t))1_{t≤β(x_{z,ω}(t))≤t+v}`，则
\[
\int_0^B\Gamma_v(t)^2dt
=W^{-2}\iint d\mu(z)d\mu(z')\iint d\varsigma_n(\omega)d\varsigma_n(\omega')
\int_0^B\chi_v(t,z,\omega)\chi_v(t,z',\omega')dt.
\tag{35}
\]
两接收点分别是 `x=z+ρ(t)ω` 与 `x'=z'+ρ(t)ω'`，只共享径向 `ρ(t)`，一般 `x≠x'`。现有同输出的近来源对/透镜支付不能直接代付 (35)。也不要求每个 `(z,z')` 的交叉壳列独立 polylog；需要的是同一 `μ⊗μ/W²` 的平均，允许实际质量补偿。

实际最大赢家尚有以下可用但未转换成 (34) 的确定性共同约束：设
`A_s(x)=μ(Q(x,a e^{s/n}))`。对 `x∈E`，第一条上限只在原允许的尺度 `a e^{s/n}∈\mathcal J` 上成立：
\[
A_s(x)\le u(x)a^ne^s,\qquad
A_{\beta(x)}(x)=u(x)a^ne^{\beta(x)},\qquad 2\lambda<u(x)\le4\lambda.
\tag{36}
\]
只有 `R` 在整个连续 `[a,b]` 上最大时，第一条才可对全部 `0≤s≤B` 使用。接触等号与硬带来自实际 `R`，不需要连续最大假设。若原硬半径只是 witness 而非最大赢家，只能保留实际接口确实给出的部分，不能免费使用第一条尺度上限。需完成的几何论证是由对应尺度集上的 (36) 把近接触壳 `χ_v` 的**不同接收点**交叉项，转成输入质量加权的可积 overlap/capture 预算；现有 Jacobian、`g≤1` 和径向正层蛋糕没有给出该转换。

当前没有证明或反驳 (34)。本稿未执行 root 另行组织的 hardband 压力数据，不能以“本次未证”当反证，也不能将只保留 hard 门的数值样本叫作 actualFIRST 交通。

## 9. 已付项、最小剩余桥与回代范围

当前严格已付结果是：核心 `≤H`；整个重/深交通 source-once 粗费 `(1+B)H`；每源输出无关的宽 `L` 尺度窗费 `(1+L)H`；保留 hardband 的单原子检查；有限硬原子数条件费用 (30)。原浅父组费用仍为 `R_{≤D}≤(J+√n)W`，`D=√n/J`。这些是替代或交通子支费用，不能把全窗粗费、浅费与 (30) 当成各自独立已付项相加。

要推进一般实际深补支，需要利用 (6) 中尚未放松的共同 FIRST/未来、真实 winner/hardband、共同重捕获、出生/CP/GP 与早第一跳信息，对其实际 profile 证明 (18)，或证明一份覆盖全部来源的固定窄窗/捕获份额分解且有可计算的补支费用。单凭 `g≤1`、固定出口质量 `≤2W`、后验条件 Bessel 或每层源质量预算都没有完成这一步。

若最终证明深交通 `≤C√n n^{o(1)}W+ελ|E|`，再将原浅支一次回代最新提供的账；`ε≤49/8192` 留有固定半吸收余量。尚须核查全尺度局部化与更新上游接口，方能得到一般中心立方体结论。

**最终状态：** 实际门控 source-once 合并与 cube-cone 身份通过；抽象长壳只否定过松 `R,g` 外类，并明确缺少 hardband；原真实 hardband 下统一 polylog 径向能量及一般深交通平方根预算仍未证。本轮未启动长数值，未修改主账。

# 低诊断 count、平均 future 权与 source-capacity 接口

2026-10-07。本轮结论：平均 future 的精确 Bernstein 标签权 w_k 单调增加；它与已付 tail 后的低 K 方向一致，没有产生矛盾。可将 `w_K≥ε_n` 标签并入既有 S 来源支付，**同价替换**该来源判据，不另加一份 fee。严格补集获得一个独立于 p_C 极小值的诊断 count 上限；source-only `log(1/p_C)` 价格没有一般 polylog 质量和。

## 1. 原交通、诊断标签与查重

完整原 μ、q、实际 σ,Ls,Rh、whole future cap、FIRST、出生、唯一 LCA、strict CP/GP、共同重捕获、far、nonconcentration、短壳、实际 early jump/full continuation、lowS、lowcoin 保持。沿 `actual_future_capacity_bridge_20261007.md`，每条实际历史仅按原核 `P^A/P` 正分配，K=|A|。已付 `K≥κ_C(x)`，κ_C=nσ+√(2nσ log(1/p_C))+2log(1/p_C)/3；剩余诊断标签 `K<κ_C`。K 不是实际跳跃、continuation 活跃数或下界源标签占用。

查重包括原 FIRST-MGF/FACE/MIXED、`future_softness_moment_budget_20261007.md`、`low_future_replacement_budget_20261007.md` 的虚拟 replacement K_*、`joint_future_scale_budget_20261007.md` 的最强 S、`actual_future_capacity_bridge_20261007.md` 的容量 tail。本稿的 w_K 是原 σ mixture 的标签 Radon–Nikodym 权；它不等于 replacement 谱权 ω_k，更不重跑既有 Beta/D/M 包络。

预先固定 α=⌈√n⌉、理论 ε_n=ln^-4(n+2)，定义既有正核

\[
\overline P_{s,L}=\alpha\int_0^1v^{\alpha-1}P_{s+(1-s)v,L}\,dv,
\quad S(u)=\sup_{s\in[0,1],L\in[a,b]}\overline P_{s,L}(u),
\quad\int S\le K_{n,\alpha}.
\tag{1}
\]

其可测性、closed faces/b端点和 `K_{n,α}=O(√nlog(n+2))` 沿已核原工具；本轮不重证该列定理。当前 σ>1/n，分母 Pσ,Ls>0。

## 2. 原 mask 标签的精确权与边界

固定0<σ≤1。对0≤k≤n，令

\[
w_k(\sigma)=\alpha\int_0^1v^{\alpha-1}
\left(1+\frac{1-\sigma}{\sigma}v\right)^k(1-v)^{n-k}\,dv.
\tag{2}
\]

σ<1 时每个原分量准确满足

\[
\overline P_{\sigma,L}^A=P_{\sigma,L}^A w_{|A|}(\sigma),
\qquad \boxed{\overline P_{\sigma,L}/P_{\sigma,L}
=\mathbb E[w_K(\sigma)\mid x,y]}.\tag{3}
\]

条件期望是原核的 mask 接受权；没有归一化原实际 history、没有声称门与自然路径 mask 独立。正装饰可在每条历史上直接执行，原 G 与所有门保留。

设ρ=(1−σ)/σ。α整数时完全正的有限式为

\[
\boxed{w_k=\alpha\sum_{j=0}^k\binom kj\rho^j
B(\alpha+j,n-k+1)},\quad
w_0=\binom{n+\alpha}{\alpha}^{-1},\quad w_n\ge1.
\tag{4}
\]

相邻 integrand 的比率 `(1+ρv)/(1−v)>1` 于0<v<1，所以 w_k 严格随 k增加。固定k>0时它严格随σ增加而减小；k=0恒定。完整计数 prior 权

\[
\eta_k=\binom nk\sigma^k(1-\sigma)^{n-k}w_k
\tag{5}
\]

是先抽 V~Beta(α,1)、再抽 Bin(n,σ+(1−σ)V)的分布，Ση=1，均值为 `n[1−(1−σ)/(α+1)]`。这份 prior 不是 y 后验或实际接受历史的独立 Bernoulli 律。

σ=1时原核只有K=n，w_n=1。形式上的其余值 `w_k(1)=1/binom(n−k+α,α)` 不改变零质量标签。σ↓0时w_0保持上述值，k>0有 `w_k∼σ^-k αB(α+k,n−k+1)`；原σ=0只有K=0，却存在新future soft分量。因此σ=0不能以有限w_0把全部future核写成原核的RN乘子，连硬盒内部也不能忽略新分量。只在σ>0使用(3)；这是原实际余项的合法范围。

root提出、作者解析复核的计算补充：令d_k=w_k−w_{k−1}，则

\[
d_1=\alpha w_0/(n\sigma),\qquad
d_{k+1}=\frac{(\alpha/\sigma)w_k+k\rho d_k}{n-k},\qquad
w_{k+1}=w_k+d_{k+1}\quad(0\le k<n),
\]

k=0时第二项置零。证明是对 `v^α(1+ρv)^k(1−v)^(n−k)` 积分分部，两端因α>0、n−k≥1消失；再用 `d_{k+1}=(α/σ)∫v^α(1+ρv)^k(1−v)^(n−k−1)dv`。此正递推只是(2)的计算身份，不添加空间结论。作者守卫已在收到递推之前終態，使用独立Beta和正项递推；不因新增计算写法重跑冻结数据。

## 3. lowS 与低K没有代数冲突；新的同价支付

最强 lowS 使每个保留来源

\[
S(x-y)/P_{\sigma,L_s}(x-y)<\varepsilon_n,
\quad \mathbb E[w_K\mid x,y]<\varepsilon_n.
\tag{6}
\]

这给 weighted upper tail信息；低K恰使w_K小，不能反向推出(6)失败。例如w_0与p_C无关，且n≥512时w_0<ε_n：α≥3、binom(n+α,α)≥(n+1)(n+2)(n+3)/6>(n+2)^2≥ln^4(n+2)。末个粗界由ln x≤√x可得。这只核函数范围，不认证实际 lowS/CP-GP/history 余项样本。

然而既有 S 来源费用可**同价扩充标签**。令 H={S/P≥ε_n}，这是原已付来源标签。当前增加H补集里的 `{w_K≥ε_n}`，只取尚未付实际 lowcoin/lowK 部分也成立。联合 kernel 边际满足

\[
\begin{aligned}
M_{paid}(x,y)
&\le P_{\sigma,L_s}1_H
+1_{H^c}\sum_{A:w_{|A|}\ge\varepsilon_n}P^A_{\sigma,L_s}\\
&\le P_{\sigma,L_s}1_H+(\overline P_{\sigma,L_s}/\varepsilon_n)1_{H^c}
\le\boxed{S(x-y)/\varepsilon_n}.
\end{aligned}\tag{7}
\]

全部actual门先保持；再正支配、按原完整hard行≤C_hq积分。原 y 边际一次支付

\[
\boxed{R_{H\,\mathrm{or}\,(H^c,\,w_K\ge\varepsilon_n)}
\le C_hK_{n,\alpha}W/\varepsilon_n.}\tag{8}
\]

这**替换**已登记的S来源分支费用，不把两个同阶 fee累加。先前高完整average输出费用、coin子支和容量tail费用各沿原互补分支保留；本稿没有用(8)重领这些历史整账。最新条件账在移项后乘8192/49一次。理论ε_n不改变，因而费用仍O(√nlog^5(n+2))W。

定义唯一整数 cutoff

\[
k_\varepsilon(\sigma)=\min\{0\le k\le n:w_k(\sigma)\ge\varepsilon_n\}.
\tag{9}
\]

由w_0<ε_n、w_n≥1，1≤kε≤n。paid包含等号；新严格标签补集是

\[
\boxed{K<\min\{\kappa_C(x),k_\varepsilon(\sigma(x))\},
\quad S/P<\varepsilon_n,\quad\ell<\delta_n,\quad\text{全部原门}.}\tag{10}
\]

kε不因p_C→0而发散；但它不限制early活跃数，也未控制接收空间或证明(10)空。

另有准确实际来源几何弱式。令 `O_s(x,y)=#{i:|x_i−y_i|>Ls/2}`，closed hard面不计outside。原每个非零分量必须K≥O_s，故(3)与单调性给

\[
\overline P_{\sigma,L_s}/P_{\sigma,L_s}\ge w_{O_s}(\sigma),
\qquad S/P<\varepsilon_n\Longrightarrow O_s<k_\varepsilon(\sigma).
\tag{11}
\]

这是y对Ls的outside数，不是y对Rh、不是y−z的far坐标或cone面，也不是早跳坐标。完整S包含其它物理参数，但(3)的分母只在原σ,Ls匹配，不能把(11)免费推到每个尺度。

## 4. 不能免费把log(1/p_C)价当polylog来源账

原dyadic父P孩子质量m_C、总M=Σm_C，born opposite质量Hopp≤M。对m_C>0，原clipped p_C给

\[
u_C=\log(1/p_C)\le\tfrac12\log n+\log(M/m_C).
\tag{12}
\]

若P至多2^n个孩子，以π_C=m_C/M计，则

\[
\boxed{\sum_{C\subset P}m_Cu_C
\le M[\tfrac12\log n+H(\pi)]
\le M[\tfrac12\log n+n\log2].}\tag{13}
\]

逐原J层不交来源求和给 `≤J[.5ln n+nln2]W_hi`，没有atom count，但有明确维数与深度费用。2^n个等质量、full-born孩子时p_C=1/[√n(2^n−1)]，Σm_Cu_C=M[.5ln n+ln(2^n−1)]≥(n−1)ln2·M，说明不能仅从来源树质量期待polylog。这里是source-only sharp范围，不认证整个actual输入门。

若把价乘回forward容量，Σp_CHopp u_C≤Σm_Cu_C/√n，仍可为Θ(√n)M；再接N_h=Θ(n)hard列会失去目标阶。无限层binary full-born树每层Σm_Cu_C=.5lnn·W，因此没有删去J的无条件全树总和；本项目原J有限可明收，不能在泛称任意树时隐去。新kε截断诊断K不改变这些完整prior质量，未证明实际接受流继承一份熵折扣。

已只读下界构造总表A/B：非均匀固定占用/Gibbs/XOR/隐变量以及增长Cantor允许任意正联合标签权，不能额外限source atom/label数。那里的源标签占用k与本稿Bernstein K不同。固定峰型或Cantor指定最细短窗的上限不能套其整个增长族；本轮未重跑旧下界模型。

## 5. 剩余空间接口

(8)是来源一次、同价增强；(11)是逐真实来源的弱几何约束。(13)说明按u_C收费的基本范围。它们未证明固定hard z对不同x的selected占用小；K=0的诊断分量具有小w_0，且K=0不排除真实路径多次跳跃。原fullfuture、strict CP/GP、birth、far与强S合同仍须在同一输入/真实接受流上取得空间交叉估计。一般√n目标尚未完成，不把有限系数压力当actual反例。

## 6. 预登记与三轮精确守卫收据

新 `low_count_future_capacity_guard_registration_20261007.json` 在执行前固定三个维数、α、rational checker ε和系数范围；理论ε_n=ln^-4仍不改变。只核w权、joint paid边际及source entropy的精确代数；不重跑旧Beta/D/M、原φ或actual样本。

`low_count_future_capacity_exact_guard_20261007.py` 已exit0，`low_count_future_capacity_exact_guard_results_20261007.json` 保存三轮 n=1024、4096、16384，α=32、64、128，共 **37,321 项 exact integer/Fraction 检查 PASS_EXACT**。没有随机种子、没有live handle。注册、脚本与结果之间的SHA收据保存在JSON。

注册的数值阈值仅是 **ε_checker=1/4096**；另用精确有理ε_eq=w_1核包含等号的分配。两者均不是理论ln^-4阈值的数值认证。注册区间0≤k≤4α+64内，精确单调搜索的代表cutoff如下：

| n | α | σ=2/n 的k_checker | σ=1/(4√n) 的k_checker | σ=1/2 的认证范围 |
|---:|---:|---:|---:|:---|
|1024|32|41|63|192<k_checker≤1024|
|4096|64|81|147|320<k_checker≤4096|
|16384|128|160|331|576<k_checker≤16384|

前两列是该注册有理阈值下的精确整数cutoff，分别查w_k≥ε_checker与w_(k−1)<ε_checker；半软度没有超出预定区间搜索，只报严格下界，≤n来自解析w_n≥1。表中增长不拟合渐近阶，也不改变理论或主账参数。

各轮σ取2/n、1/(4√n)、1/2，w样点包含0、1、2、floor(α/2)、α、2α、4α+64，正Beta直接式与项递推交叉核。独立RN身份采用全n维、分别8/16/32个可支持soft坐标的非负系数模型：原mask期望与直接future replacement正积分完全相等。这些模型只核标签接口；hardpersistent零系数并非原φ核或实际输入，不用它们声称完整future/CP-GP残余可发生。

joint marginal守卫分别检查source-high阈值等号、source-low含正high-label质量及label等号分配。联合上包只有一份S/ε，严格补集正分割全部原质量。source entropy用对称压缩的2^n等质量born孩子精确核p、完整质量、(n−1)ln2下界参数和binary深度1/2/8的显式层质量；这些是source-only范围证书，不构成actual压力模型。

作者只读复核root转录的 `阶段研究与数值检验笔记.tex.txt` 末节“同价增强来源包络与低诊断计数补集”及 `current_joint_budget_20261007.md` 末节：权公式、σ0限制、来源费用替换、严格新补集、仅Ls的outside数、熵与J范围均准确。未修改这两个root文件，未重复执行其数据。

本轮新增记账仅同价增强(8)，没有新增总费。一般目标依然需原全部门下的真实selected空间占用估计；(10)尚未支付。

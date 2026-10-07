# 全 future 软度的正矩预算与 Beta 平均空间费

2026-10-07。本稿只处理当前未付交通的互补子集。已付近全软带、逐源退出桥不重加。采用已读取的 [B03](/Users/zhengzhihao/.codex/skills/math-b03-generalized-moments/SKILL.md) 方法中的逐点 majorant／正积分弱对偶；不调用矩锥强对偶、有限原子归约或后验独立性。新结果是一个完整输入的 Beta 平均最大列预算及其明确可付输出支；低平均未来响应的 actual 补支仍开放。root 已全文核 (1)–(16)，tensor 已全文只读独审通过；三轮新精确系数守卫也已落盘，尚不据此认证实际 FIRST 样本。

## 1. 原合同和本轮工具变量

完整正来源 `μ`，`W=μ(ℝ^n)>0`，`n≥512`，`a≤b≤2a`。冻结当前实际 FIRST 的 `σ(x),L_s(x),q` 与硬赢家 `R_h(x)`。原合同为
\[
P_{\sigma,L_s}*\mu(x)=q,\qquad
P_{s,L}*\mu(x)\le q\quad(\sigma\le s\le1,\ L\in[a,b]),
\]
\[
2\lambda<u(x)=h_{R_h(x)}*\mu(x)\le4\lambda,
\quad 3\lambda/4\le q\le\lambda,\quad C_h=16/3.
\tag{1}
\]
完整 GOOD、score、owntrace、出生、森林、唯一 LCA、重捕获、分数、删除门以及前轮取得的 nonconcentration、固定短壳和 actual far 门均保留。这里只向当前 residual 切输出子集，不为它重新选 FIRST。

`P_{s,L}=[(1−s)h+sφ]_L^{⊗n}`。原 tex 730–865 已有 FIRST-MGF、FACE、MIXED；本稿不把这些计数尾界再次登记。新对象是以一个**预先固定** `α=⌈√n⌉` 平均整个未来区间，不按输出优化平均权。

对 `0≤v≤1` 写 `s(v)=σ+(1−σ)v`，并定义
\[
\overline P_{\sigma,L}^{(\alpha)}
=\alpha\int_0^1v^{\alpha-1}P_{s(v),L}\,dv,
\quad
\overline F_\alpha(x)=\overline P_{\sigma(x),L_s(x)}^{(\alpha)}*\mu(x).
\tag{2}
\]
`σ=1` 的端点是 `P_1`，不是硬核。完整 future cap 立即给 `0≤F̄_α≤q`，但不会自动给它相对 `q` 的目标阶可用下界；已有 tiny 下界见 (6a)。

## 2. 完整正 Bernstein 矩锥和弱对偶

固定实际 `x,σ,L`。令 `j_σ=(1−σ)h+σφ`。逐坐标有
\[
\psi_{s(v)}=(1-v)j_\sigma+v\phi.
\]
对 `A⊂[n]` 定义完整输入系数
\[
D_A(x)=\int\prod_{i\in A}\phi_L(x_i-y_i)
\prod_{i\notin A}(j_\sigma)_L(x_i-y_i)d\mu(y),
\quad D_k=\sum_{|A|=k}D_A.
\]
于是精确地
\[
F(v)=P_{s(v),L}*\mu(x)
=\sum_{k=0}^nv^k(1-v)^{n-k}D_k(x),
\quad D_k\ge0,\quad D_0=q,\quad F(v)\le q.
\tag{3}
\]
这里 `k` 是在实际 `j_σ` 基础上的 replacement 标签数，**不是**原 FIRST 后验的 soft 坐标数。全部来源相关性留在 `D_A`，不将其变成 Bernoulli 后验。

任何固定正测试测度 `η(dv)` 都给
\[
\sum_kD_k\int v^k(1-v)^{n-k}d\eta(v)\le q\|\eta\|.
\tag{4}
\]
若目标系数 `a_k≤∫v^k(1−v)^{n−k}dη` 对全部 `k` 成立，则 `Σa_kD_k≤q‖η‖`。这是直接逐点 majorant 后正积分的弱证书；无需强对偶或达到性。特别地，Beta(a,b) 概率权给
\[
\sum_k\frac{(a)_k(b)_{n-k}}{(a+b)_n}D_k\le q
\quad(a,b>0).
\tag{5}
\]
均匀权是 `Σ_k D_k/[(n+1)binom(n,k)]≤q`；本稿权则是
\[
\overline F_\alpha
=\sum_k\alpha B(k+\alpha,n-k+1)D_k\le q.
\tag{6}
\]
更准确地，它有一个无需猜测的 tiny 下界：`D_0=q` 和正性给
\[
\overline F_\alpha/q\ge\alpha B(\alpha,n+1)
=\binom{n+\alpha}{\alpha}^{-1}.
\tag{6a}
\]
这只纠正“完全没有下界”的说法；该倒数远大于允许的 subpower 费用，不足以把全输出免费归入 (15)。
没有在此把 `D_0=q` 再除一次 `q`。真正的空间费如下，不由单行约束直接推出。

## 3. inactive 参数的平均与峰值 majorant

先证明统一 continuation 工具，随后以 `t=0` 直接用于完整输入。固定 `0≤t<T_n<1/16`、`c=1−t`，`b_t=h*G_c`，真实 inactive 参数是
\[
r=(1-\sigma)/(1-t),\quad
S_{\sigma,t,L}=[rh+(1-r)b_t]_L^{\otimes n}.
\]
前轮已解析核验 `I(b_t)≤10`，并保留全部 hard 面得到固定 inactive 数 `m` 的形状列费
\[
\int\sup_{a\le L\le b}M_{A,t,L}
\le1+\log(b/a)[m+\tfrac12\sqrt{10(n-m)}].
\tag{7}
\]
`t=0` 时 `G_1=w`、`b_0=φ`；因此 (7) 同样适用于原完整 `P`，无需 first-exit 参考。

在 future 软度 `s(v)` 上，`r′=r(1−v)`。对 `r>0`，置 `ρ=r′`，平均参数的密度是
\[
f_{\alpha,r}(\rho)=\frac\alpha r(1-\rho/r)^{\alpha-1}
\mathbf1_{0<\rho<r}.
\tag{8}
\]
`r=0` 单独读作 `δ_0(dρ)`，对应全软 `b_t^{⊗n}`；不能误写成 `h^{⊗n}`。

最大化 (8)：对固定 `ρ`，令 `y=ρ/r`，要最大化 `y(1−y)^{α−1}`。临界点 `y=1/α`，因此对 `α≥2`
\[
\sup_{0<r\le1}f_{\alpha,r}(\rho)
\le\begin{cases}
1/\rho,&0<\rho\le1/\alpha,\\
\alpha(1-\rho)^{\alpha-1},&1/\alpha\le\rho<1.
\end{cases}
\tag{9}
\]
第一段精确峰值为 `(1−1/α)^{α−1}/ρ`，此处只作安全放宽；第二段最大值在 `r=1`。

令 `q_m(ρ)=binom(n,m)ρ^m(1−ρ)^{n−m}`，平均系数的归一峰值为
\[
e_m=\binom nm\sup_{0\le r\le1}
\int \rho^m(1-\rho)^{n-m}f_{\alpha,r}(\rho)d\rho.
\]
`m=0` **直接**用概率归一化得 `e_0≤1`，不把它送入发散的 `∫dρ/ρ`。对 `m≥1`，(9) 给
\[
e_m\le\int_0^{1/\alpha}\frac{q_m(\rho)}\rho d\rho
+\int_0^1q_m(\rho)\alpha(1-\rho)^{\alpha-1}d\rho.
\tag{10}
\]
第二段扩大至完整 `[0,1]`，保持正性。定义 `D̄=Σe_m`、`M̄=Σme_m`，交换有限求和及正积分，有
\[
\begin{aligned}
\overline D
&\le2+\int_0^{1/\alpha}\frac{1-(1-\rho)^n}\rho d\rho\\
&\le3+\log(n/\alpha),\qquad (2\le\alpha\le n),\\
\overline M
&\le\int_0^{1/\alpha}n\,d\rho
+\int_0^1n\rho\,\alpha(1-\rho)^{\alpha-1}d\rho\\
&=n/\alpha+n/(\alpha+1)\le2n/\alpha.
\end{aligned}
\tag{11}
\]
`D̄` 最后一行把 `[0,1/α]` 在 `1/n` 分开，分别用 `1−(1−ρ)^n≤nρ` 与 `≤1`。故 `ρ→0` 没有隐藏无限质量；单列 `m=0` 是必须步骤。

## 4. 真正的共同参数最大列和来源一次空间预算

按固定 mask 逐项用 (7)，但只使用**平均之后**的系数峰值 (10)，可得
\[
\int\sup_{\substack{0\le r\le1\\a\le L\le b}}
\alpha\int_0^1v^{\alpha-1}
[(r(1-v))h+(1-r(1-v))b_t]_L^{\otimes n}dv
\le K_{n,\alpha},
\]
\[
\boxed{K_{n,\alpha}
=D_*+\log(b/a)[M_*+\tfrac12\sqrt{10n}D_*],
\quad D_*=3+\log(n/\alpha),\quad M_*=2n/\alpha.}
\tag{12}
\]
没有乘旧 `D_n≈√n`，也没有忽略 `m` 的费用；它由 `M̄` 完整支付。这里只有固定 `t`；不取不同 `t` 核的免费共同 supremum。

对预先固定 `α=⌈√n⌉`，`2≤α≤n`，`M_*≤2√n`、`D_*≤3+(log n)/2`，所以
\[
K_{n,\alpha}\le3+\tfrac12\log n
+\log2\,[2\sqrt n+\tfrac12\sqrt{10n}(3+\tfrac12\log n)]
=O(\sqrt n\log n).
\tag{13}
\]
这对每份输入成立，未使用原子数、source 分区族或 tensor 输入假设。

直接在 `t=0` 取完整 `μ`，Tonelli 给
\[
\boxed{\int_E\overline F_\alpha(x)dx
\le\int\sup_{\sigma,L}\overline P_{\sigma,L}^{(\alpha)}*\mu(x)dx
\le K_{n,\alpha}W.}
\tag{14}
\]
只有一份 `W`；无 first-exit 的额外因子二。可测 supremum 可采用正核的可测规范及 dense 参数、物理 hard 上跳的单侧极限，包含窗口端点；原有限物理集合的版本也被这个连续**上界**覆盖。它没有把有限集合 winner 扩大成连续 winner。

因此固定 `ε>0` 的一般输出支
\[
E_{\rm av}(\varepsilon)
=E_{\rm current}\cap\{\overline F_\alpha\ge\varepsilon q\}
\]
满足
\[
\boxed{\lambda|E_{\rm av}|\le\frac4{3\varepsilon}K_{n,\alpha}W,
\qquad R_{\rm actual}|_{E_{\rm av}}
\le\frac{C_h}{\varepsilon}K_{n,\alpha}W.}
\tag{15}
\]
第二式只用实际交通 `≤∫u≤4λ|E_av|`。例如预先选 `ε=1/16` 给目标阶 `√n polylog(n)` 的来源一次支。若改为 `ε=n^{-o(1)}`，其倒数必须公开付入目标；不能把指数小 `ε` 隐藏在“常数”内。

接主账推荐采用预先指定、与数据无关的
\[
\boxed{\varepsilon_n=\log^{-4}(n+2).}
\tag{15a}
\]
`n≥512` 时 `log(n+2)>2`，故 `ε_n<1/16`。这给 actual 费用 `C_h log^4(n+2)K_{n,α}W=O(√nlog^5(n+2))W`，可放入既有 `√nlog^6` 预算。剩余则取得比固定 `1/16` 更强的完整输入条件 `F̄_α/q<log^{-4}(n+2)`；这次删支没有新增吸收项。

低支仍保留全部原资格，并得到一个新的完整输入筛选条件
\[
\boxed{\overline F_\alpha(x)<\varepsilon q.}
\tag{16}
\]
不是仅筛选实际 gate 的未来平均，也未删去任何原来源。

## 5. actual continuation 门的合法冻结与不同分母

若需要仅付某份门控交通，先在实际参数下，**条件在首跳标签 `ζ=(C,i,t,r,y,w)` 以及最终接收点 `x`**，把原 continuation 细路径门平均成 `a_x(ζ)∈[0,1]`。它是选中 endpoint 子核相对原完整 endpoint kernel 的 Radon–Nikodym 接受密度，保证乘回原完整 `S_{σ,t,L_s}(x−w)` 正好恢复 actual `H_x(0)`。`ζ` 只列首跳历史，不把它冒称全部 terminal 路径标签；hard source 与原 sourcepair 门也在这份原条件平均中积分。只投影原真实门，不假定它与 mask 独立。

随后在同一 receiver 的 future 测试中**冻结**该接受密度、实际 `L_s`、历史截止 `t<min(σ,T_n)`。这是一个 positive test family，不声称更改软度后又是同一 FIRST 或同一条件路径律。因为 `a≤1` 且 `σ′≥σ`，它仍是完整 future 正响应的子核：
\[
H_x(v)\le P_{s(v),L_s}*\mu(x)\le q.
\tag{17}
\]
每个固定 `t` 有 `j_{s(v),t}=(1−v)j_{σ,t}+v b_t`，故 `H_x(v)` 也有正 Bernstein 系数及 (4)–(6)。`H_x(0)=d_x≤q`；这里不可宣称它等于 `q`。实际交通是 `R=∫(u/q)d_x dx≤C_h∫d_x dx`。

平均 `H̄_α=∫αv^{α−1}H_x(v)dv≤F̄_α`。若使用原已核固定首跳参考 `G_{c_t,L_s}≤2G_{c_t,b}`，它的全标签总质量 `≤2W`，逐 `t` 使用 (12) 给
\[
\int\overline H_\alpha dx\le2K_{n,\alpha}W.
\tag{18}
\]
参考放松之后**不**继续冒称该参考的每行仍 `≤q`；(17) 只属于原 actual first-exit 核。

若 `d_x>0`，令 `Π_x` 为该实际门控子核除以 `d_x` 的后验，并把 future 与实际 continuation 的 likelihood ratio 的 Beta 平均记作 `Z_α(ζ)`，则
\[
\overline H_\alpha/d_x=\mathbb E_{\Pi_x}Z_\alpha.
\tag{19}
\]
于是 `EΠZ_α≥ε` 的交通支可付 `2C_hK_{n,α}W/ε`。这是相对 **d_x** 的门控支，和完整输出判据 `F̄_α/q≥ε` 不同。仅 `H̄_α≥εq` 才被 (15) 直接覆盖；低质量 `d_x` 下不能把 `H̄_α/d_x≥ε` 换成 `H̄_α/q≥ε`。优先支付 (15) 的完整输出，再对其补集使用 (19) 时才不重复。

## 6. hard-persistent 范围守卫和最小缺口

所有 future cap 是上界，所以不能仅由正矩锥推出 (16) 为空。一个真实**核组件**守卫说明 likelihood 下界可能很小：取 `σ=1/2`，完整 continuation 的端点偏移每坐标为零。由 resolvent `b_t≤φ/c`、`φ(0)<3/4`、`c>15/16`，有 `b_t(0)<4/5`；同时 `r=(1/2)/c≥1/2`。因此
\[
b_t(0)/j_{σ,t}(0)\le8/9,
\quad S_{s(v),t,L}(0)/S_{σ,t,L}(0)
\le(1-v/9)^n.
\tag{20}
\]
这份 tensor 是核组件，不是一般 FIRST 合法输入，未认证 actual hardband、nonconcentration、far 与完整 future cap 的共同实现。不得作为实际 geom 反例。它只禁止在 (19) 中免费宣称一个全域常数下界。

对任意 `δ∈(0,1)`，其平均 likelihood 有解析上界
\[
\int_0^1\alpha v^{\alpha-1}(1-v/9)^ndv
\le\delta^\alpha+(1-\delta/9)^n.
\tag{21}
\]
取 `δ=1/α`、`α=⌈√n⌉`，右侧趋于零。原已知 hard-persistent 张量所警告的尺度选择损失并未因“平均 future”被消除。

原 tex 730–865 是 FIRST-MGF／FACE／MIXED；3225–3253 已付 near-full；5522–5582 是未平均的 `N_hD_n`；前轮 `far_large_softness_dilation` 的共同 `r` FTC 对 near-full 参数带收费。本轮核验这些位置未发现 (8)–(12) 的 Beta 平均系数峰值或 (15) 的完整平均未来响应删支。原已付 allsoft 输出若与 `E_av` 重叠，只取当前 residual 的交集，不能再领一份完整来源。最新上游三个尺度原件仍缺，是否已涵盖具体该支须另核，不能据本稿声称全新历史整账。

当前可以登记的是 (15) 的一般已付条件支与 (16) 的低平均未来完整输入筛选；余下需利用 hardband、nonconcentration、固定短壳、actual FIRST／future cap，把 low-average 支付或吸收。完整 Bernstein 系数上界没有自行完成这个几何桥，也未证明当前 residual 的低支实际非空。

## 7. 三轮精确守卫收据

新独占脚本 `future_softness_moment_exact_guard_20261007.py` 和同名 `_results.json` 已实际运行，`n=512,1024,4096`、`α=⌈√n⌉` 三轮全部 PASS_EXACT。无随机种子、无 Monte Carlo、无 `G_c` 求积或阶数拟合。

实现用两个独立正表示核对平均系数：一是 (10) 的 incomplete-Beta 峰值 majorant；另一是 `W0~Beta(1,α)`、`J|W0~Bin(n,W0)`、`K|J~Bin(J,r)` 的 thinning 表示。Beta-binomial 的 common-denominator 整数权为 `binom(n−m+α−1,α−1)/binom(n+α,α)`，归一化与均值 `n/(α+1)` 精确核验。每轮对 `m=0,1,α,2α,n`、`r=0,1/α,1/2,1` 做 20 个正系数检查，共 60 个；特别 `r=0` 留全软 `m=0` 点，不送入 `1/ρ`。

首段系数采用精确恒等式
\[
\int_0^{1/\alpha}\frac{q_m(\rho)}\rho d\rho
=\frac1m\Pr\{\operatorname{Bin}(n,1/\alpha)\ge m\}\quad(m\ge1),
\]
因此 `Σm×首段=n/α` 无浮点误差。所有 `m` 都参与 D̄/M̄ 总预算；log 上界由 dyadic 缩放后的正级数及显式有理余项包络认证。

hard-persistent 的 `(1−v/9)^n` 平均也用 incomplete-Beta／Binomial tail 算出精确 Fraction，并与 (21) 的 `δ=1/α` 上界比较：三轮都严格 `<1/8`，数值 log10 约见 JSON。这只是 §6 的真实核组件守卫，不是满足 FIRST、hardband、共同 future cap、nonconcentration 和 far 的输入反例，更不能由有限模型推断低支全局覆盖。

大整数仅保存 binary SHA256、bit length 和诊断 log10，PASS 判定均来自精确整数／Fraction 比较；log10 不参与证书判断。

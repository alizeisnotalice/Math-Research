# 低平均 future 来源的 replacement 谱与高对比补支

2026-10-07。本轮完整读取 `checkpoint_20261004_general_flux_continuation.md` 最后“Beta平均future”节。对象是 `R_far,lowav`；不重新证明 Beta 高平均输出费、旧 FIRST-MGF 或逐源退出桥。沿用 B03 已核的正系数／正积分弱证书，不调用强对偶。只编辑本稿及同前缀守卫。

## 1. 当前完整合同与 root 新来源筛选

完整输入 `μ,W`、`n≥512`、`a≤b≤2a`、`q∈[3λ/4,λ]` 和所有 actual FIRST/future cap、GOOD/score/owntrace、出生、森林、唯一 LCA、分数、共同重捕获、删除门保持。当前硬赢家及 `2λ<u≤4λ`、nonconcentration 微盒合同、固定短壳、早期 firstjump exit／大 mark、选定 vertices 的 far 门均保持，不把条件化停止来源当新的 `μ`。

预先固定 `α=⌈√n⌉`、`ε_n=log^{-4}(n+2)`。每个当前 receiver 还满足
\[
\overline F_\alpha/q<\varepsilon_n,
\quad \overline P_{\sigma,L}=\alpha\int_0^1v^{\alpha-1}P_{\sigma+(1-\sigma)v,L}dv.
\tag{1}
\]
原完整 soft 来源概率为 `ρ_x(dy)=P_{σ,L_s}(x−y)μ(dy)/q`。定义
\[
Z_x(y)=\overline P_{\sigma,L_s}(x-y)/P_{\sigma,L_s}(x-y)>0.
\tag{2}
\]
当前 `σ>1/n`，故分母全域正；`σ=1` 时 `Z=1`，不能留在下面低 Z 分支。

root 已提出并负责来源标签直接收费：`Z≥τ_n=ε_n` 的 actual 交通用 `P_σ≤P̄/τ_n` 来源一次支付 `C_hK_{n,α}W/τ_n=O(√nlog^5n)W`，无新增吸收。以下只处理它的互补来源
\[
\boxed{Z_x(y)<\varepsilon_n.}
\tag{3}
\]
不重复收费 (1) 的高平均输出支或 root 的高 Z 来源支。`ρ_x` 未在 (3) 上重新归一；实际 gate 始终是原子概率下的 `[0,1]` 接受权。

## 2. 正 replacement 谱的逐来源精确身份

固定原 `x,y,σ,L_s`，令 `u_i=(x_i−y_i)/L_s`，
\[
a_i=\phi(u_i)/\psi_\sigma(u_i)>0,\qquad
B_k(a)=\binom nk^{-1}\sum_{|A|=k}\prod_{i\in A}a_i.
\]
`B_0=1`，`B_n=P_{1,L_s}(x−y)/P_{σ,L_s}(x−y)`。逐坐标代数给
\[
\frac{P_{\sigma+(1-\sigma)v,L_s}(x-y)}{P_{\sigma,L_s}(x-y)}
=\prod_i(1-v+va_i)
=\sum_{k=0}^n\binom nk v^k(1-v)^{n-k}B_k(a).
\tag{4}
\]
由正 Beta 积分，
\[
\boxed{Z_x(y)=\sum_{k=0}^n\omega_{n,\alpha,k}B_k(a),
\quad \omega_{n,\alpha,k}
=\frac{\binom{k+\alpha-1}{\alpha-1}}{\binom{n+\alpha}{\alpha}},
\quad\sum_k\omega_k=1.}
\tag{5}
\]
这是一份**预先固定的谱权**，不是实际来源 posterior 的新概率先验，也不是原 FIRST 的 soft 坐标数。`k` 代表虚拟 future 测试中把原 `ψ_σ` 换成 `φ` 的 replacement 数。

对遗漏数 `m=n−k`，谱权为 Beta-binomial：先 `V~Beta(α,1)`，再 `K_*|V~Bin(n,V)`。因此
\[
\mathbb E(n-K_*)=n/(\alpha+1),
\]
\[
\boxed{\Pr_\omega\{m>M\}
=\frac{\binom{n-M+\alpha-1}{\alpha}}{\binom{n+\alpha}{\alpha}}
\le\left(1-\frac{M+1}{n+\alpha}\right)^\alpha,\quad0\le M<n.}
\tag{6}
\]
等式是 hockey-stick；不等式由 `α` 个分式的直接乘积。取 `M≈[(n+α)/α]log(1/δ)` 即使谱权的尾 `≤δ`；其尺度为 `√nlog(1/δ)`。

在剩余逐来源 (3)，每个谱层和每个非负预先固定序列 `c_k≤ω_k` 都严格有
\[
\omega_kB_k<\varepsilon_n,
\quad\sum_kc_kB_k<\varepsilon_n,
\quad\Pr_\omega\{B_{K_*}\ge t\}<\varepsilon_n/t\quad(t>0).
\tag{7}
\]
故近全 replacement 谱中至少 `1−δ−ε_n/t` 的预先权重同时满足 `m≤M` 和 `B_k<t`。这是逐 source 的正谱预算；没有推断任意具体 coordinate mask 的概率或每坐标 `1/n`。

## 3. 与原 FIRST soft 坐标数的准确区别

在固定来源 `y` 后，**原完整核、尚未被 actual 门进一步条件化**的 latent soft 概率是 `p_i=σa_i`；外部 hard 因子为零的坐标 `p_i=1`。原 conditional mask count `K` 的 generating function是 `∏(1−p_i+p_i z)`；完整来源 posterior 混合后通常相关。

对 `v<1`，逐来源有
\[
\prod_i(1-v+va_i)
=(1-v)^n\mathbb E\left[\left(1+\frac{v}{\sigma(1-v)}\right)^K\mid Y=y\right].
\tag{8}
\]
它只是同一 kernel 的代数关系，不重新证明旧全行 FIRST-MGF。`K_*` 的均值 `nα/(α+1)≈n−√n` 来自 (5) 的虚拟 Beta 权；原 `K|y` 均值 `Σσa_i` 可以完全不同。不可把 replacement 的近全谱当成实际 continuation 的高 active mask，更不可据此删掉实际 first-exit continuation。

## 4. 小 Z 推高全软对比：AMGM 与 Jensen

定义逐来源全软对比
\[
A_x(y)=\log\frac{P_{\sigma,L_s}(x-y)}{P_{1,L_s}(x-y)}
=-\sum_i\log a_i.
\tag{9}
\]
对全部正 `a_i` 和 `v∈[0,1]`，weighted AMGM 给 `1−v+va_i≥a_i^v`。再对 `V~Beta(α,1)` 用 Jensen，得到
\[
\boxed{Z_x(y)\ge\alpha\int_0^1v^{\alpha-1}e^{-vA}dv
\ge\exp\left[-\frac\alpha{\alpha+1}A_x(y)\right].}
\tag{10}
\]
`A` 可为负；此时下界超过一，故 (3) 自动排除它。剩余每个来源满足
\[
\boxed{A_x(y)>\frac{\alpha+1}\alpha\log(1/\varepsilon_n)
=4\frac{\alpha+1}\alpha\log\log(n+2).}
\tag{11}
\]
也可直接对 `B_k` 用 AMGM 得 `B_k≥exp(−kA/n)`；它与 (5) 再 Jensen 恢复同一下界，没有假定坐标后验独立。

反方向的 useful 身份保留遗漏修正：令 `b_i=1/a_i=ψ_σ(u_i)/φ(u_i)`，`C_m(b)=e_m(b)/binom(n,m)`，则
\[
B_{n-m}(a)=e^{-A}C_m(b),\qquad
Z=e^{-A}\sum_m\omega_{n-m}C_m(b).
\tag{12}
\]
因此高全软对比不是唯一统计量；被遗漏的少数坐标仍有正乘积修正。不能把 `Z` 直接换成 `e^{-A}`。

## 5. 原全软列工具的有偿来源分层

全软 fixed-shape column 已由前轮 Fisher 核验：
\[
N_{\rm soft}=\int\sup_{a\le L\le b}P_{1,L}(x)dx
\le1+\tfrac12\sqrt{10n}\log(b/a)=O(\sqrt n).
\tag{13}
\]
这项 kernel 工具已经存在；这里仅给它在当前 actual 来源标签上的直接应用，不宣称新全软最大定理。

对任意预先确定 `B_n≥0`，只在当前未付交通取 `A_x(y)≤B_n` 的来源交集。原 actual 子核质量 `≤ρ_x(dy)`，故
\[
\begin{aligned}
R_{A\le B_n}
&\le C_h\int\!\int P_{\sigma(x),L_s(x)}(x-y)
\mathbf1_{A_x(y)\le B_n}\,dx\,d\mu(y)\\
&\le\boxed{C_h e^{B_n}N_{\rm soft}W.}
\end{aligned}
\tag{14}
\]
只用原完整 `μ` 一次，所有门和 firstjump 标签在起点保持；之后正支配上界。没有让 source 分配依输出重领质量，没有在参考 first-exit 放大后继续使用 row cap。

优先选择 `B_n=8loglog(n+2)`，则费 `O(√nlog^8(n+2))W`，补集获得
\[
\boxed{Z_x(y)<\varepsilon_n,\quad A_x(y)>8\log\log(n+2).}
\tag{15}
\]
它在 (11) 的低界之上（`α≥23`），两个阈值间有正宽度；不据此认证实际分层中有非零交通。费用仍是目标允许的 `√n n^{o(1)}`；若主账坚持已有 `log^6` 固定幂，则不得隐去这次 `log^8` 扩展。root 本轮主账尚未采用这份可选全软分层，主 residual 先保留 (3)、(11)；只有选择支付 (14) 后才能登记 (15)。

可选 subpower 版本取 `B_n=√log(n+2)`，费 `O(√n exp√log(n+2))W=√n n^{o(1)}W`。它只有在比当前 `8loglog` 更大时才加强来源筛选；不能宣称有限小维数已经加强。可统一取两者最大，但倒数／指数费用须明记。这里不按数据选择阈值。

## 6. 对比的真实坐标几何与仍欠的空间量

记 `I={i:|u_i|≤1/2}`、`O=n−|I|`。闭 hard 面属于 `I`。逐来源真实恒等式是
\[
A=O\log\sigma+\sum_{i\in I}\log\left[\sigma+\frac{1-\sigma}{\phi(u_i)}\right].
\tag{16}
\]
原 `φ(1/2)>3/8`、`φ(0)<3/4` 给
\[
\sigma\le b_i\le B_\sigma:=\sigma+\tfrac83(1-\sigma),
\quad A\le n\log B_\sigma-O\log(B_\sigma/\sigma).
\]
于是 `A>B_n` 强制
\[
\boxed{O<\frac{n\log B_\sigma-B_n}{\log(B_\sigma/\sigma)}}
\quad(0<\sigma<1).
\tag{17}
\]
分子非正即该 source 支为空；`σ=1` 已被 (3) 排除。它是来源对 `L_s` 的坐标数限制，不是来源相对硬赢家 `R_h` 的同一计数，更不是 hard cone 面标签的限制。

遗漏修正还有完全正的粗界 `C_m(b)≤(8/3)^m`。但 (6) 有效遗漏宽度为 `M≈√nlog(1/δ)`，把这项粗界直接代回会付 `exp[O(√nlog(1/δ))]`，超出目标。既不能删除遗漏修正，也不能把虚拟谱权 `ω` 当 actual mask 概率来躲开它。

真正未付量是 (15) 补集上的
\[
\int u(x)\,\rho_x(dy)\,r_x(dz)\,
A_{\rm actual}(x,y,z,\theta)\,
\mathbf1_{\rm original\ residual}\,
\mathbf1_{\{Z<\varepsilon_n,\ A>B_n\}}dx.
\tag{18}
\]
原 nonconcentration 给每个 `h=2a/n` 盒的 hard 捕获质量 `≤η_n m(x)`；(17) 给 `y` 相对另一 cube `L_s` 的内部坐标多。二者之间尚无已证 source-once 重数、角度或管道预算。固定短壳控制 `z` 对 `R_h` 的 radial depth，也未把这些 `L_s` 内坐标自动变成可付的 `y−z` 几何。不能将受限 `y` 或 replacement endpoint 作为新 `μ` 测试 complete future cap。

当前一般弱进展是 (5)–(7) 的逐来源 positive spectrum、(10)–(12) 的精确高对比桥，以及 (14) 的公开费用来源分层。没有证明 (18) 空、没有证明 actual counterexample；主目标仍开放。

## 7. 新三轮精确代数守卫和范围

`low_future_replacement_exact_guard_20261007.py` 首次执行的 session `34701` 已权威 exit0，同名 `_results.json` 保存三轮 `n=512,1024,4096`、`α=23,32,64`，全数 PASS_EXACT。无随机种子、无模拟或原核求积。未重跑旧 Beta 高支数据。

每轮直接用 integer common denominator 重构全部 replacement 权，核归一化、遗漏均值 `n/(α+1)` 及三个遗漏尾的 hockey-stick／乘积上界。共九个 uniform 正 ratio 模型 `a_i=3/4,1,2`：以 (5) 正和算 `Z`，在 `3/4` 模型另以 incomplete-Beta／Binomial tail 独立核同一值；weighted AMGM 的有理 `v` 检查以及 Jensen 的无根形式 `Z^{α+1}≥B_n^α` 均精确通过。

三轮另用 `n/8` 个 `a_i=2`、其余 `a_i=3/4` 的 mixed 正模型，在 `k=0,1,n−α,n−1,n` 共十五个谱系数上核 `B_k^n≥B_n^k` 与 (12) 的遗漏 reciprocal 身份。原完整核的 conditional count 关系 (8) 以三个有理 `v` 独立生成函数算出，与虚拟 replacement 的不同均值并列保存。

这些是正 ratio／计数代数有限模型。`a_i=1` 的 neutral 和 `a_i=2` 的 uniform 模型也只用于身份范围检查；不声称它们都符合原当前 cap。mixed 模型初始 likelihood 导数为负，加 log-concavity 可核一个 scalar future 单行下降条件；它仍不是共同输入、真实物理赢家、hardband、nonconcentration 与完整 actual FIRST 的联合实例。数值不提供原 residual 非空反例，不证明几何计数预算。

大数只记 binary SHA256、bit length 和诊断 log10；PASS 来自 Fraction／整数比较。root 新 `low_future_source_payment_20261007.md` §1–2 已由本 agent 只读独审：原 `y` 边际 domination 与高 Z 直接空间费通过，无需重选 FIRST、无额外 `q` 或 first-exit `2W`，无新增吸收。此独审没有改动 root 文件。

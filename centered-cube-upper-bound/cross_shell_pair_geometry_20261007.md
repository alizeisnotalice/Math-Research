# 固定单位短壳的双来源几何：位移身份、互捕核与未付非互捕项

2026-10-07。只新增本稿；没有启动数值，没有修改其它稿件。对象是完整输入下的 hard 放松 `Γ_1`，不是 actual FIRST 交通。一般 `∫Γ_1²≤polylog(n)` 本轮仍未证明或否定。以下给出共同径向参数的精确双源换元、真实 winner/hardband 的质量约束以及一个可直接核验的互捕几何核；同时给出解析障碍，说明不能仅凭近来源透镜、同面薄片或核的 Lebesgue 体积完成来源平均预算。

## 1. 对象和技能范围

令 `μ=f dx`，`f≥0`、`W=∫f∈(0,∞)`。允许尺度集为原 `𝒥⊂[a,b]`；连续窗口版本另取 `𝒥=[a,b]`。选定同一输入的可测实际最大赢家 `R(x)∈𝒥`，
\[
u(x)=R(x)^{-n}\mu(Q(x,R(x))),\quad
E=\{2\lambda<u(x)\le4\lambda\},\quad
\beta(x)=n\log(R(x)/a),\quad B=n\log(b/a).
\tag{1}
\]
其中 `Q(x,r)=x+[-r/2,r/2]^n`，`B≤nlog2`。沿用已审稿的 cube-cone 概率测度 `ς_n`：先均匀选面坐标 `j∈{1,…,n}` 与符号 `ε∈{−1,1}`，令 `ω_j=ε`，其余坐标独立均匀于 `[-1,1]`。这份独立性是硬方体 cone 测度的定义；不是 actual firstjump 后验的独立性。

记
\[
r_t=ae^{t/n},\quad \rho_t=r_t/2,\quad
x=z+\rho_t\omega,\quad
\chi(t,z,\omega)=\mathbf1_E(x)
\mathbf1_{\{t\le\beta(x)\le t+1\}},
\]
\[
\Gamma(t)=W^{-1}\int f(z)dz\int\chi(t,z,\omega)d\varsigma_n(\omega),
\quad \mathcal E=\int_0^B\Gamma(t)^2dt.
\tag{2}
\]
所有 `E,R,β` 在虚拟 `t` 改变时保持原定义，没有为每条射线重新选赢家。`0≤Γ≤1`、`𝓔≤B` 是现有粗界。

已读 G03/F03 的 `SKILL.md`、provenance、method 与 cube-interface；使用其“完整展开交叉项、不得用对角或 HS 术语自动支付交叉”的范围守卫。本稿不调用未定义方向捕获树、公开树嵌入、径向 Schrödinger 或维数统一 Schatten 结论；下面均是直接 Tonelli、换元和区间概率计算。

## 2. 双来源位移与两 cone 方向：精确身份

令 `Δ=z'−z`。原四重积分换元的 Jacobian 为一，严格得到
\[
\boxed{
\mathcal E=W^{-2}\int_{\mathbb R^n}d\Delta\int_{\mathbb R^n}
f(z)f(z+\Delta)dz\iint d\varsigma_n(\omega)d\varsigma_n(\omega')
\int_0^B\chi(t,z,\omega)\chi(t,z+\Delta,\omega')dt.}
\tag{3}
\]
两个接收点满足
\[
x'=x+\Delta+\rho_t(\omega'-\omega).
\tag{4}
\]
因此共同径向 `t` 并不等于同一接收点，也不等于两个来源与同一个 cube 的透镜。

为明确“源对位移分布”及不能删除的来源位置，置
\[
P(\Delta)=\int f(z)f(z+\Delta)dz,\quad
d\tau_\Delta(z)=P(\Delta)^{-1}f(z)f(z+\Delta)dz
\tag{5}
\]
于 `0<P(Δ)<∞`。Tonelli 给 `∫P(Δ)dΔ=W²`，所以 `P<∞` 对 `dΔ` 几乎处处，例外和 `P=0` 不影响积分。(3) 等价于先按概率 `P(Δ)dΔ/W²` 选位移，再按 `τ_Δ` 选共同来源位置，再选两个 cone 方向和共同 `t`。

**不能仅保留 `P(Δ)`。** `χ` 还依赖 `z` 经同一 `E,R` 的位置。若删除 `τ_Δ` 下的原位置条件，就已经扩大了模型。对于有限 Borel 来源，也可用源对映射 `(z,z')↦(z,z'−z)` 的推前及测度分解写同一身份；本稿的原 `L¹` 版本没有离散原子数费用。

另一个等价接收点表达是
\[
\begin{aligned}
\mathcal E=W^{-2}\int_0^Bdt\iint d\varsigma_n(\omega)d\varsigma_n(\omega')
\iint dx\,dh\;&f(x-\rho_t\omega)
f(x+h-\rho_t\omega')\\
&\mathbf1_{E_t}(x)\mathbf1_{E_t}(x+h),
\end{aligned}
\tag{6}
\]
其中 `E_t={x∈E:t≤β(x)≤t+1}`，`h=x'−x`。这里 `Δ=h−ρ_t(ω'−ω)`。固定方向的两次平移均保留 Lebesgue 测度；**没有**从 (3) 到 (6) 额外插入 `r_t^n` 或 `e^t`。径向 cone Jacobian 已在定义 (2) 的 surface 概率里处理。

## 3. 保留 hardband 与赢家时的实质质量约束

在 `χ=1` 的行，令 `q=R(x)/r_t=e^{(β-t)/n}`，则
\[
1\le q\le Q:=e^{1/n},\qquad
2\lambda r_t^n q^n<\mu(Q(x,r_tq))
\le4\lambda r_t^nq^n.
\tag{7}
\]
这里的原完整来源 `μ` 正是 (3) 的 `f`；不能另加独立财富源。

若 `r_t∈𝒥`，最大条件进一步给
\[
\mu(Q(x,r_t))\le q^{-n}\mu(Q(x,r_tq)),
\]
\[
\mu(Q(x,r_tq)\setminus Q(x,r_t))
\ge(1-q^{-n})\mu(Q(x,r_tq)).
\tag{8}
\]
后一条是外壳质量的**下界**，不是可免费删除的薄壳上界。即使 `q−1=O(1/n)`，原输入也可把非忽略质量放进这份壳。一般有限 `𝒥` 不包含 `r_t` 时，(8) 不可使用。

连续赢家在任意允许的外扩半径 `R e^{h/n}≤b` 上还有
\[
\mu(Q(x,R e^{h/n})\setminus Q(x,R))
\le(e^h-1)\mu(Q(x,R)).
\tag{9}
\]
这是真实增长限制，但只在允许窗口中成立，而且仍是每个接收点的来源质量上界。它没有直接约束两接收点的共同 `t` 交叉。

可以精确将每条占用射线挂回自己的完整捕获来源：令
`M(x)=μ(Q(x,R(x)))`。则
\[
\chi(t,z,\omega)
=\chi(t,z,\omega)\int
\frac{\mathbf1_{Q(x,R(x))}(y)}{M(x)}d\mu(y),
\tag{10}
\]
且 `M(x)>2λr_t^n`。这构造了合法的共同输入 capture coupling；它不会让 (3) 的两个原来源自动互捕。若用 (10) 替换 (3) 两份 `χ`，得到四来源积分；中间费用 `1/M(x)` 依赖原接收点，不可无偿改为固定源对核。

## 4. 精确互捕／非互捕分解

对 (3) 中每个原双射线，定义
\[
\mathsf M=\mathbf1_{\{z'\in Q(x,R(x)),\ z\in Q(x',R(x'))\}}.
\tag{11}
\]
按 `χχ'𝕄` 与 `χχ'(1−𝕄)` 精确分成
`𝓔=𝓔_mut+𝓔_non`；无需新门或 J。令 `δ=Δ/ρ_t`、`q'=R(x')/r_t`。互捕等价于
\[
\boxed{\|\delta-\omega\|_\infty\le q,
\qquad\|\delta+\omega'\|_\infty\le q'.}
\tag{12}
\]
尤其 `||Δ||∞≤ρ_t(1+Q)`；远于这个阈值的源对一定属于非互捕，**但并不会退出 (3)**。

若两 cone 面坐标相同 `j` 且符号相同 `ε`，(12) 在第 `j` 坐标同时要求
\[
|\delta_j-\varepsilon|\le q,\quad
|\delta_j+\varepsilon|\le q'.
\]
从而
\[
|\Delta_j|\le\rho_t(Q-1)
=\tfrac12r_t(e^{1/n}-1).
\tag{13}
\]
这是实际单位短壳下的真实 `O(r_t/n)` 薄片。异面或同面异号不具有这条双边薄片。又因为 `f` 可以集中在该薄片，薄片的 Lebesgue 宽度本身不等于 `μ⊗μ` 的 `1/n` 概率。

## 5. 两 cube-cone 的显式互捕核

仅为上界，在 (12) 将实际 `q,q'` 同时放宽为 `Q=e^{1/n}`；此处已删除 winner/hardband，不能把以下自由几何核当作原资格的等价表达。定义
\[
p_Q(s)=\Pr_{U\sim\mathrm{Unif}[-1,1]}(|s-U|\le Q)
=\begin{cases}
1,&|s|\le Q-1,\\
(Q+1-|s|)/2,&Q-1<|s|\le Q+1,\\
0,&|s|>Q+1,
\end{cases}
\]
\[
b_Q(s)=\tfrac12\left(\mathbf1_{\{|s-1|\le Q\}}
+\mathbf1_{\{|s+1|\le Q\}}\right),
\]
\[
L_Q(\delta)=\frac1n\sum_{j=1}^n
b_Q(\delta_j)\prod_{k\ne j}p_Q(\delta_k).
\tag{14}
\]
`L_Q(δ)` 正是 `Pr_ω(||δ−ω||∞≤Q)`。cone 测度对称、两方向独立，所以互捕角概率严格为
\[
\boxed{K_Q(\delta)=L_Q(\delta)^2.}
\tag{15}
\]
因此一份已明确扩大范围的几何上界是
\[
\mathcal E_{\rm mut}\le W^{-2}\iint d\mu(z)d\mu(z')
\int_0^B K_Q\left(\frac{2(z'-z)}{r_t}\right)dt.
\tag{16}
\]
角平均中同面、异面交叉已经全部包含；不能只留对角 `j=j'`。

这份核还具有可核验的 Lebesgue 体积。记
\[
A_Q=\int p_Q(s)^2ds=2Q-\tfrac23,\quad
C_Q=\int b_Q(s)^2ds=\int b_Q(s)p_Q(s)ds=2Q-1.
\]
展开 (14) 的平方，得到
\[
\int_{\mathbb R^n}K_Q(\delta)d\delta
=\frac{C_QA_Q^{n-1}}n
+\left(1-\frac1n\right)C_Q^2A_Q^{n-2}.
\tag{17}
\]
`n=1` 时第二项为零，第一项为 `C_Q`。`Q=1` 时化简为
\[
\int K_1=(4/3)^{n-2}(1+1/(3n)),
\quad
\int K_1(2\Delta/r)d\Delta
=\frac9{16}(1+1/(3n))(2/3)^n r^n.
\tag{18}
\]
`Q=e^{1/n}` 相对 (18) 只有维数无关常数损失：
`A_Q/(4/3)=1+(3/2)(e^{1/n}−1)`，其 `n` 次幂有界，而 `C_Q≤2e−1`。故单位短壳互捕的自由角核确有指数小的 Lebesgue 体积。

**为何 (17) 尚未付来源费用。** (16) 以 `μ⊗μ` 积分，而非 Lebesgue 位移积分。原 `L¹` 来源没有统一 `L∞` 密度，不能乘一个免费密度常数。`K_Q(0)=1`，且来源可集中在很小方体中；对这样的来源，(16) 删除资格后的右侧可接近 `B`。这不是 hardband 反例：该右侧已经删去原资格，在单原子来源检查中真实 hardband 只允许宽 `log2` 的赢家窗。它显示单独几何 overlap 上界丢掉了 (7) 中必须恢复的质量窗补偿；不以未经验证的窄峰转移断言全部原资格。

## 6. 解析阻碍：非互捕项不能用近源透镜代换

下面只用于反驳“完整平方项都可归入近源互捕”的推导，不作为研究特殊输入族的终点，也不反驳 `𝓔≤polylog(n)`。

取已经核验的唯一赢家守卫输入
\[
f(y)=c\left(1+\frac{\|y\|_2^2}{nM^2}\right)
\mathbf1_{[-M,M]^n}(y),\quad M=n^2b,\quad \lambda=c/2,
\tag{19}
\]
并要求允许尺度含 `b`。质量 `W=(4c/3)(2M)^n`。若接收点位于 `[-M+b/2,M-b/2]^n`，全部 `R≤b` 的 cube 都在支撑内，
\[
h_R*\mu(x)=c\left[1+\frac{\|x\|_2^2}{nM^2}
+\frac{R^2}{12M^2}\right].
\tag{20}
\]
它随 `R` 严格增加，故 `b` 是唯一赢家；响应在 `(c,2c)`，属于同一 hardband。在来源内盒 `Z_0=[-M+b,M-b]^n`，所有 cone 输出 `z+ρ_tω` 均满足上述条件。记
\[
q_0=\mu(Z_0)/W
=(1-b/M)^n\frac{3+(1-b/M)^2}{4}
\ge\tfrac34(1-1/n).
\tag{21}
\]
因此在整个 `t∈[B−1,B]∩[0,B]` 上，这些来源及全部方向都被占用。

但只要 `||z'−z||∞>b`，接收 cube `Q(z+ρ_tω,b)` 就不可能捕获 `z'`：三角不等式给其任意被捕获来源与 `z` 的距离至多 `b`。故这批源对全部属于非互捕。又 `f≤2c`，所以
\[
\frac1{W^2}\iint\mathbf1_{\{\|z'-z\|_\infty\le b\}}
d\mu(z)d\mu(z')
\le\frac{2c(2b)^n}{W}
=\tfrac32(b/M)^n.
\tag{22}
\]
对 `B≥1`，得到严格的下界
\[
\boxed{\mathcal E_{\rm non}
\ge q_0^2-\tfrac32(b/M)^n.}
\tag{23}
\]
它是常数级、可由远源对贡献的合法 hardband 与唯一赢家能量。它说明即使互捕核有 (18) 的小体积，完整平方也不能被它单独正支配。一般路线必须允许远来源的同尺度相干贡献，并另外支付它；这个例子没有超 polylog 增长。

## 7. 固定微格支配的输出：一般已付能量子支

root 在本轮提出固定格覆盖分支，以下直接核验并接到本稿能量；它适用于任意有限正 Borel 来源，没有原子数条件。固定半开方格分割，边长 `a/n`，每格中心 `c_C`、质量 `M_C=μ(C)`。给定 `η∈(0,1]`，记
\[
E_{\rm dom}=\{x\in E:\exists C,
\mu(C\cap Q(x,R(x)))\ge\eta\mu(Q(x,R(x)))\}.
\tag{24}
\]
若该格支配行，则 `M_C≥ηu(x)R(x)^n>2ηλR(x)^n`；捕获交集非空给
`||x−c_C||∞≤R(x)/2+a/(2n)≤(1+1/n)R(x)/2`。于是
\[
E_{\rm dom}\subseteq\bigcup_{C:M_C>0}
Q\left(c_C,(1+1/n)(M_C/(2\eta\lambda))^{1/n}\right).
\]
正质量格至多可数；直接体积求和和 `ΣM_C=W` 得
\[
\lambda|E_{\rm dom}|\le\frac{(1+1/n)^n}{2\eta}W
\le\frac e{2\eta}W,\qquad
\int_{E_{\rm dom}}u(x)dx\le\frac{2e}{\eta}W.
\tag{25}
\]
这只用原 hardband、`R≥a` 与固定微格几何，不要求 winner 在连续或有限尺度集中最大。

令 `Γ_dom` 为 (2) 将 `E` 换成 `E_dom` 的来源全质量 profile。cone Jacobian 给一个无需连续赢家的精确源／输出换元：若
`s(x,z)=nlog(2||x−z||∞/a)`，则
\[
\int_0^B\Gamma_{\rm dom}(t)dt
=W^{-1}\iint\mathbf1_{E_{\rm dom}}(x)
\mathbf1_{\{0\le s(x,z)\le B,\ 0\le\beta(x)-s(x,z)\le1\}}
e^{\beta(x)-s(x,z)}h_{R(x)}(x-z)dx\,d\mu(z).
\tag{26}
\]
原子对角 `x=z` 对 `dx` 为零；有限 Borel 来源可先逐来源再 Tonelli，原 `L¹` 情形无额外问题。由 `e^(β−s)≤e`、`Γ_dom≤1`，
\[
\boxed{\int_0^B\Gamma_{\rm dom}(t)^2dt
\le\int_0^B\Gamma_{\rm dom}(t)dt
\le\frac{2e^2}{\eta}.}
\tag{27}
\]
实际残余交通在此输出交集的费用也至多 (25) 的 `2eW/η`；只在当前残余中切分，不能与已删交通重复收费。

令 `E_coop=E\setminus E_dom`，仍用**完整原来源 `μ/W`**定义 `Γ_coop`。有精确 `Γ=Γ_dom+Γ_coop`，因此
\[
\sqrt{\mathcal E}\le\sqrt{2e^2/\eta}
+\left(\int_0^B\Gamma_{\rm coop}^2dt\right)^{1/2}.
\tag{28}
\]
所以固定 `η`，或 `η^{-1}=polylog(n)`，已把目标缩减为每个 winner cube 内所有固定格捕获份额 `<η` 的协同残余。若为实际平方根交通选择 `η=n^{-1/2}`，(25) 支付平方根，但 (27) 只支付平方根能量，不能据此称已付 polylog 能量。

协同约束还准确作用在 (10) 的 companion source：固定原被捕获来源 `z∈C` 后，按 `μ|Q(x,R)/M(x)` 选 `y`，则
\[
\Pr(y\in C\mid x,z)<\eta.
\tag{29}
\]
这可保证伴随来源通常不在同一个固定格，却不能直接给 `||y−z||∞≥a/n`：不同半开格可以相邻、点可以任意接近。也不能把 (29) 变成两条不同接收点原来源 `z,z'` 的独立 small-ball 界。下一几何论证必须使用这个明确 capture conditional，不能免费换成无条件 `μ⊗μ`。

## 8. 当前真正剩余的平均几何估计

令 `K_μ(z,z')=∫_0^B∫∫χχ' dς_n dς_n' dt`，则 (3) 的目标是
\[
W^{-2}\iint K_\mu(z,z')d\mu(z)d\mu(z')\le C\log^A(n+2).
\tag{30}
\]
这份核依赖同一输入 `μ` 经 `E,R`，不是可先独立选择再做普通 Schur 的固定核。第 7 节已在原输出上合法删除固定微格支配行；剩余 (30) 将 `χ` 换成 `χ_coop`，保留完整 `μ⊗μ/W²`。已获得的具体结构是：

- 每条占用射线在厚度一的 log-volume 壳中，接收质量满足 (7)，且通过 (10) 可挂回原来源。
- 互捕子核具有精确公式 (14)–(17)；同面同号产生 (13) 的薄片。要将其付成来源平均界，仍需保留 (7) 而控制浓集源对，不能借位移 Lebesgue 体积代替。
- 非互捕子核不能删除或归入同输出透镜，(23) 已给严格解析障碍。需要对不同接收点的原 capture coupling 建立全质量预算，或以支付补支证明其不导致长径向相干。

将 (10) 直接插入 (30) 后，需要控制的明确量是
\[
\begin{aligned}
W^{-2}\int_0^Bdt\iint d\mu(z)d\mu(z')\iint d\varsigma_n(\omega)d\varsigma_n(\omega')
\chi\chi'\iint
\frac{\mathbf1_{Q(x,R)}(y)\mathbf1_{Q(x',R')}(y')}
{M(x)M(x')}d\mu(y)d\mu(y').
\end{aligned}
\tag{31}
\]
它与原能量精确相等，而硬带允许 `M(x)M(x')>4λ²r_t^{2n}`。但是这条下界的费用含 `λ`，并且四来源 capture incidence 还没有预算；不把它写成已证 Carleson 或 HS 条件。

向压力 agent 已发送可执行机制建议：在现有自适应候选中记录互捕／非互捕比例及同面同号来源位移薄片占比，检查高能量是否来自不能归入透镜的非互捕，或来自集中于 `r/n` 薄片的源对。没有要求重跑旧数据，也没有构造尚未核验的反例。

**状态。** 本轮建立精确源对位移身份及可解析 cube 互捕核，并将固定微格支配输出以 (27) 支付；剩余为同一输入、hardband 和赢家下的微格协同双源交叉。排除了“共享径向等于共享接收点”“薄片宽度自动给源概率”“小核体积支付全部平方”三种不合法推进。固定 `v=1` 的一般 hard 放松 polylog 能量仍开放；所有失败机制的认证范围止于上面明确的推导，不否定 actual geom 目标。

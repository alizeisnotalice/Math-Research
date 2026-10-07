# 原捕获后验的径向 age：完整参数、熵方向与临界矩边界

2026-10-07。父任务指定 6.1-sol high。解析独审，不运行数值、不改他人文件或主账。

**结论：连续全窗口真赢家给 age 后验指数尾，方向正确，且可保留全部 floor/端点/原子。** 它产生真实 weak joint 的深部可吸收子支；但在包络熵身份里只改善 E_P A，方向不支付未知弱比。浅部的共同正包络仍恰为整个 k_*，而临界矩 E_πe^A 在完整 L¹、唯一连续赢家的正体积区域仍可为 Θ(n)。本文另证明任意实际接受门的后验 age 约束，保留接受率的真实费用，不假角独立或 free posterior。一般 source-once/geom 主核未付。

使用已读 [D04 SKILL](</Users/zhengzhihao/.codex/skills/math-d04-entropy-chain-rule/SKILL.md>)、method/cube-interface/provenance，固定联合律、KL方向与条件权。下文直接证明 Radon–Nikodym/层饼身份，不调用未读一般熵定理。查重并全文读 [winner radial transport §1](winner_radial_transport_20261007.md)、[posterior depth transfer](posterior_depth_source_once_transfer_20261007.md)、[winner envelope entropy](winner_envelope_entropy_bridge_20261007.md)、[cross-shell](cross_shell_budget_20261007.md)：中间尺度帽、当前 root 的 age/deep 子支已有；本文是独审及熵/临界矩连接，不重新领这些已知费用。

## 1. 原完整后验：准确尾身份与全参数

完整有限正 Borel μ，0<W<∞，全边长闭 cube Q(x,r)，全部连续允许尺度 [a,b]，0<a≤b<∞。令

\[
 U_r(x)=r^{-n}\mu(Q(x,r)),\qquad
 M(x)=\sup_{a\le r\le b}U_r(x).
\]

固定某个 x 及候选 L∈[a,b]，m_L=μ(Q(x,L))>0；M有限且M≥U_L>0。定义近赢家 gap（不与旧接受 g 混用）

\[
 \rho=\log(M/U_L)\ge0,\quad
 \pi_L(dy)=\mathbf1_{Q(x,L)}(y)\mu(dy)/m_L,
 \quad B_L=n\log(L/a).
\]

对真实 π_L 下的 Y，D=2||x−Y||∞，

\[
 A=n\log\frac L{\max(a,D)},\qquad0\le A\le B_L.
\]

0≤t≤B_L 时 r=L e^{-t/n}∈[a,L] 为合法原尺度，含边界原子准确有

\[
 \boxed{\pi_L(A\ge t)=\frac{\mu(Q(x,r))}{m_L}
 =e^{\rho-t-\rho_r}\le\min(1,e^{\rho-t}),}
 \quad\rho_r=\log(M/U_r). \tag{1}
\]

若 U_r=0，置ρ_r=+∞，中式按0解释。t=0为1，t=B_L为原 Q_a 捕获比例，均未删除端点。t>B_L 左端0；t<0左端1，不使用 e^{ρ−t} 的有意义尾段。L=a 或 b=a 时 B_L=A=0，保留整个 floor 原子。L=b 的端点赢家同样合法，不要求尺度导数为0。D=0由 max(a,D)处理，未使用 log0。

精确赢家 L=R 有ρ=0，但允许原 ties和端点；本论证只用完整所有中间 U_r≤M，不求 winner 导数。有限目录赢家通常不能在未列入的 r 上使用 (1)：须改为上邻尺度 r⁺ 的帽 e^ρ(r⁺/L)^n，或已核其最大值等于全连续最大。本文不替原 FIRST 的 hard witness 免费补全连续赢家资格。

式 (1) 的等号含全部原尺度 gap ρ_r：不是任选满足指数帽的后验。它比单条平均 age 更明确保留完整输入/各尺度一致性，却尚无空间来源预算。

## 2. 锐 floor 截断矩与 σ=1

由层饼，令 B=B_L，

\[
 \mathbb E_{\pi_L} A\le
 \begin{cases}
 B,&\rho\ge B,\\
 \rho+1-e^{\rho-B},&0\le\rho<B,
 \end{cases}
 \le\min(B,\rho+1). \tag{2}
\]

对 **0≤σ<1**，

\[
 \mathbb E e^{\sigma A}\le
 \begin{cases}
 e^{\sigma B},&\rho\ge B,\\
 \displaystyle\frac{e^{\sigma\rho}
       -\sigma e^{\rho-(1-\sigma)B}}{1-\sigma},
       &0\le\rho<B,
 \end{cases}
 \le\frac{e^{\sigma\rho}}{1-\sigma}. \tag{3}
\]

证明使用 E e^{σA}=1+σ∫₀ᴮe^{σt}P(A>t)dt 并在ρ处分两段。闭端点的 atom不影响积分；σ=0为准确1。σ<0时最后一条上界一般不成立，例如 L=a、ρ=0、A=0，1>1/(1−σ)。σ≥1不能沿用其带正分母的形式。

σ=1 的合法临界上界单独是

\[
 \mathbb E e^A\le
 \begin{cases}
 e^B,&\rho\ge B,\\
 e^\rho(1+B-\rho),&0\le\rho<B.
 \end{cases} \tag{4}
\]

对真赢家ρ=0，简化为

\[
 \mathbb EA\le1-e^{-B},\quad
 \mathbb Ee^{\sigma A}\le
 \frac{1-\sigma e^{-(1-\sigma)B}}{1-\sigma},
 \quad\mathbb Ee^A\le1+B. \tag{5}
\]

B=0时临界式也为1；0<B<∞ 时σ→1− 的锐中式趋于1+B，而粗式1/(1−σ)发散。均值至多1并不使临界矩有维数无关上界。

## 3. 真正联合约束及接受门的明确代价

对 E={M>τ}，原真赢家 R(x)，I=τ|E|，取原完整 weak joint

\[
 dJ=\tau\mathbf1_E(x)dx\,\pi_{R,x}(dy),
 \quad P=J/I,
 \quad P_X=dx|_E/|E|,
 \quad P_Y=\frac{S(y)}{I/W}\frac{\mu(dy)}W.
\]

source column S=τ∫E1_Q/m dx 是实际列，不把 P_Y 换成 μ/W。由 (1) 对每个非负可测 receiver 测试 ψ（可积或非负扩展值），

\[
 \mathbb E_P[\psi(X)e^{\sigma A}]
 \le\frac1{1-\sigma}\mathbb E_{P_X}\psi(X),
 \quad0\le\sigma<1,
\]

以及 ∫1_{A>u}dJ≤e^{-u}I，u≥0。这些是已核 root 深部子支，不新增 W 费用；没有给 fixed source 的多 receiver 次数上界。

现在取任何同一 joint 上的实际可测门 q(x,y)∈[0,1]，p(x)=∫q(x,y)π_R,x(dy)。p>0 的 row 下接受后验为 qπ_R/p；不假 q独立于 source/角/历史。它准确满足

\[
 P(A\ge t\mid x,\mathrm{accepted})
 \le\min\{1,e^{-t}/p(x)\},\qquad0\le t\le B_R(x).
\]

故 (2)–(3) 以有效 gap log(1/p) 合法应用，特别有

\[
 \boxed{\int q e^{\sigma A}d\pi_R
  \le\frac{p(x)^{1-\sigma}}{1-\sigma},\qquad
  \mathbb E[A\mid x,\mathrm{accepted}]
     \le1+\log(1/p(x)).} \tag{6}
\]

p=0 row 的未归一化式按0，不能定义接受后验。近赢家ρ时，第一式右边乘 e^{σρ}，第二式加ρ；全部仍由原π_L与原门产生。

更强的全 joint 合同：若 r=∫q dP∈(0,1]，P_acc=qP/r，则对真赢家

\[
 \boxed{\mathbb E_{P_{\rm acc}}e^{\sigma A}
       \le\frac{r^{-\sigma}}{1-\sigma},\qquad
       \mathbb E_{P_{\rm acc}}A\le1+\log(1/r).} \tag{7}
\]

第一式由 (6) 及 p^{1−σ} 凹性，第二式也可直接从 P_acc(A>t)≤min(1,e^{-t}/r)层饼得到。这适用于任意实际 source subset B 取 q=1_B，以及任意已积分进[0,1]的 actual gates；这里 r 是 **真实 Palm 接受质量** ∫B S dμ/I，不能换成 μ(B)/W。接受稀少会真实消耗 log(1/r)，不能免费把 accepted posterior 保持为完整π_R。

这些约束比无条件 EA≤1保留更多原 joint 信息，但仍没有固定来源列 S 的截断尾。完整后验径向帽不能经过 conditioning on Y后变成同一帽；(7) 恰记录了这种条件化的代价。

## 4. 代入精确熵身份：方向与临界 normalizer

沿旧 envelope稿令 m=I/W、α(x)=τ/M(x)在E上，k_*(z)=max(a,2||z||∞)^−n1_{2||z||∞≤b}、Z=1+nlog(b/a)，Q=(μ/W)(dy)k_*(x−y)dx/Z。原 d=A，准确有

\[
 \boxed{D(P\|Q)+\mathbb E_P[A+\log(M/\tau)]
                      =\log(Z/m).} \tag{8}
\]

因此 EA≤1只给

\[
 D(P\|Q)\ge\log(Z/m)-\mathbb E_P\log(M/\tau)-1.
\]

这不是 m 的改进上界：未知 joint KL 是非负而没有独立的小上界/几何下界；age 为正的项若有下界才会压低 m。即使把 EA精确算成常数，也没有一个原角 score的P均值支配 m。Q 的 source law仍为 arbitrary μ，Q 的宽径向混合不因 P 的 age帽变成独立角/坐标后验。

可进一步在每个 receiver x 清楚识别参考改变。令 K_*(x)=k_* *μ(x)，原 Q 的 receiver 条件后验是
q_x(dy)=k_*(x−y)μ(dy)/K_*(x)，不等于π_R。π_R≪q_x，在原查询上有

\[
 \frac{d\pi_R}{dq_x}=\frac{K_*(x)}{M(x)}e^{-A},
 \quad
 D(\pi_R\|q_x)=\log[K_*(x)/M(x)]-\mathbb E_{\pi_R}A.
 \tag{9}
\]

并且恰有

\[
 \frac1{M(x)}\int_{Q(x,R)}k_*(x-y)\mu(dy)
             =\mathbb E_{\pi_R}e^A. \tag{10}
\]

所以把原π的σ<1矩用于这个 change-of-reference 时，所需的是 **σ=1 临界矩**，其合法上界是1+B_R，仍可O(n)。查询外的 envelope来源又是额外正项；当 R=b 时它为0。用σ<1矩免费替代 (10) 会漏掉正质量 normalizer。

## 5. 完整 L¹、唯一连续赢家的临界矩 Θ(n) 实例

### 5.1 正体积 endpoint winner：准确临界积分

固定 n≥1、b>a>0，T=b(n²+1)、κ=1/4、τ=3/4。完整输入为

\[
 f(y)=\left(1+\kappa\frac{\|y\|_2^2}{nT^2}\right)
                       \mathbf1_{[-T,T]^n}(y),
 \quad W=(2T)^n(1+\kappa/3).
\]

正体积 receiver盒 E₀=[−T+b/2,T−b/2]^n上，所有尺度 query都在来源盒内。准确

\[
 U_L(x)=C_x+\frac{\kappa L^2}{12T^2},
 \quad C_x=1+\kappa\frac{\|x\|_2^2}{nT^2}\ge1.
\]

随 L 严格增大，故原全部连续窗口的唯一 winner是 R=b，M>1>τ。又 f≤1+κ=5/4<2τ，故 E₀ 属于同一真实幅度带 τ<M≤2τ；没有自由构造 posterior、极值/nearmax 认证或补加无限背景。W 是完整来源质量。E₀ 体积占来源盒比例 (1−b/(2T))^n≥1−n/[2(n²+1)]≥1/2，保持实际正体积范围。

令 B=nlog(b/a)、δ=κb²/(12T²)。对 x∈E₀，真实径向CDF为

\[
 \pi_b(D\le r)=\left(\frac rb\right)^n
                  \frac{C_x+\delta(r/b)^2}{C_x+\delta},
 \quad a\le r\le b.
\]

于是准确临界矩

\[
 \boxed{\mathbb E_{\pi_b}e^A
  =1+\frac{C_xB+(n\delta/2)(1-(a/b)^2)}{C_x+\delta}.} \tag{11}
\]

证明为1+∫₀ᴮeᵗP(A>t)dt，代入CDF；μ≪dy没有端点原子问题。由于 C_x≥1、δ≤1/48，

\[
 1+\frac{48}{49}B\le\mathbb E e^A\le1+B,
 \qquad\mathbb EA\le1-e^{-B}\le1. \tag{12}
\]

因此固定 b/a>1 下，σ=1临界 normalizer为Θ(n)，同时根提出的全部 age均值/σ<1矩帽成立。R=b使 (10) 恰为 K_*/M；(9) 可有约 log B 的条件KL而age均值是常数。这严格排除“用σ<1矩或EA≤1免费将 critical normalizer 改成常数”的步骤。

此输入密度有上界，可属于已付 bounded-density L² 支；它不反驳一般 weak平方根界或真实geom。它只验证当前熵改变参考的必要维数成本，不能被重新包装为 nearflat/高源极值压力样本。

### 5.2 严格内部 winner：相同临界障碍不依赖尺度端点

可进一步取 a=1、b=2、ε₀=1/100，完整 L¹ 输入

\[
 f(y)=\left[1+\frac{\epsilon_0}{n}
            \sum_{i=1}^n(y_i^2-y_i^4)\right]
                  \mathbf1_{[-2,2]^n}(y).
\]

在 [-2,2] 上 t²−t⁴∈[−12,1/4]，所以全来源密度范围为[22/25,401/400]，完整质量

\[
 W=4^n[1+\epsilon_0(4/3-16/5)]
       =\frac{368}{375}4^n.
\]

对 E₀=[−1/10,1/10]^n，所有尺度查询内接完整来源盒。令 v_x=n⁻¹Σx_i²，

\[
 U_s(x)=C_x+\alpha_x s^2-\beta s^4,
 \quad C_x=1+\frac{\epsilon_0}{n}\sum_i(x_i^2-x_i^4),
 \quad\alpha_x=\epsilon_0(1/12-v_x/2),
 \quad\beta=\epsilon_0/80.
\]

这是逐坐标准确积分，quartic 平均中的6x_i²E u_i²给 x_i²s²/2，E u_i⁴=s⁴/80，未漏交叉系数。U作为 s² 的函数严格凹，唯一最大值

\[
 R(x)^2=\frac{\alpha_x}{2\beta}
             =10/3-20v_x\in[47/15,10/3]\subset(1,4).
\]

故整个正 Lebesgue E₀ 有原唯一严格内部连续 winner，gap为0。τ=1/2 时 M≥22/25>τ，E₀⊂E；没有 nearmax 认证，也不声称 E₀ 是完整E。posterior 密度相对 Uniform(Q_R) 下界是
κ₀=(22/25)/(401/400)=352/401。因此从真实CDF层饼直接有

\[
 1+\frac{352}{401}n\log R(x)
       \le\mathbb E_{\pi_R}e^A\le1+n\log R(x),
 \qquad\mathbb E_{\pi_R}A\le1. \tag{12a}
\]

R²≥47/15>1，临界矩仍为Ω(n)；没有 endpoint导数或自由后验的弱点。此输入同样 bounded-density，不能将该临界障碍称作原weak型或geom反例。

为独审 root 的新注册组件算式，固定 x、q=R²、u_R=U_R，全部 moment来自相同实际CDF。准确

\[
 \mathbb EA=\frac n{R^n u_R}
 \left[\frac{C_x(R^n-1)}n
 +\frac{\alpha_x(R^{n+2}-1)}{n+2}
 -\frac{\beta(R^{n+4}-1)}{n+4}\right],
\]
\[
 \mathbb Ee^{A/2}=1+\frac{n/2}{R^{n/2}u_R}
 \left[\frac{C_x(R^{n/2}-1)}{n/2}
 +\frac{\alpha_x(R^{n/2+2}-1)}{n/2+2}
 -\frac{\beta(R^{n/2+4}-1)}{n/2+4}\right],
\]
\[
 \mathbb Ee^A=1+\frac n{u_R}
 \left[C_x\log R+\frac{\alpha_x(q-1)}2
                         -\frac{\beta(q^2-1)}4\right]. \tag{12b}
\]

证明分别是 n/(Rⁿu_R)∫₁ᴿrⁿ⁻¹U_rdr、1+(n/2)/(Rⁿ⸍²u_R)∫₁ᴿrⁿ⸍²⁻¹U_rdr、1+(n/u_R)∫₁ᴿU_rdr/r。n∈4ℕ时前两式仅需 q 的整数幂，临界式仅 logR=logq/2 需要有理区间。root负责新三轮登记与执行；本稿只独审公式，不重复运行。

**root 新证书的只读独审。** 已读 [guard](interior_winner_critical_age_guard_20261007.py)、[registration](interior_winner_critical_age_guard_20261007_registration.json)、[results](interior_winner_critical_age_guard_20261007_results.json)、[receipt](interior_winner_critical_age_guard_20261007_receipt.json)，独核三文件 SHA256 全匹配收据、72 records、每轮24、792 saved checks 全为真。n为4/8/12、16/64/256、32/128/512；ε₀=1/100与1/200、四种有理 receiver patterns。脚本的 integral_power 与 (12b) 完全一致，n∈4ℕ保证所有非 log 项为 Fraction；logq 的120项正 atanh展开及几何尾给上下区间，critical 使用两端计算，浮点仅作显示。来源密度与质量公式、内部 q范围、mean≤1、half MGF≤2和临界线性下界均有明确 saved flags。

执行者是 root，终态0；本代理没有重跑 guard、MC或空间积分。critical/n 显示范围约[0.572861,0.851897]是这些已保存记录的描述，不拟合阶数，也不声称有限receiver点检查证明整个 E₀。连续正体积结论由本节解析负责；没有 nearmax、original FIRST/history 或 weak常数反例资格。

## 6. 浅 age 共同包络恰为全 k_*，floor 不能删除

对任意 u≥0，定义保持原 floor 的浅核共同包络

\[
 k_{\rm shallow,u}(z)=
 \sup_{a\le L\le b}h_L(z)
 \mathbf1_{\{n\log[L/\max(a,2\|z\|_\infty)]\le u\}},
\]

age只在 h_L非零的 incidence上定义。令 D=2||z||∞。D≤a时选择 L=a，age=0且h_L=a⁻ⁿ=k_*；a<D≤b时选择L=D，age同为0且h_L=D⁻ⁿ=k_*；D>b时全部为0。因此全参数准确

\[
 \boxed{k_{\rm shallow,u}=k_*,\qquad
                       \int k_{\rm shallow,u}=Z.} \tag{13}
\]

closed-cube face选择在这个逐z supremum中合法；不能用 fixedL面零测删除 adaptive L=D 的图像。若改用 raw age nlog(L/D)，core会出现1−e⁻ᵘ之类因子，但 raw变量不受合法尺度 a 的后验帽同样控制，不能混用。

与 (13) 相容，root 的∫S_deep,u dμ≤e⁻ᵘI是实际 receiver后验 traffic可吸收支；它与 (13) 衡量不同对象。cross-shell已另有不要求 truewinner的长迟延 source包络支付，不能将两个重叠支相加重领。剩余浅部即使只宽u/n仍能覆盖整个source-fixed外径向包络，需真实 selected occupation 而非再次取sup。

## 7. 对 source-once 与 actual geom 的精确距离

本稿通过：完整原后验的所有中间尺度 gap身份、锐端点矩、实际门的接受率倾斜联合约束、熵身份方向，以及完整L¹ critical-moment实例。它们不假source/角坐标独立，没有 receiver 后验抽样代替真实 Lebesgue E。

当前一般可用子支是原 continuous truewinner weak joint的深部吸收。若 actual交通在同一 joint有已核密度0≤q≤C，则其深部≤Ce⁻ᵘI；若只有选定 hard幅度 U_R≤Cτ和 nearwinner gap≤G，则对应未归一化 hard捕获深部≤Ceᴳe⁻ᵘτ|E|。这些是明确条件，不是已核当前 R_angle：原软来源可能在硬cube外，原firstjump/fullfuture/LCA门可能偏置后验，original hard witness未必为连续max；缺mapping或q的上界时不得登记paid。

还需一个真正的跨receiver几何连接：让 actual shallow/fullfuture 流量在同一source上的重复占用可付，或证明一个保留原posterior及选择依赖的几何score，其P均值强迫未知m而参考MGF可核。仅有age常数、critical normalizer 线性基线、角MGF与sourcePalm身份不提供它。一般 sqrt(n)n^{o(1)} 与原geom主预算未解决。本轮没有数值执行或新阶数声称。

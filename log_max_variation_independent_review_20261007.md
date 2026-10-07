# Log-max 乘性变分：独立解析审计及有限窗口正宽化

2026-10-07，gated_radial_energy_audit，指定 6.1-sol high。只读审根 [log_max_variation_20261007.md](log_max_variation_20261007.md) §1–5，被审 SHA-256 80b051af1eb23e22797f1f357481e939931dfe211d061b4d5f987452a41a38ca。不改根稿、不执行新/旧数值。已读 M03 SKILL、method、cube-interface、provenance入口，仅应用优化方向/亏损/量词审计流程，不引用未验外部稳定性定理。

**结论：核心不等式方向、常数、可积性、tie/阈值处理及量词均通过。** 来源势平坦只是保持质量乘性扰动的必要条件，不证明极大值存在、不约束平坦常数的维数阶、不转移 actual history 门。附节给出所有有限正测度与所有非负 L1 输入之间，在固定有限连续窗口上的 constant-preserving 正宽化；该结论不延伸到任意 structured source class 或全尺度。

## 1. joint usc、真实赢家与积分域

设 μ有限正 Borel、W>0，0<a≤b<∞，原闭 cube 核 h_R=R^{-n}1_{[-R/2,R/2]^n}。若 (x_j,R_j)→(x,R)，对每个来源 y 的闭 cube indicator 有 limsup≤原 indicator。R≥a 使归一化有统一有限界，反向 Fatou 给 h_R*μ 联合上半连续。因此每个 x 的 compact 参数最大值存在，最大化参数集合非空闭，最小赢家存在。

对任意 R<b，从右侧用有理 R_j↓R，闭 cube质量从上连续收敛，响应也收敛；补上 b 足够将 M写成可数上确界。最小赢家的 {R_*≤t}={M_t=M}，其中 M_t 是原允许族在≤t参数的最大值，确实 Borel。连续窗口用 [a,t]；有限族版必须将 M_t 限定到原有限族，不能另加未允许的连续 query。M=0处取a（有限族则取最小原参数）合法。

kernel包络的精确积分为
\[
 Z=1+n\log(b/a),\qquad
 \int M(x)dx\le ZW.
\tag{1}
\]
故 E={M>τ}有限、I=τ|E|≤ZW，Φτ=τ∫log_+(M/τ)有限，且 S(y)≤Z。对任意 |t|<1 的乘性扰动，
(1-|t|)μ≤μ_t≤(1+|t|)μ，故 M_t≤(1+|t|)M、Φτ(μ_t)有限。不涉及无限量减无限量。

正 Tonelli 给
\[
 \int S\,d\mu
 =\tau\int_E\frac{h_{R(x)}*\mu(x)}{M(x)}dx
 =I.
\tag{2}
\]
这是完整来源积分恢复 receiver Lebesgue 量，不把 μ 当作 dx。v 只需 μ-可测，可选等价 Borel代表定义有符号 vμ；零集不影响扰动。

## 2. 冻结旧赢家的下界与 Taylor 方向

固定允许 v，|v|≤1、∫v dμ=0；μ_t=(1+tv)μ，正性及质量保持。E上
\[
 p(x)=\frac{h_{R(x)}*(v\mu)(x)}{M(x)},\qquad |p|\le1.
\]
原 query 对新源的响应准确是 M(x)(1+tp(x))，所以新最大值≥该值。即使该值跌到τ以下，仍可用 log_+(M_t/τ)≥log[M(1+tp)/τ]；右侧允许负值。E外原 log_+为0、新项非负。因此
\[
 \Phi_\tau(\mu_t)-\Phi_\tau(\mu)
 \ge\tau\int_E\log(1+tp)dx.
\tag{3}
\]
这一步不要求唯一赢家、不要求最大值的导数、不要求阈值 plateau零测。1+tp≥1-|t|>0，使右侧有限。

沿 s∈[-|t|,|t|]，log(1+sp) 的二阶导数为
-p²/(1+sp)²≥-p²/(1-|t|)²。Taylor 下界准确给
\[
 \Phi_\tau(\mu_t)-\Phi_\tau(\mu)
 \ge tD-\frac{t^2Q}{2(1-|t|)^2},
 \quad D=\int vS\,d\mu,\quad
 Q=\tau\int_Ep^2dx\le I.
\tag{4}
\]
有符号 Fubini 的绝对值预算≤∫|v|S dμ≤I；不存在条件收敛交换。根常数 1/2 及平方分母正确，Q无需正下界。

## 3. local max 与 near max 的精确量词

local max假设仅要求：对每个允许 v，存在相应小邻域，使 ±t 两侧都不增 Φ；邻域可以依赖 v。固定 v 令 t→0±，(4)强制 D=0。这不是假设所有扰动共享一个未证邻域。

记 c=I/W，m_s=W^{-1}∫sgn(S-c)dμ，选
\[
 v=\frac{\operatorname{sgn}(S-c)-m_s}{2}.
\]
它属于允许域，并准确有
\[
 D=\frac12\int|S-I/W|d\mu.
\tag{5}
\]
因此每个已选合法真赢家的 S 都等于 I/W，μ-a.e.。v依当前未扰动输入和该赢家选取，local max对全部允许v的量词包含它；不存在后验越权。结论仅在原来源支撑上，不限制 S 于 μ零测的receiver/source位置。

定量版固定0<t<1，必须对这个 v 的两方向都满足
Φτ(μ_{±t})≤Φτ(μ)+εW。则
\[
 |D|\le\frac{\epsilon W}{t}
       +\frac{tI}{2(1-t)^2}.
\tag{6}
\]
若有全体允许 v 的这个上界，代入 (5) 给
\[
 \int|S-I/W|d\mu\le
 \frac{2\epsilon W}{t}+\frac{tI}{(1-t)^2}.
\tag{7}
\]
全局 ε-near-max 在同一固定τ、同质量完整输入类上会提供上述扰动上界；一般输入没有这种假设。local ε上界若只在 |t|≤t0 有效，优化t不得越过t0。

未知 I 明确保留；用已有 I≤ZW也只能恢复窗口基线，不变成√n。平坦 S=c给 ∫S²dμ=cI 是恒等式，Cauchy在此正好取等，不能用它循环证明c小。没有证明平坦势反推 local/global max、极大值存在或输入接近张量/格点/特定极值形状。

## 4. Φ 与 weak 常数：层蛋糕的两个方向

对同一输入类及同一原尺度族，非负 Tonelli准确给
\[
 \Phi_\tau=\tau\int_\tau^\infty
          |\{M>u\}|\,\frac{du}{u}.
\]
若 C_w 是 strict weak 常数，用 |{M>u}|≤C_wW/u，积分给 Φτ≤C_wW。因此 KΦ≤C_w。

任意 q>1，Φτ≥τlog q·|{M>qτ}|，令 λ=qτ 后得
C_w≤(q/log q)KΦ。q/logq 最小值为e，于q=e取得。所以
\[
 K_\Phi\le C_w\le eK_\Phi
\tag{8}
\]
正确。strict level边界不改变方向，不必存在任何极大输入。全尺度版若先有有限C_w，层蛋糕仍可付Φ；但不能据有限window的 ∫M 强可积，自动传全尺度。

## 5. 常数保持的 mollification：全有限 μ 与全 L1 类

以下是独立自含附节，域严格为固定0<a≤b<∞、完整连续窗口、所有有限正Borel测度与所有非负L1函数。取非负unit mollifier ηε，support⊂[-ε,ε]^n，形成approximate identity；可选平滑unit bump的尺度族。令 με=μ*ηε、Wε=W、Mε=M_{[a,b]}με。

### 5.1 同源上界与闭面的逐点下界

正卷积及sup移到积分外给
\[
 M_\varepsilon(x)\le(M*\eta_\varepsilon)(x).
\tag{9}
\]
M∈L1，由approximate identity，右端→M 于Lebesgue points且在L1中收敛。

对任何固定x、原R<b，当 ε充分小使R+2ε≤b，source y∈Q_R(x) 的全部 mollified mass都落在Q_{R+2ε}(x)。所以
\[
 M_\varepsilon(x)\ge
 \left(\frac R{R+2\varepsilon}\right)^n h_R*\mu(x).
\tag{10}
\]
这是原闭面上的真实atom也完整捕获，不用错误的“原R kernel点态平滑保持闭面质量”。若原真赢家R(x)<b，(10)对这个固定x的原R直接给liminf Mε≥M。

固定 b端不能越出参数窗。由非负Tonelli
\[
 \int dx\,\mu(\partial Q_b(x))
 =\int\mu(dy)\,|\{x:y\in\partial Q_b(x)\}|=0,
\]
故 μ(∂Q_b(x))=0 对dx-a.e. x。对这些x，h_b*με(x)→h_b*μ(x)，由固定核边界外的indicator收敛及μ有限的dominated convergence即可。因而赢家b的点也有liminf≥M。a=b时同一个固定端证明适用。

结合 (9) 的Lebesgue point上界，
\[
 M_\varepsilon\longrightarrow M\quad dx\text{-a.e.}
\tag{11}
\]
没有要求唯一赢家或统一winner gap，也没有在无数R上取一个错误的共同exceptional null set；内部(10)完全逐点，唯一需要fixed核null set的是b。

### 5.2 L1收敛及两个常数的互传

Fatou和(9)给
\[
 \int M\le\liminf\int M_\varepsilon
 \le\limsup\int M_\varepsilon
 \le\int M,
\]
因为正卷积右包络质量恒为∫M。故积分收敛；正函数的Scheffé/min argument给 ||Mε-M||1→0。自含可写
∫|Mε-M|=∫Mε+∫M-2∫min(Mε,M)，最后min≤M、a.e.收敛，dominated convergence。

函数 Fτ(u)=τlog_+(u/τ) 在u≥0为1-Lipschitz，故
\[
 |\Phi_\tau(\mu_\varepsilon)-\Phi_\tau(\mu)|
 \le\|M_\varepsilon-M\|_1\longrightarrow0.
\tag{12}
\]
且严格水平集满足1_{M>λ}≤liminf1_{Mε>λ}，因此
λ|{Mμ>λ}|≤liminf λ|{Mμε>λ}|。
W保持不变，任何有限μ的Φ/weak ratio可由L1输入以所需方向逼近；L1又是有限μ的子类。结论是这个**完整输入类、固定window**的 KΦ 和 Cw 在两种输入类中相等，不引入常数损失。

### 5.3 finite family、structured class 与未被证明的迁移

若原允许尺度为共同有限集合，不能在(10)使用新R+2ε。此时每个固定R的 h_R*μ∈L1，正宽化响应=h_R*μ*ηε，于固定kernel的Lebesgue points收敛；有限交集即使所有原responses同时收敛。有限max再配(9)与上面L1/min argument，仍得constant-preserving transfer。

以上不保持某个固定atom-count、格点分辨率、packet形状、张量、码或其它structured sourceclass；不能对这些类称无代价迁移。它不证明极值存在、near-max序列紧性、原FIRST门随正宽化稳定或全尺度 strong L1 收敛。尤其来源势S本身的pointwise正宽化/winner选择稳定未在这里证明，仅Φ/weak常数得到互传。

## 6. 不需极值存在或紧化的全局 near-max reduction

根随后提出的加强版本也严格成立，必须与§3的局部、固定质量假设区分。令完整输入类的 K=KΦ≤Z。任取全局 ε-near-optimal pair (μ,τ)，满足
\[
 \Phi_\tau(\mu)\ge(K-\epsilon)W.
\]
现在允许任意 μ-可测 |v|≤1，不要求 ∫vμ=0。μ_t=(1+tv)μ仍非负有限，|t|<1，质量准确为
\[
 W_t=W+t\int v\,d\mu.
\]
§2的下界不使用mean-zero，仍成立；全局常数定义给
Φτ(μ_t)≤K W_t。相减准确得到
\[
 t\int v(S-K)d\mu
 \le\epsilon W+\frac{t^2Q}{2(1-|t|)^2}.
\tag{13}
\]
这里不能继续把新mass写W，但无需重新归一化μ_t；K乘Wt正好补上这一质量变动。

固定t>0选v=sgn(S-K)，无mean subtraction，得到
\[
 \boxed{\int|S-K|d\mu
 \le\frac{\epsilon W}{t}
       +\frac{tI}{2(1-t)^2}
 \le\left[\frac\epsilon t+
          \frac{tZ}{2(1-t)^2}\right]W.}
\tag{14}
\]
它是真正全局 near-optimal quotient 的source偏差证书，不是仅假设某μ是local max。由 ∫Sμ=I，
\[
 |I/W-K|\le W^{-1}\int|S-K|d\mu.
\tag{15}
\]

对固定n、固定window，Z固定有限。仅由sup定义选ε_j↓0的near-optimal pairs(μ_j,τ_j)，每个pair内使用自身原truewinner和fixedτ_j的扰动；τ_j可以不同，sources不需要收敛。取t_j=√ε_j≤1/2即使(14)右侧趋于0。因此存在source势在自身μ_j下归一化L1趋近K的near-optimal序列，并有I_j/W_j→K。不需要极大值存在、紧性、形状极限或winner一致性。ε=0若真global maximizer存在，令t↓0则S=K、μ-a.e.。

若需显式跨n预算，而不是暗用fixedn极限，在0<ε≤Z/4时可取t=√(ε/Z)≤1/2，从(14)得到
\[
 W^{-1}\int|S-K|d\mu\le3\sqrt{\epsilon Z}.
\tag{16}
\]
只有登记ε_nZ_n→0才保证这一跨n偏差趋于0；ε固定不自动给uniform flatness。式(16)只是已证一般亏损/偏差联系，仍没有控制未知K_n。

§5的constant-preserving正宽化允许将sup定义限定到所有L1输入，再选这类near-optimal序列；无需先把atomic势S的极限转移到平滑输入。不过这仍不保持structured family或actual依赖μ的history门。

## 7. 本路线的剩余问题

完整 Φ 的局部变分证书与§6的全局near-max reduction均是正确的一般接口。后者已经不需要紧化来得到 near-optimal sources 的归一化势偏差，但要付√n，仍需证明这些输入的真实共同winner/source几何能控制未知K、I或Φ；不能从接近未知常数K的flatness直接推阶数。门依赖μ的 actual FIRST/CPGP/LCA/high-source/anchor交通不等于本Φ的导数；需重新计门变化，不能把完整S flatness移到受限交通。

本审计未执行 root计划的变分probes、未审其数值结果，也不注册额外试验。解析审计通过范围只限上述根稿快照及自含finite-window transfer。

# 完整 log-max 的加性来源障碍：精确 KKT、弱空间异常与点态容量限制

2026-10-07；gated_radial_energy_audit，用户指定 6.1-sol high。只写本稿，不改总稿/主账，不运行数值。应用已读 M03 的优化方向、亏损、归一化及量词审计；已读 SKILL、method、cube-interface、provenance，不引用未读外部稳定性定理。查重对象为 log_max_variation_20261007.md、其 independent_review，以及 source_spike_background_tail_audit：乘性来源平坦性和有限窗口正宽化沿用已核结论；本稿自证加性条件及其稳定范围。这里没有证明 K 的平方根维数阶，也没有证明极值存在。

## 1. 完整输入、真赢家与来源势

固定维数 n、全边长窗口 0<a≤b<∞。闭方体 Q_R=[−R/2,R/2]^n，h_R=R^{-n}1_{Q_R}。完整有限正 Borel 测度 μ，质量 W>0；令
\[
 M(x)=\max_{a\le R\le b}h_R*\mu(x),\quad
 E=\{M>\tau\},\quad I=\tau|E|,\quad
 \Phi_\tau(\mu)=\tau\int\log_+(M/\tau)\,dx,
\]
\[
 S(y)=\tau\int_E\frac{h_{R(x)}(x-y)}{M(x)}\,dx,\qquad
 Z=1+n\log(b/a).
\tag{1}
\]
R(x) 是任一可测真正赢家；联合 usc、紧参数最大值和最小赢家可测性已经在前述根稿/独审稿自证。闭面的真实质量保留，不假设唯一赢家或零测 threshold plateau。固定窗口包络给
\[
 0\le S\le Z,\quad \int M\le ZW,\quad
 I\le ZW,\quad \int S\,d\mu=I.
\tag{2}
\]
所有观察点积分始终是 Lebesgue 测度。令 K=sup_{\mu,\tau}\Phi_\tau(\mu)/W，sup 取所有完整非零有限正来源和 τ>0。本稿扰动均在这个完整类内；有限测度与所有非负 L1 类的同一常数相等，使用前述 independent_review §5 的 constant-preserving 正宽化。此等价不保持任意 structured source class。

还可把旧 K≤Z 改成绝对常数更好的
\[
 \boxed{K\le Z/e.}
\tag{3}
\]
因为对 u≥0，τlog_+(u/τ)≤u/e。它只是基线常数改善。K>0：取越来越大的盒内常数密度 1、τ=1/e，完全内接收盒上 M=1，Φ/W 的下极限≥1/e。因此后文可在固定维数取 0<ε<K/4。

## 2. 任意正向添加来源的精确下界

令 ν≥0 有限、V=ν(R^n)>0，t≥0；添加后的完整输入是 μ+tν，质量 W+tV。在原 E 上冻结原真正赢家，定义
\[
 p_\nu(x)=\frac{h_{R(x)}*\nu(x)}{M(x)}\ge0,\qquad
 D_\nu=\tau\int_Ep_\nu\,dx=\int S\,d\nu,\qquad
 Q_\nu=\tau\int_Ep_\nu^2\,dx.
\tag{4}
\]
正 Tonelli 给 D 的等式；pν≤V/(τa^n)，故
\[
 Q_\nu\le\frac{V}{\tau a^n}D_\nu
 \le\frac{ZV^2}{\tau a^n}<\infty.
\tag{5}
\]
新最大值≥M+t h_R*ν，原 E 外新的 log 项非负。标量 log(1+z)≥z−z²/2 对全部 z≥0 成立：二者差的导数是 z²/(1+z)≥0。因此
\[
 \boxed{\Phi_\tau(\mu+t\nu)-\Phi_\tau(\mu)
 \ge tD_\nu-\frac{t^2Q_\nu}{2},\qquad t\ge0.}
\tag{6}
\]
无需 t<1、无需赢家导数、无需新赢家等于旧赢家。两边均有限；正扰动无负部问题。

若 Φτ(μ)=KW 是精确全局极值，则上界≤KtV；固定任意 ν，在 (6) 中令 t↓0，得 Dν≤KV。特别允许 ν=δ_y，故
\[
 \boxed{S(y)\le K\quad\text{每个 }y\in\mathbb R^n.}
\tag{7}
\]
上一轮任意有界乘性扰动、同时保留质量变化的一阶式还给 S=K，μ-a.e.；不能只从 (7) 与 ∫Sμ=I 自动推 I=KW，需这个双向乘性必要条件。

S 是 usc。若 y_j→y，对固定 x 的闭方体指示 limsup≤原指示；积分权 τ/[M(x)R(x)^n]≤a^{-n}，且 |E|<∞，反向 Fatou 给 limsup S(y_j)≤S(y)。所以 {S<K} 开，(7) 与 μ-a.e. 接触共同推出
\[
 S=K\quad\text{在完整 }\operatorname{supp}\mu\text{ 上}.
\tag{8}
\]
另外 S(y)≤a^{-n}|E∩(y+Q_b)|→0 当 ||y||∞→∞，由有限 |E| 的空间尾即可；故如果精确全局极值存在，其来源支撑有界。这里是必要条件，既不构造极值，也不推出其几何形状。

## 3. 近极值：测试来源的二阶容量

假定 Φτ(μ)≥(K−ε)W。全局上界对 μ+tν 给
\[
 t(D_\nu-KV)\le\epsilon W+\frac{t^2Q_\nu}{2}.
\tag{9}
\]
若 ε>0、Qν>0，任意 t>0 可取，优化 t=√(2εW/Qν)，得到
\[
 \frac{D_\nu}{V}\le K+\sqrt{2\epsilon\chi_\nu},
 \qquad \chi_\nu=\frac{WQ_\nu}{V^2}
 \le Zd,\quad d=\frac{W}{\tau a^n}.
\tag{10}
\]
Qν=0 时 Dν=0，单独处理；ε=0 用 t↓0。保留 (5) 的第一步，记 xν=Dν/V，则
\[
 x_\nu\le K+\sqrt{2\epsilon d x_\nu}
 \le K+\epsilon d+\sqrt{\epsilon d(2K+\epsilon d)}.
\tag{11}
\]
还可与 xν≤Z 取最小。ν=δ_y 给每个 y 的同一界。因此 ε_jd_j→0 是全空间一致上障碍误差趋零的一项充分条件；只说 ε_j→0、而允许完整输入和 τ_j 改变，不足以保证它。

更一般地，令 qν=ess sup_E pν，定义 dν=Wqν/V，则 (11) 的 d 可以改成 dν，因为 Qν≤qνDν。若 ν/V≤Cμ/W 作为测度，pν≤CV/W，故 dν≤C；这给受限测试类的输入尺度无关误差。或者 ν 有密度≤Cτ 且 V≥θW，则 dν≤C/θ。两种条件均不能免费赋给任意点来源 δ_y。

## 4. 不含 d 的弱 Lebesgue 障碍

在来源位置上用 dy 而不是 dμ 积分，正 Tonelli 精确给
\[
 \int S(y)\,dy
 =\tau\int_E\frac{dx}{M(x)}
 \le |E|\le\frac{ZW}{\tau}.
\tag{12}
\]
任意 u>0，令 A={y:S(y)>K+u}。S 为 Borel，(12) 给 |A|<∞。若 |A|=0 结论直接成立；否则选 ν=τ1_A dy，V=τ|A|。记 xν=Dν/V=|A|^{-1}∫_A S，严格有 xν>K+u。因为 h_R 是保质量正核，
\[
 h_R*\nu\le\tau,\qquad
 0\le p_\nu\le\tau/M<1,\qquad
 Q_\nu\le D_\nu=x_\nu V.
\tag{13}
\]
将 (13) 代入 (9)，取 t=(xν−K)/xν>0，得到
\[
 V\le\frac{2\epsilon W x_\nu}{(x_\nu-K)^2}.
\]
函数 x/(x−K)^2 在 x>K 上递减，其导数为 −(x+K)/(x−K)^3，故
\[
 \boxed{\tau\,|\{y:S(y)>K+u\}|
 \le \frac{2\epsilon W(K+u)}{u^2},\qquad u>0.}
\tag{14}
\]
这是完整近极值来源势的一般正确弱空间约束，右侧没有 a^{-n} 或 d。测试 ν 依赖原 S 是合法的：全局 K 对所有正 ν 的扰动均给上界，原赢家仍冻结。它控制来源位置的 Lebesgue 异常体积，不能改成 μ(A) 或原 receiver 体积，也不能将 K 当成已知小常数。u≥Z−K 时异常集为空；ε=0 的 Lebesgue 零异常还弱于 (7) 的逐点结论。

## 5. 真实 hot-point 组件：完整背景不能省略

下面仅审计能否从 ε-near-max 无条件推出 sup_y(S(y)−K)→0；不是极大函数弱型反例，也不是 actual FIRST/CPGP/history 剩余模型。

固定 n,a,b；取 η=a^n/10，阈值 τ>0，完整来源
\[
 \mu_{\rm bad}=\tau\bigl(1_{Q_{2b}}\,dy+\eta\delta_0\bigr),
 \qquad W_{\rm bad}=\tau[(2b)^n+\eta].
\tag{15}
\]
当 x∈Q_b，每个允许 Q_R(x)⊂Q_{2b}，背景响应准确为 τ；尖峰可捕获当且仅当 R≥2||x||∞。严格超阈集合恰为 Q_b（边界对 dx 无碍），其真实赢家唯一为
\[
 R_0(x)=\max(a,2\|x\|_\infty),\qquad
 M_{\rm bad}(x)=\tau(1+\eta/R_0(x)^n).
\tag{16}
\]
因为捕获前响应仅 τ，捕获后信号 η/R^n 严格下降。Q_b 外任何 query 捕获不到尖峰，而背景响应≤τ，故没有遗漏完整 E。

使用径向体积 d|Q_R|=nR^{n-1}dR，中心来源势准确为
\[
 S_{\rm bad}(0)
 =\frac{a^n}{a^n+\eta}
  +\log\frac{b^n+\eta}{a^n+\eta}
 \ge Z-\frac{2\eta}{a^n}=Z-\frac15.
\tag{17}
\]
最后一步用 a^n/(a^n+η)≥1−η/a^n、log(1+η/a^n)≤η/a^n，保留非负 log(1+η/b^n)。由 e>2 和 (3)，
\[
 S_{\rm bad}(0)-K>Z/2-1/5\ge3/10.
\tag{18}
\]
所有背景质量 τ(2b)^n 都进入 Wbad。实际 Φbad 也是完整的 τ∫_{Q_b}log(1+η/R_0^n)dx；不把背景免费化，不要求这个 bad 组件本身近优。

### 5.1 自含的正宽度 L1 版本

可把 ηδ0 换成 η(2ℓ)^{-n}1_{Q_{2ℓ}}dy，保持相同完整背景和质量，明确取 ℓ=a/(32n)。不援用一个未证的 S 正宽化连续性。

对 x∈Q_b，背景仍恒定。因为 R≥a>2ℓ，每个轴的 packet 捕获比例准确为
\[
 F_i(R)=\left[\frac{R-2|x_i|+2\ell}{4\ell}\right]_{[0,1]},
 \qquad m_\ell(R)=\eta\prod_i F_i(R).
\]
在正质量而未完整捕获的每个可微参数段，至少一个 F_i∈(0,1)，因而
\[
 \frac{m_\ell'(R)}{m_\ell(R)}
 =\sum_{0<F_i<1}\frac1{R-2|x_i|+2\ell}
 \ge\frac1{4\ell},\qquad
 \frac{d}{dR}\log\frac{m_\ell(R)}{R^n}
 \ge\frac1{4\ell}-\frac nR>0.
\tag{19}
\]
有限分段的连续性连接这些严格上升段；完整捕获后信号严格下降。故对 r=2||x||∞≤b−2ℓ，真正赢家唯一为 Rℓ=max(a,r+2ℓ)，它捕获整个 packet。此论证也允许多个轴同时 partial，不假设最大坐标唯一。每个 packet 来源 y∈Q_{2ℓ} 都由这些同一真实赢家完整捕获。因此对所有这样的 y，
\[
 S_\ell(y)\ge
 \frac{(a-2\ell)^n}{a^n+\eta}
 +n\int_a^b\frac{(R-2\ell)^{n-1}}{R^n+\eta}\,dR
 \ge(1-2\ell/a)^n S_{\rm bad}(0)
 \ge\frac{15}{16}S_{\rm bad}(0).
\tag{19a}
\]
第二步用 R≥a、(1−2ℓ/a)^{n−1}≥(1−2ℓ/a)^n；最后是 Bernoulli (1−1/(16n))^n≥15/16。a=b 时积分为空而同一 core 仍给该下界。由 (18) 的前两步，
\[
 S_\ell(y)-K>
 \frac{15}{16}(Z-1/5)-\frac Z2
 =\frac{7Z-3}{16}\ge\frac14 .
\tag{19b}
\]
这是固定显式正宽度、完整 L1 组件的真实赢家证书，不只是 atomic limit。

完整 L1 输入的 E 可能因 packet 宽度在 Q_b 外略扩展；本论证没有假称它等于原 Q_b，也没有删除这些输出，只以非负性得到足够的 score 下界。背景完整质量仍不变。该 L1 版本仅是点态稳定性反例组件的解析升级，不认证 actual 门。

## 6. 同阈值 compact 近优来源与有限稀释

先证明无需极值存在即可选 compact 近优组件。由 sup 定义，取完整来源 μ0、阈值 τ，使 Φτ(μ0)≥(K−ε/4)W0。令 μA 为其大闭盒截断，尾来源 λ=μ0−μA。正性和 sup 的次可加性给
\[
 0\le M\mu_0-M\mu_A\le M\lambda,\qquad
 \|M\lambda\|_1\le Z\lambda(\mathbb R^n).
\]
函数 τlog_+(t/τ) 是 1-Lipschitz，故 Φ 的损失≤Zλmass。选截断使 λmass≤εW0/(4Z)，则
\[
 \Phi_\tau(\mu_A)\ge(K-\epsilon/2)W_0
 \ge(K-\epsilon/2)W_A,\qquad W_A>0.
\tag{20}
\]
若要求所有输入 L1，先按前述同常数正宽化在 L1 类选择 μ0，再截断；不从近优推出任何紧性或形状。

取 N 个 μA 的平移副本和一个 (15) 的 bad 组件；同一 τ 保持不变。沿第一坐标将所有来源支撑盒之间距离取严格大于 b，确保其 b-query receiver halos 两两不交。任何允许 query 都不能同时捕获两个组件；因而新的完整最大值、严格 E 和 Φ 是真实的分离组件之和。它不是抽象直和，也没有按观察点重启 FIRST。bad 来源点 0 的 S 不变，因为其 halo 不与其它组件相交。

完整质量 WN=NWA+Wbad。由 (20) 和 Φbad≥0，
\[
 \frac{\Phi_\tau(\mu_N)}{W_N}
 \ge K-\epsilon/2-K\,\frac{W_{\rm bad}}{W_N}.
\tag{21}
\]
例如取有限正整数
\[
 N\ge\max\!\left(1,\left\lceil
 \frac{2ZW_{\rm bad}}{\epsilon W_A}\right\rceil\right)
\tag{22}
\]
即可使 (21)≥K−ε，且 S_N(0)>K+3/10（使用正宽度版本则 >K+1/4）。ε 可以任意趋零，每次重新按 sup 定义选 compact 近优组件和有限 N。这里没有给未知近优 μA 指定数值形状、没有假称某个明确原子输入已近优。

这个家庭的完整容量参数准确为
\[
 d_N=\frac{N W_A}{\tau a^n}
      +\left(\frac{2b}{a}\right)^n+\frac1{10}.
\tag{23}
\]
(22) 的大副本数正是输入依赖不能免费消掉的机制。所有副本和完整背景质量都计入 WN；异常组件的来源质量占比趋小，与乘性 near-flatness 无冲突。也与 (14) 无冲突：异常来源位置的绝对体积可固定，但归一化 τ|A|/WN 随质量稀释变小。

因此，即使固定 n、固定窗口，所有完整 ε-near-max 输入也不能无条件有 sup_y S(y)≤K+o_ε(1)。它不否定 (14) 的弱 Lebesgue 障碍，不否定 source L1 平坦性，更不否定任何一般 √n 弱上界。

## 7. 适用边界与执行状态

加性 (6)、精确 (7)、近极值 (9) 和弱空间 (14) 同样适用于共同有限尺度族；此时 K 与 Z 可取该族的对应常数/包络，仍允许使用连续窗口 Z 的上界。但 (16)–(17) 的精确径向公式属于完整连续 [a,b]，不能把按 x 的 arrival 尺度目录冒充共同有限族。若另外验证真实有限网的 hot 组件，需要单独计算其真正 finite winner score；它是独立证书，不是这里的连续极限或误差认证。

上述扰动改变完整输入，原 FIRST、严格 CP/GP、出生、LCA、history 等门也可能随 μ 改变。本稿没有把完整 Φ 的导数直接指认为受门交通，也没有推出 actual R_angle 的费用。真正未解的是：如何用 (7)/(14) 与完整真赢家几何控制未知 K，或给实际门变化另立有效预算；S≈未知 K 不是 K 小的证明。

作者未执行新数值。radial_hard_band_pressure 另负责 log_max_additive_obstacle_probe_20261007 前缀，计划对显式 bad 组件常数、有限网独立 score 与加性测试接口做注册有理/区间守卫；未知 compact 近优组件只允许符号稀释验证，不冒称数值构造或 actualFIRST 样本。其注册、执行与结果以独立收据为准，不在本稿预先计通过。

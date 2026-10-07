# 真实赢家径向测试的全局化：合法拼接与线性维数量词损失

2026-10-07；gated_radial_energy_audit，用户指定 6.1-sol high。仅新增本稿，纯解析，不跑新/旧数值，不改主稿或他人文件。读取 nearflat_cube_geometry、source_fragment_localization、winner_envelope_entropy_bridge、log-max 两稿及独审；继续 M03 的近优量词审计，并读取 D01/F05 的 SKILL、method、cube-interface、provenance入口。仅使用静态同一来源概率、合法测试下界与自含 Fourier/Plancherel 估计；不调用未读高级输运或谱定理。

**结论。** 一个 receiver 局部径向测试确实给 N/atom-count 无关的物理 W1 下界，但不能把随 x 变化的测试直接送进 log-max 的全局来源扰动。固定 source blocks 上的共同测试可以无块数费拼接，需显式控制几何近似误差和捕获支撑。归一化到 shell 厚度的 bounded 径向测试还存在真实连续 cube、完整 L1 来源、唯一真赢家的量词损失：局部后验差≥3/40，而每个全局有界来源测试的 receiver 平均差≤36/n。因此一般转换至少损失线性 n，不能用 √n polylog(n) 代价替换。该例未认证 near-max，故不排除依赖完整近极值信息的专门转换，也不否定未归一化物理 W1 的更弱合同。

## 1. 原后验、近赢家和确实合法的径向下界

固定完整有限正 μ、共同连续窗口 [a,b]，h_R=R^{-n}1_{Q_R}，M=max h_R*μ，E={M>τ}。R(x) 是原可测真赢家；另取同一原窗口的可测 L(x)≥R(x)，正响应，写
\[
 q(x)=\frac{h_L*\mu(x)}{M(x)}=e^{-g(x)},\qquad
 \pi_{R,x}=\frac{\mu|_{Q_R(x)}}{\mu(Q_R(x))},\quad
 \pi_{L,x}=\frac{\mu|_{Q_L(x)}}{\mu(Q_L(x))}.
\tag{1}
\]
后验均来自同一完整 μ，来源坐标可以任意依赖。对每个原允许中间边长 R≤s≤L，真赢家给 μ(Q_s(x))≤Ms^n，因此
\[
 \pi_{L,x}\{2\|y-x\|_\infty\le s\}
 \le e^{g(x)}(s/L)^n.
\tag{2}
\]
连续窗口可用于所有这样的 s；finite family 只能对原目录中的 s 使用 (2)，不能免费积分整个中间区间。

令 d_x(y)=dist∞(y,Q_R(x))=(||y−x||∞−R/2)_+。它对 y 是 1-Lipschitz，π_R d_x=0。在 π_L 支撑上 0≤d_x≤(L−R)/2，层饼给
\[
 \pi_{L,x}d_x-\pi_{R,x}d_x
 \ge\frac12\int_R^L[1-e^{g(x)}(s/L)^n]_+\,ds.
\tag{3}
\]
任何 π_R、π_L 正耦合的平均 ∞-距离≥这个差，因为 d_x 是1-Lipschitz；无需假定输运势达到或一般强对偶。这是原后验物理 W1 的合法下界，与来源 atom 数无关。当 g<n log(L/R) 时右侧严格正。

若将它归一化成 admissible bounded test，
\[
 v_x(y)=\left[\frac{2d_x(y)}{L(x)-R(x)}\right]_{[0,1]},
\tag{4}
\]
其幅度≤1，但 Lipschitz 常数为 2/(L−R)。壳厚小则该常数很大。式(3)与式(4)不能不保留长度因子便互称同一 W1 合同。

## 2. log-max 变分的原量词不能交换

nearflat_cube_geometry 已证：若完整 Φτ 是 ε-near-global-max，则对每个**同一 source function** |v(y)|≤1，及每个合法候选 L(x)，
\[
 \tau\int_E[\,t|\pi_{L,x}v-\pi_{R,x}v|-g(x)\,]_+dx
 \le2\epsilon W+\frac{t^2I}{(1-t)^2},\quad
 0<t<1,\ I=\tau|E|.
\tag{5}
\]
v 可以依赖未扰动 μ 和所选 R,L，但必须是 y 的一个函数；对应扰动仍是 (1±tv(y))μ(dy)。不能改成 (1±tv_x(y))μ(dy) 并沿各 x 保留同一个全局来源。

尤其 (5) 本身不包含
\[
 \tau\int_E[\,t\sup_{|v|\le1}|\pi_{L,x}v-\pi_{R,x}v|-g(x)\,]_+dx.
\tag{6}
\]
sup 放在 receiver 积分内是另一个算子范数。将逐 row 输运势或(4)的局部 normalized distance 代入，须有另证的 global-test budget。自由 posterior 反例不是必要的：§4–6给真实原 cube 后验的失败证书。

## 3. 可支付的合法拼接：支撑块与近似误差

这是一项一般正确的转换，不许诺所有输入都有其所需支撑结构。取固定至多可数的两两不交 μ-可测 source sets F_j，和 receiver 可测不交子集 E_j⊂E。要求 E_j 上两后验都支撑于 F_j，即 π_R(F_j)=π_L(F_j)=1，dx-a.e.。为每个 j 选择一个共同 |v_j|≤1 的 source test。若局部测试满足
\[
 \sup_{y\in F_j\ {\rm relevant}}|v_x(y)-v_j(y)|\le\omega(x),
 \qquad x\in E_j,
\tag{7}
\]
则 v=Σ_j1_{F_j}v_j 是一个全局 |v|≤1 测试（其它来源处取0），且
\[
 |\pi_Lv-\pi_Rv|\ge|\pi_Lv_x-\pi_Rv_x|-2\omega(x).
\tag{8}
\]
将 (8) 送入原(5)，在 union E_j 上仍得同一右侧；不付 j 个数、不重领 W。可数性和有界性保证正 Tonelli/可积性，不要求每块各解新的极值问题。

对共同宽度 d>0 的 clipped distance test v_x=min(1,d_x/d)，若块中用中心 x_j、半径 R_j，则
\[
 \omega(x)\le
 \frac{\|x-x_j\|_\infty+|R(x)-R_j|/2}{d}.
\tag{9}
\]
“小 receiver block”必须相对真实 shell 宽度 d 小，而不是仅相对 a 小。若宽度也变化，还须另扣宽度误差。

若 F_j 不互斥，且每个 source 最多属于 D 个 F_j，辅助独立 signs ξ_j 可给有偿版本：vξ=D^{-1}Σ_jξ_j1_{F_j}v_j，仍 |vξ|≤1。对 x∈E_j，二/四矩与 Hölder 给
\[
 \mathbb E_\xi|\pi_Lv_\xi-\pi_Rv_\xi|
 \ge\frac1{\sqrt3D}
 |\pi_L(1_{F_j}v_j)-\pi_R(1_{F_j}v_j)|.
\tag{10}
\]
hinge 凸性把该下界送入(5)的辅助平均。每个 ξ 实现都是合法全局扰动，不是 receiver 重选符号。明确复杂度是 D，不是免费 √D。有限系统可直接展开；可数且点态有限重叠时用截断、dominated convergence 得同式。若还满足 E_j 上 posterior 支撑 F_j，则(10)右侧就是共同 v_j 的差。

旧 source_fragment_localization 的 b-邻接组件确有每 row 唯一标签，故支撑拼接可无组件数费；但组件可很长，不能因此推出(7)的中心/半径近似误差小。这是与旧局部化的准确差异。

## 4. 一个 dimension-free L2 的 hard dilation 差估计

反例用到的谱估计完整自证。置 a(t)=sin(t)/t、a(0)=1，b(t)=t a'(t)=cos(t)−a(t)。对全部实 t，
\[
 |b(t)|\le13(1-|a(t)|).
\tag{11}
\]
证明：对 |t|≤1，a≥0，交错 Taylor 给 1−a≥19t²/120；而 sin t−t cos t=∫_0^t s sin s ds 给 |b|≤t²/3，比值≤40/19。对 1≤|t|≤2，sinc 在 (0,2) 上递减，a(t)≤a(1)≤101/120；对 |t|≥2，|a(t)|≤1/2。所以 |t|≥1 时 1−|a|≥19/120，且 |b|≤2，比值≤240/19<13。t=0 两侧为0。

cube 的 Fourier multiplier 为 m_R(ξ)=Π_i a(Rξ_i/2)，故
\[
 |\partial_{\log R}m_R(\xi)|
 \le13\sum_i(1-|a_i|)\prod_{j\ne i}|a_j|\le13.
\tag{12}
\]
最后的有限非负和是独立 Bernoulli 恰一次 failure 的概率；它≤1仅为标量乘积恒等式，不把真实来源坐标认作 independent。积分参数并用 Plancherel，
\[
 \boxed{\|h_L-h_R\|_{L^2(dx)\to L^2(dx)}
 \le13|\log(L/R)|\quad\text{对全部 }n.}
\tag{13}
\]
没有逐轴累加 n，没有将 Schur 范数混作逐项乘子。它只估单个 dilation 差，不是 hard 最大函数 L2/弱型付款。

## 5. 唯一真赢家、完整 L1 来源的局部强分离

固定 n≥2，原完整连续窗口 [1,2]。取 H=n²+2、ε0=1/4，完整输入
\[
 f_H(y)=\left(1-\frac{\epsilon_0\|y\|_2^2}{nH^2}\right)
 1_{[-H,H]^n}(y),\qquad
 \mu=f_Hdy,\quad W=\frac{11}{12}(2H)^n.
\tag{14}
\]
0≤f≤1，盒内 f≥3/4。取 τ=1/2，完整 E={M>τ}；每个平均≤1，因此完整 E 符合原幅度带 τ<M≤2τ。下面使用真实 receiver 子集
\[
 E_0=[-H+1,H-1]^n\subset E,\qquad V_0=|E_0|=(2H-2)^n.
\]
对 x∈E0，所有原允许 Q_R(x) 都完全在来源盒内，完整平均准确为
\[
 u_R(x)=1-\frac{\epsilon_0\|x\|_2^2}{nH^2}
                 -\frac{\epsilon_0R^2}{12H^2}.
\tag{15}
\]
严格随 R 下降，所以原唯一真赢家是 R=1，无 tie 人为选择。候选 L=exp(1/n)∈[1,2]，真实 g=log(u1/uL)>0，并有
\[
 0<g\le\frac{L^2-1}{36H^2}.
\tag{16}
\]
因为 uL≥3/4，log(u1/uL)≤(u1−uL)/uL。E0 外保留完整 M 和合法原赢家；候选设为原 R、局部测试设0即可。完整 E 没有被删成 E0，但此处量词转换只在 E0 上受检验。边界响应无须忽略：(15)只对 E0 使用，外部不套内盒公式。

取(4)的局部测试
\[
 v_x(y)=\left[\frac{2\|y-x\|_\infty-1}{L-1}\right]_{[0,1]}.
\tag{17}
\]
π1 v_x=0。Q_L(x) 上 f≥3/4，捕获质量≤L^n，因此 πL v_x≥(3/4) times uniform-cube expectation。令 ρ=2||Y−x||∞，uniform Q_L 的 CDF 是 (s/L)^n，故
\[
 \pi_Lv_x\ge\frac{3}{4(L-1)}
 \int_1^L[1-(s/L)^n]\,ds\ge\frac3{40}.
\tag{18}
\]
在 s≤(1+L)/2 的前半段，(s/L)^n≤[(1+e^{-1/n})/2]^n≤(1−1/(4n))^n≤e^{-1/4}≤4/5。
使用 e^{-z}≤1−z/2 对0≤z≤1（Taylor 二阶上界），以及 e^{1/4}≥5/4。积分前半段长度为 (L−1)/2，给 uniform expectation≥1/10。没有数值拟合。

所以每个 E0 receiver 的 normalized radial test 差至少3/40；全部原(2)真实中间尺度约束也成立。这不是 arbitrary posterior 模型。

## 6. 每个全局来源测试却只有 O(1/n) 平均差

任意 μ-可测 |v|≤1，令 F=f_Hv 并在来源盒外取0；||F||2≤(2H)^{n/2}。E0 上精确分解
\[
 \pi_Lv-\pi_1v
 =\frac{(h_L-h_1)*F}{u_L}
   +(h_1*F)\frac{u_1-u_L}{u_Lu_1}.
\]
uL≥3/4、|h1*F|≤u1；由(13)和 Cauchy–Schwarz，
\[
 \frac1{V_0}\int_{E_0}|\pi_Lv-\pi_1v|dx
 \le\frac{52}{3n}\sqrt{\frac{(2H)^n}{V_0}}
      +\frac{L^2-1}{36H^2}\le\frac{36}{n}.
\tag{19}
\]
最后 (H/(H−1))^{n/2}≤exp[n/(2(H−1))]<2，L²≤e<3，H=n²+2≥n+1。36 是统一粗常数。第一个比值完整保留了 source-box/receiver-inner-box 的边界体积损失，未先用无限常数背景或 torus 替代；此界对全部 global v 同时成立，不是某一个 Fourier mode 的失败。

若要求对所有真实输入/赢家/候选存在统一转换常数 Cn，使
\[
 \int_{E_0}|\pi_Lv_x-\pi_Rv_x|dx
 \le C_n\sup_{|v|\le1}
       \int_{E_0}|\pi_Lv-\pi_Rv|dx,
\tag{20}
\]
那么(18)–(19)强制
\[
 \boxed{C_n\ge n/480.}
\tag{21}
\]
√n polylog(n) 的普遍转换不成立。例子是完整原连续中心方体、唯一真赢家、完整 L1、真实幅度带和真正捕获后验，无 source-count/有限网漏洞。

hinge 也有相同损失。取 t=1/4，local hinge 的 E0 平均≥3/160−sup_E0 g；任意 global v 的 hinge 平均≤t·36/n=9/n。n充分大时(16)的 g 很小，因此 local 常数级与 global O(1/n) 不能交换。该例没有认证 f_H 是 ε-near-max；它排除真实赢家/小 gap 本身自动全局化的步骤，不违反附带完整近极值前提的专门定理。原(5)若应用于某个实际 nearmax source，则还须实质利用该额外前提。

### 6.1 本例没有否定的范围

(17)测试的 Lipschitz 常数2/(L−1)≈2n；改回1-Lipschitz物理距离时，差须乘 (L−1)/2=O(1/n)。本例不证明未归一化 W1 有常数阶局部分离却 global bounded tests 完全不可见，也不排除 O(√n) 的某个更弱物理运输转换。它精确排除的是：把 N无关物理下界按薄壳厚度归一化，再无代价当成一个全局 |v|≤1 扰动。

例子的 FIRST/CPGP/LCA/history 门未认证，不是 actual geom 余项反例。完整 μ、M、连续窗口和 E0 真赢家已经认证；nearmax 与其它依 μ 的原门仍是 unknown，不将 unknown 当通过。

## 7. 可行接口与最小缺口

局部径向(3)是正确一般事实，(8)是合法共同 source 拼接，(10)给明确重叠复杂度。解析反例(21)说明，仅有真实 winner/intermediate-scale dominance/很小 response gap/完整幅度带，不保证 normalized local test 能以 √n polylog 复杂度全球化。

继续使用 log-max 原(5)需要实质证明近极大输入还有一种低复杂度结构：足够多 receiver 的局部测试能在同一 source block 近似一致，且误差/重叠预算可付；或者另一真正共同 source-test family 能捕捉所需后验差。不能通过逐 receiver 单选一个势、细化 atom partition 或删除量词获得。

已读 winner_envelope_entropy_bridge 的完整 KL/source Palm预算；它只有同一总 log(Z/m)，没有替所有 local 输运势付选择复杂度。旧 component localization 消除了远离重复组件的数量费，却允许组件任意长，没有(7)的中心一致性。当前仍无一般 √n 空间付款，不给主账新增 paid支。

本稿只有解析证明，无新数值登记/执行，未重跑 nearflat、log-max 或 localization probes。


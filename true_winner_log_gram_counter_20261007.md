# 真实连续赢家下 volume-only log Gram proxy 的反例

本稿独立核验根提供的输入。结论：即使 n=1、共同输入为非负 L1 函数、尺度连续取 [1,2]、各计入 receiver 区间的赢家严格唯一，以下放大后的 log Gram proxy 也没有仅依赖维数的统一费用。反例不否定精确 source occupation 的平方预算；精确 occupation 在同一输入上有常数界。它不认证原 actual FIRST、CP/GP、出生、唯一 LCA 或 history 余项资格，也不作为完整 geom 反例。

仅使用 L03 的预登记、有理残差和独立证书方法；本稿解析证明自含，无外部定理迁移。连续尺度是主合同，未宣称共同有限尺度集版本。

## 1. 完整原输入及正质量

记 Q_L(x)=[x-L/2,x+L/2]、h_L=L^{-1}1_{Q_L(0)}。固定整数 J≥2，
\[
 h=\frac1{100J},\quad s_i=ih,\quad
 m_i=\frac1{J(1+s_i)},\quad
 z_i=\frac{(9/10)s_i-15s_i^2}{1+s_i}
 \quad(0\le i\le J).
\tag{1}
\]
取 α_i=m_i-m_{i+1}>0 (i<J)，α_J=m_J>0。于是
\[
 \rho=\sum_{i=0}^J\alpha_i\delta_{z_i},\qquad
 \rho(\mathbb R)=m_0=1/J,\qquad
 \mu=\rho+1_{[-3,3]}dx,\quad W=6+1/J.
\tag{2}
\]
背景不是免费来源，全部 W 包括其质量 6。阈值 τ=β=1，实际响应及完整水平集为
\[
 M(x)=\max_{1\le L\le2}(h_L*\mu)(x),\qquad E=\{M>1\}.
\tag{3}
\]
在 0≤s≤1/100 上，
\[
 z'(s)=\frac{9/10-30s-15s^2}{(1+s)^2}\ge29/50.
\tag{4}
\]
理由：分子≥1197/2000，而 (29/50)(101/100)^2<1197/2000。故位置严格递增，z_0=0，z_J=3/404<1/100。不能使用下界 .59：端点导数约为 .5867；.58 已充分。

## 2. 真实连续最大值的线包络

对 x∈[3/5,9/10]，所有位置在 x 左侧，各到达尺度
\[
 L_i(x)=2(x-z_i)\in(1,2).
\]
增大 L 时捕获集合依次为后缀 A_i={i,…,J}，质量正好 m_i。所有 query 均在背景 D=[-3,3] 内，故在固定捕获集合的尺度间隙中响应为 1+m_i/L，严格递减。L=1 未捕获 spike，响应为1；L=2 在所有 spike 捕获之后且严格劣于 L_0。因而整个连续最大值只需比较这些到达点，不能用任意选择器代替：
\[
 M(x)=1+\max_i\frac{m_i}{2(x-z_i)}.
\tag{5}
\]
signal 的倒数准确等于
\[
 \frac{2(x-z_i)}{m_i}
 =2J\{x+(x-9/10)s_i+15s_i^2\}.
\tag{6}
\]
这是关于 s_i 的严格凸二次式，最小值在最接近 s_*(x)=(9/10-x)/30 的节点取得。每个内部 i=1,…,J-1 的严格唯一赢家区间
\[
 I_i=\bigl(9/10-30(s_i+h/2),\
             9/10-30(s_i-h/2)\bigr)
\tag{7}
\]
长度为 v_i^0=30h=3/(10J)。区间端点仅相邻两个候选等值，其他候选严格较差；端点为零 Lebesgue 测度。没有凭人工 tie 规则制造正体积赢家。

完整 E 也须保留。一个 spike 可在某 L≤2 中被捕获，当且仅当 x∈∪_i[z_i-1,z_i+1]=[-1,1+z_J]。在此区间所有 query 的背景完整；有 spike 被捕获就使 M>1。区间外没有 spike，背景平均≤1。因此原子模型精确有
\[
 E=[-1,1+3/404],\qquad |E|=2+3/404<201/100.
\tag{8}
\]
局部 I_i 并非任意指定响应，而是这个完整 E 中由同一输入实际连续赢家选出的不交部分。由于 L≥1、spike 总质量1/J，整个 E 还满足 1<M≤1+1/J<2；若只要求 hard 幅度带，可取 λ=1/2。此处仍不认证原 soft/FIRST 资格。

## 3. log proxy 的下界与精确系数的区别

令 E_{A_i} 为完整 E 中实际赢家完整捕获后缀 A_i、未捕获其他 spike 的 receiver 集。可以只使用 I_i 作为它的已核子集；其他捕获类型保持在完整输入中。于是 v_{A_i}=|E_{A_i}|≥3/(10J)。选取实际有限 capture-block 中这些内部后缀，记
\[
 \ell_A=\log(1+v_A/m_A),\quad
 b_A=\sqrt{m_A}\ell_A,\quad
 \Gamma_{AB}=\frac{m_{A\cap B}}{\sqrt{m_Am_B}}.
\tag{9}
\]
因为 v_{A_i}/m_i≥(3/10)(1+s_i)≥3/10，
\[
 \ell_{A_i}\ge\log(13/10)\ge3/13.
\tag{10}
\]
最后一步自含：log(1+t)=∫_0^t(1+u)^{-1}du≥t/(1+t)，取 t=3/10。没有以浮点 log 认证端点。

所有内部后缀共有最后一个原子 z_J，其完整质量 α_J=100/(101J)。因此非负 proxy 的任意其他 block 只会增大它，且
\[
 Q_{\log}:=\sum_{A,B}b_Ab_B\Gamma_{AB}
 =\sum_{A,B}m_{A\cap B}\ell_A\ell_B
 \ge\frac{900(J-1)^2}{17069J}.
\tag{11}
\]
右侧线性趋于无穷，而完整 W=6+1/J、τ|E|<2.01 均有界。这排除仅凭真实 capture-overlap Gram、block 质量与体积定义 (9) 后要求其统一 n-only 上界的合同。不是忽略 log 非线性后产生的 fake block 复制；这些 block 为严格唯一连续赢家的真实不交输出，且每个都有体积≥.3/J。

精确 occupation 则使用
\[
 s_A=\int_{E_A}\frac{dx}{R(x)+m_A},\qquad
 c_A=\sqrt{m_A}s_A.
\tag{12}
\]
在完整捕获 block 上，每个来源 i∈A 的系数都是同一个 s_A；它不依赖 i，不能保留无必要的 i 下标后误估。精确 spike 平方是 Σ_{A,B}c_Ac_BΓ_{AB}。log 重排给 s_A≤ℓ_A 仍成立，但 (11) 说明这一次放大不能整体保留统一 Gram 费用。

尤其每个 source y 的真实 profile
\[
 S(y)=\int_E\frac{1_{Q_{R(x)}(x)}(y)}{M(x)R(x)}dx
\]
满足 M>1，故
\[
 S(y)\le\int_{\mathbb R}\sup_{1\le L\le2}h_L(x-y)dx
 =1+\log2.
\tag{13}
\]
后一个积分在 |x-y|≤1/2 为1，在 1/2<|x-y|≤1 为1/(2|x-y|)，其他处为0。Tonelli 还给 ∫S dμ=|E|，因此 ∫S²dμ≤(1+log2)|E|。它针对完整来源，包括背景。

另有实际下界 R≥1，直接使 s_A≤v_A/(1+m_A)。在已核局部 block I_i 上精确系数≤.3/J，而 proxy 的 ℓ_A 为常数阶。这清楚定位了失败：体积重排放大丢掉了这个真实尺度门，并非真实 S 自身线性增长。

## 4. 精确 L1 输入：完整 packet 的严格赢家

无需仅靠原子极限。取
\[
 \varepsilon=\frac1{10^6J^2},\qquad
 \rho_\varepsilon(y)=\sum_{i=0}^J
       \frac{\alpha_i}{2\varepsilon}
       1_{[z_i-\varepsilon,z_i+\varepsilon]}(y),\qquad
 f_\varepsilon=1_D+\rho_\varepsilon.
\tag{14}
\]
这是非负 L1 函数，质量仍为6+1/J。由 (4)，相邻中心间距≥(29/50)h>2ε，packet 不重叠。对 x∈[3/5-ε,9/10-ε]，所有 packet 在 x 左侧，且各捕获带
\[
 [\,2(x-z_i-\varepsilon),\,2(x-z_i+\varepsilon)\,]
\tag{15}
\]
互不相交并在 (1,2) 内。背景响应为1。

对 i<J，
\[
 \alpha_i=\frac{h}{J(1+s_i)(1+s_{i+1})}
 \ge\frac1{103J^2},
\tag{16}
\]
末 packet 的 α_J 更大。当 L 在 packet i 的部分捕获带内，捕获 spike 质量 m(L) 的导数为 α_i/(4ε)，且 m(L)≤1/J。因此
\[
 \frac{d}{dL}\frac{m(L)}L
 =\frac{L\alpha_i/(4\varepsilon)-m(L)}{L^2}>0,
\quad
 \frac{\alpha_i}{4\varepsilon}
 \ge\frac{10^6}{412}>\frac1J.
\tag{17}
\]
间隙中的 m(L)>0 固定，响应严格递减；首次捕获之前响应为1。故所有可能的最大值均在完整捕获端
\[
 L_i^\varepsilon(x)=2(x-z_i+\varepsilon).
\tag{18}
\]
这些端点完整捕获 packet i,…,J，不捕获 i 之前任何 packet。倒数式 (6) 仅把 x 替换为 x+ε。所以 I_i^ε=I_i-ε 仍是严格唯一的整个连续最大值区间，长度仍为3/(10J)。所有 (9)–(11) 中来源集合可替换为这些真正固定 packet，质量不变，重叠公共来源为最后一个完整 packet，proxy 仍线性发散。

全域严格水平集精确为
\[
 E_\varepsilon=(-1-\varepsilon,\ 1+3/404+\varepsilon),
 \quad |E_\varepsilon|=2+3/404+2\varepsilon.
\tag{19}
\]
极端 receiver 仅触 packet 的端点，捕获质量为0，故不在严格水平集；内部捕获正质量。区间内所有尺度 query 都在 D，区间外没有 spike 且背景平均≤1。因此完整背景与全 receiver 域无缺失。任意全局最大值选择都必须在 I_i^ε 取 (18)，没有 tie 选择依赖。完整 L1 profile 仍满足 (13)，∫S²fε≤(1+log2)|Eε|。

## 5. 范围与缺口

本构造否定的是 (9) 定义的 volume-only log proxy 的统一费用，而不是 (12) 的 exact-coefficient Gram。单 block 的 log 上界合法；把所有合法上界同时放大后，公共来源使其平方失控。原来源数量可任意大，不能限制 J 或将最后一大 packet 丢弃。

本稿主结论严格在连续 [1,2] 成立。依 receiver 变化的到达目录不是共同有限尺度集；本次未给有限网格右侧逼近/捕获标签稳定证明，不能据此宣称原任意 finite family 的同结论。完整原 input/真实 centered h_L/winner/阈值/L1 已认证，原 actual FIRST 与 CP/GP/LCA/history 未认证。无原 actual geom 反例、无一般弱界反例、无主账新增费用。

一般输入的 dominance 重排可用阈值本身给 kernel 第二下界，不一定需要真实背景；该独立推广由另一作者负责，本稿不重做。

## 6. 新守卫预登记

三个确定性精确轮次 J=16,64,256；无随机种子。所有有理输入和残差使用 Fraction。检查：来源正性/质量 telescope、.58 导数残差、位置排序、整个连续极大值线包络的内部点及相邻 tie 边界、真实区间长度、全 E 的背景支持、后缀 overlap、(11) 公共来源下界，以及 L1 packet 间距/部分捕获严格升/完整捕获端的 x+ε 同一包络。

这里没有 receiver 时间网格或启发式优化；连续极大值由 (5)–(7) 和 (17)–(18) 解析减少到有限精确候选。脚本用 log1.3 的解析有理下界3/13，不以浮点 log 作证。数值守卫核有限实例的解析组件，不代替所有 J 证明，更不认证 actual FIRST。

## 7. 终态执行收据

新脚本只执行一次，1.308 秒完成；无 live handle、无随机采样，无旧数值重跑。共 702,607 项 Fraction 检查 PASS（含预登记一致性1项）。其中多数为真实候选的逐对线包络及后缀来源重叠组件检查，数量不是不同定理数。

| J | 本轮检查数 | 公共末 packet 的 proxy 有理下界 |
|---|---:|---:|
| 16 | 2,746 | 50625/68276 ≈ .74148 |
| 64 | 41,674 | 893025/273104 ≈ 3.26991 |
| 256 | 658,186 | 14630625/1092416 ≈ 13.39291 |

完整 E 的原子体积三个轮次均为811/404；L1 体积再加2ε，完整来源质量为6+1/J。上述小数仅便于阅读，下界认证是有理数与 (10) 的解析积分不等式。结果 JSON 同时保存全部内部后缀 overlap 的有理 proxy 下界摘要及分数 SHA，不输出巨大分数。

- 脚本：[true_winner_log_gram_counter_20261007_exact_guard.py](true_winner_log_gram_counter_20261007_exact_guard.py)，SHA-256：4c07385c7c0294d6f2ba79bebf3334d484cf4b46236ac762ec1e77bcee4f99e2。
- 预登记：[true_winner_log_gram_counter_20261007_registration.json](true_winner_log_gram_counter_20261007_registration.json)，SHA-256：9eacb9220f4bccb55616952f74817bf3c40b9c3354e023e20c190922662debf9。
- 终态：[true_winner_log_gram_counter_20261007_results.json](true_winner_log_gram_counter_20261007_results.json)，SHA-256：b81f36379b3fdb67bcd649f0741e819b116d902536e6b616602c8a624af52ad7。

此收据是原 hard 连续赢家构造的精确有限组件与 L1 完整 packet 资格，不是抽象 incidence toy；但仍不能补全原 actual FIRST/history 合同。先前 single_capture_source_rearrangement 稿中 log 单 block 上界和正 Gram 放大恒等式保持正确；需要撤回的是把这个放大 proxy 作为一般可支付统一合同的候选，精确 c_A 仍必须保留。

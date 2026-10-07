# Gamma 占用到原接收点赢家：持续窗口、加权峰及剩余费用

2026-10-07。本稿将已证 Gamma log-\(r\) 占用接到原 \(K_r\nu_b\) 的**完整参数全局赢家**，严格剔除 fixed \(H_r\nu_b\) 已付高值支，并保留同一个 \(\sigma=Su=\nu_b-\mu_b\)、\(\Omega^c\) 与来源质量。得到小赢家分支的显式 \(n^{1/3}\) 体积费，以及更强的维数无关加权峰预算；overshoot 子支因此可付常数费用。没有证明整个 polylog 弱端点，没有把 \(D_r\sigma\) 当正多项式，不报最优阶数，不转入 actual geom 主账。

## 1. 固定输入与精确分支

沿用 [原域外差接口](ordered_frozen_exterior_difference_20261007.md)、[Gamma score 稿](resolvent_gamma_score_20261007.md)：
\[
S=\sum_i(I-G_{c,i}),\quad
K_r=\prod_i[(1-r)I+rG_{c,i}],\quad H_r=e^{-rS},\quad D_r=K_r-H_r,
\]
\[
u\ge0,\quad \Omega=\{u>0\},\quad
\sigma=Su=\nu_b-\mu_b,\quad
\nu_b=\nu\mathbf1_\Omega,\quad
0\le\mu_b\le\kappa,\quad \int\nu_b=\int\mu_b=W_b.
\tag{1}
\]
每段固定同一个 \(c>0\)、物理尺度与输入；常数 uniform \(c\)，不是跨 \(c\) 或跨尺度换源的 maximal。原 obstacle 的 \(u\in L^1\cap L^2\)、\(\|\sigma\|_1\le2W_b\)、\(|\Omega|\le W_b/\kappa\) 保留，势能 cap 不使用。

定义真实非负响应与全参数最大：
\[
F_x(r)=K_r\nu_b(x),\quad
M(x)=\max_{0\le r\le1}F_x(r),\quad
p(x)=\min\operatorname*{argmax}_{0\le r\le1}F_x(r).
\tag{2}
\]
在除去有限个 mask 系数的共同零集后，\(F_x\) 是有限非负 Bernstein polynomial、全 \(r\) 连续且有限。不重判来源，不把截断最大点当完整赢家。

令
\[
\mathcal G=\{x\in\Omega^c:\sup_{t\ge0}H_t\nu_b(x)\le\kappa\}.
\tag{3}
\]
它剔除整个 fixed-semigroup 已付支：
\[
\kappa|\Omega^c\setminus\mathcal G|
\le A_{\rm fix}(n)W_b,\qquad A_{\rm fix}(n)=O(\log(n+2)).
\tag{4}
\]
本稿引用既有 fixed endpoint，不重跑其数值。另由正 contraction，
\[
-\kappa\le D_r\mu_b\le\kappa\quad\text{全空间、全 }r.
\tag{5}
\]
因此在 \(\mathcal G\)，
\[
D_r\sigma=F_x(r)-H_r\nu_b(x)-D_r\mu_b(x)
\ge F_x(r)-2\kappa.
\tag{6}
\]
这是从正多项式回到原 signed 差的合法回代。没有假设 \(D_r\sigma\ge0\) 全参数、log-concave 或 Bernstein 系数非负。

对 \(0<b\le1/2\)，定义
\[
E_b=\{x\in\mathcal G:M(x)>4\kappa,\ p(x)\le b\},\qquad
G_b=\{x\in\mathcal G:M(x)>8\kappa,\ p(x)\le b\}.
\tag{7}
\]
所有结论保持 full-global winner 条件。若只在 \([0,b]\) 取最大而最大恰在 \(b\)，导数不必为0，下文窗口论证不能照用。

可测性初等：\(M\) 是 countable-rational supremum；最小 compact argmax 存在，且 \(\{p\le a\}=\{\max_{[0,a]}F=M\}\)，右侧 Borel。按整数紧时间区间的 bounded-generator exp 正 majorant，\(H_t\nu_b\) 可取共同满测集的全 \(t\) 连续版本，(3)用有理时间检测即可。无需不可测 winner 选择。

## 2. 后验驻点与左侧持续窗口

固定 \(x\in E_b\)。写
\[
F_x(r)=\sum_{k=0}^n a_k(x)r^k(1-r)^{n-k},\qquad a_k\ge0.
\]
在 \(\Omega^c\)，\(a_0=\nu_b(x)=0\)。\(M>0\) 且 \(p\le b<1\)，所以 \(p\in(0,1)\) 是完整全局最大处的内点驻点。令
\[
\pi_k=\frac{a_kp^k(1-p)^{n-k}}M.
\]
求导给 \(E_\pi k=np\)，而 \(\pi_0=0\) 给
\[
np\ge1.
\tag{8}
\]
这里后验只是固定响应的代数权；不假设来源坐标独立。由 [已审正 Bernstein 恒等式](positive_bernstein_rational_mesh_20261007.md) 或直接 Jensen，
\[
F_x(r)\ge M e^{-nD_{\rm KL}(p\Vert r)},\qquad
D_{\rm KL}(p\Vert r)\le\frac{(p-r)^2}{r(1-r)}.
\tag{9}
\]
只用驻点与非负系数，不要求 \(F_x\) log-concave，允许多峰及真 ties。

取
\[
\epsilon=\frac1{4\sqrt{np}}\le\frac14,\qquad
I_x=[(1-\epsilon)p,p]\subset(0,b].
\tag{10}
\]
窗口仅向左，故即使 \(p=b\) 也完整落在 Gamma 占用带中，不需扩大 \(b\)。对 \(r\in I_x\)，
\[
nD_{\rm KL}(p\Vert r)
\le\frac{np\epsilon^2}{(1-\epsilon)(1-p)}
\le\frac16.
\]
于是
\[
F_x(r)\ge e^{-1/6}M\ge\frac56M.
\tag{11}
\]
由(6)，在 \(E_b\) 窗口内 \(D_r\sigma>4\kappa/3>\kappa\)；在 \(G_b\) 内更有
\[
D_r\sigma\ge\frac56M-2\kappa>\frac12M.
\tag{12}
\]
端点等号无害；严格水平集由 \(M>4\kappa\) 或 \(M>8\kappa\) 给出严格剩余裕量。

两个重要宽度是
\[
\int_{I_x}\frac{dr}{r}=-\log(1-\epsilon)\ge
\frac1{4\sqrt{np}},
\]
\[
\int_{I_x}\sqrt{nr}\,\frac{dr}{r}
=2\sqrt{np}[1-\sqrt{1-\epsilon}]
\ge\sqrt{np}\epsilon=\frac14.
\tag{13}
\]
上穿正性只在此已证窗口使用；并未将全 \(D\sigma\) 换成正 polynomial。

## 3. 原 signed 来源的两种 Gamma 占用

已证真实单参数取消，对 \(r\le1/2\)，
\[
\|D_r\|_{\rm TV}\le\sqrt{(1+\tfrac43r^3)^n-1}
\le\sqrt{\tfrac43n}\,r^{3/2}e^{(2/3)nr^3}.
\tag{14}
\]
保留同一个 \(\sigma\)、\(\|\sigma\|_1\le2W_b\)。定义
\[
L_b=\int_0^b\int_{\Omega^c}[D_r\sigma]_+\,dx\,\frac{dr}{r},
\quad
T_b=\int_0^b\sqrt{nr}\int_{\Omega^c}[D_r\sigma]_+\,dx\,\frac{dr}{r}.
\]
正 Tonelli、(14)给
\[
\boxed{L_b\le\frac8{3\sqrt3}\sqrt n\,b^{3/2}
e^{(2/3)nb^3}W_b,\qquad
T_b\le\frac2{\sqrt3}nb^2e^{(2/3)nb^3}W_b.}
\tag{15}
\]
这是时空 signed 来源一次占用，不按输出或窗口再领取 \(W_b\)。

当 \(n\ge8\)、\(b=n^{-1/3}\le1/2\) 时，
\[
L_b<\frac{16}{5}W_b,\qquad
T_b\le\frac2{\sqrt3}e^{2/3}n^{1/3}W_b.
\tag{16}
\]
用 \(e^{2/3}<2,\sqrt3>5/3\) 可给完全有理的安全常数。

## 4. 真实加权峰预算与 overshoot 已付子支

由(12)(13)，每个 \(x\in G_b\) 满足
\[
\int_0^b[D_r\sigma(x)]_+\frac{dr}{r}
\ge\frac{M(x)}{8\sqrt{np(x)}}.
\]
因此
\[
\boxed{\int_{G_b}\frac{M(x)}{\sqrt{np(x)}}dx
\le8L_b
\le\frac{64}{3\sqrt3}\sqrt n\,b^{3/2}
e^{(2/3)nb^3}W_b.}
\tag{17}
\]
canonical \(b=n^{-1/3}\) 下
\[
\int_{G_b}\frac{M}{\sqrt{np}}dx<\frac{128}{5}W_b.
\tag{18}
\]
这比仅控制体积保留更多实际输出信息：对任何预先固定 \(a>0\)，
\[
G_{b,a}=\{x\in G_b:M(x)\ge a\kappa\sqrt{np(x)}\}
\]
有
\[
\boxed{\kappa|G_{b,a}|<\frac{128}{5a}W_b
\quad(b=n^{-1/3}).}
\tag{19}
\]
故真正 overshoot 子支可以常数/\(a\) 支付；若 \(a^{-1}\) 取 polylog，费用也是 polylog。保留 \(M>8\kappa\)、完整赢家与 H-good 条件，不能删掉它们后扩成所有 signed-source 最大值。剩余 moderate-height 子支可能满足 \(8\kappa<M<a\kappa\sqrt{np}\)，不由(19)付清。

## 5. 丢权重得到的小赢家体积费及其范围

由(11)(13)，\(x\in E_b\) 的带权占用至少 \(\kappa/4\)，故
\[
\boxed{\kappa|E_b|\le4T_b
\le\frac8{\sqrt3}nb^2e^{(2/3)nb^3}W_b.}
\tag{20}
\]
canonical 带有
\[
\kappa|E_b|<10n^{1/3}W_b.
\tag{21}
\]
这是明确、可回代的弱分支费，不是整个 polylog 端点，也不是新的最优阶数。\(n^{1/3}\) 来自舍弃(17)的峰权重后取最坏 \(np\)；不能据此推断真实最优维数阶。

若 \(b=O(n^{-1/2})\)，(20)为常数级，但旧窄窗冻结已有支付，本稿不重报它作新的全主账成果。较大赢家 \(p>b\) 没有进入 Gamma 小带，仍未付；H-good 条件外由(4)已付，不混同两种费用。

原 \(D_r\sigma\) 若在小参数处有大值，但 \(K_r\nu_b\) 的完整赢家在带外，亦不能强套(20)。此输出只说明 \(M\) 大，不能凭小带中的 signed 峰自行指定一个驻点。使用截断 winner 或忽略端点导数会失去(9)。

## 6. 直接 \(D_r\nu_b\) 的解析简化，来源质量减半

Root/gated 在守卫运行后提出更直接的接口，本节只解析附注，不修改已运行 guard。保留同一 \(M,p,\mathcal G\)。在
\[
G_b^{\rm dir}=\{x\in\mathcal G:M>4\kappa,\ p\le b\},
\]
窗口中
\[
D_r\nu_b=F_x(r)-H_r\nu_b\ge\frac56M-\kappa>\frac12M.
\]
来源质量是 \(W_b\)，不需要 \(\sigma\) 的 \(2W_b\) 因子。因此(17)右侧减半：
\[
\boxed{\int_{G_b^{\rm dir}}\frac{M}{\sqrt{np}}dx
\le\frac{32}{3\sqrt3}\sqrt n\,b^{3/2}
e^{(2/3)nb^3}W_b
<\frac{64}{5}W_b\quad(b=n^{-1/3}).}
\tag{22}
\]
对应 overshoot 子支 \(\kappa\) 体积费 \(<64W_b/(5a)\)。也有直接体积费
\[
\kappa|G_b^{\rm dir}|
\le\frac2{\sqrt3}nb^2e^{(2/3)nb^3}W_b
<\frac52n^{1/3}W_b
\quad(b=n^{-1/3}),
\tag{23}
\]
因为窗口中 \(D_r\nu_b>2\kappa\)，带权占用至少 \(\kappa/2\)。

完整输入拆为 \(\nu=\nu_b+\nu_g\)，\(\nu_g\le\kappa\)。例如固定 \(\kappa=\lambda/12\) 后，\(\sup K_r\nu>\lambda\) 在域外蕴含 \(M>11\kappa\)；\(\Omega\) 先付 \(12W_b\)，H-bad 支付 \(12A_{\rm fix}W_b\)。余下小赢家 overshoot 支可按(22)乘12回代；大赢家或 moderate overshoot 仍待证。这里使用同一 bad 来源的完整赢家作外覆盖，不冒充原 cube fullfuture/FIRST/birth/交通门；不同原约束的实际回代仍须另审。

## 7. 新三轮纯有理守卫与限制

先写 [注册](gamma_occupancy_winner_bridge_registration_20261007.json)，首次运行新 [guard](gamma_occupancy_winner_bridge_guard_20261007.py)，终态 exit0。[结果](gamma_occupancy_winner_bridge_results_20261007.json) 保存52829条 Fraction/有理区间断言，三轮 \(n=8,32,128\) 为1966/8018/42845。各用129/257/513个窗口节点、Exp 项数64/128/256；\(b=1/\lceil n^{1/3}\rceil\)，全部 \(p=k/n\)、\(1\le k\le\lfloor nb\rfloor\)。

标量响应为单-mask polynomial（确有全局峰 \(p\)），以及 \(k\ge2\) 的 \(k-1,k,k+1\) 后验权 \(1/4,1/2,1/4\)（只核驻点，**不宣称该混合的全局峰为 \(p\)**）。63份保存 profile 核验后验质量/均值、零常数项、窗口 KL/chi 上界、实际有理节点响应持续值，以及原 H/cap 扣除阈值。三轮扩大节点，所有决定性运算精确有理；节点断言只核实现，连续窗口由(9)–(11)解析负责，不以节点一致代替连续证明。

Guard 为避免非平方根用更窄有理窗口 \(\epsilon=1/(4\lceil\sqrt k\rceil)\)，比解析窗口至多损失2倍宽度。它核验的是同一持续/扣除机制，解析主常数(17)–(23)由本文的精确平方根窗口证明负责；不称更窄窗口数值独立认证全部空间预算。\(e^{2/3}<2\) 由正 Taylor 和/几何尾有理区间负责，\(\sqrt3>5/3\) 平方比较精确。

这些 scalar 系数是普适 Bernstein 引理的条件守卫，不声称任意系数可由原 \(G_c\) 来源实现，不是 actual input 反例。原 \(\sigma,\Omega^c,\nu_b,\mu_b\) 合同和 full-global winner 在解析回代中完整保留；没有数值构造新的 obstacle、实际 fullfuture/CP/GP/FIRST 样本。单项 Beta 积分只核有限 factorial identity，不据它宣称原源 sharpness、弱下界或最优 \(n^{1/3}\)。

注册 SHA256：66b00b79d8a99e660b3e3af2f66cd2ef5ccaee10e0d26f6e3cfa67a00d8d7372。
终态脚本 SHA256：13800a865478387bd8784b80174459bf53fbb0078ffbbe507d6b2404488752a1。
旧6425/6030、root mesh/face/原 profiles 均未重跑。§6是运行后的解析 corollary，收据不宣称它另有新 numerical pressure。

**准确未付接口**是同一来源下的大完整赢家，或小赢家中未达到 overshoot 比例的 moderate-height 分支。占用与赢家可以通过已证持续窗口接通，但不能无代价抹去 \(1/\sqrt{np}\) 权重，亦不能把本固定 ordered 工具直接加入原 actual geom 主账。

# High-m packet 桥：任意细化的量词障碍与同源稀疏化近恒等极限

2026-10-07。仅负责本文件。未执行新数值，未修改主 TXT 或其他工作者文件。使用 A02、A03、M04 已读 SKILL/method/cube-interface 的误差、极限和可测参数核验流程；以下结果均自含证明，不引用未核外部课题定理。

结论范围：本稿严格排除“任意足够细的空间分割，或对这种分割免费随机平移，必然产生固定强度 posterior 检测／固定时间稀疏化简化”这一量词。它不排除在 \(d\asymp a/n\) 选择合适空间尺度，不排除检测信号与曲率同时趋零但阶数不同的策略，也不排除真实几何停止规则。尚未得到 \(\sqrt n\,n^{o(1)}\) 普适支付。

## 1. 已有合同与查重

完整有限正 Borel 来源 \(\mu\)，\(W>0\)，原连续全边长窗口 \([a,b]\)，\(0<a<b\)，闭 cube \(Q(x,r)\)。令 \(U_r=\mu(Q(x,r))/r^n\)，\(M=\max_{[a,b]}U_r\)，\(E=\{M>\tau\}\)，\(I=\tau|E|=mW>0\)，原可测真赢家 \(R(x)\)。接收端为 Lebesgue，来源坐标不独立。写 \(Z=1+n\log(b/a)\)，则 \(\int M\le ZW\)。

已读最新 `source_cell_thinning_budget_20261007.md` 和 `nearmax_conditional_covariance_20261007.md`：对一个固定、所有接收点共同使用的半开格子分割，边长 \(d\le\min\{a/(8n),(b-a)/4\}\)，原 posterior 碰撞 \(\chi_R=\sum_C\pi_R(C)^2\) 满足
\[
 \tau\int_E\chi_R\le6W,\qquad
 \tau\int_E\chi_R\log(M/\tau)<\tfrac52W.
\tag{1}
\]
对合法候选 \(L(x)\in[a,b]\)，若 \(g=\log(M/U_L)\in[0,\eta]\)，则 \(\tau\int_E\chi_L\le6e^{2\eta}W\)。same-global-checker 的 variance-sensitive hinge 合同已付 \(W\) 级曲率。固定格子、独立指数寿命的同源删除／补偿增长给
\[
 \mu_d^{(u)}=e^u\sum_C\mathbf1_{\{T_C>u\}}\mu|_C,
 \quad \mathbb EW_d^{(u)}=W,\quad
 \mathbb E\Phi_\tau(\mu_d^{(u)})\ge\Phi_\tau(\mu)-16uW.
\tag{2}
\]
本稿不重证或重新记为新成果。

旧 `nearflat_cube_geometry_20261007.md` 已提醒非原子分割的 \(\ell^2\) posterior 差可趋零；`general_packet_share_band_square_20261007.md` 已提醒质量深度随微细格子扩大。`log_max_variation_independent_review_20261007.md` §5 与 `nearmax_smoothing_obstacle_20261007.md` 已证明保持质量的平滑近优归约。本文新增的是：对所有平移和所有 bounded-gap 真实候选一致的定量退化，以及同一来源稀疏化的全连续最大值 \(L^1\) 近恒等极限及其显式有限 mesh 误差。停止费用由另一工作者负责，不在这里重复推导。

## 2. 非原子来源的最大小块质量一致趋零

固定 \(n\) 和一个有限、非原子 \(\mu\)。定义
\[
 q_d=\sup_{y\in\mathbb R^n}\mu(y+[-d,d]^n).
\tag{3}
\]
则 \(q_d\to0\) 当 \(d\downarrow0\)。这是对任意格子平移一致的量，不仅是某个嵌套网格中的最大格质量。

证明：若存在 \(d_j\downarrow0\)、\(y_j\) 与 \(c>0\)，使这些 cube 质量至少 \(c\)，取紧集 \(K\) 使 \(\mu(K^c)<c/2\)。相应 cube 均与 \(K\) 相交，故 \(y_j\) 有界。提子列 \(y_j\to y\)。任意 \(r>0\) 时，充分大 \(j\) 的 cube 包含于 \(y+[-r,r]^n\)，所以后者质量至少 \(c\)。有限测度的从上连续性给 \(\mu(\{y\})\ge c\)，矛盾。

任意 \(\ell^\infty\) 直径至多 \(d\) 的非空格 \(C\) 都包含于以其中一点为中心的边长 \(2d\) cube，因此 \(w_C=\mu(C)\le q_d\)。此处不假设来源密度有界，亦适用于非原子奇异来源。

## 3. 对全部真实候选、网格平移的一致小碰撞与小差

取上述任意可数 Borel 分割。令 \(m_R(x)=\mu(Q(x,R))>\tau a^n\)。对 \(x\in E\)，
\[
 \chi_R(x)=\frac{\sum_C\mu(C\cap Q_R)^2}{m_R(x)^2}
 \le\frac{q_d}{m_R(x)}\le\frac{q_d}{\tau a^n}.
\tag{4}
\]
对任意合法 \(L(x)\) 满足 \(0\le g\le\eta\)，有
\(m_L=e^{-g}ML^n>e^{-\eta}\tau a^n\)，故
\[
 \chi_L(x)\le\frac{e^\eta q_d}{\tau a^n},\qquad
 B_d(x)^2:=\sum_C[\pi_R(C)-\pi_L(C)]^2
 \le\frac{2(1+e^\eta)q_d}{\tau a^n}.
\tag{5}
\]
右边不依赖 \(x\)、候选的选择或格子平移。真实赢家与候选可在所有原尺度中选择；无需将它们替换成自由 posterior。

若每格赋同一个全局独立符号 \(\xi_C\)，\(v_\xi|_C=\xi_C\)，则
\[
 \mathbb E_\xi|\pi_Rv_\xi-\pi_Lv_\xi|
 \le B_d(x)\le\left[\frac{2(1+e^\eta)q_d}{\tau a^n}\right]^{1/2}.
\tag{6}
\]
对所有 \(0<t\le1\)，其 averaged hinge 也有上界
\[
 \frac\tau W\int_E\mathbb E_\xi
 [t|\pi_Rv_\xi-\pi_Lv_\xi|-g]_+dx
 \le tm\left[\frac{2(1+e^\eta)q_d}{\tau a^n}\right]^{1/2}\longrightarrow0.
\tag{7}
\]
同样，使用已证 deterministic \(B_d/\sqrt3\) 的 hinge，左边至多 (7) 的右端。若 \(g\ge g_0>0\)，当 \(tB_d/\sqrt3<g_0\) 时该 deterministic hinge 恰为零。

**精确 no-go 量词。** 固定任何实际非原子输入及 \(m,\tau,a,b,\eta\)，不存在一个只依赖这些参数的正数 \(F>0\)，使 (7) 的左端对每个足够细的分割均至少 \(F\)。high \(m\)、来源 nearflat 和径向窄峰若原先成立，在这一步没有改变输入，仍成立；故不能把这个障碍归因于换成低 \(m\) 示例。

已知保常数平滑归约使全局任意精度 nearmax 可取为光滑密度，因此此量词障碍实际进入 nearmax 类。这里仅引用已有归约，未另报一个新的平滑定理。

**不能过度解释。** (4)--(7) 没有否定方差敏感合同。真实平方曲率亦可能趋零得更快，或候选 gap 随 \(d\) 更快趋零；此时小绝对信号仍可能大于曲率。特别是光滑非退化内部赢家的 checker 位移路线不被 (7) 排除。也没有证明在具体 \(d\asymp a/n\) 的网格上 \(B_d\) 一定小；(3) 没有关于 \(n\) 的统一收敛速率。

## 4. 随机平移与全都微细的跨尺度混合不能修复这一量词

随机平移不改变 (3)--(7)，所以先取任意平移，再对平移求平均，也不能产生统一的正检测下限。不能用一个接收点专属平移作为同一个全局来源测试。

更一般，先抽一个测试层 \(j\)，概率为 \(\lambda_j\ge0\)、\(\sum_j\lambda_j=1\)，每层有自己的固定共同分割和独立 checker。若所有层的格直径均至多 \(d_*\)，则对整个测试混合仍有 (7)，把 \(q_d\) 换成 \(q_{d_*}\)。当整体最粗尺度 \(d_*\downarrow0\) 时，有限或可数的这种混合均退化。

若以共同测试 \(v=\sum_j c_jv_j\) 合并层，\(\sum_j|c_j|\le1\) 足以保持 \(|v|\le1\)。三角不等式再次给同一个 \(q_{d_*}^{1/2}\) 上界。无界地重新归一化微细 checker 则必须同时缩小允许的扰动幅度，以保持正来源 \((1\pm tv)\mu\)；原对称变分不提供免费的放大。

因此跨尺度策略若要躲开本 no-go，需要保留某个非消失的物理尺度、利用精细信号／曲率比，或证明其他真正共用测试。这里只说明必要的量词改变；不宣称这些方式已足够，更不免费给每层领取一次 \(W\)。

## 5. 固定时间、同源微细稀疏化的真实连续 maximal 近恒等极限

**命题。** 固定任意有限非原子 \(\mu\)、\(a,b,n\)、\(u<\infty\)，取任意平移、边长 \(d\downarrow0\) 的半开网格，并使用同一格子的 survival 硬币。记 \(p=e^{-u}\)，\(\mu_d=p^{-1}\sum_C\eta_C\mu|_C\)，\(\eta_C\) 独立 Bernoulli\((p)\)。则
\[
 \mathbb E\|M_{[a,b]}\mu_d-M_{[a,b]}\mu\|_{L^1(dx)}\longrightarrow0,
 \qquad
 \mathbb E|\Phi_\tau(\mu_d)-\Phi_\tau(\mu)|\longrightarrow0,
\tag{8}
\]
并且
\[
 \mathbb E|W_d-W|^2=(p^{-1}-1)\sum_Cw_C^2
 \le(e^u-1)q_dW\longrightarrow0.
\tag{9}
\]
这些是期望下的真实全连续尺度陈述，不是有限 \(J\) 的 exact-winner 冒充连续赢家。

**第一步：固定查询的同源方差。** 对每个原固定 \(s\in[a,b]\)，令 \(V_{s,d}=U_s(\mu_d)-U_s(\mu)\)。独立性作用在完整空间格硬币上，因而逐点
\[
 \mathbb EV_{s,d}=0,\quad
 \mathbb E|V_{s,d}(x)|^2
 =\frac{p^{-1}-1}{s^{2n}}\sum_C\mu(C\cap Q(x,s))^2
 \le\frac{(p^{-1}-1)q_d}{s^{2n}}\mu(Q(x,s)).
\tag{10}
\]
可数硬币的求和以有限截断收敛；\(\sum w_C=W\) 保证定义和二阶界有效。Tonelli 得
\[
 \int\mathbb E|V_{s,d}|^2dx\le(p^{-1}-1)q_dW/s^n.
\tag{11}
\]
在 \(A_H=[-H,H]^n\)、\(H>b/2\) 上 Cauchy--Schwarz 给
\[
 \mathbb E\int_{A_H}|V_{s,d}|dx
 \le\left[(2H)^n(e^u-1)q_dW/a^n\right]^{1/2}.
\tag{12}
\]
对 \(x\notin A_H\) 的 cube 查询，其来源点至少有一个坐标在 \([-H+b/2,H-b/2]\) 外。核质量为 1，\(\mathbb E\mu_d=\mu\)，所以
\[
 \mathbb E\int_{A_H^c}|V_{s,d}|dx
 \le2\mu\bigl(\mathbb R^n\setminus[-H+b/2,H-b/2]^n\bigr)=2T_H.
\tag{13}
\]
这一步明确支付来源尾部；不由总质量有限直接断言全空间 \(L^1\) 的 Cauchy--Schwarz。

**第二步：有误差的有限原尺度 mesh。** 取 \(0<\xi\le1\)，令 \(\Delta=a\xi/(2n)\)、\(J=\lceil(b-a)/\Delta\rceil\)，
\[
 s_j=\min\{b,a+j\Delta\},\quad j=0,\ldots,J,
 \qquad M_J(\lambda)=\max_{0\le j\le J}U_{s_j}(\lambda).
\tag{14}
\]
相邻尺度满足 \(s_{j+1}/s_j\le1+\xi/(2n)\)，而
\(n\log(1+\xi/(2n))\le\xi/2\le\log(1+\xi)\)。任意 \(r\in[a,b]\) 向上取最近 \(s_j\)，由 cube 包含关系，对每个有限正来源 \(\lambda\) 准确有
\[
 M_J(\lambda)\le M_{[a,b]}(\lambda)\le(1+\xi)M_J(\lambda).
\tag{15}
\]
这里原连续窗口并未改变，\(b\) 是准确端点，不允许向窗外 padding；(15) 是显式误差包络，绝非有限 \(J\) 与连续 winner 相同。

令 \(A=M(\mu_d), B=M(\mu), a_J=M_J(\mu_d),b_J=M_J(\mu)\)。两个方向分别用 (15)，给
\[
 |A-B|\le(1+\xi)|a_J-b_J|+\xi(A+B),\qquad
 |a_J-b_J|\le\sum_{j=0}^J|V_{s_j,d}|.
\tag{16}
\]
所有随机来源的窗口强包络仍给 \(\mathbb E\int A\le ZW\)。结合 (12)--(16)，得到可登记、对实际输入有效的全连续误差界
\[
 \boxed{\quad
 \mathbb E\|M(\mu_d)-M(\mu)\|_1
 \le(1+\xi)(J+1)
 \left[\sqrt{(2H)^n(e^u-1)q_dW/a^n}+2T_H\right]
 +2\xi ZW.
 \quad}
\tag{17}
\]
固定 \(\xi,H\) 后令 \(d\downarrow0\)，再令 \(H\to\infty\)，最后 \(\xi\downarrow0\)，即得 (8) 第一式。目标顺序清楚；没有关于 \(n\) 的免费一致速率。\(t\mapsto\tau\log_+(t/\tau)\) 的 1-Lipschitz 性给 (8) 第二式。(9) 独立地由质量方差得出。

**两个合法加强。** 第一，(8)(9) 给 \(\Phi_\tau(\mu_d)/W_d\to\Phi_\tau(\mu)/W\) 于概率中，零质量实现的概率趋零。因此若原输入的分式得分至少 \(K-\varepsilon\)，则对每个固定 \(\delta>0\)，以趋于 1 的概率新输入的得分仍至少 \(K-\varepsilon-\delta\)。已有全局 nearflat 合同适用于这些真正新输入及其真正新赢家；在其条件范围内，\(|m_d-K|\le3\sqrt{(\varepsilon+\delta)Z}\)。故这是实际 nearmax 类上的近恒等结论，不仅是低值来源示例。

第二，全部证明只用 \((e^u-1)q_d\)。因此对确定的时间安排 \(u_d\ge0\)，只要 \(e^{u_d}q_d\to0\)，(8)(9) 同样成立，允许 \(u_d\to\infty\)。它仍不允许来源依赖的随机停止免费代入。这里来源、维数和原窗口固定，未给跨维数的一致速率。

## 6. 数量与幅度：固定时间不强制进入少格类

若 \(\mu\) 紧支撑，网格中非空格数 \(N_d<\infty\) 且 \(N_d\ge W/q_d\to\infty\)。存活非空格数为真正的 Binomial\((N_d,p)\)，故每个固定 \(K<\infty\) 满足
\[
 \mathbb P\{N_d^{\rm alive}\le K\}\longrightarrow0.
\tag{18}
\]
例如用其方差 \(N_dp(1-p)\) 与均值 \(N_dp\to\infty\) 的 Chebyshev 界即可。非紧支撑时若非空格可数无穷，固定 \(p>0\) 后仍有无穷多个存活格 a.s.。所有实现同时有
\[
 \sup_C\mu_d(C)\le e^u q_d\longrightarrow0.
\tag{19}
\]
原来源非原子，则 \(\mu_d\) 也非原子。因而固定 \(u\) 下先把格子细化，既没有生成少量物理原子，也没有生成少量存活小格；其原尺度响应和 \(\Phi\) 还按 (8) 逼近原输入。

对上述确定时间 \(u_d\) 的加强，若 \(e^{u_d}q_d\to0\)，亦有 \(N_dp_d\ge W/(e^{u_d}q_d)\to\infty\)，而 Binomial 方差至多其均值，所以 (18) 仍成立；(19) 的右边仍趋零。

这是对“固定时间并对任意微细网格自动简化”的严格障碍，不是对全部稀疏化策略的否定。几何上选 \(d\asymp a/n\)，让 \(u\) 随实际几何停止，或以另外的可支付类作为目标，仍需独立证明。尤其 (2) 只直接给确定时间 \(u\) 的期望合同；任意来源依赖的停止时间不能免费代入 \(16uW\)。

## 7. 径向窄峰提供的合法接口与仍未检测的量

固定混合 \(H=\ell^{-1}\int_a^bU_s ds/s\)，\(\ell=\log(b/a)\)，已有 \(\mathbb E_X H/M\le1/m\)。令真实 nearwinner 对数尺度宽度
\[
 A_\eta(x)=\int_a^b\mathbf1_{\{U_s(x)\ge e^{-\eta}M(x)\}}\frac{ds}s.
\tag{20}
\]
直接由 \(\mathbf1_{\{U_s/M\ge e^{-\eta}\}}\le e^\eta U_s/M\)，得到
\[
 \mathbb E_X A_\eta\le e^\eta\ell/m,
 \qquad \tau|\{x\in E:A_\eta(x)\ge w\}|
 \le\frac{e^\eta\ell}{w}W\quad(w>0).
\tag{21}
\]
这是 fixed-mixture 门的透明推论，不报为独立强上界。它控制 nearwinner 集的尺度宽度；没有给 (5) 中空间格 posterior 差的下界。

真实 mass 的单调性保证：对 \(L\ge R\)，\(U_L/M\ge(R/L)^n\)，所以 \([R,\min\{b,Re^{\eta/n}\}]\) 是合法 bounded-gap 候选区间。但质量若没有进入外层，则 \(\pi_L=\pi_R\)，仅因体积归一化产生 gap。因此“存在近赢家候选”与“候选具有可检测的新来源”是不同条件。此观察不被用于构造 low-m 反例冒充 high-m 反例。

在真正 high-m nearmax 来源上，(21)、小碰撞预算和来源 nearflat 共同成立；目前尚缺利用适当物理尺度或信号／曲率比，由这些共同条件强迫 (variance-sensitive) switch 增益的普适定理。本稿明确堵住任意细化及免费 shift 的强量词，但没有关闭这一一般空间桥。

## 8. 若登记后处理，最精确的连续尺度 mesh 守卫

不启动新实例。可登记 \(\xi\in(0,1]\)、\(\Delta=a\xi/(2n)\)、\(J=\lceil2n(b-a)/(a\xi)\rceil\)，以及准确 \(s_j=\min(b,a+j\Delta)\)。若 \(a,b,\xi\) 为有理数，所有 \(s_j\) 与守卫
\[
 a=s_0<s_1<\cdots<s_J=b,\qquad
 (s_{j+1}/s_j)^n\le1+\xi
\tag{22}
\]
均可用 Fraction 精确核验。任意其他有限原尺度 mesh 若满足 (22)，同样得到 (15)--(17)，查询数改为它的实际数量。

数值可核 fixed-query 的真实 full-source posterior 方差、\(W_d\) 方差和网格质量 \(q_d\) 的有效上包；全 \(dx\) 积分及来源尾部仍须解析承担。只保存 finite \(J\) 响应时，必须保留 (15) 的 \(1+\xi\) 费用，不能将其称为 exact continuous winner。三轮实例或 saved-posterior 守卫应另行预登记，不由本稿的极限证明假定已经完成。

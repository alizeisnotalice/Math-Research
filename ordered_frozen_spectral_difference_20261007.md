# 原 ordered 核与固定生成元：谱差最大估计、障碍能量反例及 signed 来源接口

2026-10-07。本轮暂停小 contrast 重排，独立核查新 spectral-difference 策略。本文证明一个原 joint-symbol 的 dimension-free L2 工具，严格排除误用的障碍能量上界，并保存可回代的未付接口。没有改进 \(A_{\rm ord}\)，没有付原 cube 时变空间账。旧 Abel/obstacle/原核实验不重跑。

## 1. 原固定入口与比较对象

固定 \(c>0\) 与物理尺度，令
\[
G_i=G_{c,i},\qquad A_i=I-G_i,\qquad S=\sum_i A_i,
\]
\[
K_r=\prod_i(I-rA_i)=T_r^c,\qquad H_r=e^{-rS},
\qquad0\le r\le1.
\tag{1}
\]
各 \(G_i\) 是原轴向正对称保质量 convolution operator，其乘子 \(g_c(\xi_i^2)=[1+cB(\xi_i^2)]^{-1}\)，\(B(v)=v/\log(1+v)-1\)。各 \(A_i\) 交换、自伴，谱位于 \([0,1]\)。以下证明更一般地适用于任意有限个 commuting positive self-adjoint contractions \(G_i\)。

原 joint spectrum 写为
\[
a_i=1-g_c(\xi_i^2)\in[0,1],\quad
s=\sum_i a_i,\quad q=\sum_i a_i^2\le\min(s,s^2),
\]
\[
k_r=\prod_i(1-ra_i),\quad h_r=e^{-rs},
\quad m_r=\sqrt{s}(k_r-h_r).
\tag{2}
\]
这里小写 \(s\) 是谱值，大写 \(S\) 是算子；不是 softness 变量。\(H_r\) 是现有 fixed 工具的同一个 generator-family 成员：若旧冻结生成元 \(L_{\rm fr}=c^{-1}S\)，则 \(H_r=e^{-crL_{\rm fr}}\)。uniform \(c\) 来自证明只使用 \(0\le a_i\le1\)，没有把 \(K_r\) 当作 semigroup。

## 2. 准确谱积分界

对所有 joint spectrum，有
\[
\boxed{\int_0^1\left(|m_r|^2+|r\partial_rm_r|^2\right)\frac{dr}{r}
\le40\min(s,1/s).}
\tag{3}
\]
\(s=0\) 时右边定义为0，左边也为0。

先取 \(0<r\le1/2\)。令
\[
d_r=-\sum_i\log(1-ra_i)-rs\ge0.
\]
由 \(-\log(1-x)-x\le x^2\) 对 \(x\le1/2\) 成立，
\[
d_r\le r^2q,\qquad
r d'_r=\sum_i\frac{r^2a_i^2}{1-ra_i}\le2r^2q.
\tag{4}
\]
写 \(u=rs\)。于是
\[
|m_r|\le\sqrt{s}e^{-u}r^2q,\qquad
|rm'_r|\le\sqrt{s}e^{-u}r^2q(u+2).
\]
把 \(u\) 积分范围扩大到 \((0,\infty)\)，
\[
\begin{split}
I_{\rm small}
&\le\frac{q^2}{s^3}\int_0^\infty u^3e^{-2u}[1+(u+2)^2]du\\
&=\frac{27}{4}\frac{q^2}{s^3}
\le\frac{27}{4}\min(s,1/s).
\end{split}
\tag{5}
\]
积分使用 \(\int_0^\infty u^je^{-2u}du=j!/2^{j+1}\)；三项分别 \(15/8,3,15/8\)。这保留 product 与 exponential 的 cancellation，没有用二者 triangle 代替 small-\(r\) 差。

对 \(1/2\le r\le1\)，\(0\le k_r\le e^{-rs}\)，且不除以可能为0的因子而直接求 polynomial 导数，得
\[
|k'_r|\le\sum_i a_i e^{-r(s-a_i)}
\le e\,s e^{-rs}.
\]
故 \(e<3\) 给
\[
|m_r|\le\sqrt{s}e^{-u},\qquad
|rm'_r|\le4\sqrt{s}\,u e^{-u}.
\tag{6}
\]
若 \(s\le1\)，直接在 \(r\in[1/2,1]\) 上积分给 \(I_{\rm late}\le17s\log2<17s\)。若 \(s\ge1\)，把 \(u\in[s/2,s]\) 扩大到 \([s/2,\infty)\)：
\[
I_{\rm late}
\le s\int_{s/2}^\infty e^{-2u}(1+16u^2)\frac{du}{u}
\le[1+4s(s+1)]e^{-s}
\le33/s.
\tag{7}
\]
最后一步用 \(s^je^{-s}\le j!\) 对 \(j=1,2,3\) 成立。式(5),(7)总系数 \(27/4+33=159/4<40\)；\(s\le1\) 总系数 \(27/4+17<40\)。\(a_i=1,r=1\) 的端点由 polynomial 导数界处理，未误用发散的 log 表达式。

## 3. log-\(r\) 最大 L2：没有额外因子2

对任意 \(g\in L^2\)，由(3)及共同谱分解，
\[
\boxed{
\left\|\sup_{0\le r\le1}|(K_r-H_r)S^{1/2}g|\right\|_2^2
\le40\left\langle g,\min(S,S^{-1})g\right\rangle
\le40\|g\|_2^2.}
\tag{8}
\]
\(\min(S,S^{-1})\) 在零谱处定义为0。

证明：设 \(F_\ell=(K_{e^\ell}-H_{e^\ell})S^{1/2}g\)，\(\ell\le0\)。每个 receiver 上 \(F_\ell\to0\) 当 \(\ell\to-\infty\)，且
\[
\sup_{\ell\le0}|F_\ell|^2
\le2\int_{-\infty}^0|F_\ell||\partial_\ell F_\ell|d\ell
\le\int_{-\infty}^0(|F_\ell|^2+|\partial_\ell F_\ell|^2)d\ell.
\]
空间积分后 Plancherel/共同谱定理给(8)。因此 norm 常数是 \(\sqrt{40}\)，不是 \(\sqrt{80}\)。

all-\(r\) 的版本并非只由 strong L2 连续推断。这里 \(K_r\) 是有限 operator polynomial；\(S=nI-\sum_iG_i\) 有界。\(H_r\) 的 exp 幂级数被正算子 \(\sum_k(nI+\sum_iG_i)^k|S^{1/2}g|/k!\) 支配，此 majorant 在 L2 中有限，故在 \([0,1]\) 几乎处处绝对一致收敛。给定同一满测集，轨道及导数连续，以上 pointwise log-\(r\) 基本微积分合法。

对一个 \(u\in L^2\)，令 \(g=S^{1/2}u\)，(8)还给更准确的
\[
\left\|\sup_r|(K_r-H_r)Su|\right\|_2^2
\le40\|\min(S,1)u\|_2^2
\le40\langle u,Su\rangle.
\tag{9}
\]
这只是谱权改进；下一节说明两种右侧均不能由 cap 与来源质量支付。

## 4. 任意输入障碍没有 \(\langle u,Su\rangle\le\kappa W\)

现有 [fixed jump obstacle 接口](fixed_jump_obstacle_interface_20261007.md) 对原冻结轴核的非负 \(L^1\cap L^2\) 密度构造
\[
u\ge0,\quad \Omega=\{u>0\},\quad
0\le\mu_{\rm tot}=\nu-Su\le\kappa,\quad
\mu_{\rm tot}=\kappa\ \text{于 }\Omega.
\tag{10}
\]
按常数缩放 potential，这正是同一 \(S\)：若旧构造为 \(\nu-L_{\rm fr}U=\mu\)，则取 \(u=U/c\) 即有 \(Su=\nu-\mu\)。不是新存在性猜想。保留真实拆分
\[
\nu_b=\nu\mathbf1_\Omega,\quad\nu_g=\nu\mathbf1_{\Omega^c},\quad
\mu_b=\nu_b-Su,\quad
\mu_{\rm tot}=\nu_g+\mu_b.
\]
有 \(\nu_g\le\kappa\)、\(0\le\mu_b\le\kappa\)、\(\int\mu_b=\int\nu_b=W_b\le W\)，\(|\Omega|\le W_b/\kappa\)。fixed killed operator 是 \(S_\Omega+\kappa/u\)，源是 \(\nu_b\)，不能用含域外 good 部分的总 \(\nu\) 写 killed 方程。unbounded \(\kappa/u\) 已由旧 graph-domain 正迭代处理，不由本稿重新假设为有界。

现在用同一个原 \(S\) 给严格能量障碍。取任意固定 \(n,c,\kappa,W\)，令 \(\nu\) 在体积 \(W/H\) 的有限盒上等于 \(H>\kappa\)，别处0。它是非负 \(L^1\cap L^2\) 密度，不是原子、有限链或 reset proxy。盒内必属于 \(\Omega\)：在 \(u=0\) 处 \(Su=-\sum_iG_i u\le0\)，故 \(\mu_{\rm tot}\ge\nu\)，与 \(H>\kappa\) 矛盾。于是盒内
\[
Su=H-\kappa.
\]
因 \(0\le S\le nI\)，
\[
\|Su\|_2^2\le n\langle u,Su\rangle,\qquad
\min(S,1)\ge S/n.
\]
所以
\[
\boxed{\langle u,Su\rangle
\ge\frac{(H-\kappa)^2}{nH}W,\qquad
\|\min(S,1)u\|_2^2
\ge\frac{(H-\kappa)^2}{n^2H}W.}
\tag{11}
\]
固定 \(n\) 后令 \(H/\kappa\to\infty\)，两者相对 \(\kappa W\) 均无界。没有用错误的 \(\int_\Omega u(\nu-\kappa)\) 正项估计去忽略其余区域的负项；(11)使用正确非负谱能量。不能据(8)–(9)宣称 ordered 弱端点已闭合。

旧 Abel 漏泄估计是 source/exit-layer 的 \(\sqrt N\kappa W\) square-function，而不是(11)的 potential energy 上界；两者不能混称。

## 5. 可回代的域外 joint-symbol 缺口

令 \(\sigma=Su=\nu_b-\mu_b\)。完整来源仍是 \(\nu=\mu_{\rm tot}+\sigma\)，且
\[
\sigma|_{\Omega^c}=-\mu_b,\qquad \|\sigma\|_1\le2W_b.
\tag{12}
\]
对 cap 部分，positive contraction 给 \(|(K_r-H_r)\mu_{\rm tot}|\le\kappa\)。fixed \(H\) 最大弱工具已付 \(A_{\rm fix}(n)=O(\log(n+2))\)；没有理由把(8)的全空间 potential 能量当其额外 source coupon。

一个明确但**未证**的域外 joint-difference 接口是
\[
2\kappa\left|\left\{x\in\Omega^c:
\sup_{0\le r\le1}[(K_r-H_r)\sigma(x)]_+>2\kappa\right\}\right|
\le B_nW_b.
\tag{13}
\]
它只要求原 obstacle 来源、原同一个 \(u,\Omega,\nu_b,\mu_b\)，不要求一般 signed 输入。阈值 \(2\kappa\) 是可调整的常数，故不是强 L2 合同；需要 \(B_n=\mathrm{polylog}(n)\) 才能接目标。

例如选 \(\kappa=\lambda/6\)。\(\Omega\) 先付 \(\lambda|\Omega|\le6W_b\)。域外
\[
K_r\nu=H_r\nu+(K_r-H_r)\mu_{\rm tot}+(K_r-H_r)\sigma.
\]
若 \(H\) max≤\(\lambda/2\)，正差 max≤\(2\kappa\)，则总响应≤\(\lambda\)。故(13)若成立，严格水平集有合法覆盖，总费用
\[
\lambda|\{\sup_rK_r\nu>\lambda\}|
\le2A_{\rm fix}(n)W+(6+3B_n)W_b.
\tag{14}
\]
这只是充分接口及常数回代，不是已付结果。

为什么旧 free−killed Abel theorem 不能直接给(13)：\(K_r\) 和差乘子依赖各个 \(a_i\)，不是 \(S=\sum a_i\) 的单变量函数；源 \(u\) 的 killed operator 是 \(S_\Omega+\kappa/u\)，其中压缩坐标不再交换。旧 \(\rho_j=T s^j(s+K)^{-j-1}\nu_b\) 的正序/source-once/pointwise \(\mu_b\) 支配，只覆盖 fixed resolvent 的对应展开。把 joint multiplier 直接放在 killed operator 或 exit measure 上需要新 intertwining/泄漏预算，不能由 free 端 Fourier 界推得。

真实有序 killed 过程亦可在固定 \(\Omega\) 定义，Duhamel 的正 exit source 总质量≤\(W_b\)。但再对其未来传播取 receiver max，立即返回未证 \(A_{\rm ord}\)；这是现有循环，不是(13)的证明。尤其原过程的 \(G_{c(1-r),i}\) 随时间变化，不能把它们替换成同一个 fixed \(G_c\) 后仍称为原 exit flux。

## 6. 不依赖 potential 能量的 signed \(L^1\) 来源表示

另一个准确入口是
\[
h=(I-e^{-S})u=\int_0^1e^{-tS}\sigma\,dt,
\qquad \boxed{\|h\|_1\le2W_b}.
\tag{15}
\]
Bochner/Fubini 合法，因为 \(\sigma\in L^1\)、\(H_t\) 是 L1 contraction；旧 obstacle 的 \(u\in L^1\) 还给 \(\int h=0\)。在 \(\Omega^c\)，\(u=0\) 故 \(h=-e^{-S}u\le0\)。**\(h\) 是 signed 新来源，不是非负输入，也没有 cap-L2 上界。**

定义 joint multiplier/operator
\[
\mathcal M_r=(K_r-H_r)\frac{S}{1-e^{-S}},
\tag{16}
\]
其中 ratio 在零谱处连续取1、整体 \(\mathcal M_r=0\) 于零谱。则准确
\[
(K_r-H_r)\sigma=\mathcal M_rh.
\tag{17}
\]
来源 \(h\) 对整段参数是同一份，式(15)只支付一次 \(2W_b\)。

还有一个可复核 L2 后果：
\[
\int_0^1\left(|\mathcal M_r|^2+|r\partial_r\mathcal M_r|^2\right)\frac{dr}{r}
\le40\frac{\min(s^2,1)}{(1-e^{-s})^2}\le160.
\]
因为 \(1-e^{-s}\ge s/2\) 对 \(s\le1\)，且≥\(1/2\) 对 \(s\ge1\)。故
\[
\|\sup_r|\mathcal M_rh|\|_2\le\sqrt{160}\|h\|_2.
\tag{18}
\]
这不提供 \(L^1\to L^{1,\infty}\)。若未来能证明该差核对一般 signed \(h\) 的 weak constant \(B_n\)，由(15)、\(\kappa=\lambda/4\) 及阈值 \(\lambda/4\) 可给 \(A_{\rm ord}\le2A_{\rm fix}+8B_n\)；目前这只是条件回代。更弱地保留 \(h\) 的 obstacle 来源及域外符号，可能足以给(13)，不应先强求全空间一般 signed 定理。

**本稿真正未付项**是 joint 差 kernel 的 weak 源预算，或(13)的原域外 exit-layer 版本。没有把 \(L^1\) 来源表示、dimension-free L2 或固定 semigroup 的 weak 常数彼此混换。

## 7. 小型三轮原 symbol 守卫收据

先写 [注册](ordered_frozen_spectral_difference_registration_20261007.json)，再运行新 [脚本](ordered_frozen_spectral_difference_guard_20261007.py) 一次，[结果](ordered_frozen_spectral_difference_results_20261007.json) 保存全部45组 joint symbol。三轮 \(n=8,32,128\)，\(c=1,1/2,1/1000\)，各有 zero、统一低频 \(10^{-3}\)、统一1、统一64、周期频率 \(0,1/8,1,8,64\) 五种预设。不改模型，不拟合阶数。

终态 exit0：55项 Fraction 常数/真实 spike 必要能量下界，加225项原 symbol 浮点 screen，共280项。在 log-\(r\in[-40,0]\) 上16个 panel 的 positive GL16/32/64，遗漏的小 \(r\) 端有式(4)的解析上界
\[
\frac{s q^2r_0^4}{4}[1+(r_0s+2)^2],\qquad r_0=e^{-40}.
\]
最大嵌套求积差 \(1.135\times10^{-11}\)，直接 product/导数与 cancellation 表示的最大残差 \(6.453\times10^{-14}\)。浮点 quadrature、低频 Taylor 与 symbol 计算不是 interval；它们只审查实现，dimension-free 定理由§2–3解析负责。没有解新障碍矩阵或把必要下界当存在性数值。

注册 SHA256 \(\texttt{84ca5bc878c63ab5942f244513a63c2afe82c2b86e9b1ae4359488ac1a2caafb}\)，脚本 SHA256 \(\texttt{55b79520130b75bf7cb59a2509550af609db3f24e9eb27e41fdfa32e7def473e}\)。旧三轮原核/Abel/face 收据未重跑；无 live 进程、无 main cube 费用。

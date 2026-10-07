# 原 resolvent Gamma score：共同潜变量取消与截断差核预算

2026-10-07。核查 [root 候选](resolvent_score_next_interface_20261007.md)，接 [原域外差接口](ordered_frozen_exterior_difference_20261007.md)。两个真实原核 score 上界成立；共同 Gamma 时间比较给 \(r^3\) 单轴 \(\chi^2\) 取消；截断差核 maximal 与 log-参数占用各有明确来源一次预算。较大参数的域外正来源弱费仍未证。本文不将单参数 TV 或占用当 maximal，不重复 root 原 face-TV 下界或旧6030项。

采用 [F04](/Users/zhengzhihao/.codex/skills/math-f04-poisson-gamma-abel-square-function/SKILL.md) 的参数/迁移审计规范，已读其 provenance、method、cube-interface；技能不是本课题端点定理。下述推导直接针对原 resolvent 与既有 obstacle。

## 1. 真实核及来源合同

固定 \(c>0\) 和物理尺度。令
\[
G_i=(I+cB_i)^{-1},\quad A_i=I-G_i,\quad S=\sum_iA_i,\quad
K_r=\prod_i(I-rA_i),\quad H_r=e^{-rS}.
\tag{1}
\]
\(B_i=B(-\partial_{ii})\)，\(B(v)=v/\log(1+v)-1\)。原 \(B\) 是 Bernstein function，可由 \(B(v)=\int_0^1[(1+v)^t-1]dt\) 核对；原 \(\Pi_t^c=e^{-ctB_i}\) 正、对称、保质量。因此
\[
G_i=\int_0^\infty e^{-t}\Pi_t^c\,dt,\qquad
G_i^j=\int_0^\infty\frac{t^{j-1}e^{-t}}{\Gamma(j)}
\Pi_t^c\,dt\quad(j\ge1).
\tag{2}
\]
Gamma 取 shape \(j\)、rate 1，不是 reset 替代模型；\(j=0\) 是时间0的 identity。送至真实位移的 Markov map 可含奇异部分，RN/条件 Jensen 都以测度表述。所有上界 uniform \(c>0\) 与物理尺度，因为 data processing 不增加潜变量上界。每个结论仍固定同一个 \(c\) 与来源；uniform 常数不等于 \(\sup_c\) 或不同尺度共同换源的 maximal 结论。

保留已构造的原 fixed obstacle：
\[
\sigma=Su=\nu_b-\mu_b,\quad 0\le\mu_b\le\kappa,\quad
\int\nu_b=\int\mu_b=W_b,\quad
\Omega=\{u>0\},\quad u\in L^1\cap L^2,\quad u\ge0.
\tag{3}
\]
\(\nu_b\) 支撑 \(\Omega\)。共同源与目标差算子为
\[
R=(I+S)^{-1},\quad h=R\sigma=h_+-h_-,\quad
N_r=(K_r-H_r)(I+S).
\tag{4}
\]
旧稿已证 \(h_+\subset\Omega\)、\(h_-\le\kappa\)、各质量至多 \(W_b\)，整个负源费用 \(320W_b\)，域外 holding 项作用于 \(h_+\) 为0。不按参数或 receiver 重选来源。

## 2. 原 \(SK_r\) 的 Gamma score

单轴 \(p_r=(1-r)\delta_0+rG\)。基准核用 holding 标签概率 \(1-r\)、Exp(1) 时间标签概率 \(r\) 表示。目标 \(Gp_r=(1-r)G+rG^2\) 相对这份带标记基准的密度为
\[
\ell_p=
\begin{cases}0,&\text{holding},\\
(1-r)/r+T,&\text{Exp 标签， }T\sim\mathrm{Exp}(1).
\end{cases}
\tag{5}
\]
有 \(E\ell_p=1\)、\(\operatorname{Var}\ell_p=1/r+r-1\)。把标签/时间送到原位移，条件 Jensen 也给 \(\chi^2(G^2\Vert G)\le E(T-1)^2=1\)。坐标乘积的内部 jump 标签独立，\(SK_r\) 是 \(n-\sum_i\ell_{p,i}\) 的真实 Markov pushforward。故
\[
\boxed{\|SK_r\|_{\rm TV}
\le\min(2n,\sqrt{n(1/r+r-1)}),\quad0<r\le1.}
\tag{6}
\]
TV 是 signed measure 总变差范数，无概率距离的 \(1/2\) 因子。独立性属于生成元内部噪声，不要求任意来源坐标独立，不声称实际后验独立。

更精确的可计算上界是
\[
\|SK_r\|_{\rm TV}\le
\sum_{k=0}^n\binom nk r^k(1-r)^{n-k}
E|n-k(1-r)/r-T_k|,\qquad T_k\sim\Gamma(k,1).
\tag{7}
\]
\(T_0=0\)。令 \(b=n-k(1-r)/r\)。当 \(k>0,b\le0\)，期望为 \(k-b\)；当 \(b>0\)，期望为
\[
k-b+2[bF_k(b)-kF_{k+1}(b)],\quad
F_k(b)=1-e^{-b}\sum_{j=0}^{k-1}b^j/j!.
\tag{8}
\]
这先合并了完整 Gamma score 的正负部分；仍是空间 TV 的上界，不声称相等。\(r=0\) 不用含 \(1/r\) 的公式，粗界 \(2n\) 有效。

## 3. 原 \(SH_r\) 的 Poisson score

原 \(H_r\) 的坐标 jump 数 \(J_i\sim\mathrm{Pois}(r)\)，每跳时间 Exp(1)，总数 \(J\sim\mathrm{Pois}(nr)\)。插入 \(G_i\) 的 size-bias 密度是 \(J_i/r\)，所以
\[
SH_r=\operatorname{push}[(n-J/r)\mathbb P],
\]
\[
\boxed{\|SH_r\|_{\rm TV}\le E|n-J/r|
=2n\,\mathbb P\{\mathrm{Pois}(nr)=\lfloor nr\rfloor\}
\le\min(2n,\sqrt{n/r}).}
\tag{9}
\]
等式由 \(E(J-m)_+=m\,\mathbb P\{\mathrm{Pois}(m)=\lfloor m\rfloor\}\) 得到；整数 \(m\) 的等号事件自身贡献0，不需更改端点。空间位移再合并可以减少 TV。

先将 Poisson 数合并到总 Gamma 时间，还给一轴条件 score 方差
\[
\frac1r-\frac{E\operatorname{Var}(J_i\mid T_i)}{r^2}\le\frac1r.
\tag{10}
\]
不能把重复 jump 标签的绝对值先取掉；但本文不从(10)宣称新的维数阶数，也不把空间投影 score 范数当完整核上界。

## 4. 共同 Gamma 时间的 \(r^3\) 取消

单轴时间空间上的 Bernoulli law 为 \(P_r=(1-r)\delta_0+r e^{-t}dt\)，compound-Poisson law 为
\[
Q_r=e^{-r}\delta_0+r e^{-r-t}F(rt)dt,\qquad
F(z)=\sum_{j=0}^\infty\frac{z^j}{(j+1)!j!}.
\tag{11}
\]
这是全部 Poisson/Gamma 标签先合并为同一真实时间，再送到原位移。准确有
\[
\chi^2(P_r\Vert Q_r)
=e^r\left[(1-r)^2+
r\int_0^\infty\frac{e^{-t}}{F(rt)}dt\right]-1.
\tag{12}
\]
Root 在首次 guard 执行前提供了更简短的常数链；独立核验：\(F(z)\ge1+z/2\)，所以
\[
1/F(z)\le1/(1+z/2)\le1-z/2+z^2/4,
\]
最后一步等价于 \((1-z/2+z^2/4)(1+z/2)-1=z^3/8\ge0\)。积分至多 \(1-r/2+r^2/2\)，方括号至多
\[
1-r+r^2/2+r^3/2\le e^{-r}+\tfrac23r^3,
\]
因为 \(e^{-r}\ge1-r+r^2/2-r^3/6\)。因此
\[
\boxed{\chi^2(P_r\Vert Q_r)\le\tfrac23 e^r r^3
\le
\begin{cases}\frac43r^3,&r\le1/2,\\2r^3,&r\le1.
\end{cases}}
\tag{13}
\]
\(r=0\) 连续取0，\(e^{1/2}<2,e<3\) 足够，常数不称最优。

潜变量 product 满足 \(1+\chi^2(P^{\otimes n}\Vert Q^{\otimes n})=(1+\chi^2(P\Vert Q))^n\)。Markov data processing 与 Cauchy 给真实 \(D_r=K_r-H_r\)
\[
\boxed{\|D_r\|_{\rm TV}\le
\min\left(2,\sqrt{(1+C(r)r^3)^n-1},\,2nr^2\right),}
\tag{14}
\]
其中 \(C=4/3\) 对 \(r\le1/2\)，否则2。最后一项来自原 Bernoulli(r)/Poisson(r) 标签 TV \(2r(1-e^{-r})\le2r^2\) 与 product coupling。共同 Gamma 上界保留了更多 repeated-jump 取消；\(nr^3\) 小时为 \(O(\sqrt n\,r^{3/2})\)。这是单参数非渐近上界，不是 maximal 阶数拟合。

## 5. 完整差 score 的共同源表示

在(11)中，positive time 的
\(m_r(t)=E[J_i\mid T_i=t]=1+rtF'(rt)/F(rt)\)，holding 取0。Gamma 系列另给 \(E[J_i(J_i-1)\mid T_i=t]=rt\)。定义 holding 时 \(\ell_{p,i}=\ell_{q,i}=0\)，否则
\[
\ell_{p,i}=(1-r)/r+t_i,\quad \ell_{q,i}=m_r(t_i)/r,
\quad U_p=n-\sum_i\ell_{p,i},\quad U_q=n-\sum_i\ell_{q,i}.
\]
共同 \(Q_r^{\otimes n}\) 基准下 likelihood 为
\[
w_r(t)=\frac{dP_r^{\otimes n}}{dQ_r^{\otimes n}}(t)
=\frac{e^{nr}(1-r)^{n-k}}{\prod_{i:t_i>0}F(rt_i)},
\quad k=|\{i:t_i>0\}|.
\tag{15}
\]
于是准确
\[
\boxed{N_r=\operatorname{push}
\left([w_r(1+U_p)-(1+U_q)]Q_r^{\otimes n}\right).}
\tag{16}
\]
这在一次 pushforward 前保留完整跨标签取消，不先 triangle 两个 score。原位移条件核 \(\bigotimes_i\Pi_{t_i}^c\) 作用于同一任意相关 \(h_+\)。

(13)控制 \(w_r-1\)，不自动支付(16)的 score 权。Root 的 [独立真实 face-TV 结果](original_face_tv_obstruction_20261007.md) 沿 \(n=m^3,r=1/m\) 给 \(\|N_r\|_{\rm TV}\gtrsim n^{2/3}\)。这与 \(nr^3=1\) 时(13)完全相容，排除以全 TV 给 polylog 费用的捷径，不排除原域外弱合同。本文只读该证明，不重复其数值。

固定 \(r\) 时，(14)可付 \(\kappa|\{[D_rh_+]_+>\kappa\}|\le\|D_r\|_{\rm TV}W_b\)；改为 receiver winner 后不能免费平均该费用。

## 6. 截断连续 maximal：来源只收一次

此为严格差核截断预算，不是新完整 \(A_{\rm ord}\) 上界。写单轴 \(p_i=I-rA_i\)、\(h_i=e^{-rA_i}\)、\(\delta_i=p_i-h_i\)。正 Markov contraction 给
\[
\|\delta_i'\|_{\rm TV}=\|A_i(I-h_i)\|_{\rm TV}\le4r,\qquad
\|\delta_i\|_{\rm TV}\le2r^2.
\]
保留取消的 product telescoping 是
\[
D_r=\sum_i\left(\prod_{j<i}p_j\right)\delta_i
\left(\prod_{j>i}h_j\right).
\]
其他 factors 范数1、导数范数至多2，故 \(\|D_r'\|_{\rm TV}\le4nr+4n(n-1)r^2\)。用 \(\|I+S\|_{\rm TV}\le2n+1,N_0=0\)，对 \(0\le b\le1\) 与同一 \(g\in L^1\)，
\[
\boxed{\|\sup_{0\le r\le b}|N_rg|\|_1
\le C(n,b)\|g\|_1,\quad
C(n,b)=(2n+1)[2nb^2+\tfrac43n(n-1)b^3].}
\tag{17}
\]
证明用 receiver 上 \(N_rg=\int_0^rN_t'g\,dt\)，再 Tonelli，不是 conditional weak 平均。Bounded-generator exp 正 majorant 保证同一满测集上的全参数连续版本。若 \(b=L/n\le1\)，
\[
C(n,L/n)\le6L^2+4L^3.
\]
因此原正源支有
\[
\kappa|\{x\in\Omega^c:\sup_{r\le L/n}[N_rh_+]_+>\kappa\}|
\le(6L^2+4L^3)W_b.
\tag{18}
\]
这仅付该段。若 \(r\ge L/n\) 的同一正部未知费用为 \(C_n^{\rm late}\)，则旧回代
\[
B_n\le2[C(n,L/n)+C_n^{\rm late}]+320.
\tag{19}
\]
没有逐窗另收 \(W_b\)。旧窄窗正核冻结已付某些原 \(K_r\) 分支，故(18)只作差核辅助截断接口，不报为新完整主账成果，也不声称优于旧原核窄窗覆盖。

## 7. 同源 log-参数占用：真实但不是选时端点

Root 在 guard 终态后指出的解析 corollary，可直接核验，无需追加实验。令 \(b=n^{-1/3}\)。在 \(r\le b\)，
\[
(1+2r^3)^n-1\le e^{2nr^3}-1\le2e^2nr^3.
\]
由(14)，对任何 signed \(f\in L^1\)，
\[
\boxed{\int_0^{n^{-1/3}}\|D_rf\|_1\,\frac{dr}{r}
\le\frac{2\sqrt2 e}{3}\|f\|_1.}
\tag{20}
\]
特别取同一原 \(\sigma=\nu_b-\mu_b\)，\(\|\sigma\|_1\le2W_b\)，有
\[
\int_0^{n^{-1/3}}\int_{\Omega^c}[D_r\sigma(x)]_+\,dx\,\frac{dr}{r}
\le\frac{4\sqrt2 e}{3}W_b.
\tag{21}
\]
这是真正 source-once 时空 \(L^1\) 占用费用。Receiver 选时图不由 log-time 平均控制；还缺参数 persistence/宽度或 signed 交叉控制，不能把(21)写成 maximal 弱端点。此 corollary 解析负责，不宣称包含在已冻结的 guard 实验项中。

## 8. 三轮新潜变量精确证书

先写 [注册](resolvent_gamma_score_registration_20261007.json) 再启动新 [guard](resolvent_gamma_score_guard_20261007.py)。Root 的(13)简化在首次执行前进入注册及代码，明确记录 preliminary 常数6改为2/4/3；没有常数6结果或重跑。

首次执行真实 session39898，终态 exit0。[结果](resolvent_gamma_score_results_20261007.json) 保存6425条 Fraction/有理区间断言，各轮 \(n=8,32,128\) 为439/1207/4779。\(r\) 为去重后的 \(1/(4n),1/n,1/\lceil\sqrt n\rceil,1/\lceil n^{1/3}\rceil,1/4,1/2,3/4,1\)，截断 \(L=1,2,4\)；Exp 项数64/160/512。覆盖完整 binomial-Gamma score 均值/方差、(8)(9)绝对值区间与 Cauchy/2n 界、Gamma-Laplace 代数一致性、(13)常数链、(18)费用。理论积分与原 Markov 映射由本文解析负责，潜变量 screen 不冒充实际空间 TV 数值。

Gamma 整数 shape CDF 用精确有限 polynomial 与有理 \(e^{-b}\) 外包；正 Exp Taylor 和作下界，尾上界为首遗漏项除 \(1-b/(M+2)\)，reciprocal 反向给负 Exp 区间。无浮点 quadrature、Monte Carlo 或 rare-event 漏尾。三轮 \(r=1\) 的 latent SK/SH 上界同为约2.2333845、4.5017784、9.0211583，均严格小于根号方差；认证字段保存 Fraction 端点，十进制只供阅读。

注册 SHA256：efa89ae57dea44ff0d86a2d47dde6e5939a3ca9aaa9701a840878a0a8c3bb9f3。
终态脚本 SHA256：6fab306bbbb503db65d0b2b79c0587362e8b2fd65cdf81fd911ec469861aa2a7。
旧6030项、原 profiles、root face 下界均未重跑，无实际 obstacle 样本或 FIRST/nonconc/geom 门认证，不拟合阶数。

**剩余最小费用**：在原 \(h_+\subset\Omega,h_+\le R\nu_b,\mu_b\le\kappa\) 合同下，支付较大 \(r\) 的 \(\Omega^c\) 正差 maximal。共同 Gamma likelihood 把真实取消写成(16)，但尚未给 receiver 选时下的正弱预算；score TV、\(\chi^2\)、截断费与时空占用均不能代替它。本轮不入 cube 主账。

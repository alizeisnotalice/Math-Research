# 任意来源 packets 的整组时钟、真实合并与残余标签熵

2026-10-07。父任务指定 6.1-sol high。只新增本稿及同前缀 LaTeX 片段；不改旧冻结稿或他人文件。解析研究，无数值执行。

**所得与范围。** 任意固定有限来源分割都可合法整组增长／删除，并在已接受历史刷新时合并当前 packets。合并不改变当前来源、原始细标签或质量倾斜密度。预算必须保留有符号
\[
 R_a(\nu,\mathcal B)=\sum_{B\in\mathcal B}a_BD_B(\nu)-L_a(\nu),
\]
不能沿用空间小格的 \(16W_{\rm active}\)。粗标签熵的合并下降准确转移到块内残余熵；原始标签熵仍按真实整组死亡强度耗散。另独核 root 的自身包络提案：\(R_B\le K_BW_B\)，系数为一，但当前 packet 的 \(K_B\) 不免费为常数。本稿给多 packet 终端的可检验充分接口，尚未构造一般输入满足目标费用的控制。

读并使用 [J01](</Users/zhengzhihao/.codex/skills/math-j01-inhomogeneous-jump-generator/SKILL.md>)、[D04](</Users/zhengzhihao/.codex/skills/math-d04-entropy-chain-rule/SKILL.md>)、[E03](</Users/zhengzhihao/.codex/skills/math-e03-tree-bellman-cross-layer-budget/SKILL.md>) 及必要 method/cube-interface/provenance。下面只有有限律、初等积分、真实路径似然和有限望远镜，不借外部 Bellman 或熵定理认证 cube 接口。

查重：[LC](local_clock_source_budget_20261007.md) 已证同一小空间格的加权 \(16W\)；[entropy barrier](local_clock_entropy_barrier_20261007.md) 已证固定细格控制的熵耗散，并明确中途整组删除不属于其原模型。本稿不重复这些小格定理。新增的是合法 coalescing-group 模型、原始标签的准确熵账及有符号组费；几何目标本身仍未解。

## 1. 任意物理来源分割和原完整 posterior

完整有限正 Borel 来源 \(\mu\)，原连续尺度 \([a,b]\)、\(a>0\)、原阈值 \(\tau>0\)。固定有限个不交 Borel 来源集合 \(E_i\)，令
\[
 \mu_i=\mu|_{E_i},\quad w_i=\mu_i(\mathbb R^n)>0,\quad
 \mu=\sum_{i=1}^N\mu_i,\quad W_0=\sum_iw_i.
\]
它们不要求空间直径小。重合物理原子不拆成两个独立标签。当前来源允许
\[
 \nu=\sum_{i\in A}c_i\mu_i,\quad c_i>0,
\]
其当前 groups \(\mathcal B\) 是存活原始细标签 \(A\) 的分割。每组
\[
 \nu_B=\sum_{i\in B}c_i\mu_i,\quad W_B=\nu_B(\mathbb R^n),\quad W=\sum_BW_B.
\]
原始来源 \(\mu\) 的谱系固定；合并只改变 bookkeeping，不替换当前 \(\nu\)。增长和死亡本身当然是完整输入扰动。

写原完整 \(M_\nu=\max_{a\le r\le b}r^{-n}\nu(Q(x,r))\)、\(E=\{M_\nu>\tau\}\)、\(I=\tau|E|\)、\(\Phi(\nu)=\tau\int\log_+(M_\nu/\tau)dx\)。取当前共同可测真赢家 \(R(x)\)，令 \(m=M_\nu R^n\)、\(z=\log(M_\nu/\tau)\)，
\[
 \pi_B(x)=\nu_B(Q(x,R(x)))/m(x),\qquad\sum_B\pi_B=1\quad(x\in E).
\]
所有分母保持当前完整来源。共同可测赢家可取最小最大化尺度：合法子窗最大值是从右侧逼近的有理尺度与端点的可数上确界，最小赢家的子水平集等于子窗最大值与完整 \(M\) 相等。选择服务于本功能积分，不继承另一个原历史选择器。

## 2. 正帽之外保留有符号删除减增长费

对固定于接收点之前的 \(a_B\in[0,1]\)，定义
\[
 D_B(\nu)=\Phi(\nu)-\Phi(\nu-\nu_B)\ge0,\qquad
 L_a(\nu)=\tau\int_E\sum_Ba_B\pi_Bdx.
\]
原赢家用于删除后的完整输入给
\[
 D_B\le\tau\int_E\min\{z,-\log(1-\pi_B)\}dx.
\]
标量 \(\min\{z,-\log(1-p)\}\le p+(1+4z)p^2\) 因而给
\[
 \boxed{R_a:=\sum_Ba_BD_B-L_a\le C_a,\qquad
 C_a=\tau\int_E(1+4z)\sum_Ba_B\pi_B^2dx.}
\tag{1}
\]
\(R_a\) **可以为负**；其 deletion 和完整 denominator 都是真实的，但增长只使用当前冻结赢家下界。不得把 \(R_a\) 等同于所有状态下的精确 \(\Phi\) 生成元。无跳增长为 \(\nu_h=\sum_Be^{a_Bh}\nu_B\)，且
\[
 \log[U_R(\nu_h)/M_\nu]=\log\sum_B\pi_Be^{a_Bh}
 \ge h\sum_Ba_B\pi_B.
\]
局部 TV-Lipschitz 和 \(\Phi\) 的 TV 稳定性使实际漂移 \(dt\)-几乎处处至少为 \(L_a\)，足够使用 \(R_a\)。winner switches 和阈值平台不需被忽略。

任意 packet 没有空间小格碰撞常数。仅有通用粗界，记 \(Z=1+n\log(b/a)\)：
\[
 -ZW\le R_a\le C_a\le I+4\Phi
 \le Z(1+4/e)W.
\tag{2}
\]
全部来源合成一组时 \(\pi=1\)、\(C_a=a(I+4\Phi)\)、\(R_a=a(\Phi-I)\)。若当前 threshold-stationary 输入满足 \(\Phi=I\)，真实有符号组费为零，正帽仍为 \(5aI\)。此为必须保留的差异，不是沿后续过程免费维持 stationarity 的证明。

## 3. 有限合并／死亡控制的合法构造

政策只观察已接受的整组死亡、当前完整来源、现有分组及独立政策种子，不直接观察未来辅助 marks 或永生标签。在有限个确定刷新点及已接受死亡后，可合并若干当前 groups，并给新 groups 选择 \(a_B\in[0,1]\)。相邻刷新间 rates 固定；合并只做 union，从不重新分裂或复活已删来源。

合并瞬间 \(\nu,\Phi,W,c_i\) 全不变。新 packet 是合并当时真实 \(\sum\nu_B\)，不是将成员恢复成初始共同系数。此后共同乘法增长与整组死亡才改变输入。非平凡合并至多 \(N-1\) 次，整组死亡至多 \(N\) 次。

每段条件于刷新历史，各存活组独立取 \(\mathrm{Exp}(a_B)\) 等待时间，零 rate 为无穷；到下一确定刷新点截断。未死亡组用指数记忆无关重启，已接受死亡按左状态 rate，删除后才合并和刷新。有限递归构造保持真正适应性。

条件段首历史，\(\Lambda=\sum_Ba_B\)，首死亡标签 \(B\)、延迟 \(s\) 的密度是 \(a_Be^{-\Lambda s}ds\)，无死亡到 \(s\) 的概率 \(e^{-\Lambda s}\)。故该段实际删除费期望与 \(\int\sum_Ba_BD_B(\nu_s)ds\) 的存活占用期望相同；零 \(\Lambda\) 时均为零。对质量用同一计算，增长 \(\sum_Ba_BW_B\) 与删除质量补偿相消。有限拼接、条件塔式期望及 (2) 的绝对可积性给
\[
 \mathbb E_PW_T=W_0,\qquad
 \boxed{\mathbb E_P\Phi(\nu_T)\ge\Phi(\mu)
             -\mathbb E_P\int_0^TR_a(\nu_s,\mathcal B_s)ds}
\tag{3}
\]
对任意有界 accepted-history 停时 \(T\)。证明可先给确定时刻条件合同，再对 \(\Phi+\int R_a\) 向上有限网格停止；有限 horizon 下 \(W_t\le e^hW_0\)、\(|R_a|\le Z(1+4/e)W_t\)，故有符号积分无无穷抵消。合并本身既不是 \(\Phi\) 跳损失，也不是可以自动记成负费用的“熵奖励”。

## 4. 原细标签的 Q law，而不是合并后重抽标签

初始按 \(w_i/W_0\) 抽原细标签 \(J=i\)。在带标签 \(Q\) 下，保护其当前所在**整组**不死亡，其余 groups 执行同一 accepted-history 政策和原 rate；所有 groups 仍按各自 rate 增长。合并后，含标签的整个大组都被保护，不能继续把其其他原成员当作独立死亡时钟。

给定有限接受死亡路径及原标签 \(i\)，相同路径上政策值一致；非标签组死亡速率因子相消，禁止标签所在组死亡的无死亡指数比值是
\[
 c_i(t)\mathbf1_{\{i\ {\rm alive}\}},\qquad
 c_i(t)=\exp\left(\int_0^ta_{B_s(i)}ds\right)
\]
在存活标签上解释。合并是已观察路径决定的无来源变化转换，不产生新的 likelihood 因子。有限路径密度逐段计算后混合给
\[
 \frac{dQ}{dP}\Big|_{\mathcal F_t}
 =\frac{\sum_iw_ic_i(t)\mathbf1_{\{i\ {\rm alive}\}}}{W_0}
 =\frac{W_t}{W_0}.
\tag{4}
\]
使用 raw 接受历史滤过；不预先加入无限未来的 \(P\)-零测事件，也不纳入拒绝 proposal 和未来 marks 后照搬该公式。

隐藏原细标签的后验和粗组后验为
\[
 p_i(t)=Q(J=i\mid\mathcal F_t)=c_i(t)w_i/W_t,\qquad
 q_B(t)=\sum_{i\in B}p_i(t)=W_B(t)/W_t.
\]
当前组 \(B\) 的真实 \(Q\) 死亡率是
\[
 a_B(1-q_B),\qquad \bar a=\sum_Ba_Bq_B.
\tag{5}
\]
终值质量密度及非负、或绝对可积 signed Fubini 给有界停止的主合同
\[
 \boxed{\frac{\Phi(\mu)}{W_0}
 \le\mathbb E_Q\frac{\Phi(\nu_T)}{W_T}
       +\mathbb E_Q\int_0^T\frac{R_a(\nu_s,\mathcal B_s)}{W_s}ds.}
\tag{6}
\]
还可用正帽 \(C_a/W_s\)，但可能丢掉重要抵销。普通标签活动量只是 \(\mathbb E_Q\int\bar a\)，不等于一般组費 (6)。

## 5. 自身响应的精确 residual 包络：独立审查

root 另提该包络交独审；这里自证关键方向，**不**凭其名字给 \(K_B\) 常数。令
\[
 H_B(x)=M_{\nu_B}(x),\quad
 G(v)=\log(1+v)-v/(1+v),\quad
 K_B=\sup_{\rho>0}\Phi_\rho(\nu_B)/W_B.
\]
在原 \(E\)，\(u=M_\nu/\tau\)、\(v=H_B/\tau\)，有 \(u\ge\max(1,v)\)、\(0\le p=\pi_B\le v/u\)。对
\[
 h(p)=\min(\log u,-\log(1-p))-p,
\]
交点 \(p=1-1/u\) 前递增、后递减。若 \(u\le v+1\)，最优值 \(\log u-1+1/u\) 随 \(u\) 增；若 \(u>v+1\)，最优 \(p=v/u\) 随 \(u\) 减。两支在 \(u=v+1,p=v/(1+v)\) 取得共同最大 \(G(v)\)。因此
\[
 R_B:=D_B-L_B\le\tau\int G(H_B/\tau)dx.
\]
原 \(E\) 可外放，因为 \(G\ge0\)。对每个 \(v\ge0\)，
\[
 G(v)=\int_0^\infty\frac{2}{(1+t)^3}
                  t\log_+(v/t)\,dt.
\]
确证：导数右边为 \(v^{-1}\int_0^v2t(1+t)^{-3}dt
=v/(1+v)^2=G'(v)\)，零端同值。非负 Tonelli 给
\[
 \boxed{R_B\le\int_0^\infty\frac{2}{(1+t)^3}
                    \Phi_{\tau t}(\nu_B)dt\le K_BW_B,\qquad
 R_a\le\sum_Ba_BK_BW_B.}
\tag{7}
\]
权重积分为一，系数一。\(K_B\) 是**当前完整 packet** 的 all-threshold 常数；其形状可含成员不同的过去 \(c_i\)，合并后必须重核。今后同组共同增长保持该 sup 的尺度齐次性，不能因此从 children 的 \(K\) 继承统一常数。唯一全通用值仍是 \(K_B\le Z/e\)。小空间 packet 的几何常数只是已核特殊证书，不代表任意 group。

自身功能随正阈值连续（被积函数由 \(H_B/e\) 支配），故 \(K_B\) 的 sup 可取正有理阈值；当前有限来源系数下它是可测量，不需要未定义的不可数 argmax。零质量组直接忽略。

精确 nearmax 也不免费给 signed 符号：若 \(\Phi(\nu)\ge(K-\varepsilon)W\)，删除后的合法来源比值至多 \(K\)，则 \(D_B\ge KW_B-\varepsilon W\)。而 \(L_a=\int S_\nu a\,d\nu\)，若已有 nearflat 误差 \(\Delta=\int|S_\nu-K|d\nu\)，则
\[
 R_a\ge-\varepsilon W\sum_Ba_B-\Delta.
\]
packet 数可能很多，且不能把原输入 nearmax 当作整条路径每个状态 nearmax。严格 \(\varepsilon=\Delta=0\) 时当前 \(R_a\ge0\)；这也不证明它沿后续控制为零或已付款。

## 6. 合并费用可能双向：真实完整来源诊断

common rate \(a\) 的 children 合并成一组时，\(L_a\) 不变；正帽精确增加
\[
 \Delta C_a=2a\tau\int_E(1+4z)
                    \sum_{C<D\ {\rm merged}}\pi_C\pi_Ddx\ge0.
\]
真实 signed 变化则是 \(a(D_{\rm union}-\sum D_C)\)，没有同样单调性。

真实一维来源可展示两种方向。全边长固定为 \(1\)，两个不同点 packets \(\nu_\pm=v\tau\delta_{\pm\epsilon}\)、\(0<\epsilon<1/2\)。重叠 receiver 区间长 \(1-2\epsilon\)，单来源足迹总长 \(4\epsilon\)，准确
\[
 \Phi(\nu_++\nu_-)
 =\tau[(1-2\epsilon)\log_+(2v)+4\epsilon\log_+v],\qquad
 \Phi(\nu_\pm)=\tau\log_+v.
\]
所以
\[
 \Delta R_a=a\tau(1-2\epsilon)
             [2\log_+v-\log_+(2v)].
\tag{8}
\]
\(v=3/4\) 为负，\(v=3\) 为正；无重合标签或自由 posterior。退化窗口例只作明确计算。正宽连续窗口 \([1,1+\delta]\) 在 \(\delta\downarrow0\) 时响应单调降到固定尺度响应，\(\Phi\) 由有限支撑及确定包络支配收敛，严格符号因而也存在于真正正宽窗口。它不是目标阶数样本或 nearmax 认证。

## 7. 熵合并账：原始标签不消失

保留细标签 \(p_i\)，粗后验 \(q_B=\sum_{i\in B}p_i\)。自然对数 Shannon 熵有逐路径有限链式
\[
 H_{\rm fine}=H(p)=H(q)+K_{\rm in},\qquad
 K_{\rm in}=\sum_Bq_BH((p_i/q_B)_{i\in B}).
\tag{9}
\]
零组质量贡献零；\(Q\) 下标签所在组保证整体非空。

若合并 children \(B_1,\ldots,B_\ell\) 为 \(B\)，source 和 \(p_i\) 不变，粗熵准确下降
\[
 \delta_{\rm merge}=q_BH((q_{B_j}/q_B)_{j=1}^{\ell})\ge0.
\]
块内残余 \(K_{\rm in}\) 同时增加这一量。这里是信息从 coarse label 移到 conditional label，不能记作原 fine entropy 被物理耗散。

在不合并的段内，组共同增长使组内细比例不变。无跳 \(\dot q_B=q_B(a_B-\bar a)\)，删除组 \(B\) 后其余组重归一化，真实强度 (5) 给
\[
 \mathcal L_QH(q)=\sum_Ba_B(1-q_B)\log(1-q_B)
 =-\mathcal D(q,a),\qquad
 0\le\mathcal D\le\bar a.
\tag{10}
\]
这是已有固定格熵计算作用于**当前真实 groups**，没有套原细格的独立删除。块内熵 \(h_B=H(p_i/q_B)\) 在该段常值；\(K_{\rm in}=\sum_Bq_Bh_B\) 的连续漂移 \(\sum_Bq_B(a_B-\bar a)h_B\)，组删除跳漂移 \(\sum_Ba_Bq_B(K_{\rm in}-h_B)\)，两者相消。因此
\[
 \mathcal L_QK_{\rm in}=0,\qquad
 \mathcal L_QH_{\rm fine}=-\mathcal D.
\]
fine entropy 在合并时也不跳，故有限 first-jump 补偿、有界停止给准确账
\[
 \begin{aligned}
 \mathbb E_QH_{\rm fine}(T)&=H_{\rm fine}(0)-\mathbb E_Q\int_0^T\mathcal Dds,\\
 \mathbb E_QH(q_T)&=H(q_0)-\mathbb E_Q\int_0^T\mathcal Dds
                              -\mathbb E_Q\sum_{\rm merges\le T}\delta_{\rm merge},\\
 \mathbb E_QK_{\rm in}(T)&=K_{\rm in}(0)
                              +\mathbb E_Q\sum_{\rm merges\le T}\delta_{\rm merge}.
 \end{aligned}
\tag{11}
\]
首次跳跃密度可由质量倾斜直接写：段初 \(q\)、\(r(s)=\sum_Bq_Be^{a_Bs}\)，
\[
 S_Q(s)=e^{-s\sum_Ba_B}r(s),\quad
 f_{Q,B}(s)=a_Be^{-s\sum a}[r(s)-q_Be^{a_Bs}]
          =S_Q(s)a_B(1-q_B(s)).
\]
它证明所用补偿，死亡后的合并熵损失另作为有限显式项扣除；不援引一般熵率定理。\(q_B=1\) 时强度零，不代入删除后的零分母。\(H\le\log N\)，有限删除和合并允许有界网格停止；无界停止需 \(Q\) 下终止与相应极限条件。

由粗熵非负和 \(\mathcal D\ge0\)，还有 \(\mathbb E_Q\sum\delta_{\rm merge}\le H(q_0)\)。这是平均 residual 信息账，不是逐路径的 \(\Phi\) 费用抵扣。

若终端粗组至多 \(k\)，工作量 \(A=\mathbb E_Q\int\bar a\)，则
\[
 \mathbb E_QK_{\rm in}(T)\ge H_{\rm fine}(0)-A-\log k.
\tag{12}
\]
合并成一个大 packet 后立即冻结可令 coarse entropy 为零、工作量为零，却保留全部原 fine entropy、完整 \(\Phi/W\)。这说明“标签少”本身不是可支付终态。新控制允许这种 large residual entropy 的 terminal，但必须真正估计其完整响应。

## 8. 可验证的多 packet 终端接口

给一个明确、强于只重写 (6) 的充分证书。终端每个原 \(Q(x,b)\) 至多捕获 \(k\) 个正质量 packets；且每个 packet 有实际已证、all-threshold 的自身 cap
\[
 \Phi_\rho(\nu_B)\le\kappa_BW_B\quad(\rho>0).
\]
同一 \(x\) 上任意合法尺度都只有这些至多 \(k\) 个 packet 贡献，所以 \(M_\nu(x)\le k\max_BM_{\nu_B}(x)\)。于是
\[
 \Phi_\tau(\nu)\le\sum_B\Phi_\tau(k\nu_B)
 \le k\sum_B\kappa_BW_B,\qquad
 \Phi_\tau(\nu)/W\le k\sum_Bq_B\kappa_B.
\tag{13}
\]
允许多个共同捕获 packets 和任意多原细标签；不是强制独立集。证明没有 source coordinates 产品化。

若同一控制在全部运行状态也有实际 own caps \(\kappa_B(t)\)，(7) 和 (6) 给
\[
 \frac{\Phi(\mu)}{W_0}
 \le\mathbb E_Q\!\left[k_T\sum_{B\in\mathcal B_T}q_B(T)\kappa_B(T)\right]
   +\mathbb E_Q\int_0^T\sum_{B\in\mathcal B_t}
                      a_B(t)q_B(t)\kappa_B(t)dt.
\tag{14}
\]
也可保留更强的 signed \(R/W\) 取代第二正帽。这里需要证明的是**实际 overlap、packet shape caps 和真实 Q 占用**三者，不是把未证平方根估计改名；取“全部来源一组”时 own cap 就返回原未知问题，毫无收益。廉价合并若使 \(\kappa_B\) 升高，必须在运行和终端都支付。

E03 允许的有限深 Bellman 组织只是充分验证流程：在实际 \(Q\) 决策历史中，以 signed 阶段期望 \(\int R/W\) 为节点费，终端势支配 (13)，逐节点核“父势至少节点费加真实条件转移期望的子势”。死亡时间可以连续，不能擅改成自由有限分支树；有限深塔式期望保留终端后得到 (6) 的上界。合并节点物理费用为零，不能凭 coarse entropy 降低直接给负节点费；状态必须保存 fine 标签或等价 residual 信息。这一 Bellman 条件本身是形式接口，不是已构造的维数预算势。上述终端合同先在有界停止范围成立；无界推广仍需真实 \(Q\) 终值和 signed 积分的相应收敛条件。

## 9. 推进与剩余

可合法更换为整组控制，原细标签／RN／signed 费没有丢失；中途合并的 entropy barrier 也已重新计算。root 自身 residual 包络独核给运行费系数一，但没有为任意多格 group 给常数 own cap。因此现在能明确“多标签、多packet合作终态”究竟需要哪份证书，而非硬压成独立集。

尚无一般几何控制同时使 (14) 的终端项和运行项达到 \(\sqrt n\,n^{o(1)}\)。原 FIRST/fullfuture/CP--GP/LCA/allhistory 也未因 grouped clock 自动继承。这里的新合法模型和精确账可供下一步研究；主一般上界仍未闭合。无新增目标阶数、无数值执行或主账修改。

# 原域外联合差：截高可付子支与退出接口审计

2026-10-07。一般冻结原核的解析子支；没有证明 ordered 弱端点，更没有支付实际移动尺度的 \(R_{\rm angle}\)。只编辑本文件及同前缀新守卫，不重跑已有数值。采用 J03 的整个域外落点、杀死与未退出质量要求；L03 用于文末有限代数证书。

## 1. 同一障碍与准确目标

沿用 [联合谱差稿](ordered_frozen_spectral_difference_20261007.md) 及 [固定跳障碍稿](fixed_jump_obstacle_interface_20261007.md)：\(G_i=G_{c,i}\) 是原固定坐标对称、保质量正算子，\(A_i=I-G_i\)，\(S=\sum_iA_i\)，
\[
 K_r=\prod_i(I-rA_i),\quad H_r=e^{-rS},\quad D_r=K_r-H_r,\qquad 0\le r\le1.
\tag{1}
\]
原障碍 \(u\ge0\) 属于 \(L^1\cap L^2\)，\(\Omega=\{u>0\}\)，
\[
 \sigma=Su=\nu_b-\mu_b,\quad
 \nu_b\ge0\text{ 支撑于 }\Omega,\quad
 0\le\mu_b\le\kappa,\quad\mu_b=\kappa\text{ 于 }\Omega,
 \quad\int\nu_b=\int\mu_b=W_b.
\tag{2}
\]
所有等式是几乎处处等式。域外条件作用于整个 \(\Omega^c\)，而非几何边界；原来源没有依 receiver 重启。

待付的是 \(\Omega^c\) 上 \(\sup_r[D_r\sigma]_+\) 的弱水平集。已核联合谱工具给任意 \(w\in L^2\)
\[
 \left\|\sup_r|D_rSw|\right\|_2^2
 \le40\|\min(S,1)w\|_2^2\le40\langle w,Sw\rangle.
\tag{3}
\]
共同满测集上的 \(r\)-连续版本沿用谱差稿，不将仅有的强 \(L^2\) 连续误作点态 sup 的定义。本稿不重证、也不重跑该谱守卫。全文的空间弱费用来自 (3)，而非 source 标签的行概率界。

## 2. 真正可付的低势源算子片

固定确定高度 \(L>0\)，在原输入上定义
\[
 w_L=\min(u,L),\qquad v_L=(u-L)_+,
 \qquad D_r\sigma=D_rSw_L+D_rSv_L.
\tag{4}
\]
这是同一原势的算子分解，不是按接收点选新来源或另造一个 FIRST。两者均在 \(L^1\cap L^2\)；\(w_L=0\) 于原 \(\Omega^c\)。

对每个对称跳核，单调且 1-Lipschitz 的截断给
\[
 (w_L(x)-w_L(y))^2
 \le(w_L(x)-w_L(y))(u(x)-u(y)).
\]
核积分、求和合法，得到
\[
 0\le\langle w_L,Sw_L\rangle
 \le\langle w_L,Su\rangle
 =\int w_L\,d\nu_b-\kappa\int w_L
 \le L W_b.
\tag{5}
\]
这里 (2) 的饱和条件仅用于原 \(\Omega\)；没有用错误的原总势能上界。由 (3) 和 Chebyshev，对任意 \(a>0\)，
\[
 a\,\big|\{\sup_r|D_rSw_L|>a\}\big|
 \le {40L\over a}\,W_b.
\tag{6}
\]
特别 \(a=\kappa\)、\(L=\kappa Q_n\) 给 \(40Q_nW_b\) 的来源一次空间付款；\(Q_n\) 可预先固定为常数或 polylog。该子支适用于任意原非负密度，未限制原子数、来源深度或形状。

准确回到原域外目标：
\[
 2\kappa\,\big|\{x\in\Omega^c:\sup_r[D_r\sigma(x)]_+>2\kappa\}\big|
 \le {80L\over\kappa}W_b
 +2\kappa\,\big|\{x\in\Omega^c:\sup_r[D_rSv_L(x)]_+>\kappa\}\big|.
\tag{7}
\]
因此低势片确实已付，高势片仍需原域外弱预算。(7) 是冻结入口的子支，不是对实际 cube 账直接增加一笔费用；原 moving \(h_L\)、赢家、LCA/history 和 fullfuture 的连接尚未解决。

## 3. 高势尾的正来源范围与不可继承的 cap

凸函数 \(v(t)=(t-L)_+\) 满足逐点正算子不等式
\[
 Sv_L\le\mathbf1_{\{u>L\}}Su.
\tag{8}
\]
在 \(u=L\) 处选次梯度 0，此时 \(v_L=0\)、\(Sv_L\le0\)，所以边端也成立。令
\[
 W_L=\int_{\{u>L\}}d\nu_b.
\]
由 (8)、保质量与 \(v_L\in L^1\)，
\[
 (Sv_L)_+\le\mathbf1_{\{u>L\}}\nu_b,
 \quad\int Sv_L=0,
 \quad\int(Sv_L)_+=\int(Sv_L)_-\le W_L,
 \quad\|Sv_L\|_1\le2W_L.
\tag{9}
\]
这是实际原势上的正源预算，不是抽象源模型。特别 \(W_L\to0\) 随 \(L\to\infty\)，由对有限来源测度的单调收敛；但该事实没有给仅依赖 \(W_b\) 的下降速度。

**负部封顶范围必须保持原域。** 在原 \(\Omega^c\)，
\[
 Sv_L=-\sum_iG_iv_L\ge-\sum_iG_iu=-\mu_b\ge-\kappa.
\tag{10}
\]
不可将 (10) 延拓到整个空间或新的 \(\{u>L\}^c\)。可无条件证明的全空间粗界仅是
\[
 Sv_L\ge-(\kappa+nL).
\tag{11}
\]
若 \(u\le L\)，用 \(Sv_L\ge-\sum_iG_iu=Su-nu\ge-\kappa-nL\)；若 \(u>L\)，用 \(Sv_L=Su-Sw_L\)、\(Sw_L\le nL\)。尤其原 \(\{u>L\}\) 上 \(w_L=L\)、\(Sw_L\ge0\)，
\[
 Sv_L=\nu_b-\kappa-Sw_L.
\tag{12}
\]
将 \(\nu_b\mathbf1_{\{u>L\}}\) 作为候选新来源时，其在新活跃域内的诱导 cap 为 \(\kappa+Sw_L\)，一般大于 \(\kappa\)。故不能从原饱和障碍直接继承一个 cap 为 \(\kappa\) 的高势新障碍，也不能免费对 \(v_L\) 重启旧 killed-Abel 付款。

## 4. 多高度没有自动可和的来源份额

这不是把抽象尾例当 actual 反例。固定原坐标卷积核及任意固定 \(n,c,\kappa,W>0\)，取实际有限盒输入
\[
 \nu=H\mathbf1_B,\qquad |B|=W/H,\qquad H>\kappa.
\]
固定障碍存在性由固定跳稿原核构造给出。\(B\subset\Omega\) 几乎处处，因为 \(u=0\) 时 \(Su=-\sum_iG_iu\le0\)，不能承受原来源 \(H>\kappa\) 且 cap \(\le\kappa\)。盒内
\[
 Su=H-\kappa\le nu,\qquad u\ge(H-\kappa)/n.
\tag{13}
\]
给任何预先指定的有限高度 \(L_1,\ldots,L_m\)，选有限 \(H>\kappa+n\max_jL_j\)，便有
\[
 W_{L_j}=W_b=W\quad(1\le j\le m),\qquad
 \sum_{j=1}^mW_{L_j}=mW.
\tag{14}
\]
输入仍是原核上的有限盒、非负 \(L^1\cap L^2\) 密度。它否定的是从嵌套尾正来源预算直接宣称跨高度份额可和的交换；不是原 ordered 弱端点的反例，更不是具备全部原 FIRST/cube 门的反例。挑一个随输入增大的 \(L\) 虽能减小 \(W_L\)，却会增大 (6) 的 \(L/\kappa\) 费用；目前没有一般权衡关闭 (7)。

## 5. 根的新正预解障碍与一个额外负源子支

独立核对 [正预解重装稿](resolvent_obstacle_repacking_20261007.md)：\(R_t=(I+tS)^{-1}\)，\(h_t=R_t\sigma=(u-R_tu)/t\) 满足全空间 \(h_t\ge-\kappa\)、原 \(\Omega^c\) 上 \(h_t\le0\)、\(\int h_t=0\)、\(\int h_t^-\le W_b\)。
\[
 \nu'_t=\mathbf1_\Omega(\kappa+h_t),\qquad
 \mu'_t=\kappa\mathbf1_\Omega+\mathbf1_{\Omega^c}R_tu/t,
\]
确实给同一 \(u,\Omega\) 对生成元 \(B_t=(I-R_t)/t\) 的合法 cap 障碍，并有
\[
 \int\nu'_t=\int\mu'_t=\kappa|\Omega|+\int_{\Omega^c}h_t^-\le2W_b.
\tag{15}
\]
常数不属于 \(\mathbb R^n\) 的 \(L^1\)，下界以正算子的保常数延拓解释；质量恒等式使用可积的 \(h_t\)，不积分一个全空间常数。

固定 \(t=1\)，\(\mathcal N_r=D_r(I+S)\)，根已核 \(\|\sup_r|\mathcal N_r h|\|_2^2\le160\|h\|_2^2\)。于是全空间封顶的负来源片可直接支付：
\[
 \|h^-\|_2^2\le\kappa W_b,\qquad
 \kappa|\{\sup_r|\mathcal N_r h^-|>\kappa\}|\le160W_b.
\tag{16}
\]
完整差项是 \(D_r\sigma=\mathcal N_rh^+-\mathcal N_rh^-\)。因此其 \(2\kappa\) 弱水平集有一个 \(320W_b\) 的负源付款及仍未付的 \(\mathcal N_rh^+\) 正源水平集。这是 (7) 的另一个替代分解，不能将两分解当成同一交通必然依次互补的两笔费用。

正部 \(h^+\) 支撑于原 \(\Omega\)，质量至多 \(W_b\)，但不封顶。\(R_1\) 有 holding 原子 \(I/(1+n)\)，因此原盒尖峰的高密度可以保留在 \(h^+\) 中。合法重装并没有把整个源变成 \(L^2\) 封顶数据。

## 6. 为什么新退出结构尚不能付联合差

\(B_t\) 的核为 \(R_t/t\)，有 holding、坐标面及多坐标落点。它确可按原整个 \(\Omega^c\) 建 killed 过程；退出落点必须保留实际外部位置，另记势杀死及未退出质量。新来源一次 \(W'_t\le2W_b\) 是合法的，不能因此宣称所有新退出测度都是概率或认同旧 FIRST 标签。

困难仍在所需算子：\(D_r\) 依赖各个 \(a_i\)，不是 \(s=\sum_i a_i\) 的单变量函数。\(B_t\) 只依赖 \(s\)，并未保存足够的坐标轮廓使 \(D_r(I+S)\) 获得正的 free/killed 交换。新 killed 压缩也破坏自由坐标的交换性。原 Abel 的正来源层只对其对应的单生成元预解展开成立；把联合 \(K_r\) 插入层间需另证泄漏/重入预算。

具体差的 Duhamel 恒等式（自由算子交换，\(K_t^{\widehat i}=\prod_{j\ne i}(I-tA_j)\)）为
\[
 D_r=-\sum_i\int_0^r t\,H_{r-t}A_i^2K_t^{\widehat i}\,dt.
\tag{17}
\]
源是有符号的 \(A_i^2\sigma\)，并非一个可直接正退出结算的 source-once 流。展开 \(-A_i^2=2G_i-I-G_i^2\) 后，删除负项只得到粗 \(n\)-费用；在 \(\Omega^c\) 原 cap 只控制首个 \(\sum_iG_i u\)，不控制插入不同坐标前缀后的再进入高势。先释放一个总质量 \(O(W_b)\) 的退出来源，再调用其 ordered continuation 极大弱费用，则使用的正是尚待证明的 ordered 端点，形成循环。

故本轮严格收获是 (6)/(7) 的通用截高空间子支、(9)/(10) 的准确尾源范围，以及 (16) 的替代负源空间子支；高势/正预解正源的原域外差仍未付。没有把新 \(R_t\) 障碍或局部行条件登记为原 geom 已闭合。

## 7. 新守卫的范围与收据

先保存 [注册](exterior_difference_exit_registration_20261007.json)，再一次执行 [脚本](exterior_difference_exit_exact_guard_20261007.py)，[结果](exterior_difference_exit_results_20261007.json) 为 PASS：三轮有限维度 \(n=2,3,4\)，每轮 211 项，合计 633 项 Fraction 精确检查，无失败；另对应核原有限盒尖峰的维度 \(8,32,128\) 系数。执行已终态，无 live handle。脚本 SHA256 为 `46db8baf38047544d9bd8a8a166cbdf22b23cf3d0ee3d3500c96902a6506ec0a`。

新守卫只检查截高正算子代数及 (13) 的尖峰系数；不重跑原联合谱或根的预解守卫。有限对称正 Markov 模型不是实际 \(\mathbb R^n\) 原核/FIRST 样本，尖峰系数检查也没有数值重建原障碍。普遍命题由 §§2–5 的证明承担，数值不认证弱端点或原 geom 费用。

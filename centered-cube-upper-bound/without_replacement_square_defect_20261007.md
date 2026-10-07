# Without-replacement 的精确平方扣项与 source-square 转换缺陷

2026-10-07。只新增本稿及同前缀守卫。读取 late_multiface_source_budget 的 \(\mathcal R_k\)/without-replacement 损失、weak_normalized_projection 的实际加权赢家与重叠、deep_gram_bridge 的原同源/不同孩子合同，以及最新 actual \(R_{\rm angle}\) 主账。本轮使用 F05/G03 的谱与来源合同，已读 SKILL、provenance 与 cube-interface；只用已核的谱/Schur基础及真假 Gram 区别，不调用未定义 TC-A14、Bergman 或外引双谱定理。

结果：平方扣项恒等式与总谱 source 预算成立；一个局部 bilinear 部分真有一份 W 预算。把它直接换成原正来源/receiver 的产品余项则不成立：本稿给原 \(G_1\)、有限 L1 来源及互斥 receiver 分割的线性下界。特殊分割不是真实弱型赢家/实际 FIRST 门；原弱型与 cube 目标仍未否定。新的实际乘积/投影剩项准确保留，未新增主账费用。

## 1. 准确算子与没有舍掉的平方

固定原 \(c,L\)，令 \(G_i=G_{c,L}^{(i)}\)。它们为 commuting、self-adjoint Markov contractions，原 Fourier symbols \(g_i\in[0,1]\)。令
\[
 Q=\sum_iG_i,\quad S=nI-Q,\quad
 e_k=\sum_{|A|=k}G_A,\quad
 \mathcal R_k=\sum_{|A|=k}\sum_{i\in A}G_i^2G_{A\setminus\{i\}},
\]
\(e_0=I,e_{n+1}=0,\mathcal R_0=0\)。旧准确计数仍为
\[
 Qe_k=(k+1)e_{k+1}+\mathcal R_k.
\]
记 \(e_j^{ij}\) 是删除坐标 i,j 后的 elementary 多项式，\(e_{-1}^{ij}=0\)。**新的精确恒等式**为
\[
 \boxed{n\mathcal R_k-kQe_k
       =\sum_{i<j}(G_i-G_j)^2e_{k-1}^{ij}},
 \qquad 0\le k\le n-1.                                  \tag{1}
\]
逐 monomial 核对：\(G_i^2G_B,\ i\notin B,|B|=k-1\) 的系数是 \(n-k\)；squarefree degree k+1 的系数是 \(-k(k+1)\)。左侧用旧计数也恰为 \((n-k)\mathcal R_k-k(k+1)e_{k+1}\)。因此 (1) 不是 \(Q^k/k!\) 上界，也没有 without-replacement 指数损失。

令 \(E_k=e_k/\binom nk\)，\(\mathsf P=Q/n\)，则
\[
 \boxed{E_{k+1}=\mathsf P E_k-\mathsf A_k},\qquad
 \mathsf A_k=
 \frac{\sum_{i<j}(G_i-G_j)^2e_{k-1}^{ij}}
      {n(n-k)\binom nk}.                                \tag{2}
\]
\(\mathsf A_0=0\)。此处 \(\mathsf A_k\) 是新 square-defect 算子，不是原 endpoint mask A、诊断 K、parent 或 source。

所有乘子 commute，且各 \(e_{k-1}^{ij}\) 谱非负，所以 \(\mathsf A_k\) 为 PSD。其准确平方因子可取
\[
 \mathsf F_{ij,k}
 =\frac{(G_i-G_j)(e_{k-1}^{ij})^{1/2}}
        {\sqrt{n(n-k)\binom nk}},\qquad
 \mathsf A_k=\sum_{i<j}\mathsf F_{ij,k}^*\mathsf F_{ij,k}.
\]
这只是 Hilbert 谱正性。**\(\mathsf A_k\) 不保持逐点正。**

## 2. 总谱预算与真实的一份 source-square 付款

求和的准确范围必须是 \(k=0,\ldots,n-1\)：
\[
 \boxed{\sum_{k=0}^{n-1}\mathsf A_k
   =I-E_n-\frac Sn\sum_{k=0}^{n-1}E_k,\qquad
   0\le\sum_k\mathsf A_k\le I.}                         \tag{3}
\]
若右侧最后的和包含 \(E_n\)，身份就不成立。上界使用 \(E_n,S,\sum E_k\) 的共同谱非负性及 commutation；不是使用它们的点态序。

任意非负密度 \(\nu_b\in L^1(dx)\)，取 \(h=\sqrt{\nu_b}\in L^2\)、\(W_b=\int\nu_b\)。不需要 \(h\) 另在 L1。由 (3)，
\[
 \sum_k\|\mathsf A_k^{1/2}h\|_2^2
 =\sum_{k,i<j}\|\mathsf F_{ij,k}h\|_2^2\le W_b.          \tag{4}
\]
固定源始终同一 \(h\)，没有每 k 再归一化来源。

另有一个真正可用的一般 bilinear 弱式。对任意有界可测 \(g_k\)，若
\(\sum_k|g_k(x)|^2\le1\)，则
\[
 \left|\sum_k\langle g_kh,\mathsf A_kh\rangle\right|
 \le\left(\sum_k\langle h,\mathsf A_kh\rangle\right)^{1/2}
    \left(\sum_k\langle g_kh,\mathsf A_kg_kh\rangle\right)^{1/2}
 \le W_b.                                               \tag{5}
\]
第二步用每个 \(0\le\mathsf A_k\le I\) 与
\(\sum_k\|g_kh\|_2^2\le W_b\)。因此该部分的一倍乘积系数2可付 \(2W_b\)，不含 J、n 或来源标签复杂度。

## 3. 从 square source 到原正来源：不能删的 signed 产品余项

\(\mathsf A_k\) 是有限 signed convolution measure \(a_k\)，\(a_k(\mathbb R^n)=0\)。对上述 h，Fubini 精确给 L1 身份
\[
 \boxed{\mathsf A_k(h^2)
     =2h\,\mathsf A_kh+\mathcal L_k(h)},\qquad
 \mathcal L_k(h)(x)
   =\int a_k(dv)\,[h(x-v)-h(x)]^2.                       \tag{6}
\]
一般 h 用 L2 近似或直接 \(\int|a_k|(dv)\|h(\cdot-v)-h\|_2^2<\infty\) 定义；右侧可积。它的全空间积分是
\[
 \int\mathcal L_k(h)=-2\langle h,\mathsf A_kh\rangle\le0.
                                                               \tag{7}
\]
**只有全空间 signed 积分有该符号**；局部 \(g_k\) 门下没有此符号。PSD 不能把 \(a_k\) 当成非正 jump measure 或把 (6) 当 carré-du-champ 的正/负版本。

实际 source \(\nu_b\) 支撑原障碍 \(\Omega\)。若 \(g_k\) 支撑 \(\Omega^c\)，那么 \(g_kh=0\)，(5) 的整个 bilinear 项为零，准确剩下
\[
 \sum_k\int_{\Omega^c}g_k\,\mathsf A_k\nu_b
   =\sum_k\int_{\Omega^c}g_k\,\mathcal L_k(h).             \tag{8}
\]
所以域外转换不是“用 (4) 剩一点小误差”：在原 source 零域上，待估正源项全部就是 product defect。

原 \(G_c\) 没有坐标原子时，(1) 的正部分支撑 k-coordinate faces，负部分支撑 k+1-coordinate faces；两类 mutually singular。由 (2)，准确 Jordan 质量为
\[
 \|a_k^+\|=\|a_k^-\|=k/n,\qquad
 \|a_k\|_{\rm TV}=2k/n,\qquad
 \sum_{k=1}^{n-1}\|a_k\|_{\rm TV}=n-1.                  \tag{9}
\]
正部来自 \((n-k)\mathcal R_k\)，负部来自 \(k(k+1)e_{k+1}\)。不同 k 的面在总体求和中还能有符号抵消，不能据 (9) 否定联合预算；但不同 receiver 门会选择不同的 face，不能免费依赖此抵消。

## 4. 原核有限 L1 证据：任意接收分割不足以清算产品余项

此节针对一个**无条件 source-square 迁移合同**，不针对 actual winner：
\[
 \sum_k\int g_k\,\mathsf A_k\nu
 \stackrel{?}{\le}\operatorname{polylog}(n)\int\nu
 \quad(0\le g_k,\ \sum g_k\le1,\ g_k\nu=0).              \tag{10}
\]
它在原核上不成立。固定真实 \(c=L=1\)，原
\[
 G_1=w,\quad w(x)=\int_0^1e^{-|x|/s}ds,\quad
 \int w=1,\quad \|w\|_\infty=1,\quad \|w*w\|_\infty\le1.
\]
这等价于原 Fourier \(g_1(\xi)=\log(1+\xi^2)/\xi^2\) 的正 Laplace 混合；没有替代空间核。

给定 n≥2，置
\[
 \delta=\frac1{64n},\qquad \epsilon=\delta^2,\qquad
 \nu_\epsilon=(2\epsilon)^{-n}
       \mathbf1_{[-\epsilon,\epsilon]^n},\quad W=1.
\]
令 \(N_\delta(x)=\#\{i:|x_i|>\delta\}\)，
\(g_k=\mathbf1_{\{N_\delta=k\}}\) 对 k=1,…,n−1。
各接收集合互斥、\(g_k\sqrt{\nu_\epsilon}=0\)，且其积分有限，因为 \(\mathsf A_k\nu_\epsilon\in L^1\)。未要求这些 receiver 集合自身有限体积。

任一正 k-face component 有 k 个活跃轴，一轴密度 w*w、其余密度 w；每个轴 density≤1。初始 source 位移每轴≤ε。union bound 给在该 component 下
\[
 \Pr\{N_\delta=k\}\ge1-2k(\delta+\epsilon).
\]
任一负 k+1-face component 若落入 \(N_\delta=k\)，至少一个活跃 displacement 的绝对值≤δ+ε，故其概率≤\(2(k+1)(\delta+\epsilon)\)。用准确 Jordan 质量 (9)，得到
\[
 \int g_k\,\mathsf A_k\nu_\epsilon
 \ge\frac kn[1-2(2k+1)(\delta+\epsilon)]
 \ge\frac78\,\frac kn.
\]
最后一式由 \(k\le n-1,\ n\ge2\) 的有理界直接成立。因此
\[
 \boxed{\sum_{k=1}^{n-1}\int g_k\,\mathcal L_k(\sqrt{\nu_\epsilon})
      =\sum_k\int g_k\,\mathsf A_k\nu_\epsilon
      \ge\frac7{16}(n-1)W.}                            \tag{11}
\]
同时 (4) 总谱预算≤W，(5) bilinear 部分为零。

这是显式有限 L1 原空间 source 和原 G 的 **新 product-defect** 下界，不是旧障碍势能 L2 尖峰反例，也不是抽象直和/有限链。source 并未被认证为原饱和障碍坏来源 \(\nu_b\)；\(g_k\) 只在 \(\nu_\epsilon\) 支撑外，不能据此称它在既定障碍 \(\Omega^c\)。它也未被认证为原 K_r 最大赢家、normalized projection 或实际 FIRST/CPGP/LCA 门；(11) 不否定任何这些有额外结构的合同，尤其不是 actual cube 反例。它只说明“总谱 source-square + 接收分割”不能自行完成所需换测度。

## 5. Replacement family 的完整 signed 重组

迭代 (2)，
\[
 E_k=\mathsf P^k-
      \sum_{\ell=0}^{k-1}\mathsf P^{k-1-\ell}\mathsf A_\ell.
\]
令 \(\beta_{n,k}(r)=\binom nk r^k(1-r)^{n-k}\)，原
\[
 K_r=\sum_k\beta_{n,k}(r)E_k,\qquad
 W_r=[(1-r)I+r\mathsf P]^n.
\]
则
\[
 \boxed{K_r=W_r-
      \sum_{\ell=0}^{n-1}\mathsf A_\ell
                       B_\ell(r,\mathsf P)},\qquad
 B_\ell=\sum_{k=\ell+1}^n
            \beta_{n,k}(r)\mathsf P^{k-1-\ell}.          \tag{12}
\]
\(B_\ell\) 是正 substochastic convolution operator，谱也在[0,1]；故整个 defect \(W_r-K_r\) 为 PSD、≤I。它不是逐点正 kernel。不能由 \(K_r\le W_r\) 的 Hilbert 形式声称原 receiver 响应逐点支配。

给原 weak_normalized_projection 的有限网赢家测试 \(g_j=\tau M^{-1}\mathbf1_{A_j}\)，\(\sum_jg_j\le1\)、支撑原Ωc。准确空间 pairing 为
\[
 \sum_j\langle g_j,K_{r_j}\nu_b\rangle
 =\sum_j\langle g_j,W_{r_j}\nu_b\rangle
      -\sum_\ell\langle b_\ell,\mathsf A_\ell\nu_b\rangle,
\quad
 b_\ell=\sum_jB_\ell(r_j,\mathsf P)g_j.                  \tag{13}
\]
源只出现一次，receiver 先选 winner，未按 source 选 r。但原互斥 \(g_j\) 经不同正 \(B_\ell\) 投影后，\(b_\ell\) **不继承**
\(\sum_\ell b_\ell^2\le1\)，也不再只支撑 Ωc。若令
\[
 f_\ell(x)=\sum_j\beta_{n,\ell+1}(r_j)g_j(x),
\]
则
\[
 b_\ell=f_\ell+\mathsf P b_{\ell+1},\qquad
 \sum_\ell f_\ell\le1.                                \tag{14}
\]
这是同一固定 \(\mathsf P\) 的未来 occupation potential；从 \(\sum f_\ell\le1\) 到 source-weighted projected \(b_\ell\) 的平方费用需要新的投影重叠/停止结构，不能将它当点态 partition。它与原停止 union-versus-count 缺陷呼应，但 (14) 自身没有清算该缺陷。

把 (6) 代入 (13)，必须保留的联合 signed 项准确为
\[
 -\,2\sum_\ell\langle b_\ell h,\mathsf A_\ell h\rangle
 -\sum_\ell\int b_\ell\,\mathcal L_\ell(h),\qquad h=\sqrt{\nu_b}.
                                                               \tag{15}
\]
若另证 \(\sum_\ell\langle b_\ell h,\mathsf A_\ell b_\ell h\rangle\le C^2W_b\)，第一部分≤2CW_b；但这不是 (3) 的结论，第二部分仍必须同符号处理。本文不把两个未知强项分别设成免费的 test-energy 引理。针对原赢家的最小目标是 (15) 的联合正收益费用；(11) 不认证该特定 projected test 是否失败。

## 6. 对 actual 同源 Gram 的边界

原 fixed-c/L 结论 (1)–(15) 不含 moving physical L 或来源时间 c_t 的共同 source 选择。实际 \(R_{\rm angle}\) 的 two-source核仍含原 distinct-child LCA、共享 seed、失败币及全部 gate；deep_gram_bridge 已核非零实际 cross-child核不是 PSD。平方因子 \(\mathsf F_{ij,k}\sqrt{\nu_b}\) 不是这些 actual source-pair 的自动 Gram 特征。

要使用它，必须给一条真实映射，将原核/全部门的 signed 空间付款写为上述 source-square 配对加一个已支付的 \(\mathcal L\)/投影/物理尺度缺陷。现在没有该映射：(6) 给出即使固定原核也不能略掉的 product defect，(13) 给出真实 weighted selector 新产生的 projected test。不能把原 \(\nu_b\) 替成某 x 的 captured subsource，再重复各层 W；也不能将谱正 \(\mathsf A_k\) 当 actual source-pair positive Gram，或重复旧能量反例已排除的强预算。

因此本轮尚未付 (15)，没有新的 actual sourceonce \(\sqrt n\) 费用，主账仍 \(R_{\rm angle}\) 与 \(\gamma_n\le147/143\)。

## 7. 新三轮证书终态

只对新合同准备 deterministic Fraction 守卫：

* 三轮形式多项式维数3/4/5，直接展开两侧 monomial、核 (1)，不沿用旧 Q^k 比较数据。
* 对压缩谱维数8/32/128，核 (2)/(3) 的准确求和范围、PSD、总budget与 (12) 的 replacement identity。
* 在3/4/5维有限正 commuting fixture 上核 (4)/(5)/(6)，只作普适代数组件，不冒称原 G。
* 原核有限 L1 证据取 n8/32/128，核 (9)/(11) 的精确 union-bound 有理常数及 source/receiver分割字段。原 w≤1、w*w≤1 和连续 face 支撑由上文解析提供；无需再对 G 做浮点积分或重跑旧模式。

已先保存 without_replacement_square_defect_registration_20261007.json，再执行同前缀 exact_guard。**3947/3947 项 Fraction 谓词 PASS**，无随机种子、无 live handle：

|形式/有限fixture维数|压缩谱及原face维数 n|exact checks|原核 selected 认证下界|粗下界 \(7(n-1)/16\)|
|---:|---:|---:|---:|---:|
|3|8|244|878003/262144|49/16|
|4|32|793|62280395/4194304|217/16|
|5|128|2910|4083485483/67108864|889/16|

形式 monomial 两侧独立展开；compressed common spectrum 取全1、交替0/1、\(1/(2+i\bmod5)\)，核所有 k 的 recurrence、sum-range 与四个 r 的完整 replacement signed 身份。finite fixture 取 \(G_i=(2I+\mathrm{flip}_i)/3\)，其 source-square 总和分别 \(113/243,\ 515/324,\ 533/135\)，不超过同一输入的 \(W=5,15,33\)。这些是谱与乘积算术，不能将 finite fixture 称为原 G。

原face证书只计算解析 (11) 的 union-bound 有理下界，source 高度/体积归一化精确；source 高度大整数用 SHA/位数保存。没有对原核做浮点积分、原 winner 模拟或 actual 门筛查。原核 positivity、w≤1 与 atom-free face 性质由解析证明负责，守卫没有冒称独立数值求积。

脚本 SHA-256：
\(\texttt{53873916b56bf2ee95dbcd59c110b9b742398d43d7431d6379f48d57bf5cb74b}\)。
注册 SHA-256：
\(\texttt{edb439dc2ca27e04b582a58f7264ac49dc14cf95bb979629eed4c2eaf97fceba}\)。
结果保存在 without_replacement_square_defect_results_20261007.json。

这些守卫不认证原 winner、original obstacle/CPGP/FIRST 或 (15) 的量级，不重跑旧523/6025/谱差/投影数值。一般成立的收获是准确平方扣项、总source-square及 (5)；general product bridge (10) 被原G有限L1来源排除，原特殊赢家的联合 signed (15) 仍未付。actual cube 主账不变。

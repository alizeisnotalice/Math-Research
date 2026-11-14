# 空间后验混合必能检测严格容量改进，但尚无中央量级

本稿只使用真实有限矩形目录及其原始均匀行核。核心结论是：完整容量严格高于原始 KL，当且仅当空间后验凸包内已有严格改进。该结论不表示凸包包含容量最优标签权，也不保证改进达到所需的中央尺度。本轮凸包证书 的凸包限制容量上证书若严格低于完整容量，与本稿完全相容。

## 1. 记号与精确损失分解

删去零先验标签与零概率原子。设原子核为 $K_{iz}$，先验为 $\pi_i>0$，$p_z=\sum_i\pi_iK_{iz}>0$，均匀空间参考为 $Q_z$。真实后验为 $\alpha_i(z)=\pi_iK_{iz}/p_z$，$P_\nu=\sum_i\nu_iK_i$。记
\[
F(\nu)=D(P_\nu\Vert Q)-D(\nu\Vert\pi),\quad
F_0=F(\pi)=D(p\Vert Q),\quad
g_i=\sum_zK_{iz}\log(p_z/Q_z).
\]
则 $\sum_i\pi_i g_i=F_0$。对任意标签权 $\nu$，直接展开两个 KL 得
\[
F(\nu)-F_0=\langle\nu-\pi,g\rangle-J(\nu),\qquad
J(\nu)=D(\nu\Vert\pi)-D(P_\nu\Vert p)\ge0. \tag{1}
\]
这里 $J$ 是原始矩形核的数据处理亏损，而非额外假定。若 $\alpha^\nu(i\mid z)=\nu_iK_{iz}/P_\nu(z)$，KL 链式法则还给
\[
J(\nu)=\sum_zP_\nu(z)D(\alpha^\nu(\cdot\mid z)\Vert\alpha(\cdot\mid z)). \tag{2}
\]
不需要改变行核或把行截成命中源。

定义空间函数
\[
B_z=\sum_i\alpha_i(z)g_i,\quad b_z=B_z-F_0,\quad
M=\max_z|b_z|,\quad \sigma_B^2=\sum_zp_zb_z^2.
\]
其均值为 $F_0$。当 $M>0$，取
\[
t=\min\{1/2,1/(2M)\},\qquad
\tau_z=p_z(1+t b_z),\qquad \nu_i^t=\sum_z\tau_z\alpha_i(z).
\tag{3}
\]
$\tau$ 是严格正的空间概率权，所以 $\nu^t$ 真正属于全部原子后验的凸包。写 $d_i=\sum_zK_{iz}b_z$，则 $\nu_i^t=\pi_i(1+t d_i)$，且
\[
\langle\nu^t-\pi,g\rangle=t\sigma_B^2,\qquad
\sum_i\pi_i d_i^2\le\sigma_B^2. \tag{4}
\]
第一式由 $\sum_i\pi_iK_{iz}g_i=p_zB_z$ 得出；第二式是各真实行内的 Jensen 不等式。由 $\log(1+x)\le x$ 及 $\sum_i\pi_i d_i=0$，
\[
D(\nu^t\Vert\pi)
\le\sum_i\pi_i(1+t d_i)t d_i
=t^2\sum_i\pi_i d_i^2\le t^2\sigma_B^2.
\]
把 $0\le J\le D(\nu^t\Vert\pi)$ 代入 (1)，得到显式下界
\[
\boxed{\quad F(\nu^t)-F_0\ge t(1-t)\sigma_B^2
\ge\tfrac12t\sigma_B^2.\quad} \tag{5}
\]
这是可直接从真实原子核计算的改进证书，不涉及非线性优化器。

## 2. 严格改进的充要条件

由于 $F$ 凹且在 $\pi$ 处的梯度为 $g$（单纯形切方向），
\[
\max_\nu F(\nu)>F_0
\quad\Longleftrightarrow\quad
\operatorname{Var}_{\pi}(g)>0. \tag{6}
\]
若 $g$ 常数，切平面上证书直接给 $F\le F_0$；若不常数，沿 $\nu_i=\pi_i[1+\varepsilon(g_i-F_0)]$ 作足够小的正扰动，一阶增益是 $\varepsilon\operatorname{Var}_\pi(g)>0$，而 (1) 的亏损为 $O(\varepsilon^2)$。

关键是空间后验不会完全抹掉这个一阶方向。令 $L_z=\log(p_z/Q_z)$，则精确有
\[
\operatorname{Cov}_p(B,L)
=\sum_i\pi_i g_i^2-F_0^2
=\operatorname{Var}_\pi(g). \tag{7}
\]
因此 $\sigma_B^2=0$ 必推出 $\operatorname{Var}_\pi(g)=0$；反向也成立。结合 (5)—(6)，
\[
\boxed{\quad
\max_\nu F(\nu)>F_0
\iff
\max_{\nu\in\operatorname{conv}\{\alpha(z)\}}F(\nu)>F_0.
\quad} \tag{8}
\]
若 $\operatorname{Var}_p(L)>0$，还可量化为
\[
\sigma_B^2\ge
\frac{\operatorname{Var}_\pi(g)^2}{\operatorname{Var}_p(L)}. \tag{9}
\]
这只给严格改进或依赖目录的下界：没有证据保证 (5) 或 (9) 在近中央目录中具有 $\log n$ 量级，更没有得到 $\Lambda_{\rm LF}\ge\mu n-O(\log n)$。

## 3. 与 NP 后验、本轮单后验分析 和 本轮凸包证书 的区别

式 (8) 使用**全部空间原子**。若只允许 NP 源内的后验，不能直接沿用 (7)：该条件空间内的行平均已改变。若每行命中率恰为 $q$，自然命中空间权 $\tau^0=P_\pi(Z\in\cdot\mid E)$ 仍满足 $\sum_z\tau^0_z\alpha(z)=\pi$。在此基准上用 $B$ 作 (3) 的倾斜，同样得到 (5)，但方差改成 $\operatorname{Var}_{\tau^0}(B)$；它可能为零，不能从完整容量严格改进自动推出其为正。若行命中率只近似 $q$，基准标签权则是 $\nu^E$，不能冒写成 $\pi$。

本轮单后验分析 的自然 NP 后验平均受凹 Jensen 控制，确实不超过 $F(\pi)$（逐行命中恰为 $q$ 时）。本稿改进来自对空间权的非平凡倾斜，再取后验混合，既不是自然平均的升级，也不是取单点后验最大值。所有单点的 $F$ 都小于 $F_0$，并不阻止内部混合因凹性达到更大的值。

本轮凸包证书 研究凸包限制最大值是否等于完整容量；(8) 只回答二者是否都严格高于 $F_0$，不回答最大值相等。有限小维实际目录也远离一般近中央所需的指数规模，不能从任何这些测试外推中央稳定性。

## 4. 本轮独立审查记录与唯一余缺

只读审查 本轮单后验分析：点态恒等式、条件互信息平均、反例必要信息费用及 Jensen 亏损式均通过。实际调用其原有有理对数区间函数，重放全部 27 例并与保存 JSON 的全部确定性字段逐项比较，完全一致；没有修改 本轮单后验分析 文件。3 例 NP 单后验提高 KL；21 例完整容量严格高于 KL，其中 18 例所有单原子候选仍严格低于 KL。因此这 21 例由 (8) 都必有空间混合改进，而单点最大值遗漏了至少 18 例的改进方向。

另实际用 Python 对同一批 27 个真实原子核计算 (3) 的单步倾斜，浮点核对 (5)、(7) 全部通过；21 例增益超过 $10^{-12}$，这些例子的最小正下界为 $0.00012429296380508934$，最大实际单步增益为 $0.0049168518523842325$。这里只作为公式实现诊断，严格正性依据 (5)—(8) 与此前精确容量证书；浮点输出不承担区间认证或近中央证据的角色。

现在留下一个明确的中央缺口：在真实负责 LF 近中央目录中，能否量化空间倾斜的线性收益与数据处理亏损 $J$ 的竞争，使所得 $F$ 距离 $\mu n$ 仅 $O(\log n)$？本稿没有该估计，也没有相应近中央反例。式 (5) 不应被称作已经具有中央量级的容量稳定性。

## 有理候选的独立认证

`tilt_gain.py/json`将B按10^12尺度舍入后精确中心化为φ，并采用严格正有理空间权。对任意E_P φ=0、|tφ|≤1/2，原证明同样给
$$F(ν)−F_0≥(ν−π)·g−t²E_P φ².$$
全部量用有理对数区间核验，未把舍入B冒作精确B。27例中21例的直接增益与此解析方向增益均严格为正；`replay_gain.py`以独立64项Fraction重放。这里只认证固定有限目录的增益，不声明中央尺度。

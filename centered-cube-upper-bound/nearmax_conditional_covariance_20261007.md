# 真实赢家联合分布、空间 packet 碰撞与方差敏感切换

日期：2026-10-07。责任文件仅本稿；不修改主 TXT 或其他工作者文件。

状态：下面的条件期望恒等式、固定尺度预算、packet 碰撞预算和切换不等式均为自含证明。它们没有证明任意输入的 \(\sqrt n\,n^{o(1)}\) 弱型上界。新增可支付的部分是：一个真正共同的空间随机测试，其原赢家和有界 gap 候选的平均平方预测可以由常数倍完整来源质量 \(W\) 支付，不再先支付 \(I=mW\)。光滑内部赢家的 checker 二阶增益仍未证明。

本稿使用 M03、M02、I02、F05 skills 的输入与底测度核验工作流；只调用下文重证的 Jensen、Tonelli、Cauchy–Schwarz 和初等概率事实，不把 skill 的迁移接口当定理。查重对象包括 `high_marked_cover/source_operator_entropy_localization_audit.md`、20261002 的 `conditional_cube_transport/band_matrix/source_operator.md`、`sparse_pair_cover/cluster_operator_compression.md`，以及本目录 `radial_test_globalization_20261007.md`、`source_fragment_localization_20261007.md`。旧 atom 对角界及径向包络是背景构件；第 4 节的固定空间 packet 推广及其候选版本才是本稿新增接口。

## 1. 完整输入与真实联合分布

令 \(\mu\) 为 \(\mathbb R^n\) 上有限正 Borel 测度，\(W=\mu(\mathbb R^n)>0\)，\(n\ge1\)，\(0<a<b<\infty\)。全段 \([a,b]\) 均为原本允许的尺度。写
\[
 Q(x,r)=x+[-r/2,r/2]^n,\qquad
 U_r(x)=r^{-n}\mu(Q(x,r)),\quad M(x)=\max_{a\le r\le b}U_r(x).
\tag{1}
\]
采用原问题已建立的可测真实赢家 \(R(x)\)，保留所有原尺度。闭 cube 边界和来源原子均允许；后面的空间格采用 half-open 分割。固定 \(\tau>0\)，记
\[
 E=\{M>\tau\},\ I=\tau|E|>0,\ m=I/W,
 \quad m_x=\mu(Q(x,R(x)))=M(x)R(x)^n.
\tag{2}
\]
接收端始终是 Lebesgue 测度。令真实赢家后验为
\[
 \pi_x(dy)=\frac{1_{Q(x,R(x))}(y)}{m_x}\mu(dy),\qquad
 J(dx,dy)=\tau1_E(x)dx\,\pi_x(dy),\quad P=J/I.
\tag{3}
\]
于是 \(P_X=1_E dx/|E|\)，而
\[
 P_Y=\rho=\frac{S(y)}I\mu(dy),\quad
 S(y)=\tau\int_E\frac{1_{Q(x,R(x))}(y)}{m_x}dx,
 \qquad \nu=\mu/W.
\tag{4}
\]
在所有公式里来源质量只通过同一个完整 \(\mu\) 使用一次。来源坐标不要求独立。

对 \(\Phi_\tau(\mu)=\tau\int\log_+(M/\tau)dx\)，使用已证全局归约
\[
 K=\sup_{\mu\ne0}\Phi_\tau(\mu)/W,
 \qquad \Phi_\tau(\mu)\ge(K-\varepsilon)W.
\tag{5}
\]
此处 \(K\) 指固定窗口的全局分式上确界，未假设存在极值。已有 nearflat 结论给
\[
 \Delta=\int|S-K|d\mu\le3\sqrt{\varepsilon Z}\,W,
 \quad Z=1+n\log(b/a),\quad |I-KW|\le\Delta.
\tag{6}
\]
所以 \(\int|d\rho-d\nu|\le2\Delta/I\)。这是总变差的未除以 2 的 \(L^1\) 约定。对有界测试 \(f\)，两边积分差至多 \(2\Delta\|f\|_\infty/I\)。准确 source-flat \(S=m\) 给 \(\rho=\nu\)。

## 2. 条件期望 channel：合法恒等式及其限制

定义 \(T:L^2(\rho)\to L^2(P_X)\)，\(Tv(x)=\pi_xv\)。Jensen 给 \(\|T\|\le1\)。它的伴随在 \(S>0\) 上为
\[
 T^*f(y)=\frac\tau{S(y)}\int_E
 \frac{1_{Q(x,R(x))}(y)}{m_x}f(x)dx.
\tag{7}
\]
在 \(S=0\) 集任意定义，不影响 \(L^2(\rho)\)。因此 \(A=T^*T\) 是 \(\rho\) 上满足细致平衡（reversible）的正半定 Markov 算子，\(A1=1\)，核为
\[
 A(y,dz)=\frac\tau{S(y)}\int_E
 \frac{1_{Q(x,R(x))}(y)1_{Q(x,R(x))}(z)}{m_x^2}dx\,\mu(dz).
\tag{8}
\]
交换 \(y,z\) 后的 \(\rho(dy)A(y,dz)\) 相同。特别地
\[
 \langle v,(1-A)v\rangle_\rho
 =\mathbb E_X\operatorname{Var}_{\pi_x}v
 =\frac12\mathbb E_X\iint(v(y)-v(z))^2\pi_x(dy)\pi_x(dz).
\tag{9}
\]
式 (7)--(9) 仅是实际联合分布的恒等式。即使 \(\rho=\nu\)，它们也没有给出 \(K\) 上界或谱隙；将这种 Markov 结构本身称为新的弱费用支付不成立。

## 3. 固定尺度／固定概率尺度混合的真正 \(1/m\) 预算

取任一预先固定的概率卷积核 \(\kappa\)，令 \(H=\kappa*\mu\)，\(\pi_H(dy)=\kappa(x-y)\mu(dy)/H(x)\)（零响应处相关项定义为零）。对有界 Borel \(v\)，Jensen、\(\tau/M<1\) 和 \(\int\kappa=1\) 给
\[
 \begin{aligned}
 I\,\mathbb E_X\left[\frac HM(\pi_Hv)^2\right]
 &\le\tau\int_E\frac1M\int\kappa(x-y)v(y)^2\mu(dy)dx\\
 &\le\int v^2d\mu.
 \end{aligned}
\tag{10}
\]
所以 \(v\mapsto\sqrt{H/M}\,\pi_Hv\) 从 \(L^2(\nu)\) 到 \(L^2(P_X)\) 的范数至多 \(m^{-1/2}\)。核固定于来源和接收点之前；自适应真赢家 \(R(x)\) 不能免费代入。

同理对预先固定概率尺度分布 \(\omega\)，
\[
 \mathbb E_X\int\frac{U_s}M(\pi_sv)^2\omega(ds)
 \le\frac1m\mathbb E_\nu v^2.
\tag{11}
\]
\(v=1\) 恢复父任务已经提出的 sharp-peak 门：
\[
 H(x)=\frac1{\ell}\int_a^b U_s(x)\frac{ds}s,
 \quad\ell=\log(b/a),\qquad
 \mathbb E_X(H/M)\le1/m.
\tag{12}
\]
因此实际高弱比输入，uniform \(E\) 上的径向响应峰平均必须窄。旧 nearuniform 平滑宽峰的 global-test 反例没有这一高 \(m\) 资格，不能据此否定带 (12) 的接口。

顺带独立核验已在父稿 `posterior_depth_source_once_transfer_20261007.md` 使用的临界年龄恒等式。令 \(D=2\|x-Y\|_\infty\)，\(Y\sim\pi_x\)，\(A_R=n\log[R/\max(a,D)]\)。有限区间层蛋糕积分（端点原子不影响 \(ds\)）给
\[
 G_R(x)=\mathbb E_{\pi_x}e^{A_R}
 =1+n\int_a^{R(x)}\frac{U_s(x)}{M(x)}\frac{ds}s,
 \qquad \mathbb E_X G_R\le1+\frac{n\ell}{m}.
\tag{13}
\]
此处 \(G_R\) 是真实前历史面积，不是假设的自由 profile。式 (13) 约束高 \(m\) 下临界深度的平均，但没有证明深度下界与高 \(m\) 必然冲突。

## 4. 固定空间 packet 的后验碰撞预算

**命题。** 取
\[
 0<d\le\min\{a/(8n),(b-a)/4\}.
\tag{14}
\]
令 \(\mathcal C\) 为一固定可数 Borel 分割，每格的 \(\ell^\infty\) 直径至多 \(d\)，每格预先取一点 \(y_C\in\overline C\)，使所有 \(y\in C\) 满足 \(\|y-y_C\|_\infty\le d\)。例如固定 half-open 边长 \(d\) 的通常格及其左下角。记
\[
 w_C=\mu(C),\quad a_C(x)=\mu(C\cap Q(x,R(x))),\quad
 \chi_R(x)=\sum_C\left(\frac{a_C(x)}{m_x}\right)^2.
\tag{15}
\]
则
\[
 \boxed{\quad\tau\int_E\chi_R(x)dx\le6W,
 \qquad \mathbb E_X\chi_R\le6/m.\quad}
\tag{16}
\]

**证明。** 设 \(\alpha=n\log(1+2d/a)\le1/4\)，把 \(E\) 分为 \(R\le b-2d\) 和 \(R>b-2d\)。在第一部分，只要 \(a_C(x)>0\)，就存在 \(y\in C\cap Q(x,R)\)；直径条件使整个 \(C\) 包含于闭 cube \(Q(x,R+2d)\)。这对 half-open 来源格与闭接收 cube 的所有边界点同样成立。由于 \(R+2d\) 是合法原尺度，真实赢家帽给
\[
 w_C\le M(x)(R+2d)^n
 =m_x(1+2d/R)^n\le e^\alpha m_x.
\tag{17}
\]
阈值又给 \(m_x>\tau R^n\ge e^{-\alpha}\tau(R+2d)^n\)。并且
\(r_C(x)=2\|x-y_C\|_\infty\le R+2d\)。从而对 \(a_C>0\)，
\[
 m_x\ge e^{-\alpha}\max\{w_C,\tau r_C(x)^n\},\qquad
 (a_C/m_x)^2\le
 \frac{e^{2\alpha}w_C^2}{\max\{w_C,\tau r_C(x)^n\}^2}.
\tag{18}
\]
\(a_C=0\) 时左边为零，无需 (17)。\(w_C=0\) 的格直接忽略。

对 \(c>0\)，以 \(r=2\|x-y_C\|_\infty\) 分层，\(|\{r\le t\}|=t^n\)，故
\[
 \int_{\mathbb R^n}\frac{c}{\max(c,r^n)^2}dx
 =\int_0^\infty\frac{c\,n r^{n-1}}{\max(c,r^n)^2}dr
 =1+1=2.
\tag{19}
\]
在 (18) 取 \(c=w_C/\tau\)，得到每格内部费用至多 \(2e^{2\alpha}w_C\)。Tonelli 求和，只出现 \(\sum_Cw_C=W\) 一次。

顶部部分用 \(\chi_R\le1\)。对任何 \(0<c_0<b\)，概率 cube 核的逐点尺度上包络满足
\[
 \int\sup_{c_0\le r\le b}r^{-n}1_{Q(0,r)}(z)dz
 =1+n\log(b/c_0).
\tag{20}
\]
这是半径 \(2\|z\|_\infty\) 从 \(c_0\) 至 \(b\) 的直接积分。故通过完整 \(\mu\) 卷积并 Tonelli，
\[
 \tau|\{x\in E:R>b-2d\}|
 \le\left[1+n\log\frac b{b-2d}\right]W=C_{\rm top}W.
\tag{21}
\]
条件 (14) 确保 \(b-2d>a\)，包括 \(b\) 任意接近 \(a\) 的窗口；且 \(2d/b\le1/(4n)\)。由 \(-\log(1-u)\le u/(1-u)\) 得 \(C_{\rm top}\le4/3\)。因此总费用
\(2e^{2\alpha}+C_{\rm top}\le2e^{1/2}+4/3<6\)。证明完毕。

**有界 gap 候选推广。** 令 \(L(x)\in[a,b]\) 为任意可测原尺度候选，\(U_L>0\)，
\[
 g(x)=\log(M/U_L),\quad0\le g\le\eta,
 \quad \pi_{L,x}=\mu|_{Q(x,L)}/\mu(Q(x,L)),\quad
 \chi_L=\sum_C\pi_{L,x}(C)^2.
\tag{22}
\]
无需 \(L\ge R\)。内部 \(L\le b-2d\) 时，查询 \(L+2d\) 仍被原 \(M\) 控制；(18) 的右边只增加 \(e^{2\eta}\)。顶部 \(L>b-2d\) 时，\(U_L=e^{-g}M>e^{-\eta}\tau\)，故同一个顶部包络以 \(e^{-\eta}\tau\) 阈值支付 \(e^\eta C_{\rm top}W\)。因此
\[
 \boxed{\quad\tau\int_E\chi_L dx
 \le[2e^{2(\eta+\alpha)}+e^\eta C_{\rm top}]W
 \le6e^{2\eta}W.\quad}
\tag{23}
\]
有界 gap 仅在要积分的接收集合成立时，可在其补集取 \(L=R\)，保持 (23)。不能把一个实际稀疏尺度族的缺失查询 \(R+2d\) 免费补入；本命题明确依赖完整原连续窗口。

旧 physical-atom 对角界取每个来源原子为一格，原子一旦捕获就有 \(a_C=w_C\)，无需扩张尺度。本命题允许真实空间格被 cube 边界切开，使用合法 \(R+2d\) 罩住整格以及顶部单独费用；不把边界混格当原子。上式没有 \(N\) 或独立坐标假设，也不是旧 support-cover 的 \(\log K_{\rm cover}\) 算子范数结论。

## 5. 真正共同的随机测试与 conditional variance

给固定可数格独立 Rademacher 标记 \(\xi_C\)，定义同一个全局来源函数 \(v_\xi(y)=\xi_C\) 对 \(y\in C\)。每个实现均为 Borel 函数且 \(|v_\xi|=1\)；不是接收点另选符号。独立性仅是外加测试随机性，绝不施于 \(\mu\) 的坐标。Tonelli 给
\[
 \mathbb E_\xi\mathbb E_X(Tv_\xi)^2=\mathbb E_X\chi_R\le6/m,
 \qquad
 \mathbb E_\xi\mathbb E_X\operatorname{Var}_{\pi_x}v_\xi
 =1-\mathbb E_X\chi_R\ge1-6/m.
\tag{24}
\]
式 (24) 对全输入成立，nearflat 并非其前提。这是一个实际共同测试的平方预算；它不宣称 \(T\) 的全算子范数小于 1。

若准确 \(\rho=\nu\)，则 \(\mathbb E_X\pi_x(C)=\nu(C)\)，Jensen 给
\[
 \sum_C\nu(C)^2\le\mathbb E_X\chi_R\le6/m.
\tag{25}
\]
令 \(\widetilde v_\xi=(v_\xi-\nu v_\xi)/2\)，则 \(|\widetilde v|\le1\)，\(\nu\widetilde v=0\)，
\[
 \mathbb E_\xi\mathbb E_X(T\widetilde v)^2\le\frac{3}{2m},\quad
 \mathbb E_\xi\mathbb E_\nu\widetilde v^2
 =\frac{1-\sum_C\nu(C)^2}{4}\ge\frac{1-6/m}{4},
\quad
 \mathbb E_\xi\mathbb E_X\operatorname{Var}_{\pi_x}\widetilde v
 \ge\frac{1-6/m}{4}.
\tag{26}
\]
nearflat 时 (25) 的无误差写法不成立。可改用 \(\rho\) 中心化，保留测试幅度和 (24) 的平方预测上界；若必须使用 \(\nu\) 中心化，需显式支付 (6) 的比较误差。下面的全局分式切换不要求单个扰动保持质量，故原 \(v_\xi\) 已足够。

## 6. 方差敏感的对称 switch 合同

令 \(v\) 为任意同一个全局 Borel 来源函数，\(|v|\le1\)，\(0<t<1\)，\(\mu_\pm=(1\pm tv)\mu\)。取候选 (22)，写
\[
 p_R=\pi_xv,\quad p_L=\pi_{L,x}v,\quad d_v=p_R-p_L,
 \quad Q_R(v)=\mathbb E_Xp_R^2,\quad Q_L(v)=\mathbb E_Xp_L^2.
\tag{27}
\]
则
\[
 \Phi_\tau(\mu_+)+\Phi_\tau(\mu_-)-2\Phi_\tau(\mu)
 \ge\tau\int_E[t|d_v|-g]_+dx
 -\frac{t^2\tau}{2(1-t)^2}\int_E(2p_R^2+p_L^2)dx.
\tag{28}
\]
**证明。** 对 \(|p|\le1\)，Taylor 积分余项给
\(\log(1\pm tp)\ge\pm tp-t^2p^2/[2(1-t)^2]\)。在每个 \(x\in E\)，两边都保留原 \(R\) 时，两项一次部分相消；若 \(t|d_v|>g\)，在有利的一个扰动里改用合法 \(L\)，一次部分增加 \(t|d_v|-g\)。不切换时平方余项是 \(2p_R^2\)，切换时是 \(p_R^2+p_L^2\)，均由 \(2p_R^2+p_L^2\) 控制。截断 \(\log_+\) 只提高这一响应下界，\(E\) 外的原基值是零而新值非负。积分即 (28)。

由 \(W_++W_-=2W\) 和 (5)，左边至多 \(2\varepsilon W\)，即
\[
 \mathbb E_X[t|d_v|-g]_+
 \le\frac{2\varepsilon}{m}
 +\frac{t^2}{2(1-t)^2}\{2Q_R(v)+Q_L(v)\}.
\tag{29}
\]
这是比旧统一 \(t^2 I\) 费用更细的合同：固定赢家曲率是实际 \(p_R^2\)，不是默认 1。

例如若 \(Tv=0\) 且 \(|v|\le1\)，令 \(D=\mathbb E_X|p_L|\)，\(G=\mathbb E_Xg\)。由 \(p_L^2\le|p_L|\) 和 \([t|p_L|-g]_+\ge t|p_L|-g\)，在 \(t=1/4\) 时得到
\[
 D\le\frac{36}{7}\left(G+\frac{2\varepsilon}{m}\right).
\tag{30}
\]
但未经证明存在可检测几何差别的真实全局 \(v\in\ker T\)，不能将 (30) 当作 \(K\) 上界。

将 (28) 用于第 5 节同一 grid 的随机测试，记
\[
 B(x)=\left(\sum_C[\pi_x(C)-\pi_{L,x}(C)]^2\right)^{1/2}.
\tag{31}
\]
独立符号和四阶矩给 \(\mathbb E_\xi|d_v|\ge B/\sqrt3\)：二阶矩为 \(B^2\)，四阶矩至多 \(3B^4\)，Hölder 插值给所述下界。可数和由 \(\ell^1\) 收敛或有限截断极限核验。再对凸 hinge 使用 Jensen，并代入 (16)、(23)，得到
\[
 \boxed{\quad
 \tau\int_E\left[\frac{tB(x)}{\sqrt3}-g(x)\right]_+dx
 \le 2\varepsilon W+
 \frac{3t^2}{(1-t)^2}(2+e^{2\eta})W.
 \quad}
\tag{32}
\]
这是真实 nearmax 对全部共同测试的一个已付 \(W\) 级 packet switch 约束；没有付来源平方 \(I_2\)，也没有把物理原子数 \(N\) 换成维数 \(n\)。它要求候选本身有实际 gap 上界，并未证明必须存在使左边很大的候选。

## 7. 已推进的桥与仍缺的一步

高 \(m\) 的真实输入同时满足：固定混合的径向面积小 (12)，以及固定小空间格的 posterior 平方碰撞小 (16)。在 nearflat 情形，共同 checker 测试仍保有常数级来源方差和条件方差，固定赢家预测只有 \(O(1/m)\)。候选有界 gap 时，nearmax 要求其 packet posterior 差服从 (32)。这把原先未付的 \(I\) 级曲率缩为 \(W\) 级，属于实质一般接口；条件方差恒等式本身不应算作同样成果。

缺口是：必须由真实连续尺度几何与 sharp-peak 门推出某些有界 gap 候选，其共同 packet 差的 hinge 积分违反 (32)。窄峰只约束响应面积，尚未给 packet 差的下界；高条件方差也只表达来源在赢家 cube 中的不可预测性。当前不能将两者相乘宣称 \(\sqrt n\) 支付。现有真实宽峰反例缺高 \(m\) 资格，亦不能直接证明这个组合不可能。

另一个待证结构候选是：对足够光滑的密度、唯一且非退化的内部赢家，极细 checker 的赢家位移二阶增益可能超过固定赢家曲率。要变成定理，至少需证明同一个固定全局分割下的边界格导数平方下界、连续 winner 的扰动正则性、\(\partial_r^2 U_R\) 的定量非退化界、局部正增益在 Lebesgue 接收集合上的可积下界，并控制该集合外的负曲率。仅猜测 \(p_R^2\sim h^n\)、\((\partial_Rp_R)^2\sim h^{n-1}\) 不足以完成这些步骤；本稿不声称排除这类局部极值。原子近极值与边界赢家也不在此候选的假设中。

本稿没有启动新数值。后续守卫应先登记真实连续 winner 的 packet 扩张查询、半开格边界和顶部费用，并复用实际保存来源与响应；稀疏有限 \(J\) 的 winner 数据不能用于认证 (17)。

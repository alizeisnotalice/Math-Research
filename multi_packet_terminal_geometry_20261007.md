# 完整 packet 的锐化空间收费与多格终端的分组合并预算

2026-10-07。本稿只新增此文件，不修改旧冻结材料或主 TXT。使用已读 E03、A02、D04 的同来源望远镜、正支配和熵链式审计规则；不引未核文献定理，不启动数值。删除余费的锐化与真实输入锐度构造由 root 提出，本稿逐项独审并给出完整证明；分组图预算及合并判据在该锐化上推导。

**结论与范围。** 任一完整 packet 的实际删除余费由其 own-logmax 费用以系数 1 支付。对任意有限完整 packet 分割、已证 own-logmax 上帽及覆盖真实共捕获的图，得到
\[
 \Phi_\tau(\mu)\le\sum_B W_B\kappa_B(1+H_{\deg(B)}).
\tag{1}
\]
它保留 packet 内的全部原来源，可以在终态保留很多细格。下面给出使 (1) 的右边确实下降的合并判据。没有证明任意输入存在平方根费用的分割；一般 packet 的 \(\kappa_B\) 不能免费置为常数。静态分组可以在运行过程前依完整输入反复选择，本文不新增随机途中合并的控制 law。

## 1. 全输入、原窗口与完整 packet

令 \(\mu\) 为 \(\mathbb R^n\) 上有限正 Borel 测度，\(W=\mu(\mathbb R^n)>0\)，\(0<a<b<\infty\)，\(\tau>0\)。沿用全部原尺度
\[
 Q(x,r)=x+[-r/2,r/2]^n,\quad U_r(\nu)(x)=r^{-n}\nu(Q(x,r)),
 \quad M_\nu(x)=\max_{a\le r\le b}U_r(\nu)(x),
\]
\[
 \Phi_\lambda(\nu)=\lambda\int_{\mathbb R^n}\log_+(M_\nu/\lambda)dx,
 \qquad Z=1+n\log(b/a).
\tag{2}
\]
接收端始终为 Lebesgue。来源可以是原子、奇异或连续测度，不假设坐标独立。已有窗口包络 \(\int M_\nu\le ZW(\nu)\) 给
\[
 0\le\Phi_\lambda(\nu)\le (Z/e)W(\nu),\qquad
 |\Phi_\lambda(\nu)-\Phi_\lambda(\sigma)|\le Z\|\nu-\sigma\|_{\rm TV}.
\tag{3}
\]
后式由 \(|M_\nu-M_\sigma|\le M_{|\nu-\sigma|}\) 及 \(t\mapsto\lambda\log_+(t/\lambda)\) 的 1-Lipschitz 性得到。

完整 packet 是一份固定正来源分量 \(\sigma\le\nu\)，不是只取当前赢家 cube 中的截断。令其质量为 \(w>0\)，定义其实际原窗口费用
\[
 K(\sigma)=\sup_{\lambda>0}\frac{\Phi_\lambda(\sigma)}w,
 \qquad A(\sigma)=\sup_{\lambda>0}\frac{\lambda|\{M_\sigma>\lambda\}|}w.
\tag{4}
\]
它们满足 \(K(\sigma)\le A(\sigma)\le Z\) 与 \(K(\sigma)\le Z/e\)，且对整体乘法缩放 \(c\sigma\) 不变。确实，\(\Phi_\lambda(c\sigma)=c\Phi_{\lambda/c}(\sigma)\)；又对数层蛋糕与 \(\lambda|\{M_\sigma>s\}|\le A(\sigma)w\lambda/s\) 给 \(\Phi_\lambda(\sigma)\le A(\sigma)w\)。这里只对同一实际 packet 改变阈值，未扩大原尺度窗口。

## 2. 实际删除余费：一个 sharp packet 引理

固定完整当前来源 \(\nu\) 及 \(0\le\sigma\le\nu\)，令 \(H=M_\sigma\)，\(M=M_\nu\)，\(E=\{M>\tau\}\)。在 \(E\) 取原真赢家 \(R(x)\)，并令
\[
 p(x)=\frac{\sigma(Q(x,R(x)))}{\nu(Q(x,R(x)))},\quad
 D=\Phi_\tau(\nu)-\Phi_\tau(\nu-\sigma),\quad
 L=\tau\int_E p(x)dx.
\tag{5}
\]
\(D\ge0\)，但 \(D-L\) 可能为负。定义非负函数
\[
 G(v)=\log(1+v)-\frac{v}{1+v},\qquad
 \Gamma_\tau(\sigma)=\tau\int_{\mathbb R^n}G(H/\tau)dx.
\tag{6}
\]
则有完整空间、完整输入的不等式
\[
 \boxed{D-L\le\Gamma_\tau(\sigma)\le K(\sigma)w\le A(\sigma)w.}
\tag{7}
\]

**冻结原真赢家。** 删除后的响应至少为 \((1-p)M\)，因为这个原尺度依然合法。故每个 \(x\in E\)，写 \(u=M/\tau>1\)、\(v=H/\tau\)，有
\[
 \log(M/\tau)-\log_+(M_{\nu-\sigma}/\tau)
 \le \min\{\log u,-\log(1-p)\}.
\tag{8}
\]
当 \(p=1\) 时把 \(-\log(1-p)\) 理解为 \(+\infty\)，右边仍为有限 \(\log u\)。当 \(p=0\) 时右边为 0。来源正性与原赢家性给
\[
 0\le p\le\min(1,v/u),\qquad u\ge v.
\tag{9}
\]
在 \(M\le\tau\) 上，删除前后两个 log-positive 均为 0；无需在那里定义后验。

**准确标量优化。** 设 \(u\ge\max(1,v)\)。对
\(h_u(p)=\min(\log u,-\log(1-p))-p\)，它在 \(0\le p\le p_*:=1-1/u\) 单调增加，在 \(p_*\le p\le1\) 单调下降。当 \(u\le v+1\) 时，\(p_*\) 可行，最大值为 \(\log u-1+1/u\)，随 \(u\ge1\) 增加。当 \(u\ge v+1\) 时，最大可行值在 \(p=v/u\)，其值 \(-\log(1-v/u)-v/u\) 随 \(u\) 减少。两支在 \(u=v+1\) 相接，得到
\[
 \sup_{u\ge\max(1,v)}\ \sup_{0\le p\le\min(1,v/u)}h_u(p)=G(v).
\tag{10}
\]
\(v=0\) 时 \(p=0\)、两边为 0；\(u=1\) 时阈值帽为 0，\(h=-p\le G(v)\)。因此阈值等号不会产生遗漏。由 (8)–(10) 积分并外放非负 \(G\) 到全空间，得到 (7) 第一项。

**弱费用。** \(G'(v)=v/(1+v)^2\)，非负层蛋糕及 Tonelli 给
\[
 \Gamma_\tau(\sigma)
 =\tau\int_0^\infty G'(t)|\{H>\tau t\}|dt
 \le A(\sigma)w\int_0^\infty\frac{G'(t)}t dt=A(\sigma)w.
\tag{11}
\]
**更强的 own-logmax 费用。** 准确恒等式为
\[
 G(v)=\int_0^\infty\frac{2t}{(1+t)^3}\log_+(v/t)dt.
\tag{12}
\]
对 \(v>0\)，右边导数为 \(v^{-1}\int_0^v2t/(1+t)^3dt=v/(1+v)^2\)，两边在 0 的极限均为 0；这证明 (12)。全部被积函数非负，故包括无界接收空间的 Tonelli 合法，且
\[
 \Gamma_\tau(\sigma)
 =\int_0^\infty\frac{2}{(1+t)^3}\Phi_{\tau t}(\sigma)dt
 \le K(\sigma)w,
 \qquad\int_0^\infty\frac{2}{(1+t)^3}dt=1.
\tag{13}
\]
有限性由 (3) 或 (11) 保证，不需先假设 \(H\) 有有限支撑。

## 3. 系数 1 的完整真实输入锐度

取完整来源
\[
 \nu=\tau\mathbf1_{[-b,b]^n}dx+w\delta_0,
 \qquad\sigma=w\delta_0.
\tag{14}
\]
背景有限且保持不动。当 \(\|x\|_\infty\le b/2\) 时，每个合法 cube 都完全落在背景盒内，背景响应恰为 \(\tau\)。能够捕获原子的最小原尺度是 \(r_* =\max(a,2\|x\|_\infty)\le b\)，故
\[
 H=w/r_*^n,\quad M=\tau+H,\quad R=r_*,\quad p=H/(\tau+H).
\tag{15}
\]
这实现 (10) 的最优 \(u=v+1\)。外域无法捕获原子，删除前后完全相同；且所有背景 cube 响应不超过 \(\tau\)。因此 (7) 第一项在此准确等号：\(D-L=\Gamma_\tau(\sigma)\)。边界由闭 cube 保留，但其 Lebesgue 质量为零。

令 \(V=w/(\tau a^n)\)、\(m=w/(\tau b^n)\)。因为 \(|\{x:2\|x\|_\infty\le r\}|=r^n\)，核心及外围的准确积分给
\[
 \frac{D-L}{w}=\frac{G(V)}V+\int_m^V\frac{G(v)}{v^2}dv
 =\frac{\log(1+m)}m-\frac1{1+V}.
\tag{16}
\]
第二等号由分部积分与 \(G'(v)/v=(1+v)^{-2}\) 得到。

同一个原子 packet 的 own-logmax 费也可自含算出。写 \(q=\lambda b^n/w\)、\(r=a^n/b^n\)。当 \(0<q\le1\)，
\[
 \frac{\Phi_\lambda(w\delta_0)}w=q\{\log(1/q)+1-r\};
\tag{17}
\]
其峰在 \(q=e^{-r}\)，值为 \(e^{-r}\)。当 \(1\le q\le1/r\)，该比值为 \(1-qr\le1-r<e^{-r}\)；更高阈值为 0。因此 \(K(w\delta_0)=e^{-a^n/b^n}\)。

取 \([a,b]=[1,2]\)、\(w=1\)、\(\tau=2^{-n/2}\)，沿偶 \(n\to\infty\) 有 \(V=2^{n/2}\)、\(m=1/V\)，(16) 趋于 1，而 \(K(\sigma)\to1\)。所以 (7) 中 own-logmax 的系数 1 不能在所有维数、所有实际完整输入上统一换成任何 \(c<1\)。背景质量虽大但有限。该构造没有 nearmax 资格，也不是原弱型常数的一般下界。

## 4. 实际碰撞与合并费用必须保留

对完整当前来源 \(\nu\) 的一个来源分割 \(\mathcal P\)，定义原真赢家后验份额 \(p_B\)，及实际空间量
\[
 C(\mathcal P;\nu)=\tau\int_E\sum_Bp_B^2dx,\qquad
 C_z(\mathcal P;\nu)=\tau\int_E\log(M/\tau)\sum_Bp_B^2dx.
\tag{18}
\]
固定同一 \(\nu,E,R\)，把两个 packet 合并的准确增量为
\[
 C(\mathcal P';\nu)-C(\mathcal P;\nu)=2\tau\int_Ep_Ap_Bdx,
\]
\[
 C_z(\mathcal P';\nu)-C_z(\mathcal P;\nu)
 =2\tau\int_E\log(M/\tau)p_Ap_Bdx.
\tag{19}
\]
它们非负。多次合并只在这一个固定当前输入上望远镜。把全部来源合成一包后，\(C=\tau|E|\)、\(C_z=\Phi_\tau(\nu)\)；一般不能用某个输入无关常数乘 \(W\) 支付。若输入随增长或删除变化，不能把不同状态的 (19) 当作一次免费的望远镜。

作为粗帽比较，(9) 与弱层蛋糕给
\[
 \tau\int_Ep_B^2dx\le2A(\nu_B)W_B,
 \quad \tau\int_E\log(M/\tau)p_B^2dx
 \le(2/\sqrt e)A(\nu_B)W_B.
\tag{20}
\]
第二式用 \(v^2\sup_{u\ge\max(1,v)}\log u/u^2\)：它在 \(v\le\sqrt e\) 为 \(v^2/(2e)\)，在 \(v\ge\sqrt e\) 为 \(\log v\)；导数除以 \(v\) 的全积分为 \(2/\sqrt e\)。这些粗帽正确但不如 (7) 的实际删除余费帽，不把它们作为最终时钟系数。

## 5. Sharp 费用进入有限控制：收费的是完整 packet

固定有限 Borel 来源分割 \(\mu=\sum_{B\in\mathcal P}\mu_B\)，\(w_B>0\)。packet 为原测度的完整限制；没有按 receiver 改变分组。设已有确定上帽 \(\kappa_B\ge K(\mu_B)\)。每包按接收点无关的 \(a_B(t)\in[0,1]\) 整体增长，以同速率整体删除；控制范围与已证 LC 相同：有限个确定刷新及至多包数次已接受死亡刷新，政策只读 raw accepted history 和已揭示的独立外部决定。无跳段速率固定，不读取未来死亡或隐藏 tag。

在一个无跳段，原真赢家冻结给 \(\Phi_\tau\) 的右导数至少为
\(\sum_Ba_BL_B\)：因为冻结的赢家响应对整包增长的对数导数是 \(\sum_Ba_Bp_B\)，阈值边界的新增贡献非负。有限时域 \(W_t\le e^hW\)，且 (3) 给路径函数绝对连续，可在几乎处处导数上使用此下界。删除包 \(B\) 的实际损失是 \(D_B\)，(7) 给
\[
 \sum_Ba_BD_B\le\sum_Ba_BL_B+\sum_Ba_B\Gamma_\tau(\nu_B)
 \le\sum_Ba_BL_B+\sum_Ba_B\kappa_BW_B(t).
\tag{21}
\]
包内自运行以来只统一乘法增长，故 \(K(\nu_B)=K(\mu_B)\)，上帽没有被改变内部相对权重。

沿 LC 已自证的有限 first-jump 密度，把无跳增益与每个已接受删除的损失拼接，有限跳补偿得到对有界停时 \(S\)
\[
 \boxed{\quad
 \mathbb E_P\Phi_\tau(\nu_S)\ge\Phi_\tau(\mu)
 -\mathbb E_P\int_0^S\sum_Ba_B(t)\Gamma_\tau(\nu_B(t))dt
 \ge\Phi_\tau(\mu)-\mathbb E_P\int_0^S\sum_Ba_B(t)\kappa_BW_B(t)dt.
 \quad}
\tag{22}
\]
所有有限时域损失和增长变化均由 (3)、\(e^hW\)、有限包数控制；因此这里只复用已证明的 finite-refresh law 和有限 first-jump 拼接，不诉诸一个未构造的任意 predictable 控制过程。质量补偿仍给 \(W_t\) 为鞅。

细格几何证书是一个例子而非额外假设：若完整 packet 位于边长 \(d\) 半开轴格，令 \(\alpha=n\log(1+2d/a)\)，则其 own envelope 为
\[
 M_{\mu_B}(x)\le\frac{e^\alpha w_B}{\max(a,2\|x-c_B\|_\infty)^n},\qquad
 \Phi_\lambda(\mu_B)\le(e^\alpha w_B-\lambda a^n)_+.
\tag{23}
\]
证明：捕获正质量的 \(Q(x,r)\) 满足 \(2\|x-c_B\|_\infty\le r+2d\)，再按 \(r_c^n\) 的全空间径向 Lebesgue 积分；核心和外围对数项准确抵消。故 \(\kappa_B=e^\alpha\) 对所有阈值合法。外放的只是 envelope 积分空间，原允许尺度仍为 \([a,b]\)。相比旧 LC 的粗系数 16，这个几何证书在 (22) 中直接只付 \(e^\alpha\)；它本身不产生平方根维数界。

## 6. 任意完整 packet 分割的调和图预算

对有限固定分割选无向图 \(G\)，要求覆盖所有真实 co-capture：若某个 \(Q(x,b)\) 同时捕获两包的正质量，必须有边。允许额外边；删除来源正质量不会创造新的危险边。无需 packet 小直径或有限支撑。定义 \(a_B=1\) 当且仅当包尚存且有尚存邻居，否则为 0；每次已接受死亡刷新。令 \(T\) 为存活图首次无边时刻。

这里复述所需停止证据，防止免费丢掉尾项。原 \(P\) 下有边状态至少有两个 active 包，总死亡率至少 2，至多 \(N-1\) 次删除。所以
\[
 T\le\sum_{j=1}^{N-1}\operatorname{Exp}_j(2),\quad
 \mathbb E_Pe^{qT}\le(2/(2-q))^{N-1}<\infty\quad(q<2).
\tag{24}
\]
\(W_{t\wedge T}\le We^T\)，取 \(1<q<2\) 得一致可积；每个 packet 质量同样是停止 UI 的非负鞅。于是 \(\mathbb E_PW_B(T)=w_B\)，而不是只假设总终态质量守恒。

质量倾斜的有限时域密度为 \(dQ_h=(W_h/W)dP\)。自含共同路径构造是：按 \(w_B/W\) 抽原 packet tag \(J=B\)，使该包永生，其余包按同一 raw-history 政策运行。给定 \(J=B\) 的有限路径似然比是 \(e^{\int_0^h a_Bdt}\mathbf1_{B\text{ 尚存}}\)；混合恰为 \(W_h/W\)，故隐藏 tag 的后验为 \(W_B(h)/W_h\)。给定 tag，每个初始邻居只要尚存就有这个永生邻居，因此始终 active，寿命为相互独立的单位指数。tag 的活动时间正好是这些邻居寿命的最大值，其期望
\[
 \mathbb E_{Q(\cdot\mid J=B)}\int_0^T a_B(t)dt
 =H_{\deg_G(B)},\qquad H_0=0,\quad H_k=\sum_{j=1}^k1/j.
\tag{25}
\]
这是 tag 的局部时钟，并非全局停止时间。\(Q\) 下有边状态仍有至少一个可死亡的非 tag 包，总死亡率至少 1，故 \(Q(T<\infty)=1\)。非负 Tonelli、有限时域密度及 tag 后验得到来源一次的精确费用
\[
 \mathbb E_P\int_0^T\sum_Ba_B\kappa_BW_B(t)dt
 =W\mathbb E_Q\int_0^T a_J(t)\kappa_Jdt
 =\sum_Bw_B\kappa_BH_{\deg_G(B)}.
\tag{26}
\]
\(\kappa_B\) 固定于完整原包形状；不按 receiver 或时刻重复使用全部 \(W\)。

终态图无边，任一原 \(Q(x,b)\) 至多捕获一个存活包。因此每点的完整 maximal response 等于各存活包 own response 的最大值；log-positive 上界可逐包求和：
\[
 \Phi_\tau(\nu_T)\le\sum_B\Phi_\tau(\nu_B(T))
 \le\sum_B\kappa_BW_B(T).
\tag{27}
\]
在 (22) 取 \(S=T\wedge h\)，再用 (3)、(24) 的支配收敛及 (26) 的单调收敛，得到 (1)。或者未停 \(Q\) 尾以 \((Z/e)WQ(T>h)\to0\) 单独控制。全部停止与质量交换只针对该有限图政策，未推广到任意无界控制停时。

有限分割仍完整覆盖任意输入；特别可取一个大 packet \(\mu\)，而不截掉无限支撑的来源尾。因此 (1) 是对每个一般完整输入都有效的分组接口，不是仅特殊来源类的定理。一个包时无图工作量，但 \(\kappa\ge K(\mu)\) 仍是未知的空间费用；把这种选择当成已证明平方根界会循环。

## 7. 可核的合并判据：同时计算终端与工作量

在过程运行之前，可以据完整输入重复合并 packet，并每次更新其完整 own-fee 证书及图。记
\[
 J(\mathcal P,G,\kappa)=\sum_Bw_B\kappa_B(1+H_{d_B}),\qquad d_B=\deg_G(B).
\tag{28}
\]
这不是仅付删除数或细格熵的目标；它准确保留已证明的终端项和局部时钟项。

合并 \(A,B\) 为完整 \(C=A\cup B\)，\(\mu_C=\mu_A+\mu_B\)、\(w_C=w_A+w_B\)。新图取原图的 quotient：\(N(C)=(N(A)\setminus\{B\})\cup(N(B)\setminus\{A\})\)，其他边保持，去掉自环及重边。这个图仍覆盖真实 co-capture；若原图就是实际共捕获图，quotient 也是新实际图，因为正质量的 \(\mu_C(Q)>0\) 当且仅当 \(\mu_A(Q)>0\) 或 \(\mu_B(Q)>0\)。

令 \(\mathcal N=N(A)\cap N(B)\)、\(e=\mathbf1_{A\sim B}\)。则
\[
 d_C=d_A+d_B-2e-|\mathcal N|.
\tag{29}
\]
共同邻居 \(D\in\mathcal N\) 的度数减少 1，其他存留包的度数不变。取新完整包的已证帽 \(\kappa_C\ge K(\mu_A+\mu_B)\)，原包的帽保持，准确目标增量为
\[
 \begin{split}
 J'-J={}&(w_A+w_B)\kappa_C(1+H_{d_C})
 -w_A\kappa_A(1+H_{d_A})-w_B\kappa_B(1+H_{d_B})\\
 &-\sum_{D\in\mathcal N}\frac{w_D\kappa_D}{d_D}.
 \end{split}
\tag{30}
\]
共同邻居的 \(d_D\ge2\)，分母无零问题。因而合法且可直接审计的合并准则是
\[
 (w_A+w_B)\kappa_C(1+H_{d_C})\le
 w_A\kappa_A(1+H_{d_A})+w_B\kappa_B(1+H_{d_B})
 +\sum_{D\in\mathcal N}\frac{w_D\kappa_D}{d_D}.
\tag{31}
\]
每次满足 (31) 都降低或保持一个对完整原输入有效的已证上界。至多 \(N-1\) 次预运行合并，无随机 law 或停时问题。若 \(C\) 可被同一小直径轴盒容纳，(23) 给完全几何的 \(\kappa_C\)；若更大，必须另证实际 own-fee 上帽，不能沿用叶包帽或只以包质量较大作为廉价证书。

原细格图的 clique 可以合成一个有内部许多细格的大包，从而把图度数工作量清零；新的终端 own-fee 则必须留下。精细来源标签的熵在这种预合并中不变，只是由包间熵转为包内条件熵。熵链式 \(H(\text{fine})=H(\text{packet})+H(\text{fine}\mid\text{packet})\) 因此解释了为什么本终端没有声称以廉价删除消灭旧 \(n\log n\) 熵。本文只证明合并的总价目判据，不证明任意输入必有一连串满足 (31) 的合并直到某个可支付的平方根终端。

## 8. 未闭合接口与可登记的实际测量字段

新的普适弱形式是 (7)、(22)、(1) 和 (31)：实际 packet profile 付删除余费，完整来源分组付终端及局部图工作量。若能对每个完整输入构造一个有限分割及可核帽，使
\(\sum_Bw_B\kappa_B(1+H_{d_B})\le\sqrt n\,n^{o(1)}W\)，则该窗口的 log-max 目标随即成立；目前没有这个存在性定理。这个展示欠项的条件没有被单独登记成新假设引理。

后续数值若由 root 授权，可严格冻结以下对象：完整原来源与来源分割；原 \([a,b]\)；每个包 own maximal profile 及阈值费用证书；原真赢家 \(E,R,M,p_B\)；\(D_B,L_B,\Gamma_\tau(\mu_B)\) 的独立误差区间；真实共捕获或覆盖它的图边证书；每次合并的 \(d_A,d_B,e,\mathcal N,d_C\)、新完整包 own-fee 帽及 (30) 区间。不得以少量 receiver 点或有限尺度列表直接认证全空间 \(K_B\)，也不得据本锐度例声明 nearmax 资格。本文未运行任何数值、未改旧数据或旧 proof hashes。

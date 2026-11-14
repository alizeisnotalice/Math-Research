# 一般实际目录的检验、几何容量与标签费用

本稿回到任意有限真实轴向矩形目录。新增严格结果是一个有限凹优化形式的几何容量及其一般上界，并计算真实共核积星对**全部标签重加权**留下的损失。没有证明一般中央目录的所需比较；低指数率例子不被当作近中央反例。r728—r729 的单相位结果只用于检查覆盖范围，不被外推为一般目录结论。

## 1. 一般中位数结算仍需要的精确不等式

由 r709，证明 $\log H_n(1/2)\le\mu n+o(\log n)$ 即足以确定一般端点的对数系数为 $c_Q$。r718 已允许每维选到稳定门槛 $q=1/2+O(1/(\sqrt n\log n))$ 及有限目录：令 $X$ 为完整并、$m$ 为真目录最小源质量、$s$ 为负责最不利检验质量，
\[
G=|X|/m\ge H_n(1/2)/2,\quad \rho=s/m\ge c(\log n)^{-4},\quad
V=\log(|X|/s),\quad \log G=V+\log\rho.
\]
取任意负责最不利混合 $P_\pi=\sum_i\pi_i\operatorname{Unif}(R_i)$，$Q=\operatorname{Unif}(X)$。r710 给某个绝对 $C_0$ 下
$D(P_\pi\Vert Q)\le B_n:=\mu n-\frac12\log n+C_0$。定义实际熵缺额与检验差
\[
\mathcal A_\pi=B_n-D(P_\pi\Vert Q)\ge0,\qquad
\mathcal T_\pi=V-D(P_\pi\Vert Q).
\]
于是精确有
\[
\log G=\mu n-\tfrac12\log n+C_0+
\mathcal T_\pi-\mathcal A_\pi+\log\rho. \tag{1}
\]
所需联合结算是 $\mathcal T_\pi-\mathcal A_\pi+\log\rho\le\frac12\log n+o(\log n)$。单独要求 $\mathcal T_\pi\le\frac12\log n+o(\log n)$ 是较强充分条件。尤其 $\mathcal T_\pi-\mathcal A_\pi=V-B_n$ 与选择哪个最不利权无关；只把 KL 最大化本身不会改变这个联合量，仍需新的几何估计。

r719 的 NP 事件 $E$ 有 $P(E)=q$、$Q(E)=\delta=s/|X|$。其条件空间熵满足
\[
\mathcal T_\pi=(1-q)V+h_B(q)+(1-q)\log(1-\delta)
-qD(P_E\Vert Q_E)-(1-q)D(P_{E^c}\Vert Q_{E^c}). \tag{2}
\]
标签后验的 $o(1)$ TV/KL 集中没有提供 (2) 中空间条件熵的下界。这就是尚未关闭的一般接口。

## 2. 新的一般有限接口：几何容量

固定任意负责行权 $\pi$，不要求最不利。对 $f\ge0$、$Qf=1$，定义
\[
\mathscr J_\pi(f)=\sum_i\pi_i\exp\left(\fint_{R_i}\log f\right),\qquad
\Lambda(\pi)=\log\sup_{f\ge0,\ Qf=1}\mathscr J_\pi(f). \tag{3}
\]
$f$ 只需在有限成员原子上常值：对每个原子以算术平均替换 $f$，保持 $Qf$ 并由 Jensen 增大每行的平均对数。零值按几何均值的连续延拓处理。原子质量记为 $q_A=|A|/|X|>0$，行核 $K_{iA}=|A|1_{A\subset R_i}/|R_i|$；则 (3) 是紧单纯形 $\sum_Aq_Af_A=1$ 上的有限加权几何均值最大化。目标连续且凹，故上确界达到。

**普适上界。** r725 已从锐强型范数证明一般几何极大函数 $G_{
m geom}$ 满足 $\int G_{
m geom}f\le e^{\mu n}\int f$。在实现 $\pi$ 的实际负责分割 $X=\bigsqcup U_i$ 上，每个负责点的几何极大值至少为其行几何均值。因此
\[
|X|\mathscr J_\pi(f)\le\int_XG_{
m geom}(f1_X)
\le e^{\mu n}|X|Qf,\qquad \boxed{\Lambda(\pi)\le\mu n.} \tag{4}
\]
这里使用真实矩形与同一个源 $f$，无需相位表示、低秩压缩或反馈动作。负责权若以分数运输给出，可在有限成员原子内分割实现。

**精确标签变分。** 对所有 $\nu\ll\pi$ 的行概率向量，记 $P_\nu=\sum_i\nu_i\operatorname{Unif}(R_i)$。则
\[
\boxed{\displaystyle
\Lambda(\pi)=\max_{\nu\ll\pi}
\{D(P_\nu\Vert Q)-D(\nu\Vert\pi)\}.} \tag{5}
\]
证明：有限 log-sum-exp 变分给
$\log\mathscr J_\pi(f)=\max_\nu[E_{P_\nu}\log f-D(\nu\Vert\pi)]$。
再对 $f$ 取上确界；两次都是上确界，可以交换。Gibbs 变分给
$\sup_f[E_{P_\nu}\log f-\log Qf]=D(P_\nu\Vert Q)$，得到 (5)。所有原子质量正、目录有限，边界通过 $0\log0=0$ 连续延拓，故最大值达到。

与 r710 的负责熵界不同，(5) 允许 $\nu$ 本身不负责，但精确计入标签费用 $D(\nu\Vert\pi)$。特别地，一般实际目录满足
\[
D(P_\nu\Vert Q)\le\mu n+D(\nu\Vert\pi). \tag{6}
\]

## 3. 凹性与可计算证书

令 $F_\pi(\nu)$ 为 (5) 的目标。在标签—成员原子的联合律 $\nu_iK_{iA}$ 下，
\[
F_\pi(\nu)=H_\nu(I\mid A)+\sum_i\nu_i\log\frac{\pi_i|X|}{|R_i|}. \tag{7}
\]
条件熵对联合律凹，故 $F_\pi$ 对 $\nu$ 凹。亦可直接在正权内部验证：对方向 $z$，
\[
d^2F_\pi(\nu)[z,z]
=\sum_A\frac{(\sum_i z_iK_{iA})^2}{\sum_i\nu_iK_{iA}}
-\sum_i\frac{z_i^2}{\nu_i}\le0, \tag{8}
\]
由每个原子上的 Cauchy–Schwarz 再用 $\sum_AK_{iA}=1$ 得到；边界靠连续性延拓。因此这里不是把一般“凸函数最大化”误称可验证凹规划。

对 $\pi$ 的正支撑内任取严格正的 $\nu$，令 $w_\nu=dP_\nu/dQ$、
\[
g_i=\fint_{R_i}\log w_\nu-\log(\nu_i/\pi_i).
\]
一阶导数与 Euler 恒等式给 $\sum_i\nu_i g_i=F_\pi(\nu)$。凹性于是给可复核的双侧证书
\[
F_\pi(\nu)\le\Lambda(\pi)\le\max_i g_i. \tag{9}
\]
若全部 $g_i$ 相同，则它们等于全局最优值。有限原子上的对数可由区间算术核验；这提供了适用于任意实际目录的具体计算入口。式 (9) 认证的是固定 $\pi$ 的内层优化，不自动认证下文最不利面上对 $\pi$ 的外层最大化。

令 $\Pi_{\rm LF}$ 为 r719 的负责最不利紧多面体，并置
\[
\Lambda_{\rm LF}=\max_{\pi\in\Pi_{\rm LF}}\Lambda(\pi),\qquad
\Gamma=V-\Lambda_{\rm LF}.
\]
(3) 的紧集最大值表示说明 $\Lambda$ 对 $\pi$ 连续，故此最大值达到。由 (4)，
\[
\log G\le\mu n+\Gamma+\log\rho. \tag{10}
\]
所以稳定中央目录若能证明 $\Gamma+\log\rho\le o(\log n)$，便关闭一般中位数归约；在 $\rho\ge c(\log n)^{-4}$ 下，$\Gamma\le o(\log n)$ 是一个更简洁的充分条件。它不被宣称为一般结论的必要条件。新的分解为
\[
V-D(P_\pi\Vert Q)=[V-\Lambda(\pi)]+[\Lambda(\pi)-D(P_\pi\Vert Q)], \tag{11}
\]
其中第二项非负，并已通过精确几何容量计入联合结算。

## 4. 真实共核积星：全部重加权仍不能普遍消除对数损失

这一节给 (10) 新接口的真实几何检验，不把它当近中央反例。取 $h$ 个 $k$ 维块，每块核心为单位立方体，$k$ 个行盒分别沿一个正坐标轴延长 $a>0$；本块每行体积 $1+a$，每个外片体积 $a$，各外片不交。目录取全部 $h$ 块的积；剩余坐标取单位区间。记
\[
v=(1+a)^h\le2,\qquad |X|=(1+ka)^h,\qquad p=\frac a{1+a}.
\]
每个行盒体积相同为 $v$，全部包含单位核心。任何行混合的 Lebesgue 密度不超过 $1/v$，核心上恰为 $1/v$。因此 $m=s=v/2$，所有负责权都是半阈值最不利权，$V=\log(2|X|/v)$。

标签为 $T=(T_1,\ldots,T_h)$。按任意行权 $\nu$ 抽行并在该盒均匀取点，令 $J$ 为落入外片的块集合；$J$ 独立于标签，各块以概率 $p$ 被选择。与 r719 相同的分原子计算给
\[
h(P_\nu)=\log v+E_J H_\nu(T_J),
\]
故对**任意**负责参考权 $\pi$，(5) 等价于
\[
V-\Lambda(\pi)=\log2+
\inf_{\nu\ll\pi}\{E_JH_\nu(T_J)+D(\nu\Vert\pi)\}. \tag{12}
\]
实际负责容量限制每个标签投影原子：对任意固定 $J$、标签元组 $t$，
\[
\pi(T_J=t)\le\left(\frac{1+a}{1+ka}\right)^{|J|}.
\]
理由是这些标签的全部负责接收点都落在指定块的指定行盒内，其余块只能使用完整并。令 $\beta=\log[(1+ka)/(1+a)]$。KL 的投影数据处理给
\[
H_\nu(T_J)+D(\nu\Vert\pi)
\ge H_\nu(T_J)+D(\nu_J\Vert\pi_J)
=-E_\nu\log\pi_J(T_J)\ge |J|\beta.
\]
对 $J$ 平均，再在 (12) 取下确界，得到下界；取 $\nu=\pi$ 并用 $H_\pi(T_J)\le|J|\log k$ 得上界。因此对全部负责权、包括在最不利面中最大化几何容量的权，一致有
\[
\boxed{\displaystyle
\log2+hp\log\frac{1+ka}{1+a}
\le V-\Lambda(\pi)\le\log2+hp\log k.} \tag{13}
\]
这已经允许 (5) 中所有不负责的 $\nu$，不能用另一组标签重加权消除该下界。

取 $h=\lfloor\sqrt{\log n}\rfloor$、$k=\lfloor n/h\rfloor$，及有理
$a=2^{-h^2}\lfloor2^{h^2}(2^{1/h}-1)\rfloor$。则 $v\le2$、$hp\to\log2$、$\log k/\log n\to1$、$|\log a|=o(\log n)$。故
\[
V-\Lambda_{\rm LF}=(\log2+o(1))\log n. \tag{14}
\]
本例真实得分却只有 $\log G=\log2+h\log[(1+ka)/(1+a)]=O((\log n)^{3/2})=o(n)$。它严格否定无近中央条件的普遍 $o(\log n)$ 比较，未否定第 3 节针对稳定中央目录的充分接口。

## 5. 单相位结果覆盖了什么，以及本轮留下的唯一实质缺口

r725 的一般几何极大函数界已在 (4) 用于任意实际目录。r728—r729 另外证明了某些实际相位输入上，指标中位数与相应几何平均的差只有 $o(\log n)$；它们没有证明一般目录的 NP 检验能由该相位变量表示。r730 更在均值最优独立增量中证明：完整切向向量可通过 $-n\min(PZ)$ 以 $O_P(1)$ 误差恢复径向和，条件 TV 趋一。因此不能用单投影条件 TV 去控制任意 NP 事件，也不能把它当实际目录的反例。

本轮可继续的具体对象是 (3)—(9) 的有限几何容量程序：它保留同一个源、真实行核与负责容量，并有双侧可计算证书。仍缺的是稳定、近 $\mu n$ 得分实际目录上的检验—几何容量比较，例如 (10) 所述充分估计。第 4 节说明近中央条件不可删除；现有低秩相位证明没有给这一比较。本稿没有关闭一般中位数或端点系数，也没有新增一个被误称近中央的低增益反例。

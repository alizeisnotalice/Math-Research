# 中央容量权衡与可核查的实际目录约束

本稿尚未证明中央权衡，也没有构造近中央反例。以下先固化精确目标与已证一般目录引理；它们不以单相位表示为前提。r732 小目录的计算不能作为中央权衡证据，父级的实验代码不在本稿修改范围。

## 1. 精确目标与其必须包含的稳定性

沿用 r731：$q$ 是 r718 的稳定中央门槛，$X$ 为有限真实矩形目录的完整并，$m$ 为真实最小源质量，$s$ 为负责最不利检验质量，$\rho=s/m\ge c(\log n)^{-4}$，$V=\log(|X|/s)$。令
\[
\Lambda=\max_{\pi\in\Pi_{\rm LF}}\Lambda(\pi),\quad
\Gamma=V-\Lambda,\quad \Delta=\mu n-\Lambda\ge0.
\]
则严格有
\[
\log G-\mu n=\Gamma-\Delta+\log\rho. \tag{1}
\]
任何 $\Gamma\le a\Delta+o(\log n)$、$a<1$ 的统一中央估计都会闭合一般中位数归约。$\Gamma\le C\log(2+\Delta)+o(\log n)$ 也足够，因为 $C\log(2+\Delta)-\Delta$ 有统一上界。这里不声称这些较强充分估计必然成立。

对于 $\log G\ge\mu n-D_0\log n$ 的稳定目录，
\[
\Gamma\ge\Delta-D_0\log n-\log\rho. \tag{2}
\]
所以若存在这种真实目录同时有 $\Delta/\log n\to\infty$，则 $\liminf\Gamma/\Delta\ge1$，立即否定任何固定 $a<1$ 的权衡，也否定上述对数型权衡。反之，证明这些权衡首先必须排除“真实得分近中央但最不利几何容量损失远大于对数”的目录。r719 的标签后验集中与 r731 的容量上界目前均未给这项下界。

## 2. 中央检验的实际重叠与标签熵必要条件

任取负责权 $\pi$，令 $v_i=|R_i|$、$Q=\operatorname{Unif}(X)$、$w=dP_\pi/dQ$。Hall 单行约束给 $\pi_i|X|\le v_i$，故
\[
w(z)=\sum_{i:z\in R_i}\frac{\pi_i|X|}{v_i}
\le N(z):=\#\{i:z\in R_i\}\le k. \tag{3}
\]
以下 (4)—(5) 取 $\pi\in\Pi_{\rm LF}$，因此其 NP 最小体积确为第 1 节的 $s$，可沿用 $V$。若采用任意其它负责权，应在这两式中改用 $s_\pi=b_q(P_\pi)$、$V_\pi=\log(|X|/s_\pi)$。对最不利权的质量为 $q$、体积为 $s$ 的 NP 事件 $E$，$\delta=Q(E)=e^{-V}$，于是
\[
q=P(E)\le\int_E N\,dQ\le k\delta,\qquad
k\ge qe^V. \tag{4}
\]
因此维数 $n$ 中接近 $V\simeq\mu n$ 的实验必须具有指数级行数，或严格的因子化/对称压缩表示；固定几个行的小目录离中央太远。

还有一个保留实际核的更精确必要条件。按 $I\sim\pi$、$Z\mid I=i\sim\operatorname{Unif}(R_i)$ 取样，置 $\tau_i=\log[v_i/(\pi_i|X|)]\ge0$。在 $w(z)>0$ 处，标签后验为 $P(I=i\mid Z=z)=e^{-\tau_i}1_{R_i}(z)/w(z)$，所以
\[
H(I\mid Z=z)=\log w(z)+E[\tau_I\mid Z=z]\ge\log w(z).
\]
又 $D(P(\cdot\mid E)\Vert Q(\cdot\mid E))\ge0$ 给
\[
\boxed{\ H(I\mid Z,\ Z\in E)\ge E[\log w(Z)\mid E]\ge V+\log q.\ } \tag{5}
\]
这里左端表示命中事件下对空间点的条件标签熵平均。稳定中央目录必须在 NP 命中区保留约 $\mu n$ 的标签歧义；这比仅计算行数更具体。但它不把条件输出通道变回均匀矩形通道，因此不能直接下推 $\Lambda$。

## 3. 真正作用于乘积目录的容量张量恒等式

取任意两个有限实际矩形目录 $\mathcal R_A,\mathcal R_B$，并取各自负责权 $\pi_A,\pi_B$。其积权由负责分割的直积实现。记各目录的行—空间核为 $K_A,K_B$、均匀并参考律为 $Q_A,Q_B$。r731 的容量目标为
\[
F_\pi(\nu)=D(P_\nu\Vert Q)-D(\nu\Vert\pi).
\]
允许积目录的输入标签 $(I_A,I_B)\sim\nu$ **任意相关**，给定标签后两个空间样本 $Z_A,Z_B$ 独立。KL 链式法则和 Markov 关系 $Z_A-I_A-I_B$ 给精确恒等式
\[
\begin{aligned}
F_{\pi_A\otimes\pi_B}(\nu)
={}&F_{\pi_A}(\nu_A)
+E_{Z_A}F_{\pi_B}(\nu_{B\mid Z_A})
-I_\nu(I_A;I_B\mid Z_A). \tag{6}
\end{aligned}
\]
事实上输出 KL 按 $Z_A$ 分解；输入 KL 按 $I_A$ 分解。两个条件标签 KL 的差等于
$H(I_B\mid Z_A)-H(I_B\mid I_A)=I(I_A;I_B\mid Z_A)$，就是 (6) 的最后一项。零权按有限空间条件律惯例处理。

因此
\[
\boxed{\Lambda(\pi_A\otimes\pi_B)=\Lambda(\pi_A)+\Lambda(\pi_B).} \tag{7}
\]
上界来自 (6)，下界由两个内层最优标签律的独立积达到。更精确地，
\[
\begin{aligned}
\Lambda_A+\Lambda_B-F(\nu)
={}&[\Lambda_A-F_A(\nu_A)]
+E_{Z_A}[\Lambda_B-F_B(\nu_{B\mid Z_A})]
+I(I_A;I_B\mid Z_A). \tag{8}
\end{aligned}
\]
三项皆非负。故任何距离积容量至多 $\epsilon$ 的相关标签倾斜，都必须满足 $I(I_A;I_B\mid Z_A)\le\epsilon$；交换两块也有 $I(I_A;I_B\mid Z_B)\le\epsilon$。这给出了相关标签试图复用两块空间信息时的精确损失，而不只是在独立输入上检查乘积公式。

这条定理对实际积矩形目录成立，但**没有**证明 $\pi_A\otimes\pi_B$ 在相同半门槛上仍为积目录的负责最不利权。也没有证明真实半质量 $m_{A\times B}(1/2)$ 等于 $m_A(1/2)m_B(1/2)$；r711 已说明共同源约束不能如此相乘。因此 (7) 本身不提供中央 $\Gamma$ 的张量反例或权衡证明。

已用本机 Wolfram 在真实区间目录 $[0,3/2],[-1/2,1]$ 的二维积上独立核验 (6)：负责参考权各为 $(1/2,1/2)$，相关输入 $\nu=\left(\begin{smallmatrix}3/8&1/8\\1/8&3/8\end{smallmatrix}\right)$。有限原子核全为有理数，链式恒等式的符号残差精确为零。相关费用约为 $0.0872080$，总容量亏损约为 $0.1162774$。源码和实际输出见 `tensor_identity_check.wl`、`tensor_identity_check_output.txt`；小维检查只核验恒等式，不检验中央权衡。

## 4. 当前可执行核查与明确未闭合点

最直接的新核查是 (8)：在任何父级已有的两个真实小目录上，固定已认证的负责参考权，独立计算积目录容量（允许全部相关 $\nu$），其值必须等于两个单块容量之和；同一个相关候选的目标亏损还必须逐项等于 (8)。这可以用真实行核验证相关项与链式分解，但小目录的此项通过仍不是中央权衡证据。若以后使用因子化大目录寻找近中央反例，(8) 可先排除仅靠遗漏相关标签输入而产生的虚假容量降低。

当前真正缺口仍是 (2) 描述的容量稳定性，或直接证明 (1) 的联合补偿。现有共核积星真实例只给 $\Gamma\asymp\log n$ 而 $\Delta\asymp n$，与两种拟议权衡都相容；重复添星块不构成新的中央机制。本稿尚无满足近中央条件、同时 $\Delta\gg\log n$ 的实际目录。

另一个待检验机制是 r730 的高秩切向统计量。它能在单个固定相位处编码径向和，却不能省去接收相位的 Haar 平均和同一轮廓的源质量。非周期的未截断函数 $h(v)=-n\min_i v_i-b$ 精确把源条件化为 $n\min_i u_i>T-b$，即一个小立方体，而不是一般近中央源。下面给出这个诊断的真实全 halo 费用。

## 5. 未周期化的径向编码精确退化为低增益立方体

令 $P=I-\boldsymbol1\boldsymbol1^T/n$、$\ell=\sum_i u_i$。若直接使用 $h(Pu)=-n\min_i(Pu)_i-b$，则
\[
\ell>T+h(Pu)
\iff\ell>T-n\min_i u_i+\ell-b
\iff\min_i u_i>(T-b)/n.
\]
当 $T>b$，真实源就是 $[-\epsilon,\epsilon]^n$，$\epsilon=e^{-(T-b)/n}$。这个例子不仅对指定角点族，而且对完整强矩形极大函数都只有次指数的半阈值增益。

事实上对一维区间 $I=[-\epsilon,\epsilon]$，$M_1 1_I(x)=1$ 于 $I$ 内，在外部等于 $2\epsilon/(|x|+\epsilon)$。乘积指标满足 $M_S1_{I^n}(x)=\prod_iM_1 1_I(x_i)$。令 $z_i=-\log M_1 1_I(x_i)$；以 $|I|$ 归一化的一维体积测度推前为 $\delta_0+2e^z1_{z>0}dz$。因此，写 $L=\log2$，半阈值 halo 的真实体积比精确为
\[
G_{\rm cube}(n)=1+\sum_{j=1}^n\binom nj\frac{2^j}{(j-1)!}\int_0^L e^s s^{j-1}\,ds.
\]
利用 $e^s\le2$、$\binom nj\le n^j/j!$，
\[
G_{\rm cube}(n)\le2\sum_{j\ge0}\frac{(2n\log2)^j}{(j!)^2}
\le2\exp(2\sqrt{2n\log2}). \tag{9}
\]
最后一步可逐项比较 $\binom{2j}{j}\le4^j$ 与指数函数的偶数项。故 $\log G_{\rm cube}(n)=O(\sqrt n)=o(n)$。径向编码本身并未产生近中央实际源。

周期化或截断该高秩轮廓以后，上述全局代数恒等式不再成立，(9) 不能替代它们的计算。对周期化候选，须先核对真实相位是否确实趋于所选环面上的 Haar 律，再计算完整接收相位平均及同一轮廓的源预算；不能只在相位零处利用编码。对纯截断但仍非周期的候选，没有自动的 Haar 平均，应直接使用真实深度层上的角分布并计算源与接收体积。当前没有执行足以支持该类新下界的实际目录实验；这一候选与中央容量权衡都保留为未闭合工作。

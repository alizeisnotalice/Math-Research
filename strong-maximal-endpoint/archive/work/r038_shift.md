# r038：从 r021 CDF 屏障到 shifted late barrier 的审计

日期：2026-11-06。只作数学审计，不修改仓库。结论分三层：

1. **不能**把 r021 单独升级为无误差的、全实轴
   \(\mathcal U_{T,N}(r)\le v_{T,N}^a(r-H)\) 且 \(H=O(\log T)\)；
   最右尾有明确指数反例。
2. r037 已有的全局 \(O_T\) 与局部 \(I_T\) 的 min-splice 可以避开该反例。
   只需在较宽自然管中证明 shifted-reference 比较；管外由 \(O_T\) 的单调性给
   shifted \(W_T\) 的上界。
3. 对固定 \(\gamma,A\)，这一初始化很可能可由 r021 + 一致倾斜 local CLT 完成，
   但其位移常数目前必然依赖于外带指数/管宽。因晚时段拼接需要
   \(A^2\gtrsim1/\gamma\)，现有 r021 的常数给出至少 \(H\gtrsim A^2\log T\)。
   所以它是严格的 \(O(f_N\log T)\) 初始化候选，**不是**令
   \(\gamma\downarrow0\) 后自动消失的 sharp 系数论证。

## 1. 两个必要检查

### 左无穷端的均值条件

对 \(r\le0\)，终端全在线性支，故
\[
 \mathcal U_{T,N}(r)=1-r/N+T\mu/N,
 \qquad
 v_{T,N}^a(r-H)=1-r/N+(H+Tm(a))/N.
\]
因此任何全局 shifted inequality（加法误差 \(\epsilon\)）至少要求
\[
 H\ge T[\mu-m(a)]-N\epsilon=T\delta-N\epsilon.\tag{1}
\]
这说明匹配动作不能只写 \(H=C\log T\)，除非已验证 \(T\delta\lesssim\log T\)。
所提预算 \(H=O(\log T)+O(T|a-\mu|)\) 覆盖 (1)，因为在 \(\mu\) 的紧邻域
\(\delta=O((a-\mu)^2)\le C|a-\mu|\)。

### 右无穷端的精确反例（无外侧补丁）

令 \(p_c=(1+e^c)^{-1}\)。对一个固定动作，Gamma--Bernoulli 和的最右端
\(e^{-r}r^T/T!\) 系数由
\[
 A(c)=E e^{c\operatorname{Bern}(p_c)}
      =\frac{2e^c}{1+e^c}\tag{2}
\]
的 \(T\) 次幂给出；有限 \(N\) 的线性尾部只产生低一阶 \(r^{T-1}e^{-r}\)
项，故同一首项仍在。取常值 \(c\to\infty\)，可使 \(A(c)\to2\)。在
\(a\to\mu\) 时，利用 \(e^{-\mu}=\mu-1\)，
\[
 A(\mu)=2/\mu.
\]
因而若不加误差/外侧屏障，\(K_c^T\phi_N(r)\le v_{T,N}^a(r-H)\) 在
\(r\to\infty\) 必须有
\[
 e^H A(a)^T\gtrsim 2^T,
 \quad\hbox{即}\quad H\ge T\log\mu+o(T).\tag{3}
\]
这否定了全局 \(O(\log T)\) shift 的说法。Wolfram 15 实际符号计算给
\(m(\mu)=v(\mu)=\mu\), \(v'(\mu)=\mu-1\), \(A(c)=2e^c/(1+e^c)\)，
并数值给 \(\mu=1.2784645427\ldots\), \(A(\mu)=1.5643765885\ldots\),
\(2/A(\mu)=\mu\).

若容许 \(\epsilon_T=T^{-P}\)，极右端最终会小于 \(\epsilon_T\)，但这不是从
普通 CLT/CDF 误差自动得到的全局比较；仍须用正 mgf 外尾界确定进入该区间的位置。

## 2. 正确的 shifted-W 初始化结构

令 r037 的 \(W_T=\min(O_T,I_T)\)。两条全局 \(O\) 屏障都支配真值，因而
\(O_T\ge\mathcal U_{T,N}\ge v_{T,N}^a\)，而 \(I_T\ge v_{T,N}^a\)，所以
\[
 W_T\ge v_{T,N}^a\quad\hbox{全局}.\tag{4}
\]
取 splice 半径 \(A\)，在更宽的 \(2A\) 管内证明
\[
 \mathcal U_{T,N}(r)\le v_{T,N}^a(r-H)\qquad
 |r-Tm(a)|\le2A\sqrt{T\log T}.\tag{5}
\]
若 \(H=o(\sqrt{T\log T})\)，则管外 \(|r-Tm(a)|\ge2A\sqrt{T\log T}\)
蕴含 \(r-H\) 位于 \(A\) 外部。r037 的外部吸收给
\(W_T(r-H)=O_T(r-H)\)，且 \(O_T\) 递减，故
\[
 W_T(r-H)=O_T(r-H)\ge O_T(r)\ge\mathcal U_{T,N}(r).\tag{6}
\]
管内由 (4)--(5) 同样成立。于是能得到真正全空间的初始化
\[
 \mathcal U_{T,N}(r)\le W_T(r-H).\tag{7}
\]
这回避了 (3)，也不要求 shifted \(v^a\) 在外侧支配任何控制。卷积核平移不变，
所以 (7) 可随晚时段 supersolution 传播。这是逻辑上正确的 splice；关键新工作仅为
(5) 的有限终端证明。

## 3. 管内比较需要的、尚未写成定理的输入

r021 给
\(F_{S_T^\pi}\ge F_{G_T}-\epsilon_T\)，其中
\[
 G_T\sim N(T\mu+H_0,\mu(T+t_0)),\qquad H_0,t_0=O(\log T).\tag{8}
\]
对有界部分，加 Exp 卷积即可把它变成 Gaussian payoff 上界。有限终端的 call
须先将 CDF 误差只积到一个 \(O(T)\) cutoff，再用正 mgf 收尾；中心区域中该误差为
\(O(T\epsilon_T/N)+O(\epsilon_T)\)。

接着需一条**一致倾斜、相对型** local CLT/Edgeworth 比较，而不只是普通 Berry--Esseen：
在 \(|r-Tm(a)|\le2A\sqrt{T\log T}\) 内，将 (8) 的 Gaussian 和
\(S_T^a\) 的共同倾斜密度及截断 payoff 逐项比较。普通绝对 CDF 误差
\(O(T^{-1/2})\) 在管边的密度 \(T^{-1/2-C A^2}\) 处不够；必须是倾斜后的相对误差。
正态均值和方差差可由
\[
 H\ \ge\ H_0+T\delta+
 C A\sqrt{\log T}\,\bigl|\sqrt{\mu(T+t_0)}-\sqrt{T\,v(a)}\bigr|
 +C_A\log T\tag{9}
\]
吸收。利用 \(v'(\mu)=\mu-1\)，(9) 确为
\(O_A(\log T)+O(T|a-\mu|)\)，并且在 r037 的匹配尺度、\(T=N^\gamma\)
下 \(H=o(\sqrt{T\log T})\)。左侧有限终端的仿射部分正由 \(T\delta\) 补偿，
而不是被 CDF 误差忽略。

若 r021 的外带误差可调为 \(T^{-P}\)，为在管边把它吸收进 shift gain，需取
\(P\gtrsim C A^2\)。这点也必须在证明中定量写出。

## 4. 常数/锐性障碍

r021 现有构造取外带参数 \(A_0\) 时，其 CDF 误差由
\(T^{-A_0^2/8}\)（右带）控制。因此达到 \(T^{-P}\) 至少需
\(A_0^2\ge8P\)。其状态移位的截断 \(\psi\) 上界为常数乘 \(A_0^2\log T\)，
故 r021 所给的 \(H_0\) 常数至少随 \(P\)（从而随 \(A^2\)）增长。

另一方面 r037 的外侧吸收要求 \(K>2/\gamma-1\)，并相应需要
\(A^2\gtrsim K\gtrsim1/\gamma\)。因此用当前 r021 直接供给 (5) 时通常有
\[
 H\gtrsim \frac{\log T}{\gamma},\qquad
 f_N H\gtrsim f_N\log N
\]
（\(\log T=\gamma\log N\) 后常数不趋零）。

所以：对每个固定 \(\gamma\) 这仍给 \(O(f_N\log T)\) 的合法 shifted initialization，
但不应宣称让 \(\gamma\to0\) 即把它变成低于已知 \(g\log N\) 的 sharp 余项。
要消除该障碍，需要一个在宽管内不带 \(A^2\) 级统一漂移的 state-dependent/Doob
比较，或一种不要求把 CDF 误差压至管边相对密度以下的 payoff 直接估计。

## 5. 独立核验：面积中的固定动作比较可改善到 \(\gamma<1\)

这不产生点态 shifted dominance，却改善 r037 的面积账本。令
\(G_\mu\sim N(T\mu,T\mu)\)、
\(G_a\sim N(Tm(a),T\operatorname{Var}Z_a)\)。对 \(a\) 位于 \(\mu\) 的
固定紧邻域，iid Wasserstein--1 CLT 的未标准化版本一致给
\[
 W_1(S_T^\mu,G_\mu)+W_1(S_T^a,G_a)\le C.\tag{10}
\]
以同一标准正态耦合两个 Gaussian，且 \(\operatorname{Var}Z_a\) 在邻域内
Lipschitz，便有
\[
 W_1(G_\mu,G_a)
 \le T\delta+C\sqrt T\,|a-\mu|.\tag{11}
\]
故三角不等式给严格、且比逐步耦合强的
\[
 W_1(S_T^\mu,S_T^a)
 \le C+T\delta+C\sqrt T\,|a-\mu|.\tag{12}
\]
这不使用不同动作的 FOSD。将 r037 的同均值 call 论证原样保留，补偿参考的
有限终端正面积可写为
\[
 C\log T+C+T\delta+C\sqrt T|a-\mu|+CT/N.\tag{13}
\]
对匹配尺度 \(\delta\asymp\log N/N\)、
\(|a-\mu|\asymp\sqrt{\log N/N}\)、\(T=N^\gamma\)，
后三个非对数项均为 \(o(1)\) 于每个固定 \(0<\gamma<1\)。Wolfram 15 对
\(N^{\gamma-1}\log N\)、\(N^{(\gamma-1)/2}\sqrt{\log N}\)、\(N^{\gamma-1}\)
在 \(0<\gamma<1\) 下的极限均输出 0。此改进仅扩展面积初始化适用的 \(\gamma\)
范围；不改变第 4 节的点态 shift 常数/锐性障碍。

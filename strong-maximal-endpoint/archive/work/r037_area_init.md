# r037：r021 全局 Gaussian 屏障可给出的面积初始化审计

日期：2026-11-06。审计对象为 r021 的式 (15)，不修改仓库。记
\(\mu=\log B\)，\(\phi_\infty(r)=\min(1,e^{-r})\)，
\[
 U_T^\infty(r)=\sup_\pi E_\pi\phi_\infty(r-S_T),\qquad
 v_T^\mu(r)=E\phi_\infty(r-S_T^\mu),
\]
其中参考策略每一步都取 \(\mu\)。结论：**是，r021 已证明的全状态
CDF 屏障，加一个标准的 iid Wasserstein--1 CLT，即严格推出**
\[
 \int_{\mathbb R}(U_T^\infty-v_T^\mu)_+\,dr\le C\log T.
 \tag{A}
\]
常数与 \(T\) 无关。这里并没有、也不需要声称所有策略的随机平移支配。

## 1. r021 实际提供的输入

r021 (15) 对每一可预测策略、所有 \(q\in\mathbb R\) 给出
\[
 F_{S_T^\pi}(q)\ge F_{G_T}(q)-\varepsilon_T,
 \quad \varepsilon_T=2T^{-4},
 \tag{1}
\]
其中
\[
 G_T\sim N(T\mu+H_T,\,\mu(T+t_0)),\quad
 0\le H_T\le C_H\log T,\quad t_0=C_0\log T.
 \tag{2}
\]
这正是所需的全实轴（不是仅中央管）输入。屏障中 \(H_T\) 的大小为
\(O(\log T)\)，而误差为多项式小；后二者都须保留。

令 \(Y\sim\operatorname{Exp}(1)\) 独立。逐点恒等式
\[
 \phi_\infty(r-s)=P\{Y\ge r-s\}=P\{s+Y\ge r\mid s\}
 \tag{3}
\]
将 (1) 与 \(Y\) 卷积，得到对每个 \(\pi\)
\[
 P_\pi\{S_T+Y\ge r\}\le P\{G_T+Y\ge r\}+\varepsilon_T.
 \tag{4}
\]
右式不依赖策略，故对上确界同样成立。特别地，写
\(w_T(r)=P\{G_T+Y\ge r\}\)，则 \(U_T^\infty(r)\le w_T(r)+\varepsilon_T\)。

## 2. 为什么不能把 \(\varepsilon_T\) 直接积分，以及正确截断

这一步是 r021 本身没有写出的桥梁。常数 \(\varepsilon_T\) 不能在
\(r\in[0,\infty)\) 上直接积分。取 \(\theta=1/2\)。由每步的条件正矩母函数
界（r033/r034 的 \(M_*(\theta)\)），对任何预测策略
\[
 E_\pi e^{\theta S_T}\le M_*(\theta)^T,
 \qquad
 U_T^\infty(r)\le 2M_*(\theta)^Te^{-\theta r}.
 \tag{5}
\]
这里 \(Ee^{\theta Y}=2\)。令
\[
 R_T=\theta^{-1}\{T\log M_*(\theta)+5\log T\}=O(T).
 \tag{6}
\]
则 \(\int_{R_T}^\infty U_T^\infty(r)dr=O(T^{-5})\)，而
\(\varepsilon_T R_T=O(T^{-3})\)。又 \(S_T,Y\ge0\)，所以当 \(r<0\) 时
\(U_T^\infty(r)=v_T^\mu(r)=1\)。于是 (4)--(6) 给出
\[
 \int (U_T^\infty-v_T^\mu)_+
 \le \int_{\mathbb R}(w_T-v_T^\mu)_+ +O(T^{-3}).
 \tag{7}
\]

## 3. Gaussian 与 iid 参考的面积距离

令 \(N_T\sim N(T\mu,\mu T)\)。一维恒等式
\[
 \int_{\mathbb R}|P\{X\ge r\}-P\{X'\ge r\}|dr=W_1(X,X')
 \tag{8}
\]
以及用同一个 \(Y\) 的耦合，给
\[
 \int(w_T-v_T^\mu)_+\le W_1(G_T+Y,S_T^\mu+Y)
 \le W_1(G_T,S_T^\mu).
 \tag{9}
\]
对右端作三角不等式。用共同标准正态耦合，
\[
 W_1(G_T,N_T)\le H_T+
 \sqrt{2\mu/\pi}\,[\sqrt{T+t_0}-\sqrt T]=O(\log T).
 \tag{10}
\]
而 \(S_T^\mu\) 是 iid、均值 \(\mu\)、方差 \(\mu\)、三阶绝对中心矩有限的和。
标准 Wasserstein--1 Berry--Esseen/Stein 界给
\[
 W_1(S_T^\mu,N_T)\le C\frac{E|Z_\mu-\mu|^3}{\mu}=O(1).
 \tag{11}
\]
（等价地，先把和标准化，\(d_{W_1}\le C\rho_3/(\mu^{3/2}\sqrt T)\)，再乘
\(\sqrt{\mu T}\)。）(7)--(11) 即证明 (A)。

这也说明该路线并未与“真实总和不可能有 \(O(\log T)\) 的全局随机平移支配”冲突：
(A) 只控制按阈值分别取上确界后的尾函数面积；它由 r021 的 CDF 下界、尾截断和
参考分布的 \(W_1\) 接近性组成，完全没有构造同一策略的路径平移耦合。

## 4. 本机 Wolfram 检查（仅代数 sanity check，不承担证明）

实际运行 WolframScript 15.0。它以符号积分核验
\[
 \int_{-\infty}^\infty\left|\Phi((x-h)/s)-\Phi(x/s)\right|dx=h
 \quad(h\ge0,s>0),
\]
并以 50 位数值在 \((h,s)=(3,2)\) 得 \(3.000\ldots\)。这与 (8)、(10) 的
平移项一致；上述严格论证仍只依赖可检验的一维 transport 恒等式和标准 CLT 界。

## 5. 有限 \(N\ge T\) 的说明：可行，但不是“直接从 (A) 自动推出”

精确分解为
\[
 \phi_N(r-s)=\phi_\infty(r-s)+\frac{(s-r)_+}{N}.
 \tag{12}
\]
因此，若写 \(C_X(r)=E(X-r)_+\)，有限终端还包含 call 项
\(N^{-1}C_{S_T^\pi}(r)\)。把 r021 的 (1) 从 \(r\) 积到 \(R_T\)，再用 (5)
收尾，可在 \(0\le r\le R_T\) 得
\[
 C_{S_T^\pi}(r)\le C_{G_T}(r)+\varepsilon_T R_T+O(T^{-5}).
 \tag{13}
\]
故屏障误差经 \(1/N\) 缩放仍为可忽略项。

但是，要将 Gaussian call 与 iid-reference call 的**面积**压到 \(O(\log T)\)，需要
一个二阶 transport / 非一致局部 CLT 输入；仅有 r021 的 CDF 屏障和 (11) 的 \(W_1\)
并不自动给出这一点。一个充分的补充引理是
\[
 W_2(S_T^\mu,N_T)=O(1).
 \tag{14}
\]
它对本项目的光滑、指数矩 iid 增量应可由一维 Edgeworth--quantile 耦合证明；在 (14) 下，
\[
 \int_0^\infty|C_X(r)-C_{X'}(r)|dr
 \le\tfrac12\|X-X'\|_2(\|X_+\|_2+\|X'_+\|_2),
 \tag{15}
\]
再加 Gaussian 平移 \(H_T\) 后、除以 \(N\ge T\)，为 \(O(\log T)\)。

同样地，所谓“仿射相减”必须保留动作代价：若
\(L_t(r)=1-r/N+t\mu/N\)，则
\[
 U_{t,N}-L_t=\sup_c\{K_c(U_{t-1,N}-L_{t-1})-[\mu-m(c)]/N\}.
 \tag{16}
\]
所以不能把无限终端的面积结论逐字移到有限 \(N\) 而忽略这项。结论是：
**无限终端的 (A) 已闭合；有限 \(N\) 版本有清晰的 (14) 这一额外、尚未从 r021/r020
现有文本中逐项证明的 bridge。**

## 6. 补充审计：有限 \(N\) 实际不必用 \(W_2\)

以上第 5 节的 \(W_2\) 是充分但并非必要的路线。下面的论证是更好的有限
\(N\ge T\) 版本；已逐项核对通过。

对任意预测策略，令
\[
 M_T=S_T-\sum_{i\le T}m(c_i),\qquad m=\mu T.
\]
则 \(M_T\) 是中心鞅，且逐步条件指数矩在原点的固定邻域一致有界（直接由
\(Z_c\) 的统一 \(\operatorname{Exp}\) 尾/ \(M_*(\theta)\) 界）。故条件
Bernstein 估计给
\[
 P_\pi\{|M_T|>x\}\le C\exp\{-c\min(x^2/T,x)\}. \tag{17}
\]
这一步允许完整的过去依赖；不能误写成独立和。又路径上
\(S_T\le m+M_T\)，因为 \(\sum_i m(c_i)\le \mu T\)。写
\[
 H_T(s)=C\int_{|s|}^\infty\exp\{-c\min(x^2/T,x)\}\,dx . \tag{18}
\]
对 call \(C_X(r)=E(X-r)_+\)，按 \(s=r-m\) 的正、负两种符号以及
\(EM_T=0\) 分别展开 \((M_T-s)_+\)，严格有
\[
 \sup_\pi C_{S_T^\pi}(r)\le(m-r)_++H_T(r-m). \tag{19}
\]
例如 \(s<0\) 时，左边减去 \(-s\) 恰为 \(E(s-M_T)_+\)，由左尾积分控制；
\(s\ge0\) 同理由右尾控制。参考 call 由凸性/Jensen 满足
\(C_{S_T^\mu}(r)\ge(m-r)_+\)。并且 Fubini 给
\[
 \int_{\mathbb R}H_T(r-m)dr
 \le C\int_0^\infty x e^{-c\min(x^2/T,x)}dx\le CT. \tag{20}
\]
因此
\[
 \int\bigl(\sup_\pi C_{S_T^\pi}-C_{S_T^\mu}\bigr)_+dr\le CT. \tag{21}
\]

将精确分解 (12) 代回；对同一策略的两个非负项取上确界后以
\(\sup(A+B)\le\sup A+\sup B\)，便得
\[
 \int\bigl(U_{T,N}-v_{T,N}^\mu\bigr)_+dr
 \le C\log T+C\,T/N. \tag{22}
\]
故只要 \(N\ge T\)，有限终端也有 \(O(\log T)\) 初始化面积。这是
initial-area module，不是整个 sharp MAX transport 或一般 \(\mathcal R_N\) 上界。

### 与固定匹配动作 \(a\) 的补偿

令 \(\delta=\mu-m(a)\ge0\)，并以
\(\widehat v_{T,N}^a=v_{T,N}^a+T\delta/N\) 作比较。对 \(a\) 位于
\(\mu\) 的固定紧邻域，用同一个 Exp 和同一个 uniform Bernoulli 变量逐步耦合，
\[
 W_1(S_T^\mu,S_T^a)\le CT|a-\mu|. \tag{23}
\]
这里虽不同 \(a\) 间没有 FOSD，\(W_1\) 耦合完全合法；每步的
\(c1_{\{U\le p_c\}}\) 在 \(L^1\) 中局部 Lipschitz 即可。因 \(\phi_\infty\)
表示为加 Exp 后的 survival，(23) 控制两固定参考无限终端的 \(L^1\) 差。

对 call，逐点 \(C_{S_T^a}(r)+T\delta\ge C_{S_T^a+T\delta}(r)\)，且后者与
\(S_T^\mu\) 同均值 \(m\)。任意同均值、有限方差的 \(X,X'\) 满足
\[
 \int_{\mathbb R}|C_X(r)-C_{X'}(r)|dr
 \le\frac{\operatorname{Var}X+\operatorname{Var}X'}2. \tag{24}
\]
理由是 \(C_X\ge C_m\)（call 的凸性）且
\(\int(C_X-C_m)=\operatorname{Var}X/2\)。两固定 iid 和方差均为 \(O(T)\)。
把 (21)、(23)、(24) 用三角不等式合并，得到
\[
 \int\bigl(U_{T,N}-\widehat v_{T,N}^a\bigr)_+dr
 \le C\log T+CT|a-\mu|+CT/N. \tag{25}
\]
这验证了所提均值补偿的方向和大小。若 \(T=N^\gamma\), \(\gamma<1/2\)，且
\(|a-\mu|\asymp\sqrt{\log N}/N\)，后两项为 \(o(1)\)。

本机 WolframScript 15.0 额外以对称 \(\pm1\) 分布符号核验
\(\int[C_X(r)-C_{EX}(r)]dr=\operatorname{Var}(X)/2\)，输出精确 \(1/2\)，
并以 50 位积分复核；这只是 (24) 的 sanity check，(24) 的证明如上。

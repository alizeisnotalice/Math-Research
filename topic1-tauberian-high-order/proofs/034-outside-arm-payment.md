# R034　三临界长尾在二轴源层之外的支付

日期 2026-11-05。固定 \(n\ge3,\ a=v^n\)。沿 R031 已审引理的六阶饱和序列，取同一末层实际源、目标和中心方体。以下 \(S=\sum_i r_i^n\)。采用 R031 的两项已证输入：

1. 若深处边恰有三个临界轴，所有者取 \(R_i\ge R_j\) 端，并按 \(p_1\le p_2\le p_3\) 写三个交长除以 \(r_i\)，则对任意固定 \(T<\infty\)，
   \(\sum_{p_2\le Tv^3} I_e=o_T(v^7S)\)，其中 \(I_e=|X_i\cap X_j|\)。
2. R031 在 R023 原有二轴缺源层之外保留的物理成本 \(Z_\sigma=(v/a)\sum_i|[(x_i+W_i')\setminus F_i]\setminus Y_{i,\sigma}|\) 满足
   \(\limsup_{v\to0}Z_\sigma/(v^6S)\le C\sigma\)。每个固定 \(\sigma>0\) 的二轴层 \(Y_{i,\sigma}\) 位于源有效端盒内，所用两坐标到对应面的深度均不超过 \(K_\sigma v^3 r_i\)，\(K_\sigma<\infty\) 先于 \(v\) 选定。

令 \(d_{i,l}(x)\) 为 \(x\) 到第 \(l\) 轴两张有效盒面的较小距离。对任意固定 \(K>0\)，取整个二轴角区
\[
 \mathcal C_i(K)=
 \{x\in x_i+W_i':\text{至少两个不同坐标 }l\text{ 满足 }d_{i,l}(x)\le Kv^3r_i\},
\qquad
 Z(K)=\frac va\sum_i|[(x_i+W_i')\setminus F_i]\setminus\mathcal C_i(K)|.
\]
零测边界任取。用 \(\mathscr E_{\ge3}\) 表示全部外部正交对中临界轴不少于三个的边，各无向对仅计一次。

**定理。** 对任意固定 \(K>0\)，
\[
 \boxed{\displaystyle
  \sum_{e\in\mathscr E_{\ge3}} I_e
       \le C_{n,A}v^{1/4} Z(K)+o_K(v^7S).}                              \tag{1}
\]
常数 \(A\) 是原末层固定定位常数。特别地，\(\mathcal C_i(K)\) 内任意合法的二轴实际缺源层仍可按完整鲜导数定价。另对任何固定 \(\tau>0\)，有无需扣除 \(Z(K)\) 的原交量结论
\[
 \boxed{\displaystyle
 \sum_{\substack{q(e)=3,\ \mathrm{deep}\\p_2p_3\ge\tau v^4}} I_e
          =o_\tau(v^7S).}                                               \tag{2}
\]
例如 \(p_2\ge c v^2\) 的所有三临界边落入 (2)。本稿不使用二维的角形刚性、七阶全局结算或一般维数的精确接触计数。

## 1　来自原见证的实际源矩形

对深处三临界边，以 \(R_i\ge R_j\) 的一端为源。反射并重排坐标使三个有向面对端同序。由两中点面对端距为 \((R_i+R_j)/2-p_l r_i\)，相应原始源、目标面对端距至多 \(R_i-p_l r_i\)。所有 \(0<p_l\le C_Av\)，两端尺度比为 \(1+O_A(v)\)。其余各轴上两个完整原始盒的任何坐标差均小于 \(R_i\)。记
\[
 A_e=p_1p_2p_3,\quad H_e=\min(p_1,A_e/v^3),\quad
 \delta_i=\Delta_i^{\rm end}/r_i,\quad
 \alpha_i=\max_{l,\pm}g_{i\to p(l,\pm)}/r_i .
\]
原估计给 \(\sum\varepsilon_i r_i^n=O(v^5S)\)、
\(\sum_{\rm deep}(\alpha_i+\delta_i)r_i^n=O(v^6S)\)。
非临界方向交长至多 \(C_Ar_i\)，所以 \(I_e\le C_AA_e r_i^n\)。

先去掉两种例外。若 \(\varepsilon_j\ge c_0 A_e/v^3\)，则
\(I_e\le C v^3\varepsilon_j r_i^n\)。若
\(\alpha_i+\delta_i\ge c_1H_e\)，则
\[
 I_e\le C v^2(\alpha_i+\delta_i)r_i^n ;
\]
因为 \(H_e=p_1\) 时 \(A_e/H_e=p_2p_3\le C_A^2v^2\)，另一分支 \(A_e/H_e=v^3\)。固定小常数 \(c_0,c_1\) 的选择见下段。有界相交度保证全部这些边的交损是 \(O(v^8 S)\)。

其余每边在所有者有效盒内取如下位置。三个有向有效面深分别在
\[
 [c_3H_e,2c_3H_e]r_i,\quad[c_2p_2,2c_2p_2]r_i,\quad
 [c_2p_3,2c_2p_3]r_i .                                                   \tag{3}
\]
对每个非临界坐标，只取到其两面均 \(>Kv^3r_i\) 的中心区。各原、有效盒宽介于 \(c_wvr_i,C_wvr_i\)。先取固定 \(c_2>0\) 足够小，使 (3) 长两轴的深度上限不超过该边宽四分之一。再取 \(c_1,c_3\) 小，使短轴深度也小于其边宽四分之一。

若 (3) 中某点属于 \(F_i\)，其同一个原指定严格见证已经包含全部 \(G_i\)，且 \(r(x)>R_i\)。从原始面对端算起，实际源位置深度比 (3) 至多多出 \(\delta_i r_i\)。由 \(\delta_i<c_1 H_e\le c_1p_1\)，每条临界轴上保证进入目标原盒的深度至少为固定倍数的 \(p_lr_i\)；截于目标真实原盒宽度后仍成立，因为 \(p_l\le C_Av\)。非临界方向的目标全宽进入。扣去目标原盒实际缺量 \(C_A\varepsilon_j a r_i^n\) 后，可保证的额外 \(F_j\) 质量至少为
\[
 c\, a r_i^n\,A_e/v^3 .
\]
这里用 \(F_j\cap G_i=\varnothing\)，并取 \(c_0\) 小于原盒穿透常数的固定倍数。另一方面，第 1 轴方向的真正面邻居给
\[
 (1-C_Av)(r(x)-R_i)
 \le (2c_3H_e+\delta_i+\alpha_i)r_i .
\]
所以见证允许的额外质量至多 \(C_A a r_i^n(2c_3+c_1)H_e\)。由 \(H_e\le A_e/v^3\) 知二者矛盾。常数可先取 \(c_2\)，用足够小的 \(c_1,c_3\) 统一保证各轴至少固定比例穿透，再缩小它们以满足容量矛盾，最后取 \(c_0\)。这几项均独立于 \(K,\sigma,\tau,v\)。于是 (3) 整块确实位于 \((x_i+W_i')\setminus F_i\)。

## 2　与所有固定二轴源层物理分离

固定 \(K\)，取 \(T>2K/c_2+1\)。只考虑 \(p_2>Tv^3\)。此时 (3) 中第 2、3 轴的有向面深均 \(>Kv^3r_i\)，相反面深至少是该轴全宽的四分之三。非临界轴也在各自中心区。小 \(v\) 下至多第 1 轴可以接近某个有效面。因此 (3) 与**整个** \(\mathcal C_i(K)\) 不交。其非临界中心区各保留至少一半原有效宽度，一旦 \(K\) 固定，这个固定比例不依赖 \(K\)。设 \(\mathcal R_e\) 为 (3) 的矩形，则
\[
 \frac va|\mathcal R_e|
 \ \ge\ c_A r_i^n\,\frac{p_2p_3H_e}{v^2}
 \ \ge\ c_A I_e\,v^{-2}\min(1,p_2p_3/v^3).                             \tag{4}
\]
每个源至多拥有 \(D_n\) 条相交边。矩形可彼此相交，但 \(\sum_{\text{非例外长尾}}(v/a)|\mathcal R_e|\le D_n Z(K)\)，并未把重复源点当成额外缺量。

把 \(A_e\le v^{57/8}\) 的边直接记入 \(O(v^{57/8}S)=o(v^7S)\)。其他边有
\[
 p_2p_3\ge A_e^{2/3}>v^{19/4},\qquad
 v^{-2}\min(1,p_2p_3/v^3)\ge v^{-1/4}.
\]
由 (4)，对固定 \(K,T\)，所有 \(p_2>T v^3\) 的三临界边遂满足
\[
 \sum I_e\le C_Av^{1/4} Z(K)+O(v^8S+v^{57/8}S).                       \tag{5}
\]
其余 \(p_2\le T v^3\) 边由 R031 固定窗结论直接给 \(o_T(v^7S)\)。这里的 \(K,T\) 都不依赖 \(v\)。

## 3　四个以上临界轴与外部边界

对于两端均深的 \(q\ge4\) 边，R023 的真实角估计结合 \(\sum\varepsilon_i r_i^n=O(v^5S)\)，逐个 \(q\) 给
\[
 O_{n,A}\big((v^{q+5}+v^{q+6q/(2q-1)}) S\big)=o(v^7S).
\]
其中 \(q+6q/(2q-1)\ge 52/7>7\)。涉及外部边界的 \(q\ge3\) 边应用 R014 原单源界 \(C_Av^q[\varepsilon_i+\varepsilon_j+v^{q/(q-1)}]r_i^n\)。边界端的邻接总权至多 \(C_n S(\partial I)=O(v^5S)\)，故这些边总计 \(O(v^8S)\)。这两类均未动用 (3) 的管外预算。合并 (5) 即证 (1)。

## 4　何时得到不扣预算的原交量小量

若 \(p_2p_3\ge\tau v^4\)，小 \(v\) 下 (4) 的倍率至少为 \(c_A\tau/v\)。先固定任意 \(0<\sigma<1/4\)，再取 \(K\ge K_\sigma\)，使 R031 的实际 \(Y_{i,\sigma}\subset\mathcal C_i(K)\)，从而 \(Z(K)\le Z_\sigma\)。取固定 \(T>2K/c_2+1\)。\(p_2\le T v^3\) 仍按 R031 处理，其他非例外边按 (4) 不需 \(A_e>v^{57/8}\) 的门槛。于是
\[
 \sum_{p_2p_3\ge\tau v^4} I_e
 \le C_A\tau^{-1}v Z_\sigma+O_A(v^8S)+o_T(v^7S).
\]
其中第一项的常数不依赖 \(\sigma,K,T\)。先令 \(v\to0\)，再令 \(\sigma\downarrow0\)，得到 (2)。

对一般长尾，(1) 的用途是**与二轴完整价格同时结算**。写 \(\lambda_i=\inf_{\mu_i\le m\le|W_i|}\frac{d}{dm}(R_i+m^{1/n})^n\ge c_A v/a\)。对任意
\(B_i\subset[(x_i+W_i')\setminus F_i]\cap\mathcal C_i(K)\)，强 CDF 界给
\[
 J_i\ge\lambda_i\left(|B_i|+|[(x_i+W_i')\setminus F_i]\setminus\mathcal C_i(K)|\right).
\]
故管内每个已定价的实际源层保留完整 \(\lambda_i|B_i|\)，(1) 只需花费管外物理缺量预算的 \(O_A(v^{1/4})\) 份额。此处没有从全体 \(J_i\) 抽取 \(v^{1/4}\) 份额而污染六阶主项。

仅从 (4) 并无全参数的 \(I_e/v\) 源成本：例如 \(p_1=p_2=p_3=v^{7/3}\) 时 \(A_e=v^7\)，(4) 的倍率为 \(v^{-1/3}\)。这只是该粗矩形倍率的代数边界例，不是 \(C_n\) 的真实反例。因此本文不将 (1) 改写成全部三临界原交量均 \(o(v^7S)\)，也不据此宣称一般维数七阶常数。

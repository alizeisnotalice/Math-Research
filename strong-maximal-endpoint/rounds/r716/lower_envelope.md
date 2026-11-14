# 一般负责混合的源端下包络对偶

设 $R_i$ 为非空有限目录，矩形均有正有限体积，$X=\bigcup_i R_i$，$Q=\mathrm{Unif}(X)$。令 $0<p<1$，任意 $0\le f\le1$ 仅定义在 $X$ 上，$a_i(f)=|R_i|^{-1}\int_{R_i}f$。外部源质量不参与任何约束，可直接延零。定义 $m_X(p)=\min\{\int_Xf:a_i(f)\ge p\ \forall i\}$。若 $P_\pi=\sum_i\pi_i\mathrm{Unif}(R_i)$，以 $b_p(P_\pi)$ 表示 $P_\pi$ 质量至少 $p$ 的 $[0,1]$ 检验的最小 Lebesgue 质量。令 $m_{\rm resp}(p)=\max_{\pi\text{ 可负责}}b_p(P_\pi)$。有精确式
\[
 m_{\rm resp}(p)=\min_{0\le f\le1}\left\{\int_X f:\ \int_X L_f\,dQ\ge p\right\},\qquad L_f(x)=\min_{i:x\in R_i}a_i(f). \tag{1}
\]

**固定预算证明。** 按矩形成员关系划分有限正体积原子 $A$。检验只需每个原子一个密度变量，可以等质量指标化。负责行权集合 $\Pi$ 是运输多面体的像：把每个 $Q(A)$ 任意分给包含该原子的行（原子无原子性保证可以真实分割），由各行取得的总质量给 $\pi_i$。因此 $\Pi$ 紧凸。它也由概率单纯形内的 Hall 条件 $\sum_{i\in J}\pi_i\le Q(\bigcup_{i\in J}R_i)$ 刻画。固定预算 $s\in[0,|X|]$，记 $\mathcal F_s=\{0\le f\le1:\int f\le s\}$ 为紧凸原子多面体。则
\[
 b_p(P_\pi)\le s\quad\Longleftrightarrow\quad \max_{f\in\mathcal F_s}\pi\cdot a(f)\ge p.
\]
对所有 $\pi$ 同时成立等价于 $\min_{\pi\in\Pi}\max_{f\in\mathcal F_s}\pi\cdot a(f)\ge p$。由有限双线性极小极大定理，这等价于 $\max_{f\in\mathcal F_s}\min_{\pi\in\Pi}\pi\cdot a(f)\ge p$。内层逐接收原子把全部质量分给覆盖它且 $a_i(f)$ 最小的一行，正得 $\int L_f dQ$。最小可行预算即 (1)。为核实最不利负责权亦可取到，$b_p$ 关于 $\pi$ 连续：若 $\pi_j\to\pi$，最优检验的极限子列给下半连续；将 $\pi$ 的最优检验以任意小正比例与 $1_X$ 混合，命中严格超过 $p$，对所有大 $j$ 仍可行，便有上半连续。证明无对偶锥归一化，零倍率不造成例外；$s=0,|X|$ 也由同一预算等价性覆盖。约束等号若需改为严格，可作同一混合后取下确界。

**核心检验的实物核对。** 在[第715轮](../r715/joint_cost.md)取 $p=1/2$，核心体积 $\epsilon=\epsilon_*=[2+na(1-a)/(1+a)]^{-1}$，主盒体积一、每条臂体积 $(1+a)\epsilon$、独占外片体积 $a\epsilon$，$D=|X|=1+na\epsilon$。对 $f=1_K$，大行平均 $\epsilon$、臂行平均 $1/(1+a)\ge\epsilon$。于是整个主盒（含核心）的下包络是 $\epsilon$，各外片上为 $1/(1+a)$。所以
\[
 \int L_f dQ=\frac{\epsilon+na\epsilon/(1+a)}D=\frac12.
\]
与该轮的最优负责混合（核心密度为 $1/(2\epsilon)$，整个空间无更高密度）给出的相反质量界合并，(1)确切达到 $m_{\rm resp}=\epsilon$。

## 删选最优时，共同检验的高侧质量

再假设每个非空子目录 $J$ 的得分均不超过完整目录，$|X_J|/m_J(p)\le |X|/m_X(p)$（每个子目录在自己的完整并上计价）。记 $m=m_X(p)>0$，$s=\int_X f$。对 $p\le\beta\le1$，取 $J_\beta=\{i:a_i(f)\ge\beta\}$。若非空，$(p/\beta)f$ 限于 $X_{J_\beta}$ 使所有这些行成功，所以 $m_{J_\beta}(p)\le (p/\beta)s$。因此对 $U_f(x)=\max_{i:x\in R_i}a_i(f)$ 有
\[
 Q(L_f\ge\beta)\le Q(U_f\ge\beta)=|X_{J_\beta}|/|X|\le \frac{p}{\beta}\frac{s}{m}. \tag{2}
\]
空子目录也成立。若 $f$ 为 (1) 的任一可行检验，则 $\int(p-L_f)_+dQ\le\int(L_f-p)_+dQ$；对 (2) 在 $[p,1]$ 积分，得
\[
 \int|L_f-p|\,dQ\le 2p\log(1/p)\,s/m. \tag{3}
\]
在中位数处常数为 $\log2$。故若 $s/m\to0$，负责源端的行平均下包络必趋于门槛（$L^1(Q)$），且达到或超过门槛的接收比例本身至多 $s/m$。结果同时保留了按Q计绝大多数接收点的行平均下包络从下侧靠近门槛的可能性；(2)–(3) 尚未给出删选最优几何中的常数比例下界。第715轮臂目录在最坏参数处，删除主盒后的诸臂并体积为 $\epsilon(1+na)$，其最小质量是 $(1+a)\epsilon/2$（在核心内即可满足所有臂）；子目录得分 $2(1+na)/(1+a)$ 增长，而全目录得分有界，故当维数充分大时它不满足本节的删选最优前提。

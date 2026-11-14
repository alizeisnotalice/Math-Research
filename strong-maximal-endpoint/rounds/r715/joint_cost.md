# 负责行权的真半质量价格可以比最优行权小线性倍

令 $m$ 为使有限目录每个矩形平均至少 $1/2$ 的最小源质量（$0\le f\le1$）。写 $b(P)$ 为 $P$ 质量至少 $1/2$ 的检验 $0\le f\le1$ 所需的最小 Lebesgue 体积，$m_{\rm resp}=\max_{\pi\ \rm 负责}b(P_\pi)$，其中 $P_\pi$ 混合各行的均匀律。有限目录下，$m=\max_{\pi} b(P_\pi)$；负责权另须可由实际接收分割实现。

**结论。** 不存在与维数无关的正常数 $c$ 使所有轴平行矩形目录满足 $m_{\rm resp}\ge c m$。下面的真正矩形在每一维 $n$ 上可以达到
\[
 \frac{m}{m_{\rm resp}}=1+\frac{3-2\sqrt2}{2}\,n . \tag{1}
\]
对应对数费用是 $\log n+O(1)$。这个反例针对全部目录，不要求中央规模的接收增益。

## 几何与容量

取 $0<\epsilon\le1/2$, $0<a\le1$。主盒为 $R_0=[0,1]^n$，核心为 $K=[0,\epsilon]\times[0,1]^{n-1}$。记核心边长 $l_1=\epsilon, l_2=\cdots=l_n=1$，第 $i$ 条臂盒 $R_i$ 将第 $i$ 轴由 $[0,l_i]$ 换为 $[-al_i,l_i]$，其余轴保持核心区间。有理 $a,\epsilon$ 给出有理端点。于是 $R_i=K\sqcup W_i$（模零），外片 $W_i$ 独占、各体积 $a\epsilon$，每臂体积 $(1+a)\epsilon$。设 $B=R_0\setminus K$，并体积为 $D=1+na\epsilon$。取实际源 $E=[0,1/2]\times[0,1]^{n-1}$，其体积为 $1/2$ 且包含核心，因为每个臂的核心占比为 $1/(1+a)\ge1/2$，故全部约束成立。大行本身给匹配下界 $m=1/2$，无论 $a,\epsilon$ 如何变化。单元内可精确指标化；严格阈值可在余部混入趋零正比例恢复同一下确界。

记 $w=\pi_0$。负责权必有 $w\le 1/D$，每个外片必归给相应的行。特别地，$w_*=1/D$ 可达，且达到时权唯一为 $\pi_0^*=1/D,\ \pi_i^*=a\epsilon/D$。它由把 $R_0$ 全部分给大行、各外片给各臂来实现。

## 全部负责混合的最优值

对 $\pi^*$，三个层 $K,B,W_i$ 上的 Lebesgue 密度从高至低为
\[
 d_K=\frac{1+na/(1+a)}D,\quad d_B=\frac1D,\quad d_W=\frac{a}{(1+a)D}.
\]
达到一半质量的最小体积检验按密度由高到低填。核心容量足够的条件 $\epsilon d_K\ge1/2$ 等价于
\[
 \epsilon\ge\epsilon_*:=\big[2+na(1-a)/(1+a)\big]^{-1}.
\]
若该式成立，取核心内体积 $s=1/(2d_K)$ 的子集即可。若 $\epsilon\le\epsilon_*$，取满核心后再在 $B$ 内填 $D(1/2-\epsilon d_K)$ 体积。总填入量是 $s=1/2-na(1-a)\epsilon/[2(1+a)]\le1/2$，且在此分支有 $s\ge\epsilon$，所以 $B$ 的真实容量 $1-\epsilon$ 足够。两分支在交界给同一值 $s=\epsilon_*$。

这两个检验也分别给**所有**负责权的上界。第一分支的核心子集在大行上的命中率为 $s$，在每臂上的命中率为 $s/[(1+a)\epsilon]\ge s$。第二分支的检验在大行上的命中率为 $s\le1/2$，每臂的命中率为 $1/(1+a)\ge1/2$。两种情形下检验的 $P_\pi$ 概率只依赖 $w$ 且随 $w$ 非增，所以任意负责 $w\le w_*$ 时，同一检验仍至少承载一半概率。与 $\pi^*$ 的密度排序下界合并，完整真值为
\[
 m_{\rm resp}=\begin{cases}
 \tfrac12-na(1-a)\epsilon/[2(1+a)], &0<\epsilon\le\epsilon_*,\\[2pt]
 (1+na\epsilon)/[2(1+na/(1+a))], &\epsilon_*\le\epsilon\le1/2.
 \end{cases} \tag{2}
\]
特别对原提议的**精确双体积臂** $a=1$，有 $\epsilon_*=1/2$，所有可行 $\epsilon$ 上 $m_{\rm resp}=m=1/2$。只按核心峰值倒推体积会漏掉 $\epsilon<1/2$ 时在 $B$ 内必须承载的剩余概率。

## 极值与尺度

对固定 $a<1$，(2) 第一支随 $\epsilon$ 递减，第二支递增。故最坏的相对价格在 $\epsilon_*$，此时 $m/m_{\rm resp}=1+na(1-a)/[2(1+a)]$。导数 $[a(1-a)/(1+a)]'=(1-2a-a^2)/(1+a)^2$，唯一最大点 $a=\sqrt2-1$ 给 (1)。若采用 $a=1/2$，取 $\epsilon=6/(12+n)$ 则精确比值为 $1+n/12$。这里直接得到有理端点的实际矩形；删去零体积边界后，计算所有行权与检验只需 $K,B,W_1,\ldots,W_n$ 共 $n+2$ 个成员原子。

在 $a$ 固定且小于一、$\epsilon=\epsilon_*$ 时，$|X|/m=2(1+na\epsilon_*)\le 4/(1-a)$ 随 $n$ 有界。所以此族足以排除只依据全部真实目录的常数比例转移，尚未证明在 $\log(|X|/m)\ge\mu n-O(\log n)$ 的中央近极值目录上，负责价格也会有同等损失。

## 同一负责混合的联合费用

取 $a\in(0,1)$ 固定且 $\epsilon=\epsilon_*$，$Q=\mathrm{Unif}(X)$，令 $P=P_{\pi^*}$ 为上面实现负责最优半质量价格的具体混合。相对于 $Q$，$P$ 在 $K,B,\bigcup W_i$ 上的密度比分别为 $C=1+na/(1+a),1,a/(1+a)$。因临界核心恰好承载半质量，$P(K)=1/2$；所有外片的总概率是 $P_W=na^2\epsilon_*/[(1+a)D]$。于是
\[
 D_{\rm KL}(P\Vert Q)=\tfrac12\log C +P_W\log[a/(1+a)],\qquad
 \log\frac{|X|}{b(P)}=\log(2C).
\]
第二式用 $b(P)=\epsilon_*$、$D=2C\epsilon_*$。因此当 $n\to\infty$，同一混合的 KL 为 $\tfrac12\log n+O_a(1)$，半质量检验的对数价格为 $\log n+O_a(1)$，后者减前者为 $\tfrac12\log n+O_a(1)$。同时 $\log(m/m_{\rm resp})=\log n+O_a(1)$。三项在实际价格中按下式结算：
\[
 \log\frac{|X|}{m}
 =D_{\rm KL}(P\Vert Q) +\left[\log\frac{|X|}{b(P)}-D_{\rm KL}(P\Vert Q)\right]
 -\log\frac{m}{m_{\rm resp}}
 =\log(2D)=O_a(1).
\]
以上 KL 计算针对这个已实现的负责混合，本证明未对 KL 另作全体行权优化。

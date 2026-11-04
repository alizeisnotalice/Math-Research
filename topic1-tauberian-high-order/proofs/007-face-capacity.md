# 面容量引理：CDF 有效端与角部体积盈余

这是第 004 轮真实二阶定理之后的受限三阶结果。本节所有常数只依赖固定的 \(n\ge2\) 和规则性常数 \(D_0\)。不对任意真实极值序列假设下述条件自动成立。

## 1　精确量词和结论

设 \(H=2^n,\ a=Ht,\ v=a^{1/n}\downarrow0\)。\(E\) 有限正测度，\(K\subset H_t(E)\setminus E\) 可测，\(M=|K|\)。有有限个严格见证 \(Q_i=Q(x_i,r_i)\)，中心距 \(>\max(r_i,r_j)\)，覆盖 \(K\)。\(S=\sum_i r_i^n\)。假设

1. 可测分配 \(K=\bigsqcup_i G_i,\ G_i\subset Q_i,\ g_i=|G_i|\)。另有 \(P_i\subset G_i\setminus\bigcup_{j\ne i}Q_j\)，\(m_i=|P_i|\)，满足
   \[
   (1-D_0v)ar_i^n\le m_i\le g_i\le ar_i^n,\qquad
   P_i\subset Q(x_i,D_0vr_i).
   \tag{A1}
   \]
2. 对 \(P_i\) 每个点 \(x\) 有一个严格的真实见证 \(Q(x,r(x))\)，其中 \(r(x)\le r_i\)，且整个 \(G_i\subset Q(x,r(x))\)（差零测亦可）。

令 \(R_i=(g_i/a)^{1/n}\)。则 \(R_i<r(x)\le r_i\)、\(R_i=(1+O(v))r_i\)。将上述第二条记为 (A2)。以下引理只以 (A1)–(A2) 为输入，不主张一般输入自动满足这些条件。后续尾部证明只对规则部分调用，并另行保留、计费例外部分。

**定理。** 对所有满足这些条件的族，充分小的 \(v\) 时，
\[
 |E|+a S\ \ge\ \frac{(1+v)^n}{a}M
       -C_{n,D_0}\left(\frac va\sum_i(g_i-m_i)+v^3 S\right).
 \tag{F1}
\]
## 2　目标尾部按实际半径增量收费

在坐标 \(k\)，先删去零测坐标异常点，并将所有投影端点解释为本质端点。记 \(P_i-x_i\) 的投影原始端点为 \(\ell,h\)，其质量分布函数为 \(F(s)=|\{x\in P_i:x_k-x_{i,k}\le s\}|\)。\(F\) 连续，\(F(\ell)=0,F(h)=m_i\)。取
\[
 \lambda_i=n 2^n a r_i^{n-1}.
\]
令 \(\ell'\) 为 \(F(s)-\lambda_i s\) 在 \([\ell,h]\) 上的任一极小点。反射坐标作同一构造得到右端 \(h'\)。于是
\[
 F(\ell')\le\lambda_i(\ell'-\ell),\qquad
 F(\ell'+d)\ge\lambda_i d\quad(0\le d\le h-\ell').
 \tag{F3}
\]
若 \(\ell'+d>h\)，在 \(d=O_{n,D_0}(vr_i)\) 范围仍成立，因为 \(F(\ell'+d)=m_i\) 而 \(\lambda_i d=O(v)m_i\)。右端反射式相同。

各方向尾部总质量至多 \(\lambda_i\Delta_i\)，其中
\[
 \Delta_i=\sum_k[(\ell'_{ik}-\ell_{ik})+(h_{ik}-h'_{ik})]
 =O(vr_i).
\]
对小 \(v\) 此量为 \(O(v)m_i\)，故每个 \(\ell'_{ik}<h'_{ik}\)。去尾后的盒 \(W_i'\) 至少容纳 \(m_i-\lambda_i\Delta_i\ge m_i/2\) 体积；各边至多 \(2D_0vr_i\)，所以各边至少 \(c_{n,D_0}vr_i\)。原始盒同样有这一下界。

用原始与有效端的中点作盒 \(W_i\)，即两端分别为 \((\ell_{ik}+\ell'_{ik})/2\)、\((h_{ik}+h'_{ik})/2\)。这是相对 \(x_i\) 的盒。\(W_i\) 比 \(W_i'\) 扩张的总边长是 \(\Delta_i/2\)。乘积展开给
\[
 |W_i|-m_i\ge
 \left[\tfrac12(c_{n,D_0}vr_i)^{n-1}-\lambda_i\right]\Delta_i
 \ge c'_{n,D_0}(vr_i)^{n-1}\Delta_i\ge0.
 \tag{F4}
\]
最后一步用 \(a=v^n\)。令 \(\epsilon_i=|W_i|/m_i-1\ge0\)，它一致有界。原始盒 \(W_i^0=\prod_k[\ell_{ik},h_{ik}]\) 还满足
\[
 |(x_i+W_i^0)\setminus P_i|\le C\epsilon_i m_i .
 \tag{F5}
\]
事实上 \(|W_i^0|-|W_i|\le C(vr_i)^{n-1}\Delta_i\)，代 (F4) 即得。第四轮单独删除固定分位会损失 \(O(v)m_i\) 的质量；(F4) 在这里保留了整个输入 \(m_i\)。

## 3　单临界方向没有正体积交叠

作 \(X_i=x_i+W_i+[-R_i/2,R_i/2]^n\)。它包含于 \(Q_i\)。若 \(X_i,X_j\) 相遇，则
\[
 \max(r_i,r_j)<\|x_i-x_j\|_\infty
 \le(1/2+D_0v)(r_i+r_j),\qquad r_j/r_i=1+O(v).
 \tag{F6}
\]
每格的相交邻居数为 \(O_n(1)\)，由不交半盒 \(Q(x_j,r_j/2)\) 的体积计数。记 \(s=r_i+r_j\)。称坐标 \(k\) 临界，如果 \(|x_{i,k}-x_{j,k}|\ge(1/2-A v)s\)，其中 \(A=5(1+D_0)\) 可用。至少一个坐标临界，任一临界方向的相交投影长至多 \(Cvr_i\)。非临界坐标上任意 \(x\in P_i,z\in P_j\) 均有 \(|x_k-z_k|<\min(R_i,R_j)\)，由 (A1)、(F6) 可得。临界坐标的所有源点都在目标点同一侧，对小 \(v\) 成立。

若仅一个临界方向 \(k\)，反射使 \(i\) 在 \(j\) 左边，设两中心该坐标差为 \(d_0>0\)。必须
\[
 d_0+\ell'_{jk}-h_{ik}\ge R_i,\qquad d_0+\ell_{jk}-h'_{ik}\ge R_j .
 \tag{F7}
\]
证明第一式。否则取 \(x\in P_i\) 使其临界坐标偏移充分接近原始右端 \(h_{ik}\)，使 \(d_*:=x_{j,k}+\ell'_{jk}-x_k<R_i\)。记 \(r=r(x)\)。取 \(0<\rho<\min\{R_i-d_*,vr_i\}\)。所有目标点中临界坐标不超过 \(x_{j,k}+\ell'_{jk}+(r-R_i)+\rho\) 的点都在 \(Q(x,r)\) 中。其余坐标用非临界估计，临界坐标用同序及 \(d_*+\rho<R_i\)。由 (F3)，这部分质量至少 \(\lambda_j[(r-R_i)+\rho]\)。它和 \(G_i\) 不交，而 (F6) 给 \(\lambda_j\ge n a r_i^{n-1}\)。因此严格见证要求的剩余容量
\[
 a r^n-g_i=a(r^n-R_i^n)\le n a r_i^{n-1}(r-R_i)
\]
被超过。矛盾。另一式用反向源点和目标右尾。对 (F7) 取平均，正是 \(W_i,W_j\) 的中点端距至少 \((R_i+R_j)/2\)，所以法向投影至多边界接触。

## 4　双坐标穿透由盒体积盈余支付

至少三个临界方向时交体积为 \(O(v^3)r_i^n\)。恰有两个时，令其轴向相交长度为 \(L_1,L_2\le Cvr_i\)，并置 \(B=L_1L_2/(v^2r_i^2)\)。以下证明
\[
 |X_i\cap X_j|\le C v^2(\epsilon_i+\epsilon_j+v)r_i^n .
 \tag{F8}
\]
可设 \(R_i\ge R_j\)，并在两个临界坐标各自反射，使源簇 \(i\) 在目标簇 \(j\) 左侧。分别选取原始包围盒 \(x_i+W_i^0,x_j+W_j^0\) 中面向对方的双面角条：在临界坐标 \(k\) 上从原始相对极端向里取厚 \(\theta L_k\)，其余坐标完全保留。\(\theta>0\) 为只依赖 \(n,D_0\) 的固定常数，取 \(\theta\le1/4\) 并充分小，保证每个厚度小于两盒相关原始边长。由原始边长上下可比 \(vr_i\)，这两个角条的体积分别至少 \(cB m_i,cB m_j\)。

每个源角条的 \(P_i\) 点与每个目标角条的 \(P_j\) 点距离不超过 \(R_i\)。非临界坐标由上一节；在临界方向，中点端之间的距离至多 \((R_i+R_j)/2-L_k\)，原始面对端距只会更小，角条向内移动合计至多 \(2\theta L_k\)。如 \(B>C(\epsilon_i+\epsilon_j+v)\) 且 \(C\) 充分大，由 (F5) 存在一个源 \(P_i\) 点，目标角条中的 \(P_j\) 质量大于 \(cBm_j/2\)。该质量超过 \(ar_i^n-g_i=O(v)m_i\)，同时与 \(G_i\) 不交，与源点的真实见证矛盾。故 \(B\le C(\epsilon_i+\epsilon_j+v)\)，乘其余坐标的交长上界即得 (F8)。

## 5　体积结算

盒乘积的几何平均不等式给
\[
 |X_i|\ge(R_i+|W_i|^{1/n})^n
 \ge (R_i+m_i^{1/n})^n+c v\epsilon_i r_i^n .
 \tag{F9}
\]
最后一步用 \(m_i^{1/n}\asymp vr_i\) 和 \(\epsilon_i=O(1)\)。由 (A1) 及中值定理，
\[
 (R_i+m_i^{1/n})^n
 \ge \frac{g_i}{a}(1+v)^n-C\frac va(g_i-m_i).
 \tag{F10}
\]
每格邻居数有界、邻居半径可比。所有两两交叠合计至多 \(C v^2\sum_i\epsilon_i r_i^n+C v^3 S\)，其第一项由 (F9) 的严格盈余吸收。二项容斥和 \(\bigcup X_i\subset\bigcup Q_i\) 给
\[
 \sum_i|X_i|-\sum_{i<j}|X_i\cap X_j|
 \le \left|\bigcup_iX_i\right|\le |E|+\sum_i |Q_i\setminus E|\le |E|+aS .
 \]
合并即得 (F1)。

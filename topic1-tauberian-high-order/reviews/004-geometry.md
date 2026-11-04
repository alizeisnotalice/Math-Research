# 锐候选 S10–S14 独立几何严审

审查对象：`r004-global/global_sharp_candidate.md`，以 S10–S14 为主，另核查 S3 的阈值如何代入。这是对新候选的重新审查，没有把旧版面柱结论作为前提。

**判定：本审查范围 proof-valid。** 在 S1–S9 的下述明确依赖成立时，分位矩形的扩大体积、一轴尾质量矛盾、至少两临界轴的交损、有界邻居及例外交损全部闭合，适用于每个固定 (n\ge2)。未发现需修复的几何缺口。这里不单独签署 S1–S9 或完整主定理；全局审查必须合并负责前半部的独立审查结果。

## 1. 必须由 S1–S9 提供的前提

使用的不是笼统的“私有簇有大见证”，而是下列具体事实。设 (c=4^n,a=2^nt,v=a^{1/n},b=\sqrt v,\delta=Lv)。对规则指标 (i)：

1. (P_i\subset K\cap Q_i\)，且 (P_i\cap Q_j=\varnothing) 对每个 (j\ne i) 成立；(K\subset E^c)。中心满足严格分离 (\|x_i-x_j\|_\infty>\max(r_i,r_j))。
2. 有可测分配 (K=\bigsqcup_iG_i)，满足 (P_i\subset G_i\subset Q_i)，(g_i=|G_i|\le ar_i^n)，并且 (m_i=|P_i|\ge(1-\delta)ar_i^n)。
3. (P_i\subset Q(x_i,c\delta r_i))，且小指定半径点总质量至多 (a(6c\delta)^n r_i^n/(1-t))。
4. 定义 (R_i=(g_i/a)^{1/n})。每个未被小半径删除的私有点 (x) 的**指定**见证 (Q(x,r(x))) 含整个 (G_i)，且

   \[
   (1-\delta)r_i\le R_i<r(x)\le r_i.
   \]

5. 指定见证严格满足

   \[
   |Q(x,r(x))\setminus E|<a r(x)^n.
   \]

第 4 条中“含整个 (G_i)”及 (r(x)\le r_i) 都不能省略：尾质量矛盾需要同时把 (g_i) 作为已有亏损、把 (ar_i^n-g_i) 作为剩余预算。仅含 (P_i) 或仅有半径下界不足。

读到的 S9 论证与这些前提一致：由 (G_i\subset Q(x_i,0.9r_i))、私有中心位移至多 (c\delta r_i)，和 (r(x)>(1-\delta)r_i)，当 ((1+c)\delta<0.1) 时确实有整个 (G_i\subset Q(x,r(x)))。指定半径上界来自降序选择及私有性，不是任取见证的通性。

以下所有“充分小”的要求可同时由

\[
 0<\delta<\delta_0:=\min\{\epsilon,1/[64(c+1)]\},\qquad
 \delta/b<2^{-n-1}
 \tag{G0}
\]

保证。对固定 (T)，因 (T\le L\le T^2)，有 (\delta\le T^2v\) 且 (\delta/b\le T^2\sqrt v\to0)，故不存在阈值循环。

## 2. S10：分位削去及端点不取到

对每一坐标，原 (P_i) 的边缘分布

\[
 F_k(u)=|\{x\in P_i:x_k-x_{i,k}\le u\}|
\]

连续，因为每个坐标超平面是 (n) 维 Lebesgue 零测集。右侧尾分布也连续。各方向删除的质量不超过 (bm_i)，所以即使各方向删除相交，联合删除仍至多 (2nbm_i)。再减去小半径点给出 S10。

令 (P_i^*) 为这两类删除后的集合。每个其中点仍满足**相对于原 (P_i)** 的两侧分位质量下界。取坐标趋于

\[
 \ell_{ik}=\inf_{x\in P_i^*}(x_k-x_{i,k}),\quad
 h_{ik}=\sup_{x\in P_i^*}(x_k-x_{i,k})
\]

的序列，再用分布连续性，可得

\[
 |\{z\in P_i:z_k-x_{i,k}\le\ell_{ik}\}|\ge bm_i,
 \quad
 |\{z\in P_i:z_k-x_{i,k}\ge h_{ik}\}|\ge bm_i.
 \tag{G1}
\]

不要求端点取到，也不要求 (P_i^*) 闭。尾集合可以包含已删除的小半径点；证明只把它们作为外部亏损质量，并不向这些尾点索取见证。这正是采用原 (P_i) 分位数的作用。

因为 (m_i\ge(1-\delta)ar_i^n)，

\[
 \frac{m_i^*}{ar_i^n}\ge
 (1-\delta)(1-2nb)-\frac{(6c\delta)^n}{1-t}
 =1-O_n(b+\delta+\delta^n).
 \tag{G2}
\]

它在固定 (T) 的整个 (L\in[T,T^2]) 上一致趋于 1，最终为正。

## 3. S11：扩大矩形体积

记 (w_{ik}=h_{ik}-\ell_{ik}\ge0)。由于端点盒包含 (P_i^*)，

\[
 m_i^*\le\prod_kw_{ik},\qquad \sum_kw_{ik}\ge n(m_i^*)^{1/n}.
\]

每个 (X_i) 的边长为 (R_i+w_{ik})。展开乘积后其余项非负，故

\[
 |X_i|\ge R_i^n+R_i^{n-1}\sum_kw_{ik}
 \ge g_i/a+nR_i^{n-1}(m_i^*)^{1/n}.
\]

由 (R_i\ge(1-\delta)r_i) 和 (G2)，后一项至少为

\[
 [n-O_n(b+\delta+\delta^n)]v r_i^n.
\]

因此 S11 的 (n-o_T(1)) 系数正确。所有端点均在 ([-c\delta r_i,c\delta r_i])，又 (R_i\le r_i)，所以

\[
 X_i\subset Q(x_i,(1/2+c\delta)r_i)\subset A_i\subset D.
\]

例外指标的 (X_i=B_i) 也满足同一个粗外包盒。此处不需要规则矩形包含自己的半盒。

## 4. S12：交叠半径可比，邻居没有漏掉巨大或微小盒

若 (X_i\cap X_j\ne\varnothing)，所有指标均有上述粗外包，因此

\[
 \max(r_i,r_j)<\|x_i-x_j\|_\infty
 \le(1/2+c\delta)(r_i+r_j).
\]

于是

\[
 |r_i-r_j|<2c\delta(r_i+r_j).
 \tag{G3}
\]

当 (\delta<1/(6c)) 时，(G3) 给出 (1/2<r_i/r_j<2)。所以只要是与 (X_i) 真正相交的邻居，就自动排除了巨大与微小尺度；无需对原族全部邻居作额外假设。

固定 (i)，每个相关 (B_j) 的边长 (r_j\ge r_i/2)。且

\[
 \|x_j-x_i\|_\infty+r_j/2
 \le(1/2+c\delta)(r_i+r_j)+r_j/2
 <3r_i
\]

由 (r_j<2r_i) 和 (\delta<1/(6c)) 可得。因此所有相关 (B_j\subset Q(x_i,3r_i))。原半盒两两不交，体积计数给邻居数至多 (12^n)。这个论证不要求规则性，也覆盖规则–例外的有序计数。

## 5. 临界坐标至少一个，至少两个时交损为二次

对规则 (i,j)，令 (s=r_i+r_j)。中心分离给出至少一个坐标满足中心差 (>\max(r_i,r_j)\ge s/2)，因此至少有一个临界坐标。

任一临界坐标上，借粗外包盒，相交投影长度至多

\[
 (1/2+c\delta)s-[1/2-(3c+3)\delta]s
 =(4c+3)\delta s.
\]

若至少两个坐标临界，取其中两个应用此界，其余投影长度用 (X_i) 边长至多 (2r_i)，再用 (s\le3r_i)，得

\[
 |X_i\cap X_j|\le
 [3(4c+3)]^2 2^{n-2}\delta^2r_i^n.
 \tag{G4}
\]

因此在任意固定维数 (n\ge2) 中都是 (O_n(\delta^2)r_i^n)。不需要“只有两个”临界坐标；有更多时仍可只选两个。

## 6. 仅一临界轴：尾集合确实进入指定见证，且与已有质量不交

设唯一临界轴为 (q)，交换指标后 (x_{j,q}>x_{i,q})。非临界坐标的中心差严格小于阈值，加入两私有簇偏移，得到任意 (x\in P_i,z\in P_j) 满足

\[
 |x_k-z_k|<[1/2-(2c+3)\delta]s.
\]

另一方面，(G3) 给出

\[
 \min(r_i,r_j)\ge(1/2-c\delta)s,
\]

故

\[
 \begin{aligned}
 \min(R_i,R_j)&\ge(1-\delta)(1/2-c\delta)s\\
 &=[1/2-(c+1/2)\delta+c\delta^2]s\\
 &\ge[1/2-(c+1)\delta]s
 >[1/2-(2c+3)\delta]s.
 \end{aligned}
 \tag{G5}
\]

这证明 S13。临界轴上的任意两私有点又有

\[
 z_q-x_q\ge[1/2-(4c+3)\delta]s>0
 \tag{G6}
\]

在统一阈值 (G0) 下成立。因此法向的距离无需再取绝对值。

设

\[
 d=x_{j,q}-x_{i,q}+\ell_{jq}-h_{iq}.
\]

若 (d<R_i)，取 (\rho=(R_i-d)/2>0)，并由上确界取 (x\in P_i^*) 满足

\[
 x_q-x_{i,q}>h_{iq}-\rho.
\]

对原 (P_j) 的左尾

\[
 T_j^-:=\{z\in P_j:z_q-x_{j,q}\le\ell_{jq}\}
\]

中的每个 (z)，(G6) 以及端点界给出

\[
 0<z_q-x_q<d+\rho<R_i.
\]

切向坐标由 (G5) 也严格小于 (R_i)。所以整个 (T_j^-\subset Q(x,R_i)\subset Q(x,r(x)))。端点是否取到不影响此严格余量。

由 (G1)，(|T_j^-|\ge bm_j)。并且

\[
 T_j^-\subset P_j\subset K\setminus Q_i,
 \qquad G_i\subset Q_i,
\]

所以尾集合与 (G_i) 严格不交。两者都在 (K\subset E^c)，且指定见证包含两者。因此

\[
 g_i+bm_j\le |Q(x,r(x))\setminus E|
 <a r(x)^n\le ar_i^n. \tag{G7}
\]

这一步不是把 (P_i^*\) 与尾集合相加，而是把 (G_i) 与原 (P_j) 的尾集合相加。候选稿对集合选择正确。

但 (g_i\ge m_i\ge(1-\delta)ar_i^n)，且半径可比给

\[
 bm_j\ge b(1-\delta)ar_j^n
 \ge b(1-\delta)2^{-n}ar_i^n
 >\delta ar_i^n
\]

在 (G0) 下成立，与 (G7) 矛盾。所以 (d\ge R_i)。对 (P_j^*) 的负向端点及原 (P_i) 的右尾对称应用同一论证，得到 (d\ge R_j)。因此

\[
 d\ge\max(R_i,R_j)\ge(R_i+R_j)/2.
\]

这正是两个 (X) 在 (q) 轴上的闭投影不发生正长度交叠的条件。投影可能端点接触；不能把本结论改写成两闭矩形完全不交，但它已给出所需 (n) 维交体积为零。

## 7. 规则–例外交叠及 S14 的求和

令 (i) 规则、(j) 例外，此时 (X_j=B_j)。中心分离提供某个坐标中心差 (>\max(r_i,r_j)\ge s/2)。该坐标的交长至多

\[
 R_i/2+c\delta r_i+r_j/2-s/2\le c\delta r_i\le2c\delta r_j.
\]

其余坐标的交长均不超过例外半盒边长 (r_j)。故

\[
 |X_i\cap X_j|\le2c\delta r_j^n.
 \tag{G8}
\]

这明确把成本计在例外半径上。每个例外指标的相交邻居至多 (12^n)，所以规则–例外所有无序对之和至多 (2c\,12^n\delta S_b)，不会漏掉被多个规则矩形碰到的例外盒。

两例外指标的矩形就是两不交半盒。规则–规则中，仅一临界轴的交体积为零，其余用 (G4) 及有界次数求和，得到

\[
 \sum_{i<j}|X_i\cap X_j|
 \le C_n\delta^2S+C_n\delta S_b.
\]

无需限制更高重叠：有限族的二项 Bonferroni 下界总成立，所以这个 pair sum 可以直接用于并集体积下界。

## 8. S3 阈值确实提供所需的对数小量

令 (F(l)=S^{-1}\sum_{z_i>l}r_i^n)。有限族的层蛋糕恒等式给

\[
 \int_0^\infty F(l)\,dl=S^{-1}\sum_i z_i r_i^n\le C.
\]

如果在整个 ([T,T^2]) 都有 (lF(l)>C/\log T)，则对 (dl/l) 积分给出上述积分 (>C)，矛盾。故存在所需 (L\in[T,T^2])，并且使用严格尾 (z_i>L) 与规则定义 (z_i\le L) 完全匹配。

由 (S_b/S=F(L)) 以及 (\delta=Lv)，

\[
 \delta S_b=LvSF(L)\le\frac{CvS}{\log T},\qquad
 \frac{S_b}{S}\le\frac{C}{T\log T}.
\]

另一方面 (\delta^2S\le T^4v^2S)。因此 S14 中

\[
 O_{n,T}(v^2S)+\frac{C_{n,C_0}}{\log T}vS
\]

的两种依赖均正确：第一项允许依赖固定 (T)，第二项的常数无需依赖 (T)。所有几何计数常数在 (G0) 后只依赖 (n)。

若进一步取候选第 6 节的 (T=v^{-1/8})，则 (L\le v^{-1/4})、(\delta\le v^{3/4})、(\delta/b\le v^{1/4}\to0)，所以 (G0) 仍最终满足；(\delta^2/v\le\sqrt v)，S14 的代入没有遗漏额外阈值。

## 9. 实际本机符号核验的范围

相关精确计算由 Terra 的本机 Mathematica 15.0 执行，我实际读取了源文件 `r004-wolfram/sharp/sharp.wl` 和对应输出。下列违反集均实际返回 `False`：

- `S12_ratio.txt`：相交外包条件下半径比越出 ((1/2,2))。
- `S12_difference.txt`：相交外包条件下 (|r-s|>2c\delta(r+s))。
- `S12_neighbor_box.txt`：相关半盒越出 (Q(x_i,3r_i))。
- `S13_min.txt`：非临界距离上界不小于 (\min(R_i,R_j))。
- `S13_critical_length.txt`：临界投影交长超过 ((4c+3)\delta(r+s))。

初次完整多变量 `S13_noncritical.txt` 输出 `MEMORY_LIMIT`，**不记为通过**。这不影响上文 (G5) 的完整手工推导。随后实际读取 `sharp2.wl` 及补充输出：`noncritical_projected.txt`、`capacity_n2.txt`、`capacity_n3.txt`、`tv_n2.txt`、`tv_n3.txt`、`tv_n4.txt` 均返回 `False`。其中二维、三维容量检查直接使用 (g\ge(1-\delta)a,m\ge(1-\delta)as^n,s\in[1/2,2],b>2^{n+1}\delta)，排除 (g+bm\le a)。完整联合几何模型的最终计算状态以计算审查者记录为准，本审查没有把未完成任务当作证据。

还应注意：仅为 S12 使用的条件 (\delta<1/(6c)) 本身不足以推出 (G6) 的法向点序；法向点序需要更小的固定阈值，例如 (G0)。候选稿写“充分小 (t)”允许这一条件，因而不是稿件漏洞。自动检验若复用过宽域，得到的点序反例也不能据此否定候选。

**最终结论：S10–S14 几何部分可通过；S3 的截断代入可通过。** 全局主定理仍须与 S1–S9 的独立计费及见证审查合并确认。

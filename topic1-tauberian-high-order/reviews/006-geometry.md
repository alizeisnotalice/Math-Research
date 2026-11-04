# R6 空间缓冲幂余项：独立几何审查

审查对象：`r005-global/spatial_buffer_power.md` 的 B2–B11，重点为缓冲图增长、不同插值参数的单临界错位、跨例外边的局部计费及增强的 B10。R4 的选族、局部化、深分配和完整指定见证作为已审定输入，本次没有重做 R4。

**结论：proof-valid。** 新空间缓冲接口全部闭合。B10 的局部平方计数也成立，因此稿件给出的增强值

\[
 \kappa=\frac1{n+2}
\]

可通过，不必退回粗备核值 (1/(2n+3))。没有发现尚缺的实质接口。以下列出每一处依赖及证明；末尾明确说明需要保留的量词和实际工具核验范围。

## 1. 已审定输入与本轮变量

R4 提供有限正半径的严格族

\[
 \|x_i-x_j\|_\infty>\max(r_i,r_j),\quad
 B_i=Q(x_i,r_i/2),\quad |B_i|=r_i^n,
\]

且所有闭 (B_i) 两两不交。对近首项饱和分支，有

\[
 d_i=e_i+1-\frac{m_i}{ar_i^n}\ge0,qquad
 \sum_i d_i r_i^n\le C_n vS.
 \tag{R1}
\]

规则指标定义为 (d_i\le\delta)，故分别有 (e_i\le d_i\le\delta)、

\[
 m_i\ge(1-d_i)ar_i^n,\quad
 P_i\subset Q(x_i,c_ne_i r_i),\quad
 P_i\subset G_i\subset Q_i,\quad
 R_i=(g_i/a)^{1/n}\in[(1-d_i)r_i,r_i].
\]

除去小半径点后，每个保留源点的指定严格见证含整个 (G_i)，半径位于 ((R_i,r_i])。原 (P_j) 的分位尾在 (Q_i) 外，与 (G_i\subset Q_i) 不交。这两项是单临界容量反证的必要输入。

本轮选 (b=2^{n+2}\delta)。删除小半径点损失 (O_n(e_i^n)ar_i^n)，分位删除损失至多 (2nbm_i)。因此 (m_i^*\ge(1-C_n\delta)ar_i^n>0)，当 (delta) 小于固定维数阈值时一致成立。

## 2. 完整矩形端点的移动必须按真实 (d_i) 控制

设 (s_i=1-m_i/(ar_i^n)in[0,1])。由 (g_i\ge m_i)，

\[
 \frac{R_i}{r_i}\ge(1-s_i)^{1/n}\ge1-s_i,
 \qquad 0\le r_i-R_i\le s_i r_i.
\]

完整矩形 (X_i) 在一个坐标上的左、右端点相对于 (B_i) 的位移分别为

\[
 \ell_{ik}+\frac{r_i-R_i}{2},\qquad
 h_{ik}-\frac{r_i-R_i}{2}.
\]

由于 (|\ell_{ik}|,|h_{ik}|\le c_ne_i r_i)，每个端点的绝对位移至多

\[
 (c_ne_i+s_i/2)r_i\le C_1d_i r_i.
 \tag{R2}
\]

因此稿件的局部 (d_i/k) 计费成立，不能在这里仅引用统一 (delta/k) 界；后者会损失本轮改进。

## 3. 图边在全部节点上强迫半径可比

图的外盒 (Q(x_i,(1/2+C_1\delta)r_i)) 对**全部**指标定义，包括 (d_i) 可能很大的例外指标。图边相交及严格分离直接给

\[
 \max(r_i,r_j)<(1/2+C_1\delta)(r_i+r_j),
\]

故

\[
 \frac{\max(r_i,r_j)}{\min(r_i,r_j)}
 \le\frac{1+2C_1\delta}{1-2C_1\delta}
 \le e^{C_2\delta},\qquad C_2=8C_1
 \tag{R3}
\]

对固定维数的充分小 (delta) 成立。这个结论不要求两个端点规则，所以路径经过例外盒不会破坏尺度控制。

图度有维数上界：一阶邻居半径与 (r_i) 可比，其原半盒不交且包含于 (Q(x_i,C_nr_i))；每个体积至少 (c_nr_i^n)，体积打包给出固定上界。这里没有遗漏巨大或微小盒，因为真正的图邻居已经被 (R3) 排除。

## 4. 多项式增长只在 (k\delta\le c_0) 的有限图球内使用

设 (k\delta\le c_0=1/(4C_2))。沿长度 (m\le k) 的路径，由 (R3)，每个半径与根半径之比在

\[
 [e^{-C_2k\delta},e^{C_2k\delta}]subset[1/2,2]
\]

内。每步中心距离至多 (3r_j)，所以到根 (j) 图距离不超过 (k) 的全部原半盒均在 (Q(x_j,4kr_j)) 内。利用所有原半盒不交，得到加权图球体积界

\[
 \sum_{\operatorname{dist}(i,j)\le k}r_i^n
 \le (8k)^n r_j^n.
 \tag{R4}
\]

这是直接的体积界，不需要先把度的指数增长转换成多项式增长。

对每个例外根应用 (R4)，再对根求和，得

\[
 S_{\rm buf}\le C_n k^nS_b
 \le C_n k^n\frac{vS}{\delta}.
 \tag{R5}
\]

不同根重复覆盖节点只会放大右侧，因此无需互斥分配缓冲节点。若某连通分支无例外节点，其中到例外集的距离为无穷，定义 (	heta_i=1)，该分支不产生缓冲损失。图边两端不可能一端距离有限而另一端无限，所以全图上仍有 (|\theta_i-\theta_j|\le1/k)。

稿件已正确写出有限半径条件。它**不能**改成对所有图尺度均有固定维数的多项式增长。反例是固定小 (x=C_1\delta>0)，令 (r_j=(1+x)^j)，中心沿一个坐标轴排列，邻接中心距为 ((1+x/4)r_{j+1})。所有中心满足严格 Besicovitch 分离；相邻外盒相交，因为

\[
 (1/2+x)(2+x)-(1+x/4)(1+x)=\frac54x+\frac34x^2>0.
\]

根的 (m) 步图球包含 (0,\ldots,m)，其加权体积随 (m) 指数增长。这个反例不触及当前稿件，因为当前只使用 (m\delta\le c_0)。

## 5. 插值体积及缓冲成本

逐坐标端点插值确等于 Minkowski 凸组合

\[
 Y_i(\theta)=(1-\theta)B_i+\theta X_i.
\]

规则节点的第 (l) 边长是

\[
 (1-\theta)r_i+\theta R_i+\theta(h_{il}-\ell_{il})
 \ge R_i+\theta(h_{il}-\ell_{il}).
\]

乘积的常数项至少为 (g_i/a)，一次项至少为

\[
 \theta nR_i^{n-1}(m_i^*)^{1/n}
 \ge\theta(n-C_n\delta)v r_i^n.
\]

例外节点的 (	heta_i=0,Y_i=B_i)，体积仍至少为 (g_i/a)。所有 (	heta_i<1) 的节点都在缓冲集合中，因此

\[
 \sum_i|Y_i|\ge M/a+nvS
 -C_n\bigl(v\delta+v^2k^n/\delta\bigr)S.
 \tag{R6}
\]

此外 (B_i,X_i\subset A_i)，而 (A_i) 凸，所以 (Y_i\subset A_i\subset D)。无需证明插值矩形包含原半盒。

## 6. 单临界分类在实际插值相交下仍有同轴分离

对规则对，假设实际 (Y_i\cap Y_j\ne\varnothing)，令 (ho=d_i+d_j\le2\delta)、(s=r_i+r_j)。由 (R2)，

\[
 \max(r_i,r_j)<\|x_i-x_j\|_\infty\le(1/2+C_1\rho)s,
 \qquad |r_i-r_j|\le2C_1\rho s.
 \tag{R7}
\]

它来自实际 (Y) 相交，未假定完整 (X) 相交。取 (C_3=4(C_1+c_n+1))，临界轴按中心差至少 ((1/2-C_3\rho)s) 定义。

至少一个轴临界。非临界轴对全部 (x\in P_i,z\in P_j) 有

\[
 |x_l-z_l|<[1/2-(C_3-c_n)\rho]s
 <[1/2-(C_1+1)\rho]s\le\min(R_i,R_j).
 \tag{R8}
\]

最后一步用 (R_i\ge(1-d_i)r_i\ge(1-\rho)r_i) 及 (R7)。在临界轴 (q) 上，适当反射后，所有 (P_j) 点均在所有 (P_i) 点的右方，只需额外要求

\[
 (C_3+c_n)\rho<1/2.
 \tag{R9}
\]

这是固定维数的小量阈值，因 (ho\le2\delta\to0) 一致满足。

若恰有一个临界轴 (q)，则其他中心坐标差均小于 (s/2\le\max(r_i,r_j))。严格中心分离只能发生在 (q)，所以原 (B_i,B_j) 在同一个 (q) 轴上按中心序分离。

完整 (X_i,X_j) 在同轴也分离。设 (D_q=x_{j,q}-x_{i,q}+\ell_{jq}-h_{iq})。若 (D_q<R_i)，取 (P_i^*) 中坐标充分接近 (h_{iq}) 的源点，其指定见证含 (G_i) 以及原 (P_j) 的整个左尾。尾的包含用 (R8)、(R9) 和 (D_q<R_i) 的严格余量；尾质量至少为 (bm_j)，并与 (G_i) 不交。因此

\[
 g_i+bm_j<ar_i^n.
\]

但 (ar_i^n-g_i\le d_i ar_i^n\le\delta ar_i^n)，而

\[
 bm_j\ge2^{n+2}\delta(1-\delta)ar_j^n
 >\delta ar_i^n
\]

由半径比 (r_j/r_i>1/2)、(delta<1/2) 成立。矛盾。交换源与目标，再用原 (P_i) 的右尾，得 (D_q\ge R_j)。因此完整矩形在 (q) 轴上的有向间隙至少为

\[
 D_q-(R_i+R_j)/2\ge0.
\]

这个证明不需要完整 (X_i,X_j) 有交点；稿件已明确这一点。分位尾取自原 (P_i)，而源点取自删去小半径后的 (P_i^*)，端点无需取到。

## 7. 不同参数的错位确实只需 (d_j/k)

对上述单临界轴，原半盒和完整矩形都按相同顺序分离。因此对任意共同参数 (	heta\in[0,1])，两个对向端点的间隙是两个非负间隙的凸组合，仍非负。

令共同参数先取为 (	heta_i)。把 (Y_j(\theta_i)) 改为实际 (Y_j(\theta_j))，该对向端点最多移动

\[
 |\theta_i-\theta_j|C_1d_jr_j\le C_1d_jr_j/k.
\]

所以单临界轴的交长至多这个量；其他轴的交长为 (O_n(r_j))，得到

\[
 |Y_i\cap Y_j|\le C_nd_jr_j^n/k.
 \tag{R10}
\]

不论对每个无序对选择哪一端记账，每个节点承受的次数不超过固定图度。因此总损失至多

\[
 \frac{C_n}{k}\sum_i d_i r_i^n\le C_n vS/k.
\]

若一端 (j) 例外，则 (	heta_j=0)，相邻规则端有 (	heta_i\le1/k)。原半盒在至少一个轴上分离，规则端相对于 (B_i) 的端点位移至多 (C_1d_i r_i/k)。故该类对同样按规则端的实际 (d_i) 计费，求和仍为 (O_n(vS/k))。没有使用不存在的“完整规则 (X_i) 与例外 (B_j) 同轴分离”。两端例外时原半盒不交。

## 8. B10 的增强平方计数成立

每一临界轴的投影交长至多 ((C_1+C_3)\rho s)。至少两轴临界时，交体积至多 (C_n\rho^2r_i^n)。对规则对，

\[
 (d_i+d_j)^2\le2\delta(d_i+d_j).
\]

又因实际相交对是图边、半径比至多 2，

\[
 (d_i+d_j)r_i^n\le d_ir_i^n+2^nd_jr_j^n.
\]

将这一不等式对无序相交对求和，每个端点出现至多固定图度次，得到

\[
 \sum_{i<j,\,\rm reg}(d_i+d_j)^2r_i^n
 \le C_n\delta\sum_i d_ir_i^n
 \le C_n\delta vS.
 \tag{R11}
\]

所以 B10 没有漏掉一个 (S\)、图球大小或例外权重因子。这里只用一阶图度，而非 (k) 步缓冲体积。

## 9. 参数平衡及所有输入的结算

有限族二项 Bonferroni 下界、(R6)、(R10)、(R11) 给出

\[
 V\ge M/a+nvS
 -C_n\left(v\delta+v^2k^n/\delta+v/k\right)S.
\]

令 (kappa=1/(n+2))、(delta=v^\kappa)、(k=\lfloor c_0/(2\delta)\rfloor)。当 (v) 足够小时，

\[
 \frac{c_0}{4\delta}\le k\le\frac{c_0}{2\delta},\qquad k\delta\le c_0/2.
\]

三项分别为

\[
 v\delta=v^{1+\kappa},\quad
 v^2k^n/\delta=O_n(v^2\delta^{-n-1})=O_n(v^{1+\kappa}),\quad
 v/k=O_n(v\delta).
\]

局部图半径可比只要求 (k\delta\) 足够小，不要求趋零。这里 (c_0) 的固定选择已足够，不能误把增强参数 (delta\asymp k^{-1}) 视为不符合图球证明。

由 (S/|E|=1+O_n(v))、(V/|E|\le1+O_n(t))、(t=O_n(v^2))，得到

\[
 M/|E|\le a-nav+O_n(av^{1+\kappa}).
\]

近饱和之外的分支已有更强上界。常数及所有起效阈值只依赖固定 (n)，故可先对紧 (K) 作内正则逼近，再对所有有限正测度 (E) 取上确界。周期下界误差 (O_n(av^2)) 小于本轮余项。因此

\[
 C_n(1-t)=1+2^nt-n2^{n+1}t^{1+1/n}
 +O_n\!\left(t^{1+1/n+1/[n(n+2)]}\right).
\]

粗备核参数也正确：(gamma=(n+2)/(2n+3),\kappa_0=1/(2n+3)) 满足

\[
 2\gamma=1+\kappa_0,quad
 2-n\kappa_0-\gamma=1+\kappa_0,quad
 \gamma-\kappa_0>0,quad \gamma>1/2.
\]

因此即便不用 B10，也可独立得到较弱的正幂余项；当前稿件的 B10 已通过，无需停在备核值。

## 10. 必须保留的条件与文字建议

没有待补的实质接口，但正式稿应保留以下具体条件：

- 图球增长只用于 (k\delta\le c_0)；所有顶点都进入同一外盒图，包含例外节点。
- 参数 (delta) 还须小于固定维数阈值，以保证原局部化、完整见证、分位保留正质量及 (R9)。这些都由 (delta=v^\kappa\to0) 同时满足，常数无输入依赖。
- 单临界分离是同一轴、同一中心顺序的分离，不能仅引用完整矩形不交。
- 错位按真实 (d_i) 计费，缓冲遗失才按图球体积 (k^nS_b) 计费。
- B10 的临界阈值按每对 (ho=d_i+d_j) 定义；若换成统一 (delta)，只能恢复较粗备核指数。

两点纯文字建议：第 4 节“质心”宜改为已定义的“中心”；第 6 节“非负的二项容斥”宜改为“有限族的二项容斥下界”，因为二项截断本身未必非负，但所用下界总成立。这两处不影响现有数学推导。

## 11. 本机 Mathematica 的实际执行记录

工具为 `/usr/local/bin/wolframscript`，版本实际输出 `15.0.0 for Mac OS X ARM (64-bit) (May 15, 2026)`。以下只作为代数子步骤核验，几何与一般维数的求和证明均已在上文给出。

已实际执行并得到：

```wl
kap = 1/(n+2);
Together[{2-(n+1)kap-(1+kap),
          1+kap-(1+1/(n+2)),
          1+(1+kap)/n-(1+1/n+1/(n(n+2)))}]
(* {0,0,0} *)

Reduce[0<x<1/4 && 4/(1-4x^2)>8, {x}, Reals]
(* False: 图边对数半径比的导数界 *)

Reduce[cc>=1 && ee>=0 && 0<=q<=uu<=1 &&
 Abs[ll]<=cc ee &&
 Abs[ll+(1-uu)/2]>(cc+1)(ee+1-q),
 {cc,ee,q,uu,ll}, Reals]
(* False: 单个左端点局部移动界 *)

Reduce[cc>=1 && c1>=1 && rho>0 &&
 c3==4(c1+cc+1) &&
 (1/2-(c3-cc)rho)>=(1/2-(c1+1)rho),
 {cc,c1,c3,rho}, Reals]
(* False: 非临界切向严格余量 *)

FullSimplify[Exp[1/4]<2]
(* True *)

Reduce[n>=2 && (kap<=0 || kap>=1 || 2<=1+kap), {n}, Reals]
(* False *)
```

还实际核验了粗备核两项幂平衡，返回 `{0,0}`，以及共同参数下两个非负间隙的凸组合不可能为负，返回 `False`。一次试图直接对符号实指数 (u^{nn}\ge q) 作 `Reduce`，系统返回 `Reduce::nsmet`，未计为证明或通过；所需结论在第 2 节用 (q^{1/n}\ge q) 的初等单调性严格证明。

**最终判定：空间平滑增强候选 B1 的新接口可通过，幂增益为 (1/[n(n+2)])。** 它仍未达到第三尺度 (av^2)，不提供第三系数。

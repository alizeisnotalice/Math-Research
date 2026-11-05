# R061：八阶竞争族中高支撑 gross 交量的 \(43/5\) 次幂率

2026-11-05。固定维数 \(n\ge3\)。新局部引理，供独立审核。仓库、论文、前轮冻结源和其他代理文件只读。

本稿把 R050 的全高支撑 \(o(v^8S)\) 定量成
\[
 \boxed{\displaystyle \sum_{e:q(e)\ge3}|X_i\cap X_j|
                =O_{n,B}(v^{43/5}S)
                =O_{n,B}(v^{8+3/5}S).}
\tag{1}
\]
它不声称完整 \(C_n\) 余项有这个幂率，不展开九阶包络，也不声称 \(43/5\) 是真实最优指数。唯一随 \(v\) 代入的量是已经证明的 Young **标量比较参数** \(\eta\)；所有源、CDF、正常集合、管、窗口及几何长尾阈值先固定不动。

## 1. 实际序列、固定构造和已审定量输入

记 \(a=v^n=2^nt\)、\(A_v=(1+v)^n\)，及
\[
 P_7(t)=\frac{A_v}{A_v-a}
       +K_nav^6-(n+26/9)K_nav^7,\qquad
 K_n=\binom n2 n^2/27.
\]
固定 \(B<\infty\)，沿任意实际有限正测度输入和严格 halo 外部紧集的序列，假设
\[
 v\downarrow0,\qquad 1+M/V\ge P_7(t)-Bav^8.
\tag{2}
\]
该范围包含当前精确八阶附近的竞争输入，但不要求某个八阶系数已被预设。

按 R049 的顺序一次构造同一 \(I\) 的新固定八层族，使用实际源 \(F_i\subset G_i\)、实际根
\(R_i=(|G_i|/a)^{1/n}\)、原始／有效／中点盒和每点原指定严格见证。不同 \(F_j\) 两两不交并与不同指标的完整 \(G_i\) 排他。原相交度、相交尺度可比、实际帽截宽及全部非临界原宽捕获保持。所有 pair 均取外部实际族，正体积无向 pair 只计一次，\(S=\sum_i r_i^n\) 仍取原全选族。

先固定一个正常集合、旧实际正常层 \(Y_{i,N}\)、R043 的固定窗口
\(R_*=100(n+1)\)、稳健面裕量 \(\tau_*=1/128\)，以及 R045 最后大管宽 \(K\)。这些量及八层统一 CDF **不随后文 \(\eta\) 改变**。记
\[
 \mathscr D_i=(x_i+W_i^0)\setminus F_i,\qquad
 R_N=\frac va\sum_i|\mathscr D_i\setminus Y_{i,N}|,\quad
 R_{\rm all}=\frac va\sum_i|\mathscr D_i\setminus Y_{i,\rm all}|,
\]
\[
 \mathfrak E=\sum_iR_i^{n-1}\Delta_i^{\rm end}.
\]
R049 已审稳定性、正常层逐源体积及 R043 坏根界给固定常数 \(C_*\)，使
\[
 R_N+R_{\rm all}+\mathfrak E+\mathcal F+\mathcal F_{\tau_*}
       \le C_*v^8S,\qquad
 D_N\le C_*v^2S,\qquad
 |Y_{i,N}|\le Ca v^5R_i^n,
\tag{3}
\]
\[
 S(\mathcal B)\le C(D_N+v^{9/2}S)\le C_*'v^2S.
\tag{4}
\]
\(\mathcal B\) 是 R043 的真正缺槽／坏根集合，不另并入一个无速率 LT 异常集。另有 \(S(\partial I)=O(v^5S)\)。

管外真实缺量
\[
 Z(K)=\frac va\sum_i
   |[(x_i+W_i')\setminus F_i]\setminus\mathcal C_i(K)|
       \le R_{\rm all}\le C_*v^8S.
\tag{5}
\]
这里 \(\mathcal C_i(K)\) 是整个固定两近面管。\(C_*\) 只使用 R049 证明中先固定的一次预算份额 \(\eta_0>0\)。下文新的 \(\eta\) 乃至 \(\eta=v^{3/5}\) 不重选这个稳定构造或其常数。

## 2. 除坏所有者长窗外，所有误差已有显式幂率

为核查 (1) 没有保留一个旧 \(o(v^8)\)，把 R050 的四部分显式列出。

### 2.1 目标 epsilon 的局部真实分解

实际源质量比较和原始盒包含给
\[
 \varepsilon_iR_i^n
 \le C|\mathscr D_i|/a
 \le Cv^5R_i^n+
                 C|\mathscr D_i\setminus Y_{i,N}|/a.
\tag{6}
\]
第一步来自
\(\varepsilon_i|F_i|=|W_i|-|F_i|\le|\mathscr D_i|\)、
\(|G_i|/|F_i|\le2\)；不倒置原粗 epsilon 上界。

对任意实际源指标集 \(A\)，相交有限度和根可比给
\[
 \sum_{\operatorname{Adj}(A)}\varepsilon_jR_j^n
 \le Cv^5S(A)+C R_N/v.
\]
于是
\[
 \sum_{\operatorname{Adj}(\mathcal B)}
                    \varepsilon_jR_j^n=O(v^7S),\qquad
 \sum_{\operatorname{Adj}(\partial I)}
                    \varepsilon_jR_j^n=O(v^7S).
\tag{7}
\]
它们分别用 \(S(\mathcal B)=O(v^2S)\) 和
\(S(\partial I)=O(v^5S)\)。大逐源残缺仍全由总 \(R_N\) 计入，不要求其分布均匀。

### 2.2 好所有者的整个固定窗口

对深处 \(q\)-临界 pair，所有者按实际最大质量根选择。
若所有者 \(i\notin\mathcal B\)，R043 的实际面对路径给
\[
 I_e\le C[v^{3q}r_i^n+
             v^{q-1}\mathcal E_i r_i^{n-1}],\quad q\ge3.
\]
固定窗口重数和 (3) 给
\(\sum_{i\notin\mathcal B}\mathcal E_i r_i^{n-1}
 \le C(\mathcal F_{\tau_*}+\mathfrak E)=O(v^8S)\)。
因此全部好所有者高支撑 gross 为
\[
 O(v^9S).
\tag{8}
\]
此处没有给 LT 增加收敛速度；原固定布局只须最终进入一个固定几何邻域，之后用精确实际端点式。

### 2.3 坏所有者的质量／面对例外和短窗

对深处坏所有者 \(i\in\mathcal B\)，按 R044 写
\[
 0<p_1\le\cdots\le p_q\le C_pv,\quad
 p_l=L_l/r_i,\quad A=\prod_{l=1}^q p_l,\quad
 P=\prod_{l=2}^q p_l,\quad H=\min(p_1,A/v^q).
\tag{9}
\]
有 \(I_e\le C_IAr_i^n\)。目标例外
\(\varepsilon_j\ge c_0A/v^q\) 的总量由 (7) 为
\(O(v^{q+7}S)\)。面对／端移例外
\(\alpha_i+\Delta_i^{\rm end}/r_i\ge c_1H\)
的总量不超过
\(Cv^{q-1}(\mathcal F+\mathfrak E)=O(v^{q+7}S)\)。
故全部 \(q\ge3\) 例外为
\[
 O(v^{10}S).
\tag{10}
\]
参数 \(c_0,c_1\) 是原几何固定常数，独立于 \(\eta\)。

固定 \(K\) 后，R044 只选一次
\(T_K>2K/c_2+1\)。若 \(p_2\le T_Kv^3\)，排序给
\(I_e\le C_Kv^{q+4}r_i^n\)。仅对坏所有者求和，使用 (4)：
\[
 \sum_{\substack{\text{坏所有者}\\p_2\le T_Kv^3,\ q\ge3}}I_e
 =O_K(v^9S).
\tag{11}
\]
特别四临界短窗为 \(v^8S(\mathcal B)=O(v^{10}S)\)，不能沿用旧全根的 \(O(v^8S)\)。

### 2.4 外部边界

原单源界对高支撑边界 pair 给
\[
 I_e\le Cv^q[
       \varepsilon_i+\varepsilon_j+v^{q/(q-1)}]r_i^n.
\]
由 (7)，epsilon 部分为 \(O(v^{q+7}S)\)；纯幂部分为
\(O(v^{q+q/(q-1)+5}S)\)。\(q\ge3\) 时最小指数为
\(19/2\)，故全部高支撑边界 gross 为
\[
 O(v^{19/2}S)+O(v^{10}S).
\tag{12}
\]
只用原 \(S(\partial I)=O(v^5S)\)，没有假设更强的标记／边界预算。

## 3. 实际长窗矩形与完全显式的 \(\eta\) 常数

对 (10) 两类例外之外、\(p_2>T_Kv^3\) 的坏所有者边，R044 在原有效源盒内给真实矩形 \(\mathcal R_e\)。它的临界面深分别为
\([c_3H,2c_3H]r_i\) 和
\([c_2p_l,2c_2p_l]r_i\) (\(2\le l\le q\))，其他坐标距两面均 \(>Kv^3r_i\)。

所有者的原指定严格见证包含完整 \(G_i\)，并在每条临界轴捕获实际目标帽
\(\min(p_lr_i/2,w^0_{jl})\ge\gamma p_lr_i\)；非临界目标原宽完整进入。扣除非例外的实际目标缺量，得到固定正倍数的额外 \(F_j\) 质量。与真正浅面的点态容量矛盾，故整条确为实际缺集，且避开整个 \(\mathcal C_i(K)\)。这组几何常数 \(c_0,\ldots,c_3,\gamma\) 及截宽比例不依赖 \(\eta\)。

令
\[
 X_e=\frac va|\mathcal R_e|.
\]
对有限的 \(3\le q\le n\)，可统一选固定 \(c_X,C_I,D>0\)，使
\[
 X_e\ge c_Xr_i^n\frac{HP}{v^{q-1}},\qquad
 I_e\le C_IAr_i^n,\qquad
 \sum_{\substack{\mathrm{owner}(e)=i\\\text{全部相关 }q}}X_e
            \le DZ_i(K).
\tag{13}
\]
最后一个 \(D\) 是跨所有支撑数的共同重数，不为每个 \(q\) 重新分一整份 \(Z_i\)。

取
\[
 c_\lambda=\min\{1,c_X/(DC_I)\}>0,\qquad
 \lambda=c_\lambda\eta,\quad 0<\eta\le1.
\tag{14}
\]
当 \(v^{q-1}<\lambda\) 时，原解析两分支给
\[
 A\le\lambda\,\frac{HP}{v^{q-1}}
          +\lambda^{-q/(q-1)}v^{\beta_q},\qquad
 \beta_q=\frac{q(2q-1)}{q-1}.
\tag{15}
\]
其证明对每个实 \(\eta>0\) 是同一代数不等式：若
\(P\ge v^q\)，第一项已超过 \(A\)；若 \(P<v^q\)，正净额要求
\(P<v^{2q-1}/\lambda\)。排序给
\(p_1^{q-1}\le P\)，从而
\(A\le P^{q/(q-1)}\le
 \lambda^{-q/(q-1)}v^{\beta_q}\)。
这里没有一个未明示的 \(C_\eta\)。

乘 \(C_Ir_i^n\)，结合 (13)–(14)，得到
\[
 I_e\le\frac{\eta}{D}X_e
       +C_Ic_\lambda^{-q/(q-1)}
                       \eta^{-q/(q-1)}v^{\beta_q}r_i^n.
\]
有限 \(3\le q\le n\)、\(c_\lambda\le1\) 允许把前面的固定常数统一成 \(C_n\)，所以全部坏所有者非例外长边满足
\[
 \boxed{\displaystyle
 \sum I_e\le
 \eta Z(K)+C_n\sum_{q=3}^n
       \eta^{-q/(q-1)}v^{\beta_q}S(\mathcal B),\quad
       v^2<c_\lambda\eta,\quad0<\eta\le1.}
\tag{16}
\]
因为 \(v\le1\) 且 \(q\ge3\)，共同条件 \(v^2<c_\lambda\eta\) 蕴含每个 (15) 的 \(v^{q-1}<\lambda\)。式 (16) 的常数、几何有效阈值和 \(T_K\) 皆独立于 \(\eta\)；仅显式保留的前提和幂因子含 \(\eta\)。

## 4. 合法代入 \(\eta=v^{3/5}\)

现在在**标量不等式 (16)** 中代入
\[
 \eta(v)=v^{3/5}.
\tag{17}
\]
不改变任何实际源或划分。原几何阈值保持一个固定 \(v_{\rm geom}>0\)；标量前提变成
\[
 v^{7/5}<c_\lambda,
\]
故只须再取
\(v<\min\{1,v_{\rm geom},(c_\lambda/2)^{5/7}\}\)。
不存在隐藏的“先对固定 \(\eta\) 取未知阈值”步骤。若还希望 \(\eta\) 小于旧预留正系数，它也在这个同一固定小 \(v\) 范围最终成立。

由 (3)、(4)，(16) 右侧预算项为
\[
 \eta Z(K)\le C_*v^{43/5}S.
\]
三临界 Young 余项为
\[
 \eta^{-3/2}v^{15/2}S(\mathcal B)
 \le Cv^{-9/10+15/2+2}S=Cv^{43/5}S.
\tag{18}
\]
对于 \(q\ge4\)，必须仍使用坏所有者权的 \(v^2\)。其指数为
\[
 \beta_q+2-\frac{3q}{5(q-1)}
   =2q+2+\frac{2q}{5(q-1)}
   \ge\frac{158}{15}>\frac{43}{5}.
\tag{19}
\]
最小值在 \(q=4\)，可由导数
\(2-2/[5(q-1)^2]>0\) 直接核查。若错误地对高支撑长尾使用全根权而漏掉这一个 \(v^2\)，(19) 就不是这里的结论。

合并 (8)、(10)–(12)、(18)–(19)，其余最低指数为 9 或 \(19/2\)，均高于 \(43/5\)。有限维数中 \(q=3,\ldots,n\) 求和不改变阶数，得到 (1)。

## 5. 实际预算和量词

本证明使用的是 R049 已证的**数值事实**
\(Z(K)\le C_*v^8S\)，再以 R044 对实际矩形的代数估计控制交量。它不在另一份源价格中重减完整 \(J_i\)，也不把 \(\eta Z\) 当作新收入。若同源价格需要保留正份额，因 \(\eta(v)\to0\)，上述 gross 上界更不要求重新分配一份固定领先源成本。

所有 \(O\) 常数只依固定维数、竞争常数 \(B\)、预先固定有限层及正常／窗口／CDF 参数，不依有限输入数目、原半径尺度跨度或随后测试的 \(\eta\)。R049 的常数始终使用固定 \(\eta_0\)；(17) 不是把稳定性证明重新执行成一个变参数构造。LT 仍只作固定身份识别，不用其收敛速度估计幂余项。

结论可按原临界族的序列口径理解：沿任意 (2) 序列，固定上述工具后对充分小 \(v\) 有同一固定常数的 (1)。原已证输入供货和正量估计允许同一常数依 \(n,B\) 及固定工具，未知的序列收敛模量只影响几何身份最终生效的阈值，不进入右侧幂常数。

二维没有 \(q\ge3\) 的 pair，空和为零；这不解决二维盒内小源的八阶或更高阶问题。本稿也没有把 (1) 直接升级成完整 \(C_n\) 的 \(O(av^{8+3/5})\) 余项：正常积分、宽度／形状价格、二轴补集的长角／近零阈值和全局账本中的其他余项仍需各自给速率。

## 6. 最小实际计算与依赖

Terra 已在本机 Mathematica 15.0 实际运行原 r059-computation/tail_exponents.wl/log，进程 PID 32510。本目录 computation/ 保存该成功脚本与日志的逐字副本，不重新运行。其最小输出核查：

* \(q=3\) 的两项指数 \(8+\kappa\)、\(19/2-3\kappa/2\) 在
  \(\kappa=3/5\) 时同为 \(43/5\)；
* 对实参数 \(q\ge4\)，(19) 低于 \(43/5\) 或 \(158/15\) 的违反集均为 False，\(q=4\) 恰为 \(158/15\)；
* 标量前提的幂差 \(2-3/5=7/5>0\)。

R044 既有 young_branches.wl/log 的实际缓存核了原两分支；未为本轮重跑旧几何或有限代数。一般 \(q\)、显式 \(\eta\) 常数及合法代入由第 3–4 节解析证明，计算不替代实际目标截宽、原指定见证、矩形重数或源构造。

数学供货：

* R049 competitive_eighth_stability.md，SHA 5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4，及其正式独审；
* R050 high_support_gross_after_stability.md，SHA 22f848879b3c9ba839edb5d307a205ae3245e5ba3a89a99f7ec2ce2974010748，及已解除 R049 条件的独审；
* R043 normal_slots_and_triple_tail.md，SHA 75d0f89e79e8de4d00a6d3926532357f08a60a1ed5cbbc20b6abfd83e82ab4ac；
* R044 high_support_net.md，SHA ef3a6df6ee4581e990fd55fbe76c2b03f43b2e5bdbd0cb27f7eba18231cf0cc4，尤其 (18)–(24) 的统一成本、共同重数和 Young 两分支。

本证明只使用第 1 节列明的定量稳定性、固定窗口实际面对关系和第 3 节的实际矩形成本接口，因此可作为将来新固定末层的供货规格。**本稿没有供应新 \(N=9\) 的对象或先验**，不能将旧第八层改名后直接套用。当前八层全局账本仍有第六前缀 \(\theta^2o(v^6S)=o(v^8S)\) 等无速率误差；它们未由 (1) 改善，须由新的同族递推／前缀结算另行处理。

状态：本稿给出 (1) 的完整实际局部证明，待独立审核；不修改任何旧冻结源，也不宣称完整 Tauberian 下一余项幂率。

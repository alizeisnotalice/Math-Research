# Actual 门与辅助退出：精确行 LP、阈值与已补足的逐点桥

2026-10-07。应用 L03 有理 primal/dual 证书工作流。独占本稿及 `far_gate_overlap_guard_20261007` 前缀脚本和数据；未改旧几何稿或总账。已执行三轮有限概率行的有理核验，没有执行 actual FIRST 实验。

**当前结论。** 只用完整行 future-cap 的平均退出信息时，最优界为 \(\int ge\,dr\ge(t-q/u)_+\)。该平均信息单独不能在 \(49/4096\) 吸收预算内取得统一正比例门耦合。但是原特定核还有严格逐点 \(e>1/5\)，已经补足全部 actual residual 的门耦合；不能继续将抽象行 LP 的零重叠解释为 actual 门桥的缺口。真正未付的是保留完整实际门的来源一次空间列预算。

## 1. 原实际门如何进入行变量

在同一 residual 输出 x，保留完整 hard posterior \(r_x\)、原 soft 子后验、全部历史/种子与实际密度门。对每个 hard z 积分原 soft realization 及种子，可写

\[
g_x(z)=\mathbb E_\vartheta\int A_f(x,\text{soft realization},z,\vartheta)
\,d\Pi_x,\qquad 0\le g_x\le1.
\]

原 \(\Pi_x\) 是子概率，门未重新归一；此定义来自原 actual 正表示。令

\[
t(x)=\int g_x\,dr_x,\qquad
j(x)=\int g_x e_x\,dr_x,\qquad \rho(x)=q/u(x).
\tag{1}
\]

同源交通分别是 \(R_g=\int u t\,dx\)、辅助退出交通 \(R_{ge}=\int u j\,dx\)。辅助 hard 噪声是在保留全部原历史之后独立附加的正核，不声称它是原 soft first-exit continuation。完整 future cap 给 \(\int e_xdr_x\ge1-\rho\)，而原 hardband 保证 \(3/16\le\rho<1/2\)。以下 LP 先只保留这些平均约束。

## 2. 最优有限行 LP 与完整证书

对任何概率 r、\(0\le g,e\le1\)、\(\int gdr=t\)、\(\int edr\ge1-\rho\)，逐点有 \(ge\ge g+e-1\) 和 \(ge\ge0\)。因此

\[
\boxed{\int ge\,dr\ge\max(0,t-\rho).}
\tag{2}
\]

该结论不需要独立性。任意此类行可以精确送入四个 binary sites \((g,e)=(0,0),(1,0),(0,1),(1,1)\)，概率分别定义为 \(\int(1-g)(1-e),\int g(1-e),\int(1-g)e,\int ge\)。于是原均值与目标完整保留，得到 LP

\[
\min p_{11},\quad p\ge0,\quad\sum p=1,\quad
p_{10}+p_{11}=t,\quad p_{01}+p_{11}\ge1-\rho.
\tag{3}
\]

令 \(J=(t-\rho)_+\)，显式最优原始行是

\[
(p_{00},p_{10},p_{01},p_{11})
=(\rho-t+J,\ t-J,\ 1-\rho-J,\ J).
\tag{4}
\]

四项非负、总和一、gate 均值为 t、exit 均值恰为 \(1-\rho\)。对偶写为 \(a+b g+c e\le ge\)、\(c\ge0\)，目标 \(a+bt+c(1-\rho)\)。当 \(t\ge\rho\) 取 \((a,b,c)=(-1,1,1)\)，四 site slacks 为 (1,0,0,0)；当 \(t\le\rho\) 取零对偶。二者与 (4) 目标完全相等，认证了 (2) 最优。

若要求严格 \(\int e>1-\rho\) 而没有定量 slack，同一值仍是最优下确界。特别是取任意 \(0<\delta<\rho\)、\(t=\rho-\delta\)，概率 \((0,t,1-\rho+\delta,0)\) 给严格平均退出，但 overlap 仍为零。以上只是平均信息外放松行，不是原核或 actual FIRST 的反例。

## 3. 只用平均 cap 时的最佳阈值与系数

### 固定 gate mass 阈值

在 \(t\le T\) 上丢掉交通，最优普适 band 费用为

\[
R_{\rm low}\le4T\lambda|E|.
\tag{5}
\]

这个系数在仅保留行约束时可由 \(u=4\lambda,t=T\) 达到。剩余 \(t>T\) 上，由 (2) 得到的最佳统一 relative overlap 为

\[
c(T)=\begin{cases}
0,&0<T\le1/2,\\
1-1/(2T),&1/2<T<1.
\end{cases}
\tag{6}
\]

\(\rho<1/2\) 的严格性不改善统一常数，因为可以从下方趋近 1/2。要 \(c>0\)，必须 \(T>1/2\)，故删除系数 \(4T>2\)。反过来，固定要求 \(0<c<1/2\)，最佳阈值 \(T=1/[2(1-c)]\)，对应 \(\varepsilon_{\rm low}=2/(1-c)\)。\(c\ge1/2\) 没有覆盖全部参数的非空固定-threshold 高支；取 T=1 只是把全部行删掉，费用 4。

在当前宽裕 ceiling \(\varepsilon_*=49/4096\) 内，(5) 只允许

\[
T\le49/16384<1/2,
\]

因此 c(T)=0。若已用现有两项 \(49/8192\) 吸收，新增吸收还须严格小于 \(49/8192\)，相应 \(T<49/32768\)，仍不能靠这条平均信息得到正比例耦合。

### 随 q/u 改变阈值

令 \(\kappa=q/\lambda\in[3/4,1]\)，要求保留的每一行有 \(j\ge c t\)、\(0<c<1\)。基于 (2)，最小可行门阈值是

\[
t\ge\frac{\rho}{1-c}.
\tag{7}
\]

阈值大于一的行只能全部删除。低支最优 coefficient 因而是

\[
\boxed{\varepsilon_{\rm low}(\kappa,c)
=\min\{4,\kappa/(1-c)\}.}
\tag{8}
\]

这里极值取 \(u=4\lambda\)，并令 t 从阈值下方趋近；(4) 给出配对见证。若阈值≤1，系数由 \(ut=q/(1-c)\) 给出。任意 c>0 都至少要 \(\varepsilon_{\rm low}>3/4\)，远超过当前 ceiling。因此 ratio-adaptive 阈值也不能修复平均信息分层。

**零耦合端点要分开。** c=0 时根本无需删除任何行，最优费用是零；\(\kappa/(1-c)\to\kappa\) 是 c↓0 的正耦合阈值成本极限，不是 c=0 的最优费用。结果 JSON 中 c=0 的三条 cutoff 记录只核了这条阈值曲线的端点代数，其 `minimal_deleted_coefficient` 标签在这个退化端点须按此修正读取；它们不参与任何正耦合或不可吸收结论。其它 c>0 记录与 (8) 一致。

### 不分层的 additive 形式

只依靠 (2)，对任意有限 \(C\ge1\)，最优普适不等式是

\[
ut\le C u j+q,
\quad\text{故}\quad R_g\le C R_{ge}+\lambda|E|.
\tag{9}
\]

最佳 additive coefficient 仍为一：选 \(q=\lambda,u=4\lambda,t=\rho=1/4,j=0\) 即达到。对 \(0\le C<1\)，最佳 coefficient 为 \(4-3C\)，由 \(t=1,j=3/4\) 达到。严格平均信息的无定量 slack 版本取相应趋近行得到同一个 sharp supremum。这些最优性仅针对被保留的平均信息。

## 4. 原核的逐点信息已经补足 actual 门耦合

root 与 gated agent 提供并在 `auxiliary_exit_column_bridge_20261007.md` 证明了更强的核级事实。独立核对如下。原

\[
w(r)=\int_0^1e^{-|r|/s}ds,\qquad \phi=h*w
\]

有质量一，偶递减性给 \(\phi\le\phi(0)\)。而

\[
1-\phi(0)=\int_0^1 2s e^{-1/(2s)}ds
\ge\frac{3}{4e}>\frac14
\]

来自 s∈[1/2,1] 的积分与 e<3。对 hard cube 中的每个 z，原张量核逐坐标得到

\[
e_{\sigma,R}(x,z)\ge1-(1-\sigma/4)^n
\ge\frac{n\sigma}{4+n\sigma}.
\tag{10}
\]

第二步由 \((1-v)^{-n}\ge1+nv\) 给出。所有 actual residual 已有 \(\sigma>1/n\)，故逐点 \(e>1/5\)，不只覆盖大 σ。于是保留任意原可测 gate 和全部 soft 历史，直接有

\[
\boxed{j\ge t/5,\qquad R_g\le5R_{ge}.}
\tag{11}
\]

这一步不必删低 t，也不花新吸收预算。LP 极值行的 e=0 支与此特定核逐点信息冲突，因此不得用于否定 (11)，也不得把第 3 节称为 actual 门桥仍未解决。这里的 1/5 是选用的简洁核常数，没有声称它是原核的最优常数。

对 root 短稿的只读独审也通过：σ=2/n 时 (10) 给严格 \(e>1/3\)；实际 n≥512、σ≤1 的范围合法。自由辅助退出核 \(h_R e_{\sigma,R}\) 的 R-sup column 夹在 \((1/5)N_h\) 与 \(N_h\) 之间，\(N_h=1+n\log(b/a)\)；σ=2/n 的下界可取 \((1/3)N_h\)。在 b=2a 是解析的 Θ(n)。这是删除 FIRST/cap/门后的核 column，不是 source/FIRST 合法坏输入，也未把 hard 来源免费替换成 soft first-exit 的共同 2W 参考来源。

因此 (11) 完成耦合，却没有付 \(R_{ge}\) 的移动空间列。必须保留实际 complete future cap、出生、森林、分数及所有原门来争取下一步来源一次预算；直接扔门做 free Schur 已无法达到 √n。

## 5. 三轮精确核验与范围

脚本 `far_gate_overlap_guard_20261007.py` 用 Fraction 检查 (3) 的 primal/dual，并独立枚举 LP 的所有 active-set vertices：两个 equality 加 Ee lower boundary/四 nonnegative faces 中的两个 active inequalities，精确有理 Gauss 消元后逐项验证可行性，再比较最小值。没有浮点 LP solver 或随机抽样。

| 轮次 | 均匀 t 网格分母 | ρ 个数 | 含 kink/预算节点的行数 |
|---|---:|---:|---:|
| 1 | 16 | 3 | 57 |
| 2 | 64 | 7 | 478 |
| 3 | 256 | 11 | 2861 |

共 3396 行，全数 primal/dual gap 与独立枚举误差精确为零。ρ 覆盖 3/16 到 511/1024，均在原 band 平均信息范围。各轮额外包括 t=ρ、附近节点及两个吸收 threshold。已有数据保存完整四点原始行、对偶、slacks 与独立枚举 minimizer，不覆盖原轮文件。

另核 (8) 的正 c 代数、(9) 的 sharp additive 行、严格平均退出仍零 overlap 的趋近行，以及三个 \(n=512,4096,32768\)、\(\sigma=(D+1)/(Dn)>1/n\) 参数的 (10) 系数严格大于 1/5。PASS 判定全部是有理断言，运行时间只是诊断浮点字段。三轮结果分别保存于同前缀 `_round1/2/3.json`，总收据在 `_results.json`。

冻结脚本中 strict-exit 注释使用了 `strict actual exit lower bound does not rescue`，措辞过强；它只指 **abstract relaxation 的 strict row-average exit lower bound** 没有定量 slack，绝不指实际原核逐点 (10)。实际 pointwise 信息正好 rescue 门耦合，已在第 4 节及结果 `information_layer` 字段明确。脚本及其收据 hash 保留，没有偷偷修改冻结版本。

本轮认证的是 **信息层级区别**：平均 cap 的 LP 界已最优，但实际原核具有更强逐点信息。没有 actual FIRST 反例、没有新 source 空间费用、没有阶数拟合，未闭合一般主目标。

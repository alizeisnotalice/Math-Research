# 九阶收支的未闭合项与三维 inside 容量常数一

本稿研究真实中心方体问题的九阶收支。新的无条件结论是 §2 的可测局部盒引理；§3 将其接到明确列出的同层真实输入合同，得到三维 inside 小源的领先价格 (2v^9S)。这并不确定真实 (C_3) 的九阶总系数。R071 尚在终审时，涉及其新 N9 的整体结论均保留该依赖；局部引理不依赖 R071 的结论。

## 0. 已实际读取的输入与范围

以下路径均相对于 `/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/`；列出的是实际文件 SHA256。

| 输入 | SHA256 | 本稿使用的内容 |
|---|---|---|
| `scratch/r071-recentered-bootstrap/fixed_tube_recentered_settlement.md` | `a014ae167a2d83a30fb2b455d827a2b13253fce5033f9223ca9c674956e6d5f2` | 新 N9 原地质量链、固定 whole-normal 管、移动极值与正剩余；不借终审状态代替这些推导 |
| `scratch/r063-eighth-rebuild/reconstruction.md` | `16401f6d8b9688d94ae115b2d5df1974dc36accbc404e102675e479de1e11a33` | 原有限值见证、全点选择、whole-chart 帽、原/有效端点区别、实际槽与分母 |
| `scratch/r060-third-feedback-lower/third_feedback_and_micro.md` | `c7bcd28338661fa782029c87c883467d7ea6dd812d3b16301a02c74598cc9ccc` | 真算子构造；三维实际微孔增益 (2K_3av^9=2av^9)；有限族比较与总系数的范围区别 |
| `scratch/r053-normal-eighth-packet/actual_normal_eighth_price.md` | `d83dc772fc2c1ca5698a39eb2bd15702395d5e8238b11209fd2b6e7765b1e642` | §3 两个 pair 排序的真实整横宽积分；非 pair 最浅轴项及 clip 的阶 |
| `scratch/r056-sharp-small-source/sharp_inside_cost.md` | `18fa60b96a6b2a2be7092d06dc59b4741a049e4214941323a33f846dad12ffc0` | 二维全点选择/锚定方法；本稿重新证明三维 edge 管情形，不跨维引用角引理 |

固定 (n)，(a=v^n)，(S=\sum r_i^n)。局部终层不等式的九阶单位是 (v^9S)；代回真实质量比才是 (av^9)。物理删除质量的价格为 ([1+O(v)]v/a)。下文所有阶数均明确区分这两种单位。

## 1. 定量稳定性尚未给出完整有限簇泛函

若 R071 的实际链通过，可用的是
\[
 \mathsf M+R_A^-+v^6D_N+L_N=O(v^9S),\qquad Z(K)=O(v^9S),
 \tag{1.1}
\]
以及新 N9 上重新证明的宽度 RMS、旧数值余项。这里 (L_N) 是围绕移动极值的亏损；旧固定中心 \(\mathfrak G_N\) 一般内禀为八阶，不能把它替换成九阶，更不能推出每条正常边偏离仅 (O(v^{3/2}))。

在九阶候选泛函中至少必须逐一处理以下项。表内“存活”是现有全输入估计未能使其成为 (o(v^9S))，不是已证明每一项都能达到正九阶极值。

| 类别 | 九阶可能来源及所需新信息 |
|---|---|
| 同一 pair 的真实反馈 | leading (v^6) 乘三次 (O(v)) 反馈；需真正算子三次代入的有界余项。单纯包络展开不能替代实际前缀帽积分。 |
| 不同 pair 的交叉退款 | R063 中另一 normal pair/角在目标帽内的占比首次足以产生九阶；需带目标并集的加权可加性，不能把每个退款独立重复相减。 |
| 最浅轴不在 pair | 八阶三近面名义成本的下一项及其首次退款为九阶；这里有真实 free-axis 图，不能只保留平面模型。 |
| 顶部 fibre clip | (v^6) 主量乘横向相对 (v^3) collar；必须计算或证明其对所需输入族为小量。n2 无 free axis 的这项为空。 |
| 根差/不对称帽 | 现有余项 (C(v^7z^2+v^8|z|+v^9)) 不能在九阶直接丢弃；需保持双向不同外帽的真实截断和端点。 |
| 槽与移动偏差 | (v^6D_N) 和 (L_N) 本身仅 (O(v^9S))；罕见非正常参数可与正常 RMS 并存。需保留其正罚，或证明相关最优松弛的严格局部比较。 |
| 非正常 packet 与 raw | \(\mathsf M,R_A^-,Z(K)\) 仅九阶；需证明任何九阶收入可由同一份费用严格压住，或纳入局部状态。不能从 (O(v^9)) 推出不贡献系数。 |
| 高支撑接触 | q≥3 的 net 目前是九阶上包；坏 owner 的小总权不能单独排除九阶集中。需真正 q3 帽与 joint jump 的局部估计。 |
| 几何密度与符号 | 三维 signed off-diagonal shear 的 triangle determinant density 是独立候选。反射能共轭某些 forest 交叉项，不会自动消掉 triangle 的符号不变量。R072 正在证明，本文不预签该候选。 |
| 小源 | n3 的 inside leading 首次为九阶，本文 §2–3 锁定其价格上界 (2K_3=2)。n2 的 leading 为八阶，其下一个微孔/宏孔反馈修正仍是九阶，需要 R073 另算。n≥4 此种基础微源收益晚于九阶。 |
| 盒外、边界与端锚 | 已有 \(CvF\)、\(Cv\mathfrak E\) 可为九阶。边界小源用同 I 的缺陷下界与旧数值正预算，在 n3 已是 \(O(v^{11}S)\)，不再列为九阶自由机制。环境 \(Ca(F+P_\partial)\) 在 n≥2 晚于九阶；盒外费用仍需保留。 |
| 前缀与真实分母 | N9 的第八前缀粗 (O(v^8)) 误差乘一次 (v) 就是九阶。精确九阶需在同一实际族升级前缀/消项，或增加固定层数并原地重放。分母和删除导数必须沿实际 (M,T) 传递。 |

这里边界指标始终是原 \(\partial I=\{j\notin I:\exists i\in I,\ Q_i\cap Q_j\ne\varnothing\}\)。实际 `paper/sections/global.tex` 的 `eq:global-boundary-deficit-floor` 给 \(d_i\ge c_nv\)，所以 \(S(\partial I)\le C_nP_\partial/v=O(v^7S)\)；`eq:global-small-source-boundary` 更直接给 \((v/a)\sum_{\partial I}\delta_i\le C_{n,A}v^nP_\partial\)。n3 因而是 \(O(v^{11}S)\)，对应真实比值 \(O(av^{11})\)。`paper/sections/sharp.tex` 在原同 I 上明确使用该边界定义和缺陷下界；R071 的 \(P_\partial=O(v^8S)\) 改善的是数值预算，不改变 I。此处仅排除一个粗估计留下的假候选，不更动深处局部引理。

宽度误差本身已有安全间隔：(v^6\sum\xi_iR_i^n=O(v^{19/2}S)=o(v^9S))。这并不使上表其他项消失。

还可严格指出一个逻辑障碍：仅给
\(\sum d_e^2R_e^n=O(v^3S)\)，允许总权 (v^3S) 的边具有固定非零偏差，其余边在移动极值上。这恰产生 (v^6\sum d_e^2R_e^n\asymp v^9S)。故从此 RMS 不能推出“全部边可统一作三阶小偏差 Taylor”。这是稳定性条件的抽象反模型，不是新的真实 halo 反例。

因此 (1.1) 尚不足以推出全输入九阶都等于独立 pair/triple 最大值之和。有限局部簇方案至少还需：统一的局部截断/加权可加引理、罕见退化参数的收费、joint 目标同一次 J、及全输入九阶前缀接口。需要的簇也不能只按“反馈次数”命名：三次有向传播在一般邻接图上可访问四个顶点，必须先证明跨 pair 传播的额外阶数，才能缩成较小簇。

下面给出可以独立闭合的一笔费用。

## 2. 三维双薄 edge 管的局部常数一引理

### 2.1 陈述

设
\[
 R=\prod_{j=1}^3[0,\ell_j],\qquad cL\le\ell_j\le CL,
\]
其中 (c,C,A>0) 固定。令 (E_h\subset R) 为 R 的十二条 edge 的双薄管之并：每条管在自由轴取整段，另外两轴距各自指定端点不超过 h。设 (T,Z\subset R) 可测，且 (Z\subset E_h)。写
\[
 D=R\setminus T,\quad z=|Z\setminus T|=|D\cap Z|,
 \quad d=|D\setminus Z|.
\]
设 (|D|\le bL^3)，其中 (b\to0)，且 (h/L\to0)。给定可测 (\mathcal S\subset T)，每个 (x\in\mathcal S) 保留一个取有限个正值的原半径 (0<r(x)\le AL)，满足真实严格不等式
\[
 |Q(x,r(x))\cap T|<a r(x)^3,\qquad Q(x,r)=x+[-r,r]^3,
 \qquad a\to0.
 \tag{2.1}
\]
则参数充分小时，存在只依赖 (c,C,A) 的常数 (C_0)，使
\[
 \boxed{\quad |\mathcal S|\le {a\over1-a}z+C_0a(d+h^3).\quad}
 \tag{2.2}
\]

不要求 (b^{1/3}L=o(h))；不要求 Z 是缺集，也不要求任何微见证或锚定半盒包含于 Z。实际输入在 R 外任意。由全输入的 (T\subset E^c) 及原严格见证即可得到 (2.1)，即使见证跨出 R。

### 2.2 全点有限选择与一个粗费用界

对有界点集及有限正值原半径，按半径从大到小，从当前尚未被已选**闭**原方体覆盖的全部点中选中心。新中心不在先选闭方体内，所以
\[
 \|x_i-x_j\|_\infty>\max(r_i,r_j).                 \tag{2.3}
\]
由于最小半径严格为正，所选中心按该最小半径 packing，选择必有限终止；因此原 Q 覆盖整个被处理点集，不是假设某个预先 finite probe 的中心已覆盖即足够。

同一符号 \(\sigma\in\{\pm1\}^3\) 的锚定盒
\[
 P(x,r;\sigma)=\prod_j[x_j,x_j+\sigma_jr]
\]
按区间无序端点理解。其体积为 (r^3)。由 (2.3)，同符号的 P 内部两两不交：实现中心距离的那条轴上，同向长度各不超过两半径最大值。

另有一个不要求 (|D|/L^3\to0) 的粗结论。若仅 (r\le AL)、(2.1) 成立，则
\[
 |\mathcal S|\le C_{c,A}a|D|.                    \tag{2.4}
\]
事实上，原 Q 的重叠数至多八：固定一点 y，将含 y 的各 Q 按中心相对 y 的三个坐标符号分类；同一符号中两个中心的每轴距离均不超过半径最大值，与 (2.3) 矛盾。又因中心在 R，
\[
 |Q\cap R|\ge \gamma r^3,\quad \gamma=\min(1,c/A)^3,
\]
从而 (a\le\gamma/2) 时 (|Q\cap D|\ge(\gamma/2)r^3)。全点 cover 给
\( |\mathcal S|\le a\sum r_i^3\le(16a/\gamma)|D|\)。

### 2.3 总缺量先排除宏观半径

若 (r\ge\ell_{\min}/2)，中心在 R 给
\(|Q\cap R|\ge(cL/2)^3\)。于是
\[
 |Q\cap T|\ge(c^3/8-b)L^3>aA^3L^3
\]
对小参数违反 (2.1)。所以全部原半径 (r<\ell_{\min}/2)。每轴向盒中心方向至少有半条边可用，故以下八种 inward-sign 分组的 P 均包含于 R。也有
\(r^3\le |D|/(1-a)\)，但不把这个粗半径估计误用为 (r\ll h)。

### 2.4 大微半径只用管外残量

取固定 (K=8)。先对 (r>Kh) 的源点按三坐标各处于盒中点的哪侧分八组，P 指向盒中心。每组独立作全点有限选择。每条 edge 管与一个 P 的交量至多 (h^2r)，故
\[
 |P\cap Z|\le12h^2r\le{3\over16}r^3.
\]
又 (P\subset Q\cap R)，(2.1) 给
\[
 |P\cap(D\setminus Z)|\ge(1-a-3/16)r^3\ge\tfrac12r^3
\]
（取 (a\le1/4)）。同组 P 不交，因此大半径全体费用至多 (16ad)。这里不是仅凭 tube 体积分数推断输入形状；真实严格见证迫使缺量出现，减去 tube 最多占量以后才将剩余部分收费到实际 (D\setminus Z)。

### 2.5 小微半径的主 edge-half 组

余下 (r\le Kh)。若中心到 (E_h) 的 \(\ell^\infty\) 距离大于 Kh，其 Q 避开 Z；八 inward-sign 组直接给费用至多 (8ad/(1-a))。

其余中心在 (E_{(K+1)h})。固定 (H=10(K+1))。对每条 edge，记自由坐标为 u；把以下中心列为例外：它在该 edge 的 ((K+1)h) 管内，且
\[
 u\in[0,Hh]\cup[\ell_u-Hh,\ell_u]
       \cup[\ell_u/2-Hh,\ell_u/2+Hh].             \tag{2.5}
\]
十二个此类集合之并记 B。其 Kh 邻域（与 R 相交）记 \(B^+\)，由有限个边长为常数倍 h 的长方盒覆盖，故
\[
 |B^+|\le C_Kh^3.                              \tag{2.6}
\]
所有这里涉及的 midplane、盒边界及平局均可按固定编号归组；取闭 collars 不损估计。

非例外的管内中心仅属于唯一一条 edge：若同时靠近不同自由轴的两条 edge，则三个坐标都近端点，已在 (2.5) 的端点段内。按该 edge 及自由轴位于中点哪侧分组，共至多二十四组。P 在两条临界轴都指向盒内，自由轴指向盒中点。

这些 P 不仅同组不交，而且跨主组也物理不交：

1. 同一 edge 的两个 half 组，其自由坐标距中点至少 ((H-K)h)，两侧区间分开。
2. 不同平行 edge 至少有一条临界轴对应相反端点；P 距该端点不超过 ((2K+1)h)，而 (h/L\to0)，所以分开。
3. 不同自由轴的两条 edge，若相交可能性未由相反面排除，则一条 edge 的自由坐标是另一条的临界坐标。前者的 P 距对应端点至少 ((H-K)h)，后者至多 ((2K+1)h)。由于 (H-K=82>17=2K+1)，仍分开。

同组使用 (2.3)。所有主 P 因而构成一个整体不交族。它们仅需包含于 R，不需落入 host 或 Z。真实见证给
\( |P\cap D|\ge(1-a)r^3\)，从全点 cover 得
\[
 |\mathcal S_{\rm main}|\le {a\over1-a}|D|
                         ={a\over1-a}(z+d).       \tag{2.7}
\]
正常 z 在此只收一次。

### 2.6 例外只付 (h^3)

中心在 B 且 (r\le Kh) 的 Q、P 全部位于 (B^+\)（P 和 Q 均只在 R 内取质量）。仍按八 inward-sign 组，从该组全部实际源点选择。于是
\[
 |\mathcal S_{\rm exc}|
 \le {8a\over1-a}|D\cap B^+|
 \le {8a\over1-a}(d+C_Kh^3).                    \tag{2.8}
\]
这同时处理三个近面的 vertex junction、自由轴中面以及相邻 edge 的潜在重叠。没有假设正常层沿 edge 具有均匀高度；哪怕高度退化或任意可测，只用其物理支撑。

将大半径、远管、主组、例外四类费用相加，得到 (2.2)。整个证明允许同一个原 Q 穿越宏盒边界；可计算的有限值半径也没有被缩小或重选。

## 3. 同层实际三维 inside 接口

### 3.1 显式输入合同

固定 n=3、\(a=v^3\)。在同一实际有限层的深处源，记原 CDF 盒 \(B_i^0=x_i+W_i^0\)，原扩展盒 \(\mathcal R_i\supset B_i^0\)，原小源 \(T_i\cap\mathcal R_i\)，以及已删除源 \(F_i\subset T_i\)。取 R071 §1 的同层实际旧正常保证层 \(Y_i^{old}\)，并固定该正常目标族的 whole-normal 名义零次层 \(B_N(0)_i\)。这里唯一使用的上包选择为
\[
 Y_i^{old}\subset B_N(0)_i,\qquad
 Y_{i,N}:=B_N(0)_i\cap\mathcal R_i.
\]
第一包含来自名义完整目标帽减去实际目标缺量、共同侵蚀及截断只会降低保证高度；不能任取一个未必包含 \(Y_i^{old}\) 的后续迭代层。小源保留原半径，\(r(x)\le Avr_i\)，且由 \(T_i\subset E^c\) 有
\[
 |Q(x,r(x))\cap(T_i\cap\mathcal R_i)|<a r(x)^3.
 \tag{3.1}
\]
实际旧正常 raw、原入向扩张及自由余量提供
\[
 D_i^0:=B_i^0\setminus F_i,\qquad
 q_i:=|D_i^0\setminus Y_i^{old}|+|\mathcal R_i\setminus B_i^0|,
 \qquad
 \sum_iq_i\le(a/v)R_*+CvF=O(v^9S).             \tag{3.2}
\]
其中仅用 (R_*,F=O(v^8S))，无需把旧 raw 宣称成九阶。 这一具体选择还给出
\[
 |(\mathcal R_i\setminus T_i)\setminus Y_{i,N}|\le q_i,
 \qquad |D_i^0|\le |Y_{i,N}|+q_i.
\]
前者使用 \(F_i\subset T_i\) 和 \(B_i^0\subset\mathcal R_i\)，后者使用 \(Y_i^{old}\subset Y_{i,N}\)。因此好盒 q 小不仅控制删除前的 T，还控制原已删除 F 的实际缺率和 CDF 端移；仅假设前一个较弱不等式不足以推出这个结论。

还用以下几何事实；它们必须在所用新层原地成立，不能只改旧层编号：

* 宏盒三边 (\asymp L_i=vr_i)，根 (R_i\asymp r_i)。
* 每源正常层体积 \(|Y_{i,N}|\le C av^5r_i^3\)。其 whole-chart 正高度只可能有两个临界近面深 \(O(v^3r_i)\)，包括最浅轴不在 target pair 的图；没有漏掉非 pair 图。
* 好盒原/有效端点的 CDF 位移及原扩张按实际质量估计，而不是视为相同端点。
* 正常参数的旧加权平方偏差 \(O(v^2S)\)、整横宽误差 \(\sum\xi_iR_i^3=O(v^{7/2}S)\)，以及实际正常槽/根桥成立。

前两条是 R071 §1、§4 的 whole-normal 物理上包；(3.2) 与原端移来自 R063 的实际源合同和入向扩张。这些条件不依赖九阶总系数已知。

### 3.2 分好坏，不要求每盒误差小

定义坏盒为
\[
 q_i>av^5r_i^3=v^5L_i^3.
\]
由 (3.2)，
\[
 \sum_{\rm bad}r_i^3=O(vS).                    \tag{3.3}
\]
对坏盒用 §2.2 的粗式，得
\[
 \sum_{\rm bad}\delta_{i,\rm in}
 \le Ca\sum_{\rm bad}(|Y_{i,N}\setminus T_i|+q_i).
\]
定价后为
\[
 {v\over a}\sum_{\rm bad}\delta_{i,\rm in}
 \le Cv\,[av^5\,O(vS)+O(v^9S)]
 =O(v^{10}S).                                  \tag{3.4}
\]

好盒上 \(|\mathcal R_i\setminus T_i|\le Cv^5L_i^3\)，所以局部 b 可取 (Cv^5\)。实际缺量给逐盒原缺率 \(\varepsilon_i\le Cv^5\)，强 CDF 端移至多 (Cv^5L_i\)。原扩张每侧也至多
\(Cq_i/L_i^2\le Cv^5L_i\)。因此将有效 source 的正常 tube 转回实际 \(\mathcal R_i\) 后，仍由宽度
\[
 h_i=C_Nv^3r_i=C_Nv^2L_i
\]
的十二条物理 edge 管覆盖。这里的常数 (C_N\) 固定，且是 whole-normal 常数。局部引理不要求粗半径 (v^{5/3}L_i\) 比 h 小；其大半径分支已经单独证明。

取 (Z=Y_{i,N}\)，必要时先与 \(\mathcal R_i\) 相交，§2 给
\[
 \delta_{i,\rm in}\le {a\over1-a}|Y_{i,N}\setminus T_i|
                +Ca(q_i+h_i^3).
\]
由于 \(\sum h_i^3=O(v^9S)\)，对好盒求和、加坏盒得到
\[
 \boxed{\quad {v\over a}\sum_i\delta_{i,\rm in}
 \le {v\over1-a}\sum_i|Y_{i,N}\setminus T_i|
       +O(v^{10}S).\quad}                      \tag{3.5}
\]
例外的物理体积 (h_i^3\sim v^9r_i^3) 先乘 a 才是 inside 删除，再乘 (v/a) 才是 ledger 价格：结果是 (v^{10}r_i^3)，不是九阶。相应真实比值误差为 (O(av^{10})=o(av^9))。

### 3.3 正常物理层的领先体积与常数 2

为防把实际软槽容量上界反向当成无权边计数，这里列出所需桥。记
\(U_N=\sum_{e\in N}R_e^3\)，(b_e\) 为外横向 overlap 长度，\(\widehat b_e=b_e/[(1+v)R_e]\)。正常缺口强制性给
\(\sum(1-\min(1,\widehat b_e))R_e^3=O(v^2S)\)。于是实际 R045 槽桥
\[
 W_N=\sum b_eR_e^2\le\mathcal C+O(v^3S),\qquad
 \mathcal C=3(1+v)\sum_iR_i^3+O(v^{7/2}S)
\]
确实推出
\[
 U_N\le3\sum_iR_i^3+O(v^2S).                    \tag{3.6}
\]
这一步用的是 normal overlap 近满的正缺口，不能仅凭软权 \(\le1\) 得到。

R053 §3 两个排序的两个有向端点，除去 \(\ell=s_vv/a\)，在同根整横宽处给每条无向边的物理体积
\[
 av^5R_e^3\,{p_eq_e(p_e+q_e)\over3},\qquad
 p_e=q_e={1\over1+v}
\]
（n=3 的移动值只差 (O(v)\)，对此 leading 比较足够）。一般实际正常参数与此值的加权 L1 差是 (O(vS)\)：这是旧平方偏差用有限度及 Cauchy 的结果，绝不逐边假定 (O(v)\)。实际整横宽误差为 \(O(av^5\sum\xi_iR_i^3)\)，根比误差更小。最浅轴不在 pair 的名义层仅 \(O(av^7S)\)，clip 与真实反馈只会减少名义上包。因而
\[
 \sum_i|Y_{i,N}|
 \le\left({2\over3}+O(v)\right)av^5U_N+O(av^6S)
 \le (2+O(v))av^5S,                            \tag{3.7}
\]
这里最后用 \(\sum R_i^3\le S\) 的实际根定义，若采用有 (1+O(v)\) 的版本也只改余项。

将 (3.7) 代入 (3.5)，因 (a=v^3)，得到
\[
 \boxed{\quad [1+O(v)]{v\over a}\sum_i\delta_{i,\rm in}
       \le (2+O(v))v^9S.\quad}                 \tag{3.8}
\]
左侧的 (1+O(v)\) 是删除质量真实基线导数相对于 (v/a\) 的比值。在实际区间 \(\mu/g\in[1-A_Nv,1]\) 上直接微分
\(f(m)=(R+m^{1/3})^3\) 得到此比值；不得把原来的泛常数 C 留在 (3.8) 前面。未删除的 (d_h^0\) 仍走原 (f_h\)/\(\mathcal D_h\) 消项，不能把 (3.8) 另算一次 J。

这锁定了 n3 inside leading 的上价格 (2K_3=2\)。R060 已有真实有限构造给微孔相对于其宏观族的同阶增益 (2av^9+o(av^9)\)。二者所说的是同一尺度的一笔微源预算，不据此声称 R060 的宏观族已经取到完整九阶最优值，也不直接与 R072 的 signed triangle 候选相加。

## 4. 对有限簇方案的具体进展与仍需供货

本稿排除了一个原先自由的九阶常数：在 §3 的真实同层合同下，三维微源 inside leading 至多 2，且 junction、端移、自由轴中面均不增添另一笔九阶费用。该结论只用 tube 支撑及真实严格见证，不要求可疑的统一非退化 host 轮廓。

其他机制仍需独立工作：

1. 对完整正常图证明九阶加权局部展开，特别是另一 pair 的首次退款、三近面 source/target 图及顶部 clip。必须对并集和同一个 J 操作。
2. 给罕见非正常/slot 亏损和 raw 的严格局部收入-费用比较。仅有 (1.1) 无法把这些状态删去。
3. 对 n3 的 signed triangle 密度证明实际有限构造及全输入上界的对应不变量；本文不替代 R072。
4. 对 n2 微孔保护的 (O(v)\) 损失计算 leading，连同宏观第三反馈的同一设计质量根重算；本文不替代 R073。
5. 在新的固定足够深实际递归上升级所需前缀，使末层误差 (o(v^9S)\)，并保持原 I、原见证、实际质量分母。没有这一项，局部有效泛函即使算出也尚非 (C_n\) 的九阶系数。

因此本稿不提供“未知 sup 形式”的九阶可计算公式，也不宣称任意阶几何已闭合。现在可交独立复核的是 §2 的无条件引理与 §3 的逐项接口。

## 5. 实际计算及其证明范围

本目录 `finite_checks.py` 已用本机 Python 的 `fractions.Fraction` 实际执行，结果在 `finite_checks.json`。它核了：大半径常数 (12/8^2=3/16\)；主组分离常数 (82>17\)；2016 个 common-orthant 有理点比较；n2..7 的 canonical pair 体积乘实际 leading 槽数等于 (2K_n\)；三维所有 residual/collar/bad-box 的原量纲幂次。所有断言通过。

这些有限核只检查常数和索引，不证明连续全点 cover 或真实 halo 几何。后两者分别由 §2 的实数论证和 §3 明列的实际同层合同承担。本文没有实际执行新的 halo 全局数值优化，也没有将数值证据冒充全输入上界。

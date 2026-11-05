# R074 非作者独立审核：三维 edge 管常数一与实际 inside 价格

审核者：Sol，独立任务 `/root/sol_r057_global_audit`。2026-11-05。
只写本独占目录；作者源、既有证明及正文只读。本文是本次完整审核，不以此前 R057/R064/R071 的裁决标签代替数学。

## 1. 最终裁决及精确范围

**PROOF-VALID，范围限于下列两项。**

1. 冻结候选 §2 的无条件局部定理：三维可比边长轴盒，任意可测双薄 edge 管内的 Z，原有限正值半径和真实严格见证，给出
   \[
   |\mathcal S|\le \frac{a}{1-a}|Z\setminus T|+C a\bigl(|(R\setminus T)\setminus Z|+h^3\bigr).
   \]
   “无条件”指完整满足其明确的局部假设；不需要额外 host 非退化、Z 等于缺集、或半径比 h 小。
2. 冻结候选 §3 在明确的同层真实合同下的接合是有效的；这些合同在已实际读取的本次 N9/R071 对象上有实际供应。因此该 N9 深处 inside 删除在真实基线中满足
   \[
   [1+O(v)]\frac va\sum_i\delta_{i,\mathrm{in}}
       \le (2+O(v))v^9S \qquad(n=3,a=v^3).
   \]
   junction、自由轴中面、端移及坏宏盒的附加价为 O(v10 S)。真实 halo 比值中相应误差为 O(a v10)。

这不是全输入九阶总系数，也不是有限簇全局变分公式。正常第三反馈、交叉 pair、clip、非正常/raw 的九阶收入费用比较及足够精确的前缀等仍须另证。不能把这里的 2 与某个 triangle 候选直接相加后宣布 C3 的九阶系数。也不自动把本合同移到未经重建的新层数。

最终签署对象为全文实读的

`scratch/r074-ninth-cluster-upper/ninth_mechanisms_and_inside.v3.md`

SHA256 **57514c9109718eb904bb47088fb995b8733ebf302200f269a4287c90b42c6424**，308 LF，21264 B。

原 v1 SHA `88d353c924a9093cb811b447ba7f79594e552e3cb6b10d9aa714b1b8930d9e69`（294 LF，19643 B）也已全文实读。v2 SHA `deb7a9943e387add628e6dc0c53c90217b61eff582cce252799296c953390968` 只收紧边界状态。已实际核两个 diff，并用作者两份 patch 的 new→old 逐块逆替换，精确恢复 v3→v2→v1 原字节。这个逆核由本审代码实跑，并非只信作者的 True。

## 2. 本审发现及 v3 的必要修复

原 v1/v2 的 §3.1 仅列
\(q_i\ge |(\mathcal R_i\setminus T_i)\setminus Y_{i,N}|\)
和总和上界。它单独不能推出下一段所需的原始已删除源 \(F_i\) 的缺率小：\(F_i\subset T_i\)，后者未控制可能位于 T 内的已删点。因而不能仅从宏盒删前缺量小推原始 CDF 端移小。

已在终审前将此缺省告知作者和 root，并给出最小修复。v3 源现明确固定
\[
 D_i^0=B_i^0\setminus F_i,\quad
 q_i=|D_i^0\setminus Y_i^{old}|+|\mathcal R_i\setminus B_i^0|,
 \quad Y_{i,N}=B_N(0)_i\cap\mathcal R_i.
\]
且 \(F_i\subset T_i\)、\(B_i^0\subset\mathcal R_i\)、\(Y_i^{old}\subset B_N(0)_i\)。旧保证层在有效源盒内，故在原始盒内，因而也在 \(Y_{i,N}\) 内。于是确实同时有
\[
 |(\mathcal R_i\setminus T_i)\setminus Y_{i,N}|\le q_i,
 \qquad |D_i^0|\le |Y_{i,N}|+q_i.
\]
旧保证高度是同一名义完整目标帽减去实际目标缺量、侵蚀并截断所得，故其包含于名义零次层。不能将未必包含它的后续迭代 Y2 任意代入。v3 明写了唯一具体选择，本审最终没有用作者口头承诺补足源量词。此项现已解决；§2 完全不受它影响。

## 3. 局部证明的独立重建

我在收到冻结全文前已从原严格见证独立推导此局部机制；之后逐段核对作者实际稿。以下给出关键的真实逻辑链。

**严格见证产生缺量的常数。** n=3 时 a=8t，Q 的体积为 8r3，但小源满足的是 \(|Q\cap T|<a r^3\)。凡 inward 锚定立方 P 满足 \(P\subset Q\cap R\)、\(|P|=r^3\)，就有
\(|P\cap D|>(1-a)r^3\)。这里没有再乘覆盖重数八或将 a 错作 t，也没有要求 P 落入 Z。

**全点选择。** 半径有限正值，按值降序，从整个尚未被所选闭 Q 覆盖的点集选择。不同已选中心满足 \(\|x_i-x_j\|_\infty>\max(r_i,r_j)\)。正的最小半径及有界 R 使 packing 有限，故选择终止并覆盖全部原点集。常数与最小半径和半径列表长度无关。固定符号的锚定 P 不交：若它们相交，每轴中心差都至多两半径最大值，矛盾。对闭边界或零坐标平局固定编号即可；没有依赖不可测的无穷选择。

**不要求小总缺量的粗界。** 对该 greedy 所选 Q，任意公共点处把中心按八个正交符号归类，同符号至多一个，因此重数至多八。中心在 R 且 r≤AL 给 \(|Q\cap R|\ge\gamma r^3\)，\(\gamma=\min(1,c/A)^3\)。当 a≤γ/2，\(|Q\cap D|\ge\gamma r^3/2\)。故 \(|\mathcal S|\le16a|D|/\gamma\)。这可用于坏宏盒，那里不需要局部 b 小或 tube 定位。

**先排除宏观半径。** 若 r≥最短边/2，\(|Q\cap R|\ge(cL/2)^3\)，减去 bL3 后仍远大于 aA3L3，违反原严格见证。故全部 r<最短边/2，按中点向内取符号时 P 真在 R。由此也可得 r3≤|D|/(1−a)，但不能推出 r≪h；作者没有偷用这一错误推论。

**大微半径。** r>8h 时，每条 edge 管与 P 相交至多 h2r，十二条总量≤3r3/16。故 a≤1/4 时 \(|P\cap(D\setminus Z)|\ge9r^3/16\ge r^3/2\)。分八 inward 符号组后只向 residual 收 16ad。这里 tube 的小占比只是从真实见证迫出的 D 中扣掉一部分，并非假定 tube 之外的输入具有某种形状。

**小微半径与例外。** r≤8h 且中心距闭几何 E_h 超过 8h 时 Q 避 Z，八组仅向 d 收费。其余中心在 E_9h。按闭几何管分类，而非按任意可测 Z 的非空投影分类，避免一个与测度无关的错误几何推断。作者取 H=90；每条 edge 的两个端点 collar 及自由轴中点 collar 长 O(h)。它们的 8h 邻域 B+ 与 R 相交的总体积为 O(h3)。例外仍经八 inward 组、原严格见证收费，故是 \(Ca(d+h^3)\)，绝不是只付 \(Ch^3\)。这一额外因子 a 对最终阶数至关重要。

**24 主组物理不交。** 非例外中心唯一属于一条 edge，按 edge 及自由轴半边分 24 组。临界两轴指向盒内，自由轴指向中点。每个 P 的两条临界坐标距指定面≤17h，自由坐标距相关端点≥82h。逐类检查：

- 同 edge 相反 half 由中点 collar 分开；
- 平行不同 edge 由某个相反面分开（h/L→0）；
- 不同自由轴的 edge，某坐标一边自由而另一边临界，82h>17h 分开。

同组由 greedy 分离不交，跨组由上述物理包络分开。因此所有主 P 的 D 质量可统一求和一次，\(|\mathcal S_{main}|\le a(z+d)/(1-a)\)。其余三类只付 Ca d 或 Ca(d+h3)，从而得常数一结论。没有 normal 的四次、十二次或二十四次重复收费。

局部参数次序合法：固定 c,C,A,K,H，之后 a,b,h/L 各自充分小；无 b^(1/3)L/h 趋零的隐含要求。外部输入任意，原 Q 允许越宏盒边界，始终只在 P⊂R 取得质量。

## 4. 实际同 N9 接合及阶数

实际小源是 T_i 中被原半径阈值删除的点，T_i⊂E^c。因此每个原严格 E 见证确给 (3.1)，不是在宏盒内重新选择输入或重新缩半径。层 h 的 F_i 是已删小源后的集合，且 F_i⊂T_i；原大源/目标、I 和有限值半径不变。

这里实际供应限定在 R071 的固定 B 竞争分支 \(1+M/V\ge P_8-Bav^9\)（其先进入 \(P_7-B_0av^8\) 分支）；常数可依固定 B 及既定源参数，不是说任意远离近极值的输入都满足该稳定性。只要另一个同层族逐项供应相同合同，局部接合也有效。本次 R071/R064 的实际源、共同 CDF、old raw 和入向扩张给
\[
 R_*=(v/a)\sum|D_i^0\setminus Y_i^{old}|=O(v^8S),
 \quad F=O(v^8S),\quad
 \sum q_i\le(a/v)R_*+CvF=O(v^9S).
\]
其中物理 raw 实为 O(v10 S)，扩张在此用较粗的 CvF=O(v9 S) 即足够。没有把 R071 新中心化 raw 与旧 R_* 混为同量，也没有宣称旧 R_* 已为 O(v9)。

以 q_i>av5r_i3=v8r_i3 为坏宏盒，则坏原权 O(vS)。每源名义正常层体积≤Ca v5r_i3。用粗容量，坏盒的 ledger 费用≤Cv[av5 O(vS)+O(v9S)]=O(v10S)。这里 O(vS) 足够，不能无依据改写成 O(v2S)。

好盒 q_i≤av5r_i3，同时控制删前宏 D 和原始 D0：\(|D0|\le C av5r_i3\)。实际质量与 a r_i3 可比，强 CDF 给 ε_i≤Cv5，原/有效端移≤Cv5L_i=Cv6r_i。因原始盒每边≥cL_i，扩张体积又控制每侧扩张≤Cq_i/L_i2≤Cv5L_i。因此 effective 正常层回到实际扩张盒的 edge 管，只加 O(v6r_i) 位移，仍被固定宽 \(h_i=C_Nv^3r_i\) 支撑。whole-normal 的 q 不在 pair 图也在其中；不是只审平面短角。

局部定理可统一取 b=Cv5、h_i/L_i=C_Nv2。它并不要求原半径 v^(5/3)L_i 比 h_i 小。求和后
\[
 (v/a)\sum\delta_{i,in}
 \le\frac v{1-a}\sum|Y_{i,N}\setminus T_i|+O(v10S).
\]
因为 v∑q_i=O(v10S)，v∑h_i3=O(v10S)。例外从物理 h3 先乘 a 得删除量，再乘 v/a 定价；不能直接把物理 h3 视为九阶删除价格。

## 5. 领先体积、软槽桥及一次 J

实际整横宽的两 pair 排序在一个有向端点给 PQ(P+Q)/(2f)，其中 n3 的平面容量 f=3v2R，源自由横宽为 vR；两端合计的物理层系数为
\[
 a v^5R_e^3\frac{pq(p+q)}3.
\]
在 p=q=1/(1+v) 的参考点，领先系数为 2/3。实际正常参数只具加权 RMS O(v2S)，但固定窗上的 Lipschitz 与有限度/Cauchy 足够给加权 L1 O(vS)，无需逐边 O(v)。width L1 为 O(v7/2 S)，根比误差更小；非 pair 最浅轴的名义物理层为 O(av7S)。union/clip 和实际退款不会增大名义上包。

不能从软权≤1 反推无权边数。这里源已给正常缺口
\(\sum(1-\min(1,\widehat b_e))R_e^3=O(v2S)\)，故
\[
 U_N\le W_N/(1+v)+O(v2S)
 \le3\sum_iR_i^3+O(v2S)\le3S+O(v2S).
\]
这是同层实际软槽/根桥与近满缺口共同推出的上界，不是模型 lattice 计数。于是 \(\sum|Y_{i,N}|\le(2+O(v))av5S\)，代入前节得 inside 领先 2v9S。

真实删除导数也须保留尖系数。写 m=zv3R3、z∈[1−Av,1]，对 f(m)=(R+m^(1/3))3 直接得到
\[
 \frac{f'(m)}{v/a}=(v+z^{-1/3})^2=1+O_A(v)
\]
一致成立。把真实积分 f(g)−f(μ) 沿原未删质量 m0 分成 g−m0 与 δ=m0−μ 两份；前者保留 d_h^0/f_h/重复量原消项，后者才使用上述 sharp inside 价格。Bonferroni 仍是同一个 \(H=\sum I-\sum J\)，J 始终基于实际 F。既不再添一笔 J，也不对同一删除质量重复付账。最终价为 (2+O(v))v9S，余项 O(v10S)。

## 6. 边界更正与 lower 比较的范围

已实际阅读并核对正文 global.tex 的 `global-small-source-boundary`（749–751）及其证明（806–818），还有 sharp.tex 219–220。原同 I 的边界指标有 d_i≥c_nv，因此 S(∂I)≤CP∂/v；更直接的小源式为
\((v/a)\sum_{∂I}\delta_i\le Cv^nP∂\)。R071 同族旧数值 P∂=O(v8S)，n3 给 O(v11S)。v2/v3 排除边界小源作为九阶自由项是正确收紧。盒外 CvF 与端锚 CvE 仍可为九阶，不能一并删除。

已全文读取 R060 的实际宏/微构造，并核对其中 §8 的保护链：H⊂U4，\(\widetilde U=U4\cup\mathcal B(U4\setminus H)\)，反单调性给真实超解，保护新增量≤Cu|H|，且固定实际总孔质量的参数变化给 n3 的相对增益 2av9+o(av9)。本审采用的对照仅为这笔实际有限族微源预算；它没有识别宏观族的绝对九阶项。R074 源相同地限制了“匹配”的含义，没有借比较定理宣布绝对系数。

## 7. 实际计算及证据边界

本审实际执行了两个独立本机 Python Fraction 窄核：

- `edge_regions.py` / `edge_regions.json`：24 主 edge-half 的独立边长仿射包络，全部 276 对在每边≥256h 的全参数域有正分离；本审使用更小的 collar 常数，最小可证间距 16h。大半径 tube 占比精确 3/16，a≤1/4 后 residual 为 9/16。作者 H=90 的 collar 更宽，完全满足本审连续分类证明。
- `contract_checks.py` / `contract_checks.json`：两次源 patch 逆恢复精确字节；两排序的有理单项式积分给每无向边 pq(p+q)/3，三槽 leading=2；真实 n3 删除导数的精确多项式系数给 (v+z^(-1/3))2。所有断言通过。当前系统 python3 无 SymPy，首次窄核导入未成功；随后仅用标准库 Fraction 完成相同有理多项式核，没有把未运行的符号工具当证据。

本审未重跑作者旧批次；已实际读取作者 finite_checks.py/json/log，脚本 SHA `ad0bf44388043b838faa5aed0945c586b1239d6869597437a7a8f22b519e4dbe`，输出/日志 SHA `eb05b1ab1a1f3a8239d95826588c1a9841eae214ccd56f272ef036e3f829947a`。其 2016 个有理点检查仅是常数/索引辅助，不能证明实数全点 cover。本文 §3 的连续几何论证和实际同层接合不由这些有限核替代。

本审计算 SHA：

| 文件 | SHA256 |
|---|---|
| edge_regions.py | ca72b2e33299cce6d3e5e21719c8735a074de47f70941489e14e09e39c7181f7 |
| edge_regions.json | b125c3005dfa33ee35465f001dc57bf1066f92c0ed1b323f136a1c3a29bfa403 |
| contract_checks.py | a28759eb7d0f84064ebbb3bea846aa031c8e4d8e662c6eac6a28a7c0502818cb |
| contract_checks.json | 613553fadb33c97f11ecbd7389c49ac3a7faf90499bd90b1cc595c716d7f4d2c |

## 8. 实际读取并重新绑定的数学源

下列 scratch 路径以 topic1-tauberian-work 为根。R063/R064/R071 的完整推导已在本独立审核链全文实际读取；本次重新检查所用实际对象、具体 q、whole 管及质量链，并重新核 SHA。不是把它们的旧审核结论作为前提。

| 源 | SHA256 |
|---|---|
| r063-foundation-rebuild/proof_record.md | 702a70cdd735883d5adbe3d22da5acbe950e665c3023457a76416dcc5dd4f572 |
| r063-sixth-rebuild/rebuild.md | 3a5d38fe5a00392af687c09ad59fcf072c552cbe418504af12fc38dc93e3cb2c |
| r063-eighth-rebuild/reconstruction.md | 16401f6d8b9688d94ae115b2d5df1974dc36accbc404e102675e479de1e11a33 |
| r064-power-remainder/power_remainder_r064.md | 80f92dc770da78284a56e56a331b82112ea5912e892a456e968078e412a25699 |
| r071-recentered-bootstrap/fixed_tube_recentered_settlement.md | a014ae167a2d83a30fb2b455d827a2b13253fce5033f9223ca9c674956e6d5f2 |
| r053-normal-eighth-packet/actual_normal_eighth_price.md | d83dc772fc2c1ca5698a39eb2bd15702395d5e8238b11209fd2b6e7765b1e642 |
| r056-sharp-small-source/sharp_inside_cost.md | 18fa60b96a6b2a2be7092d06dc59b4741a049e4214941323a33f846dad12ffc0 |
| r060-third-feedback-lower/third_feedback_and_micro.md | c7bcd28338661fa782029c87c883467d7ea6dd812d3b16301a02c74598cc9ccc |

本次边界回源正文为 `/Users/zhengzhihao/Desktop/T/Math-Research-topic1-tauberian/topic1-tauberian-high-order/paper/sections/`：global.tex SHA `7b6154f74b56e4cb8ea133c99d4c240ceec526004930c510bc29ebcb25b5f10a`；sharp.tex SHA `8eae7fa41e8bc1f9883d51330b01183001295ca0b623da12205c807c6bbe7c39`。本轮只审上述边界节点，不扩签两份正文的全部内容。

最终源 v3 上未发现尚未解决的必要缺口。任何后续文本修改应重新绑定 SHA；本报告不签署先前较弱 q 合同的抽象版本，也不将局部有效性扩大为全局九阶系数结论。

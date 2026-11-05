# R056 二维 inside 小源常数一：独立数学审查

2026-11-05。仅写本审查目录；作者源、旧源及仓库只读。

## 1. 冻结对象和裁决

已全文独立核查 scratch/r056-sharp-small-source/sharp_inside_cost.md 的数学内容。首次收到的源 SHA256 为 **f220484aa097687535ff409a8dba6a337c483ede84bee759f878fadc88c0a058**，15623 字节、285 个 LF，另有五个误写字体命令的 CR 控制字节。

作者随后只将五个 CR 字节改成字面反斜线加 r，修复其后的字体命令。已独立以 read_bytes 核验旧备份、五处修改及精确逆恢复：旧 0-based 字节位置为 [763,817,925,1230,8427]；除这五处 0d→5c72 外字节全同，逆恢复 SHA 恰为上述 f220484a…。

**最终绑定源 SHA256：18fa60b96a6b2a2be7092d06dc59b4741a049e4214941323a33f846dad12ffc0**，15628 字节、285 个 LF、零 CR。**裁决：PROOF-VALID，无必要数学补丁。** 五处修改仅修字体转写，数学内容未变。本报告不将旧 hash 当作新版本签署。

范围固定为 \(n=2\)、\(a=v^2\)、R049 同一新八层竞争输入及原有限 Borel 严格见证函数。R056 提供实际 inside 删除成本的领先 \(8/27\) 及相应原基线接口，不识别正常 packet 的八阶价格，也不完成 \(C_2\) 的八阶极限。

## 2. 原实际半径确为有限正值，而非一般可测半径的假设

实际读取正式 main 的引理 lem:finitegreedy 及其证明、global 开头与未删除源定义。对紧 halo 集，先用固定半径平均的连续性得到有限严格半径列表，再把每点赋为列表中实际可用的最大半径；可用条件为开条件，因而函数 Borel，取有限多个正值。

global 的 \(T_i\subset G_i\subset K'\subset K\) 及后续各层中心邻域都只限制这个同一函数。小源集合按原函数的半径阈值切出，未重算或缩放半径。因此任意本稿局部集合上也保留有限正值和原严格见证；无需假设一般可测正函数有正最小值。

对任意有界点集按半径值降序，从当前整个尚未覆盖点集中逐点选中心。后选点不在任一已选闭 \(Q\) 内，给
\[
 \|x_k-x_l\|_\infty>\max(r_k,r_l).
\]
若还存在未覆盖点就可以继续选；半径至少为本次列表的正最小值，选中中心的小半盒 packing 迫使有限终止。所以获得真正覆盖全部点集的原 \(Q\) 族，而非仅覆盖某个有限中心列表。

应区分两个常数：此有限终止的选中数上界可以依赖最小半径；不存在对任意微半径都统一的选中数上界。本稿只将 packing 用于有限终止，后续价格均来自不交面积和严格缺量，不含这个上界。因此所有最终价格常数确独立于最小半径、半径列表长度和尺度跨度。

## 3. 同向锚块的不交与严格缺量

对一个固定符号，若中心分离由第 \(h\) 坐标实现，则中心距大于两个半径的最大值。两个同向长度分别为 \(r,r'\) 的区间因此不交；反向符号用同一论证。故固定符号的所有
\[
 P(x,r;\sigma)=x+\sigma[0,r]^2
\]
物理不交，即使半径不相同。这不是平移不交原半盒后的推断。

原严格见证是 \(|Q(x,r)\setminus E|<4tr^2=ar^2\)。对本地 \(T\subset E^c\)，即使实际 \(Q\) 跨出本地盒或进入其他源，
\[
 |Q(x,r)\cap T|<ar^2
\]
仍直接成立。若 inward \(P\subset Q\cap\mathcal R\)，则
\[
 |P\cap D|=r^2-|P\cap T|>(1-a)r^2,\quad D=\mathcal R\setminus T.
\]
因此全点覆盖的原 \(Q\) 用于估计源面积，不交 \(P\) 用于估计盒内缺量，两者职责不同且均由同一个原半径产生。无须原 \(Q\subset D\) 或原 \(Q\subset\mathcal R\)。

## 4. 局部盒引理的近角／远角两部分

设边宽在 \([cL,CL]\)，\(|Z|,|D\setminus Z|\le bL^2\)，故 \(|D|\le2bL^2\)。对每个源点，先处理 \(r\ge\ell_{\min}\)：盒内交面积至少 \(\ell_{\min}^2\)，又 \(r\le AL\)，小 \(a,b\) 时违反严格缺量。因此 \(r<\ell_{\min}\)。

每轴裁入盒内的交长至少为 \(r\)，于是
\[
 (1-a)r^2<|D|\le2bL^2,\qquad
 r\le\sqrt{2b/(1-a)}L=o(\epsilon L).
\]
这是每一个原见证的点态结论，不是大多数点或选中点的平均结论；它也排除了原 \(Q\) 跨边界时的宏观半径分支。

近角集合用四个边长 \(2\epsilon L\) 的角块。各角单独全点 greedy，符号均向该角内部。中心在角块内、\(r=o(\epsilon L)\)，所以对应 \(P\) 在本盒内且位于该角的 \(3\epsilon L\) 块。四个块两两不交，故近角部分的所有 \(P\) 跨角也不重复向 \(D\) 收费。原 \(Q\) 覆盖全近角集合，得到
\[
 |\mathcal S\cap V|\le\frac a{1-a}|D|
 =\frac a{1-a}(z+d).
\]
\(D\) 分成 \(Z\setminus T\) 和 \(D\setminus Z\)，所以即使 \(Z\) 中有实际小源，这个分解仍正确。

远角中心到原 \(Z\) 至少为 \(\epsilon L\)，原 \(Q\) 因 \(r=o(\epsilon L)\) 与 \(Z\) 不交。按两个坐标位于盒中点哪侧分成四组，固定符号各朝盒中心；每轴可用长度至少半个盒宽，故 \(P\subset\mathcal R\)，且 \(P\cap Z=\varnothing\)。每组内部用全点 greedy 和同向不交，只向 \(D\setminus Z\) 收费。组间允许相交，但总共四组，因此
\[
 |\mathcal S\setminus V|\le\frac{4a}{1-a}d.
\]
总式为 \(az/(1-a)+5ad/(1-a)\)，得到作者 (9)。正常 \(z\) 只有领先系数 \(a/(1-a)\)；固定四倍仅落在 residual，未污染主常数。有限选中方块和原 Borel 点集保证可测性，无不可测的全局逐点选择。

## 5. 实际宏盒的 residual、坏盒和好盒

实际 \(F_i\subset T_i\)，扩展盒 \(\mathcal R_i\supset B_i^0\)。因此
\[
 |(\mathcal R_i\setminus T_i)\setminus Y_{i,N}|
 \le|\mathscr D_i\setminus Y_{i,N}|+
        |\mathcal R_i\setminus B_i^0|=:q_i.
\]
原入向面邻余量及扩张体积给扩张和至多 \(Cv\mathcal F\)（二维）。同层 R049 则给
\[
 \sum_iq_i\le(a/v)R_N+Cv\mathcal F=O(av^7S).
\]
这不是把全部已支付正常层再当作 raw。

坏盒阈值 \(q_i>av^5r_i^2\) 给坏盒原权 \(O(v^2S)\)。每盒正常层体积 \(O(av^5r_i^2)\)；在坏盒只用原粗 dyadic
\[
 \delta_{i,\rm in}\le Ca|\mathcal R_i\setminus T_i|
 \le Ca(|Y_{i,N}\setminus T_i|+q_i).
\]
坏盒求和为 \(O(a^2v^7S)\)：第一项用正常体积乘坏盒权，第二项直接用 residual 总量。其源价格为 \(O(av^8S)=O(v^{10}S)\)，确晚于八阶。

好盒上 \(|\mathscr D_i|\le|Y_{i,N}|+q_i\le Cav^5r_i^2\)。实际质量 \(\mu_i\asymp ar_i^2\) 及 \(\varepsilon_i\mu_i\le|\mathscr D_i|\) 给 \(\varepsilon_i\le Cv^5\)。global 中点盈余下界
\[
 |W_i|-\mu_i\ge c(vr_i)\Delta_i^{\rm end}
\]
遂给 \(\Delta_i^{\rm end}\le Cv^6r_i\)。每侧扩张长度也不超过扩张面积除以原盒另一轴宽 \(cvr_i\)，故为 \(O(v^6r_i)\)。

旧正常层本来在有效盒角的 \(K_Nv^3r_i\) 窗内。有效端到原端再到扩展端的总位移只有 \(O(v^6r_i)\)，因此它确在扩展盒四个固定 \(C_Nv^3r_i\) 角窗内。未把三种端点直接相等。

在局部引理取 \(L=vr_i,\ b=Cv^5,\epsilon=C_Nv^2\)，则 \(\sqrt b/\epsilon=O(v^{1/2})\to0\)。特别好盒每一个相关原半径均满足 \(r\le Cv^{7/2}r_i\)，而非额外预设微见证尺度。实际 \(\delta_{\rm in}\) 的全部点均可应用局部引理，给
\[
 \delta_{i,\rm in}\le\frac a{1-a}|Y_{i,N}\setminus T_i|+Caq_i.
\]
加坏盒后得到作者 (1)，再乘 \(v/a\) 得 (2)。没有盒外输入截断或反射网格误差。

## 6. 正常层领先体积不依赖正常八阶价格

二维对有向正常边两个排序的名义积分为
\[
 \frac{PQ^2}{4aR_i}+\frac{QP^2}{4aR_i}
 =\frac{PQ(P+Q)}{4aR_i}.
\]
实际帽缺量、共同侵蚀、正部和纤维上截断只减少旧实际层，可用此式求和上包。

在 \(p=q=2/[3(1+v)],z=0\)，无向边两端之和除以 \(av^5R_e^2\) 为
\[
 \frac8{27(1+v)^3}=\frac{2K_2}{(1+v)^3},\qquad K_2=4/27.
\]
正常紧窗内该表达一致 Lipschitz，真实根比为 \(1+O(v^3)\)。R049 的 \(\mathfrak G_N=O(v^8S)\) 与固定正常强制性给 \(\sum_{\mathcal N}d_e^2R_e^2=O(v^2S)\)，有限度和 Cauchy 给 \(\sum d_eR_e^2=O(vS)\)。故
\[
 \sum_i|Y_{i,N}|
 \le2K_2av^5\sum_{\mathcal N}R_e^2+O(av^6S).
\]
同一二维四角实际根 incidence 给 \(\sum_{\mathcal N}R_e^2=W_N\le\mathcal C+O(v^3S)\)。二维每源四个零维角槽，故 \(\mathcal C=\sum_iR_i^2\le\sum_i r_i^2=S\)，最后一步直接用 \(g_i\le ar_i^2\)。这是容量上界的正确方向，没有借未加权高维 helper。

因此作者 (3) 成立；代入 (2) 后，领先为 \(2K_2v^8S=8v^8S/27\)，误差为 \(O(v^9S)\)。这里只用正常名义领先积分、R049 RMS 和固定四槽容量，没有使用 R053 精确八阶积分或未供参数场二次变分。

## 7. 盒外／边界与真实 fresh 导数

已实际读取 R054 §2 及其正式审报。其保留原图吸收后的内部正份额，独立给同竞争族 \(P_I+P_\partial=O(v^8S)\)，不借正常八阶系数。精确边界小源式 \((v/a)\sum_{\partial I}\delta_i\le CaP_\partial\) 为十阶，环境 \(Ca(\mathcal F+P_\partial)\) 为十阶，深处盒外小源 \(Cv\mathcal F\) 为九阶。源 (246) 的三个接口确已供给，无隐藏同阶边界成本。

R054 本次核合源 SHA **bc74ac7f20b63be1387bf1f465316fb61e00513b6567002246ac32095b4f00d1**；正式审报 SHA **fc43cc5d4b3361d60daaf080e728f36ad7843a4bae866b7e6e8b5b1855ae35c2**。只需其独立正预算细化，不需微孔下界与某个特定贪心源的 inside 质量识别。

对同一实际源，\(g_i=aR_i^2,\ (1-Av)g_i\le\mu_i\le m_i^0\le g_i\)，有
\[
 \frac{f_i'(m)}{v/a}=v+\frac{vR_i}{\sqrt m}
 \le v+(1-Av)^{-1/2}\le1+(1+2A)v
 \quad(\mu_i\le m\le g_i).
\]
故一次真实积分分为
\[
 f_i(\mu_i)\ge A_vg_i/a
 -(1+O_A(v))(v/a)(g_i-m_i^0)
 -(1+O_A(v))(v/a)\delta_i.
\]
\(\sum(g_i-m_i^0)=\mathcal D_h+\mathfrak d_h^0\) 是原未删除源的精确恒等式，不因新局部估计改变。环境扣一次 \(\mathcal D_h\) 后，其正系数
\[
 [A_v-a-(1+O(v))v]/a\ge1/(2a)
\]
仍可支撑原正回收。对 \(\mathfrak d_h^0\) 继续使用原统一固定 \(b\) 和原 \(\mathfrak d_h^0\le f_{h-1}+aP_\partial\)，因此 \(f_h\) 仍是原望远镜，没有将 \(\delta_i\) 混入源缺量或用 \(b\) 放大其 sharp 系数。

二维 \(|X_i|\ge f_i(|W_i|)=f_i(\mu_i)+J_i\)，Bonferroni 仍给 \(|\bigcup X_i|\ge\sum f_i(\mu_i)-H_h\)，其中 \(H_h=\sum I_e-\sum J_i\)。所以此 refinement 与同一个 \(H_h\) 并用，未第二次消费 \(J\)，亦未在 \(H_h\) 内再扣删除质量。

## 8. 供下一总装使用的准确 H56

在上述同一新 N8 竞争族的末层，定义原实际删除量并保留同一个 \(T_i,F_i,g_i,r(x)\)。本稿供应：
\[
 (1+O_A(v))\frac va\sum_{i\ {\rm deep}}\delta_{i,\rm in}
 \le (2K_2+O(v))v^8S,\qquad 2K_2=8/27.
\]
加原深处盒外及边界删除后，全部删除费用仍至多
\((8/27+O(v))v^8S\)；环境为 \(O(v^{10}S)\)。该删除费用进入上节同一真实基线时系数为 \(1+O(v)\)，原 \(\mathcal D_h/a\) 正预算与 \(f_h\) 望远镜均保留。

这就是可供 R057 的 H56：它是实际删除成本与原账本的接口，未预设正常模型的八阶上价、任意未审有限参数优化或全输入上确界总装。前缀如需此强化，必须在其同层输入上另核；本报告不把第八层 R049 稳定性无说明搬给其他层。

## 9. 计算材料和最终签署

实际读取 r056-computation/sharp_local.wl 和成功 sharp_local.log，SHA 分别为 **e7209ffe22d3b9f70b51c189f60cc8c38ceb0c88a47fba0dc0db7578c8e388e7**、**9440677862d9eba7594760f74ec79b350729ce4f7c7fa55ea2d9a3d5d12c0aed**。

成功日志实际给同向区间正交叠违反集 False、正常双端系数残差 0、fresh 导数残差 0、平方根及质量下界违反集 False、最终费用幂次残差 0。首次幂次混淆的 first-invalid 日志不作通过证据。本次没有重跑旧批次。上述软件只核标量；一般可测局部集合的全点覆盖、不交和盒包含、所有原半径变小、实际端移／扩张及原账本对接均由本报告解析核查。

**最终签署：冻结 SHA 18fa60b96a6b2a2be7092d06dc59b4741a049e4214941323a33f846dad12ffc0 的局部可测盒引理、实际同层 (1)–(4) 及 (22)–(23) 基线接口 PROOF-VALID，无必要数学补丁。** 五处字体控制字节修复已实际独立核验并逆恢复旧稿；作者源、备份和仓库均未由本审查者改动。第 8 节给出的 H56 可供后续二维总装，但不签署正常模型八阶系数或 \(C_2\) 八阶极限。

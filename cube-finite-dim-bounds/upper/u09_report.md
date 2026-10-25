# U09有限批报告：体对角Cauchy真实势截断已获合法核接口

2026-10-25。只原clone upper/u09_*，不改总账/旧证据/bridge，不提交推送。输入、base commit及hash见u09_inputs.json；有限批已按预登记完成，独审026通过后最终锁manifest，不新增bump/成本实验。整体U3<8、U4<16仍未实现。

## 1. 两项成品及准确范围

固定V_n={(1,±1,...,±1)}、n3/4，m=2^(n−1)，L_V=−Σ_v|v·D|，所有方向共享真实时间参数、向量未单位化。完整迁移见u09_frame_transfer.md：从独立一维Cauchy和建立正保守半群；basis det4/16给Ct^-n；Gram给ψ~|ξ|、operator域H1。uniform u L1来自能量，a.e.饱和及cutoff质量回收、无限Lévy源接触、same-D全部r与一般L1均写全。没有把二维弱型常数相乘。

实质构核命题见u09_positive_part.md：G=∫_0∞p_tdt，H_c=(G−c)+，任意c>0，有

H_c≥0、偶、全空间L1，L_VH_c=μ_c−δ0，μ_c≥0、μ_c(R^n)=1、μ_c({0})=0。

这是真实全空间分布恒等式，包含全部体对角线及无穷尾源。μ_c不声称只在level面或有界域边界上；H_c不紧支撑。关键新步骤：heat mass1+Ct^-n→G weak-Lq(q=n/(n−1))→H_c L1；P_tG=G−U_t及正部给P_tH_c−H_c≥−U_t，再用L∞测试生成元极限得L H_c≥−δ0。正源质量由H_c L1大cutoff回收；原点正atom为0由有界min(G,c)和||Lχ_ε||1=O(ε^(n−1))核定。不能只检查线外点值，也不能只凭齐次说J有限。

第二文件附录还证明n3/4恰好全部frame直线上G=∞，其外有限连续：沿线one-large-jump下界p_t(av)≥c_a t^(2−n)，线外compact至少two-large-jump，任意两非平行frame方向补basis给p_t≤C_Et^(4−n)。精确log/dist渐近未声称；分类不是截断命题前提。

## 2. 完整cube支配和质量成本接口

G lower semicontinuous且严格正，cube有统一正下界，局部L1给g=essinf_cube G有限。0<c<g时

K_c=H_c/[2^n(g−c)]≥2^-n1_cube。

由此**已给出合法三四维cube核的明确构造**，而非仅愿景；它的数值常数仍由g/J决定，不能据合法性登记改善。

令J=Σ_cube_faces∫G^(n/(n−1))dS。heat weak-Lq与齐次共同给J=n|{G>1}|<∞，且

∫H_c=(n−1)/n*c^-1/(n−1)J。

归一成本B(c)=∫H_c/[2^n(g−c)]，唯一标量最优c=g/n；此时B=n^(1/(n−1))J/[2^n g^(n/(n−1))]。cube对称性等价于主设计的n^(n/(n−1))单face均值公式。沿ray齐次严格下降允许归到cube边界，不能预设corner最小；corners反而位于G∞的frame lines。

本分支没有运行主g/J数值积分，也没有全face最小认证或准确成本。主报告的32/64 Gauss约8.171962/8.171963只是主有限发现结果，未纳入本分支数值证书，未低于8；不能由此排除本族或整个frame机制。对应主results/round009_cost_probe.json读取SHA 2e7a9159bead457f85942bfd5aadcc6d074e6c32b7d402d6fbc458a995d04eb4；该记录明确sampled g是true g的上界，所以这些成本不是认证上界。

对这个具体有限源核，比较还可直接简化：记d=2^n(g−c)，则K_c*Lu=(μ_c*u−u)/d于L2，u=0处≥0。尺度r时源权为1/(dr)，正源是μ_c的概率推送。因此不必把新核的有限μ误当一般K也必有限；一般K的完整迁移仍保留无限源链。

## 3. 已登记Python有限诊断

只一个5节点对称保守path生成元Q，邻点率1、diag为负degree，t=1/4，P=I+tQ保正保质量。准确G=(0,3/4,1/2,1/4,0)的源QG=(3/4,−1,0,0,1/4)。c=3/8，H=(0,3/8,1/8,0,0)，新源QH+δ_source=(3/8,3/8,1/8,1/8,0)，非负、总质量1。P H−(P G−c)+的准确slack为(3/32,3/32,1/32,1/32,0)。

这是有限Markov凸代数诊断，不是R^n近似或线源认证；有限模型source-node仍有正质量3/8，不能替代连续μ({0})=0的小cutoff证明。

SymPy仅核n3/4层积分、radialmass、cost对数导数(nc−g)/[(n−1)c(g−c)]、c=g/n及归一公式；线/线外时间幂积分分别发散/有限。没有额外矩阵、alpha或成本搜索。最小源u09_check.py SHA54f1293b7b5f0eb7bcc3f821ee69818c6f95947c75c728c804e5adbb2f1ca16f，实际运行PASS，精简结果u09_results.json。

MATLAB只同一矩阵和合成参数质量代数，由主串行；脚本u09_matlab.m SHAffa4861a86ef6453a798b2992136814b7856a969d5cb4bf82912b01bd76685d6，upU09前缀，无clear/path。合成g=J=1只核代数抵消，绝不当真实G成本。主已串行实际MATLAB -batch exit0；rows缺陷0、source质量1、min0，μ与Jensen完全一致；n3 radial误差0/归一误差2.78e-17，n4 radial误差8.88e-16/归一误差0。无共享会话修改，完整回执及来源hash已写u09_results.json。

## 4. 恒定背景、标准依赖与状态

恒定Lebesgue用于半群密度质量、Fourier/能量、dx-Fubini、测试自伴、weak-Lq level体积、cube-face Jacobian及归一2^n(g−c)，并用于κ|D|≤∫sigma和K*sigma≤κ∫K。非恒定坏背景的分母与这些卷积不匹配，不适用。一般L1是最终单调传递，同一共同radius的全部正尺度明确处理；不是独立矩形多尺度。

标准基础为Plancherel/Rellich、正分布为Radon、Tonelli/半群与L1连续缩放；一般K迁移的LK表示仍按U08已原文核Gao§5.1使用，不声称重证。既有Lin原模型和024的来源定位只作输入，U09两份推导实际写全，未把未审锥公式当定理前提。

|类别|本批状态|
|---|---|
|严格分析/准确代数|完整frame迁移、真实正部有限概率源、全L1/J、线分类及准确有限矩阵/质量代数；独审026完整通过两锁理论和线附录；准确范围见下方|
|实际数值|Python有限诊断PASS；主串行MATLAB -batch实际exit0，最大质量代数误差8.89e-16|
|未认证|主g/J全face准确值、成本严格区间；本分支不认证主发现数值|
|未实现|无U3<8/U4<16，也无本核族/frame族普适障碍|

## 5. 下一减bump的准确条件接口（未实验）

这不是本批新增候选。若b≥0、偶、C_c∞，要使H_c−εb保持可用，必须同时付H_c−εb≥0、重新核cube下界以及全空间measure条件

ε(L_Vb)_+ dx≤μ_c（原点外）。

满足时L(H_c−εb)≥0原点外是直接准确结论。不能只支付b支撑内，因为非局部L_Vb在支撑外也可能为正，减去它会产生负远尾。

自然的定量起点是Ω={G>c}内h=min(G,c)局部恒c，故μ_c=−Lh在Ω有非负locally bounded密度

(1/pi)Σ_v∫(c−G(x+s v))_+/s² ds。

小s项因h局部恒c为0，长s由c/s²控制。此式可直接由有限方向Fubini测试得到；但Ω外的源、全域比值及cube分母仍须另轮准确结算。本批没有新增bump计算、没有选择ε，也不把内部密度下界自动当全球支付证书。下一轮是否推进由主另登记。

## 6. 白名单和停止

仅u09_inputs.json、u09_preregistration.md、u09_frame_transfer.md、u09_positive_part.md、u09_check.py、u09_results.json、u09_matlab.m、本报告和u09_manifest.json，共9文件。两理论源现锁为4ca72f62da5d4cafe3b3b2662a54e9574424c3810332b861dbf3840ad90ef265及3745b83f22961ac7977cf1b4a1e718fa77d8a2eb0691b15398bc4fb709e48912。源不因等待独审/回执继续改写；有修补需另报hash。本有限批停止，不扩计算，不改U3/U4账本。

## 7. 正式独审与最终范围

独审026全文通过u09_frame_transfer.md（SHA4ca72f62da5d4cafe3b3b2662a54e9574424c3810332b861dbf3840ad90ef265）与u09_positive_part.md（SHA3745b83f22961ac7977cf1b4a1e718fa77d8a2eb0691b15398bc4fb709e48912），含线one/two-large-jump附录，未发现补条件/缺口。已核weak-Lq/J有限、半群正部全空间源、μ质量/无原点atom、same-D frame迁移；线外连续和线∞还使pointwise g与essential inf一致。

审计报告audit/026_audit.md SHA8e4430f59682d6051081c9509e9e293b6c0b4bb8823a6978c82b588c583ecef5；其manifest SHA8896a3ad34653763160a11b737b2dbf5b738f145d39b1f1d96333fcec2ed2602，本分支已只读核实际hash。审计已披露独立推导交流的共享，不声称完全无共享证明。报告第5节后续bump条件/源密度建议未纳入026，不借此宣称修正核可行。主g/J发现值也未获本次成本认证。两锁源在审计后完全未改，没有新增运行。

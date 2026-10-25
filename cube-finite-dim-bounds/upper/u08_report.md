# U08批报告：真实坐标Cauchy的障碍与全尺度比较接口

2026-10-25。有限批完成，停止新增计算。只新增upper/u08_*；不提交推送、不改总账/旧证据/其他代理目录。目标U3<8、U4<16未实现。本文没有新低成本核，不把接口闭合当新上界，也不替代主线程正在校准的公开模型上界。

## 1. 机制、完整输入与证明内容

固定Lebesgue R^n，n3/4，L=−Σ|D_j|，真实全空间半群密度p_t(x)=Π_j t/[pi(t²+x_j²)]。跳率每坐标ds/(pi s²)、总率无限、集中坐标轴。不是Dirichlet Laplacian的谱分数幂。先对f≥0、L∞紧支撑及κ>0构造

u_t=(P_tu_t+t(f−κ))_+，sigma_t=f+(P_t−I)u_t/t。

具体证明见[u08_cauchy_obstacle.md](u08_cauchy_obstacle.md)。关键不跳步：transient上界使V_t=tΣP_jtf有界且固定t的u_t≤(V_t−κt)_+紧支撑；准确能量式给统一u_t L1；sigma_t帽和质量给生成元L2界；Fourier频率尾给强local L2/a.e.子列，证明u∈Dom(L)=H1。饱和用a.e.子列，不靠H1嵌入连续函数。极限质量先证明sigma L1，再用u L1和||Lχ_R||∞=O(R^-1)回收。结论为

f=sigma−Lu，0≤sigma≤κ，sigma=κ于D={u>0}，∫sigma=∫f，κ|D|≤||f||1。

未声称sigma/u紧支撑；非局部接触集上Lu≥0，一般不为0。真实轴跳截断的L2极限已付账。连续能量和essential-sup尾也已写明。任意Markov模型没有相应transience/频率/cutoff假设时不适用。

核比较见[u08_kernel_transfer.md](u08_kernel_transfer.md)。准确核输入为K≥0、偶、L1且LK≥0于原点外分布，不加光滑或有限源假设。原点cutoff的||L(η_εφ)||∞=O(ε)包括全部坐标轴；直接有限指数矩阵推出rho=psi Khat条件负定。标准Lévy–Khintchine表示加0≤rho≤||K||1 psi排除Gaussian项，但ν总质量允许无限。有限ν截断作用u的L2闭包给(K*L)u≥0于D^c，不作无穷大减法。全部r共用同一u/sigma/D；先有理尺度再L1缩放连续性；f_m↑|f|传一般L1。

若另有K≥2^-n1_[-1,1]^n，则中心cube弱型常数≤∫K。这个最后核构造/低成本条件本批未实现。

## 2. 公开出处与真实依赖

已原文核[Lin网页§3 Principle2](https://coolrmal.github.io/articles/two-partial-balayage-principles/)和[作者9页PDF](https://coolrmal.github.io/assets/papers/c2-less-than-4.pdf)：Lemma2 PDFpp2–5是二维坐标stable的变分构造，Lemma3 pp5–8是核比较。网页简述不能代替作者PDF的完整证明；也不把原文二维表述冒充已陈述n3/4。本批给独立离散时间连续极限推导，不作优先权声称。

标准表示定理按[Gao作者论文§5.1式(5.2)/(5.3)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7604715/)核连续实偶CND的Lévy–Khintchine形式。Plancherel、弱紧性、局部Rellich、Fourier乘子域及非负积分/支配收敛作为明确标准分析依赖。未调用Courrège原文或Delta局部H²接触引理。来源URL、实际读取bytes哈希和定位存u08_sources.json；PDF原件未在仓库内另存。

## 3. 预登记有限运行：一矩阵、一解析族

Python实际使用Fraction/SymPy/mpmath，仅5节点周期对称正率矩阵，t=1/4。P第一行为(3/8,1/4,1/16,1/16,1/4)，其余循环移位。f=(2,0,0,0,0)、κ1时准确固定点u=(2/5,0,0,0,0)，sigma=(1,2/5,1/10,1/10,2/5)，质量2、能量2/5，E+κΣu=Σfu=4/5。40步误差准确等于(2/5)(3/8)^40，约3.66e-18。

同一矩阵取f=(6,0,0,0,0)，总接收帽5小于输入质量6，故不存在有限保质量帽固定点；每步Σu至少增1/4，20步准确Σu=221479/32768。它严格反驳无条件推广到任意有限Markov背景，不是R^n Cauchy不存在分解，也不是连续近似网格。

唯一解析trial族为f=M p_s，u_T=M∫_0^T p_(s+r)dr，sigma=M p_(s+T)。Fourier准确核Lu_T=sigma−f，||u_T||1=MT，sigma保质量，T够大可满足帽；但u_T>0 everywhere、sigma不恒κ，饱和失败。它证明帽/质量/方程三项不能替代完整障碍条件。完整势V=M∫_0∞p_(s+r)dr点态有界但Tonelli给||V||1=∞；因此u≤V不能独立付清u L1。该trial输入不紧支撑，只作这些准确恒等式诊断，不用它代替定理起始输入。

固定M100,s1,T3,κ1，n3/4峰值如下（浮点展示，数学恒等式另证）：

|n|f(0)|sigma(0)|u_T(0)|u_T L1|
|---|---|---|---|---|
|3|3.225153443319949|0.0503930225518742|1.511790676556226|300|
|4|1.026598225468434|0.00401014931823607|0.336852542731830|300|

三档t=1/2,1/4,1/8的离散势峰值用mpmath40位ζ表达核右Riemann和界，均PASS。没有追加矩阵、alpha扫描、核成本搜索或优化。

## 4. MATLAB实际回执

主线程串行MATLAB -batch实际exit0；脚本u08_matlab.m SHA a0036035547d76e13929ebdd3ea3d0708fca3db13f18ffa432686a96a4a1bc96。所有变量upU08前缀，无clear/path，未改共享会话。准确同一P的行和与固定点缺陷均0；质量2、能量0.40000000000000002；sigma与Python一致。20步超帽质量6.7590026855468741，和Python准确值差在浮点舍入。

解析族峰值/质量与Python一致；另独立一维Cauchy卷积积分0.17425723696192918，准确半群密度0.17425723696192921，误差2.78e-17。它核pi/时间归一，不是障碍极限的数值认证。完整stdout及日志来源hash存u08_results.json；原日志在scratch/round008_u08_matlab.log。

最小计算证书u08_check.py SHA 42f1e86962874d5717dac0884b6dea646ded5128a444288f40a19fbeae38ac25。运行路径使用本机scratch/venv/bin/python；该脚本重新生成Python数值结果，MATLAB回执另由本报告与原日志记录。

## 5. 四类状态、恒定背景与下一接口

|类别|本批状态|
|---|---|
|严格分析/准确代数|两份完整证明稿；矩阵固定点/容量反例与解析族恒等式准确核；标准LK等依赖明列。两锁理论全文经独审024通过；范围/引用见下方|
|实际数值|Python单矩阵40/20步、一个解析族三档势峰值；MATLAB -batch实际exit0且Cauchy卷积归一交叉|
|未认证|尚无低成本K满足全部原点外分布正性及连续cube支配；没有数值连续障碍解或成本证书|
|未实现|无U3<8/U4<16，没有排除Cauchy或更广生成元的核族障碍|

Lebesgue实际用于对称Fourier/保质量、轴源的dx-Fubini、cutoff自伴配对、K*sigma≤κ∫K及κ|D|≤∫sigma。若背景换w(x)dx，分解、分母与普通卷积不匹配，不能套到已知2^n坏背景。固定Lebesgue下原子输入的平滑传递未被删除。

下一步应寻找真正非乘积、成本低且可准确核LK≥0的坐标Cauchy核或相应连续必要条件。只对昂贵(G−a)+势截断子族增加数值分辨率没有本批依据，未开展。主处理的公开模型→显式cube核成本卡与U02严格积分不在本分支重复。本有限批停止；总体突破目标仍未实现。

## 6. 整合白名单

仅11文件：u08_preregistration.md、u08_cauchy_obstacle.md、u08_kernel_transfer.md、u08_check.py、u08_matlab.m、u08_results.json、u08_sources.json、追加定点来源说明u08_lin2_source_note.md及其u08_lin2_sources.json、本报告u08_report.md、u08_manifest.json。manifest列其余文件hash/大小；不自记自身hash。独审文件属于audit分支，不复制到upper白名单。

主追加来源核验已完成：Lin二维3.879的实际公式/参数、钉版本Lean证书和合法维数迁移边界详见u08_lin2_source_note.md；只核原常数与结构，不新增成本搜索。

## 7. 正式独审定位与最终锁定

独立审计024全文通过u08_cauchy_obstacle.md（SHA c1163c7b40ec5d5723228df0959fc00f7ae91197db49fe0bfc7f132fa6ad7aec）与u08_kernel_transfer.md（SHA 98ee459c6d24ad80e4bff39822d7f27746ff046710b0cc5ce3a64183d807f405），无新增隐藏假设/未闭合环节。实际核准的链包括V_t非L1但能量付u L1、渐近谱尾/局部紧性/H1、a.e饱和、先sigma L1再cutoff回质量、原点Oeps、直接CND/无Gaussian、无限ν_delta的L2接触极限、same-D全部尺度和一般L1。

审计报告路径为audit/024_audit.md，SHA 33b5ce08d15c9c960dd43c4707c507fe5c72bd597b35d12b5baf3738f60ee606；其manifest SHA 26d81252fdc7c6cd16718cb790b32b8cba42321ca75c14278da198819ddf3737。本分支已只读核实际文件hash。独审没有重跑本批矩阵/trial；其独立无限Lévy源范围例不计为本分支的新trial。

**审计范围仅上述两锁理论。** u08_lin2_source_note.md中的下一轮finite-frame扩展接口没有独审，不能借024扩大通过范围；作者Lean证书仅读取源码，不冒称本机重建。两理论文本在审计后未改动。

# 第45轮：已接受平方根来源的参数与可计算接口

范围：只提取用户已接受的构造，不重新评判总体定理或主√n下界。`parameter_probe.py` 已在本机运行Decimal80/Fraction参数算术及一个有限来源的实际采样/CDF检查，结果为 `PARAMETER_PROBE.json`。这些参数门和成本表不是有限维C_n证书，也没有生成新的qualified GF basis或完整安全尺度族。

## 1. 实际测试对象与零模

主来源是正的**一维周期概率密度**，不是直接以原子测度测试弱型常数。论文先构造原子位置a_j，再用正宽箱平滑：

\[
\beta_\delta=\delta^{-1}1_{[-\delta/2,\delta/2]},\quad
\mu_{K,\delta}=K^{-1}\sum_{j=0}^{K-1}\beta_\delta(\cdot-a_j).
\]

每个峰质量1/K，故∫μ=1。定义有限补偿谱

\[
c_m=\sqrt{2b_m/3},\qquad D_M(x)=\sum_{m=1}^M\pi m c_m\cos(2\pi mx),
\quad \rho=1-p+p\mu_{K,\delta}+(p/K)D_M.
\]

**D_M没有0-mode，积分为0；它不是standard Dirichlet kernel。** 因而∫ρ=(1−p)+p+0=1，不是1+p/K。来源：[paper.tex:365](../accepted/paper.tex:365)、[paper.tex:879](../accepted/paper.tex:879)、[paper.tex:886](../accepted/paper.tex:886)、[paper.tex:891](../accepted/paper.tex:891)。程序JSON保留了这些精确质量值与频率范围。

周期输入F_n(y)=∏_iρ(y_i)；欧氏L1见证为f_n=1_{[-R,R]^n}F_n，R=n²。其质量(2R)^n，使用L∈I=[10/9,11/9]时内部[-R+1,R−1]^n平均与周期输入一致，边界损失(1−1/n²)^n→1。[paper.tex:1984](../accepted/paper.tex:1984)。全文使用Lebesgue/Haar查询：X_i iid Unif[0,1)，而非X_i iid ρ。把ρ sampler当query sampler会改变原概率目标。

## 2. 完整的n依赖处方

先固定0<p<1及有限w≥4。令c_p=3p/√29。对每个充分大整数n，取满足

    v=9p² n/(29K²)≥16 log n

的最大允许周期K=2^d−1。下一周期2K+1不满足该式，故16log n≤v<(2+1/K)²16log n。置u=√v，h=v/2−2，测试阈值e^h。然后

    B=不小于v的最小二幂， M=ceil√v， δ=1/(wKv).

因此v≤B<2v，K≈√(n/log n)，v=Θ(log n)，M=O(√log n)，δ随n缩小，K/v^r→∞对每个固定r成立。[paper.tex:832](../accepted/paper.tex:832)、[paper.tex:859](../accepted/paper.tex:859)。这不是固定δ、固定v或固定尺度数的toy。

有限构造还要分别核对：K+1>H(B)（论文给较强充分门K+1>3B⁴）、1≤M<K、峰在粗格内、背景正性，以及完整核/概率工具的其他条件。K+1>H(B)是存在性充分条件，门失败不等于证明不存在某个显式合格basis。B≥v是所接受增长族的参数处方，不是任意周期密度质量合法性的必要条件。

### 条件的归属

- **B≥v、B<2v**：来源分辨率处方。它使v/B²误差消失，同时局部物理差枚举中的O_T(B/v+1)保持有界；不是一个通用弱型定理的独立假设。[paper.tex:945](../accepted/paper.tex:945)、[paper.tex:1299](../accepted/paper.tex:1299)。
- **v≥16log n**：一方面定义最大K并覆盖全部整数n；另一方面明确用于远区交集总预算Kv exp(−v/25)≤C_p√log n·n^(−7/50)→0。[paper.tex:1912](../accepted/paper.tex:1912)。它不是每个固定单点平均合法性的条件。
- **B≥256、M≥32、e≤.001**：历史第79版near工具/数值报告里保留的资格门，见 [numerical_report.tex:144](../accepted/numerical_report.tex:144)。当前展开版paper核引理写的是V_min、η_D、η_N、e等明确式及最终充分大条件，没有把这三个历史门重新列为其唯一统一有限前提；不要混称所有版本都要求同一门表。

## 3. 谱和真实平均：可直接实现的公式

固定s=1/9、a=10/27、α=729/29、γ=27/29；

\[
b_m=4\alpha\left[\left(\frac{\sin(\pi ma)}{\pi m}\right)^4+
3\left(\frac{\sin(\pi ms)}{\pi m}\right)^4\right].
\]

该源保留完整H+B及其交叉项，不能把方差直接设成v。令P为每格原函数、P_δ为箱平滑，

\[
H_{L,\delta}(x)=P_\delta(x+L/2)-P_\delta(x-L/2),\quad
B_L=\sum_{m=1}^Mc_m\sin(\pi mL)\cos(2\pi mx),
\]
\[
Y_L(x)=A_L^{(1)}\rho(x)=1+\frac p{KL}(H_{L,\delta}(x)+B_L(x)).
\]

因此真实中心立方体平均是∏_iY_L(X_i)，S_L=Σ_i log Y_L(X_i)。同一组X_i及同一ρ用于全部L。[paper.tex:114](../accepted/paper.tex:114)、[paper.tex:363](../accepted/paper.tex:363)、[paper.tex:928](../accepted/paper.tex:928)。

解析幅度a_n=77p/(15K)，q_L=nE(Y_L−1)²=v+O_w(1)在安全尺度上一致；a_nv→0、v²/n→0。[paper.tex:955](../accepted/paper.tex:955)。尾计算应存h、logY、log权，不能在h>709时先用binary64算exp(h)。

## 4. GF标签的精确定义及生成成本

q_field=2^d=K+1，不要与方差q_L或前轮floor(L)符号混用。固定GF(2^d)及本原元g，B=2^b，取β_0,...,β_{b−1}的F₂独立span W，满足论文的乘法Sidon及比率邻接排除条件。定义

    y_i(j)=Tr(β_i g^j)， y_j=Σ_i2^i y_i(j)，
    U_j=(y_j+1/2)/B， a_j=(j+U_j)/K， j=0,...,K−1.

必须遍历非零field状态；全零状态未包含，标签联合计数要减去全零纤维的一次。[paper.tex:300](../accepted/paper.tex:300)、[paper.tex:350](../accepted/paper.tex:350)。不能把各j的label当成独立uniform整数替换这条确定序列。

论文给H(B)=(B−1)+T(T−1)+2R_B²+4R_B，T=(B−1)B/2，R_B=(B−1)(B−2)。实际候选验证无需O(B⁴)逐对碰撞比较：生成W非零元素后，把O(B²)无序产品排序/去重检验Sidon；生成R_B个非平凡ratio，检查数量、与g·ratio的交及g^±1,g^±2排除。全部直接用field元素完成，不必生成离散log表。primitive g和field多项式本身仍需精确验证。

这是一条可实施的有限basis验证路线，不是本轮已生成basis的声明。字典序首个合格basis的存在性处方并未给快算法；实际有限随机尝试也不能代替资格核验。B2048时ratio原始uint64载荷33,505,296 bytes，B4096为134,119,440 bytes；排序/hash及产品表还需额外内存。d>64时每field元素需多字，JSON已另记真实字数载荷，不能继续按uint64理解。

不应存整个K标签数组：给定j可即时计算g^j及b个trace。P_δ点值因箱峰留格，仅需要附近少数粗格的label；再加M项Fourier即可得到实际Y_L。浮点全局x不能无条件支持极细δ，建议用“整数粗格j＋高精度局部相位”或足够位数的固定点表示。

## 5. 已实际计算的匹配参数门与成本

以下p=.99、w=64，仅展示处方及列出的有限代数门；**不是完整见证或显式起始维数**。

|n|d，K=2^d−1|v|B|M|δ|K+1>3B⁴|
|---:|---|---:|---:|---:|---:|---|
|10^12|14，16383|1133.2560|2048|34|8.416e−10|否|
|10^30|43，8796093022207|3931.2917|4096|63|4.519e−19|否|
|10^32|47，140737488355327|1535.6608|2048|40|7.230e−20|是|
|10^34|50，1125899906842623|2399.4701|4096|49|5.784e−21|是|
|10^36|53，9007199254740991|3749.1720|4096|62|4.627e−22|是|

n=10^32时全K的uint16标签就约2.81×10^14 bytes；按bin宽≤δ/6取统一query网格至少约8.30×10^19格。n=10^36相应为1.80×10^16 bytes和1.30×10^22格。因此第38/40/41轮“先把n坐标压成Q格counts”的方法虽避免显式n维向量，仍不能直接沿用其dense Q数组。

这两行若要求query分辨率≤δ/64，至少需要70/77 bits的全局固定点精度；原53bit uniform浮点网格不足以解析这些峰边界。n也已超过int64，NumPy multinomial的普通n参数不能直接表示；不能悄悄改成Gaussian counts后称为相同模拟。

这些门通过不等于完整渐近误差已量化。安全节点、basis、完整核误差与固定辅助T/A/N/D/w的概率比较还未实例化。尤其w64只自动覆盖A≥1/64；若采用A=T^−2，只能由这条预算直接支持T≤8。不能拿该表声称已经实现系数1/35。

## 6. 安全尺度与局部窗口

主窗口I=[10/9,11/9]。缺陷集D来自GF比率，实际安全集要求dist(K(L−1) modK,D)>2。[paper.tex:512](../accepted/paper.tex:512)。取固定有限T=NA、间隔D_gap>0、P0=T+D_gap，要求w≥max(4,1/A)。原始段及节点为

    L_i=10/9+2iP0/(Kv)， L_{i,a}=L_i+2aA/(Kv)， a=0,...,N.

任何一个节点不安全就删整段；保留数J*≈Kv/(18P0)。这是n依赖、数量增长的实际共同尺度；不是233点或19点固定网格。[paper.tex:1780](../accepted/paper.tex:1780)。局部坐标L(t)=L0+2t/(Kv)，有限t节点与箱平滑ε=1/w相匹配。[paper.tex:1231](../accepted/paper.tex:1231)。

合法极限顺序先固定p,T,A,N,D_gap,w，取全部整数n极限，再对已建立的数值下界取w→∞、T→∞、p↑1。论文没有提供T,N,w同时随n增长的统一误差率或显式n0。[paper.tex:2152](../accepted/paper.tex:2152)。

## 7. 历史pilot为什么不能冒充匹配源族

38/40/41来源为K255、B4、v8、p=.5、w4、M3、n6,704,800。它们是合法有限密度；但B4<v8，且8<16log n，所以K255不符合该n的最大允许K处方。其相位更改后的多表也未自动继承baseline全部pair核性质。完整ρ的粗格质量受D_M影响，一般不是1/K，实际二类η≈.498047800889，不可偷换为1/2。[038/REPORT.md:42](../../038/REPORT.md:42)、[040/REPORT.md:7](../../040/REPORT.md:7)、[041/REPORT.md:85](../../041/REPORT.md:85)。

论文附带历史B8实例则不同：K1023/4095，p1/4、w64、v1024、M32，n=55,249,289,216 / 885,285,273,600。实际计算显示它们**通过**v≥16log n及最大K规则，却**失败**B≥v；若照处方应B1024而不是8。旧数值near门B≥256、e≤.001亦失败，e≈.217573313/.103414627。不要把这两批pilot的失败原因混在一起。[numerical_report.tex:113](../accepted/numerical_report.tex:113)、[numerical_report.tex:144](../accepted/numerical_report.tex:144)。

## 8. 正背景来源sampler：实际小例已执行

记h_bg=1−p+(p/K)D_M，H_M=Σ_{m≤M}1/m，C=20pH_M/(3K)。解析上下界为

    1−p−C ≤ h_bg ≤ 1−p+C.

处方中的正性门确保h_bg≥(1−p)/2。故ρ可作为正混合抽样：概率p取均匀j及该峰的Uniform[−δ/2,δ/2]；概率1−p从h_bg/(1−p)抽样。后者用Uniform[0,1) proposal，接受概率h_bg(x)/(1−p+C)，理想平均接受率(1−p)/(1−p+C)。这保留signed Fourier补偿，未把D_M错误当成独立概率密度。

按root要求，实际使用历史K1023、B8、v1024、p1/4、w64、M32及原NPZ标签执行了100000个ρ样本。解析h_bg上下界为.743387919199354 / .756612080800646，理想背景接受率.991260936789；共25012峰样本。对257个节点，用显式CDF

\[
F_\rho(x)=(1-p)x+\frac pK\sum_j\operatorname{clip}\left(\frac{x-(a_j-\delta/2)}\delta,0,1\right)
+\frac p{2K}\sum_{m=1}^Mc_m\sin(2\pi mx)
\]

比较，最大经验差.00255205251，端点0/1。这是有限来源实现诊断，B8<v1024明确不匹配，trig/PRNG未作区间或随机性认证；不由此宣称高维下界。

## 9. 下一步可执行接口及尚缺的压缩

先实现并冻结一个**隐式**source manifest：n,p,w,K,d,B,M,δ，field多项式、primitive g、实际验证过的β basis、defect ratio集合，及上面的非零频率系数。优先n10^32/B2048，basis验证载荷比B4096小；不存整个K标签。不具备qualified basis时只报告候选参数，不能换成独立随机labels。

随后做少量真实L0及有限local packet的source-average oracle：用整数粗格+高精度局部相位计算P_δ端点及M项补偿；对原ρ的直接窗积分交叉核对。标准query基准为Haar；若用anchor倾斜Q_L=Y_L·P，必须保留准确RN及同一组查询在全部局部节点的耦合。ρ sampler只是source抽样工具，不能替代这个query law。

要把它推进到n10^32的稀有尾，仍需新的一层计数压缩：单/少数尺度的端点labels可按同一GF state的线性trace映射之秩精确计数，非零state需作zero-fiber减1；不能逐尺度独立采labels。原始公式在 [paper.tex:321](../accepted/paper.tex:321)。与缓变Fourier补偿的**联合**律仍需解析分层/可控包络，不能在collapse时删掉它或假设它与labels独立。旧dense Q counts只解决n大而Q尚小的场景，当前Q同样巨大。

建议第一轮只产出有限packet的真实均值/协方差/倾斜矩及误差账本，利用已有解析核预算作控制变量，再决定大n的稀有尾算法；未取得的高维count sampler、有限packet概率和全增长尺度付款分别记缺口。单点Gaussian/Brownian极限可作accepted baseline的接口参考，不可冒充已运行的真实n维MC或新的有限C_n证书。此计划不保证额外log收益。

## 10. 接受常数与本提取的边界

C*=87/(220√(58π))≈.0292959535703970；论文的严格比较为C*>87/2992=1/35+53/104720>1/35。[paper.tex:2160](../accepted/paper.tex:2160)。p先固定时对应下界系数为pC*，最后才取p↑1的上确界；实际构造从不代入p=1。计算上p要大于约.97526877，pC*才超过1/35，但这只是常数目标的必要选择，不代替有限T/w及未知n0资格。总体定理按用户接受状态使用，本轮没有重新审计它。

# 031 正权body+face Cauchy frame：完整迁移条件独审

2026-10-25。结论：主 `proofs/round011_weighted_frame_transfer.md` SHA256 `78ff3d964c555d9dd857da7a4bd9bbb243d5011ebff4b00798cf217ff9239f7f` 的正权frame迁移通过。对固定n3/n4、body权1、全部face方向同一固定τ>0，026的障碍分解、任意允许无限Lévy源的核比较、正部势全空间源及一般L1/全部共同尺度结论均适用。τ=0须删除face方向后回原frame。没有发现需修改作者稿的数学缺口。

这不是cb11真实Green锥表示已审，也没有认证g≥13/2500或I≤57/100000。两项目标若另行成立，确实给条件成本约7.8986291251932403343<8；本审计完成时它们仍只是目标，没有新上界。

## 锁源、依赖和独立计算

先读011f登记，七个来源锁在 `031_source_snapshot.json`。依赖026审计SHA `8e4430f59682d6051081c9509e9e293b6c0b4bb8823a6978c82b588c583ecef5`；原两全文u09_frame_transfer和u09_positive_part分别为 `4ca72f62da5d4cafe3b3b2662a54e9574424c3810332b861dbf3840ad90ef265`、`3745b83f22961ac7977cf1b4a1e718fa77d8a2eb0691b15398bc4fb709e48912`。这些旧证明不重跑；此次重新检查新权frame是否满足每一实际前提，不是仅把标题中的生成元替换。

独立Python第一次执行通过，内部约0.235秒，无实现修补。新代码重建Gram、立方体群保权、正齐次、全部正幂成本代数，并用256bit Arb独立核目标成本；未导主程序、未计算Green值或网格。共享SymPy1.14.0、Fraction、python-flint/Arb及026既有证明明确披露。

## 半群、Gram、域和权重正常数

每方向独立Cauchy变量尺度为ci t，body的ci=1，face的ci=τ。故特征函数是exp(−tψτ)，方向跳率ci ds/(πs²)；未单位化的body长度√n、face长度√2均保留，权不应平方。这是正的自伴、保质量/常数、L2收缩卷积半群。

独立准确Gram为

\[
 \sum_b vv^T=2^{n-1}I=m_bI,\qquad
 \sum_f vv^T=2(n-1)I,\quad m_f=n(n-1).
\]

由l1≥l2得body下界√mb|ξ|，两组分别用Cauchy–Schwarz给上界

\[
\sqrt{m_b}|\xi|\le\psi_\tau(\xi)
\le\{m_b+\tau(n-1)\sqrt{2n}\}|\xi|.
\]

所以form域H1/2、L2 operator域H1，范数等价含通常L2项。τ固定有限；上界常数允许依赖τ，下界只需body。没有使用n3/n4中不存在的H1到C0嵌入。

显式body基的det绝对4/16已独立代数核。其n项和的正连续C0密度以1/(d*π^n t^n)控制；其余body及新增face变量仅作概率卷积，故该峰值上界不随τ变坏。密度在每点严格正、总质量1，并且全体变量的共同尺度t给pt(x)=t⁻n p1(x/t)，从而正时间L1连续及joint连续。这些是026实际用的密度条件。

body和face各自的无向方向集在每个坐标反射、相邻坐标置换下不变；独立代码核了cube群这些生成元。权类别没有混合，故ψ、pt、G均具cube对称性。

## cutoff、障碍分解与接触集

所有测试的有限方向二差分现在乘ci。大cutoff保持||LτχR||∞=R⁻¹||Lτχ||∞；小原点cutoff gε=η(x/ε)φ，φ≥0且φ(0)=0，sup为Oε²、Hessian O1。每方向以|s|=ε/|v|分割后，其界为Oε ci|v|，总和有限。具体总长度为n3的4√3+6√2τ、n4的16+12√2τ。

Lτφ属于L1∩C0∩L∞，用min(s²|v|²||D²φ||,4||φ||)作sup/L1控制。小cutoff的L1缩放是ε^(n−1)，大cutoff的sup缩放是R⁻¹；所需常数均有限，不会漏掉新增方向的测试项。

026固定t最小不动点构造及能量证明的输入为真实正保守半群、上述峰值Ct⁻n及ψ可比|ξ|，这里均已成立。逐项核：固定t的u_t先因V_t−tf∈C0和u_t≤(V_t−κt)+有essential紧支撑；统一L1来自能量而非V_t可积性。σ_t帽、保质量、统一L_tu_t的L2界、ψ_t频率尾、局部Rellich/a.e.子列、正集饱和与大cutoff回收质量都不依赖方向数的旧数值。所得u∈L1∩L∞∩H1、0≤σ≤κ，σ在D={u>0}上κ，∫σ=∫f。

原方向生成元截断|s|>ε的总率必须改为2(mb+τmf)/(πε)，即n3的(8+12τ)/(πε)、n4的(16+24τ)/(πε)。其符号0≤ψε≤ψτ并逐点趋ψτ，故Lεu→Lτu于L2；在u=0处非负jump-in给Lτu≥0 a.e.。没有沿用未加权的旧方向个数，也没有断言非局部零集上Lu=0。

## 任意核的无限Lévy源闭包

令K≥0、even、L1且LτK≥0于原点外分布。小cutoff把正性升级到所有非负φ、φ(0)=0；对χR|Σzi exp(iξi·x)|²，Σzi=0，有限加权二差分的统一支配和K L1允许R→∞，得到ρ=ψτ Khat连续实偶CND。由K≥0和CND，0≤ρ≤Mψτ，M=∫K。

标准实偶Lévy–Khintchine表示仍允许总质量无限的ν。ρ的线性增长排除Gaussian，ρ(0)=0排除killing，实偶排除漂移。新增方向不意味着这个ν仍只在原有限方向上，也不意味着它有限。

由u∈Dom(Lτ)、ρ≤Mψτ，得到u∈Dom(B)、Bu=K*Lτu于L2。将ν截为|z|>δ的有限νδ，其符号0≤ρδ≤ρ并趋ρ，因此Bδu→Bu于L2。在u=0处有限卷积Bδu=νδ*u−νδ(Rn)u非负；dx-Fubini和L2子列给Bu≥0 a.e.于D^c。原方向的|s|截断与本段一般ν的|z|截断是两种不同有限化，均付清；没有作∞−∞。

## 同一D、全部共同r与一般L1

ψτ的正一齐次给

\[
\rho_r(\xi)=r^{-1}\rho(r\xi)\le M\psi_\tau(\xi).
\]

因此所有K_r的接触正性都能用同一个u、σ、D，不重选接收坏集；Kr*σ≤κM。对有界紧支撑输入取κ=λ/M（M=0情形平凡），正有理r的零集可数合并，再由L1核的缩放连续及||f||∞控制扩至全部r。每个核缩放虽有order1源缩放因子，ρ_r的齐次抵消恰给相同M，没有额外r常数。

一般f用非负有界紧支撑f_j↑|f|逼近。非负卷积与尺度sup都与该单调极限相容，得到坏集和质量的单调极限；同时sup可测。没有给一般L1输入强加H1假设，也没有要求所有j共用同一个障碍解。这里的同一D指一个分解下的所有共同尺度，不是对所有逼近输入强行同D。

若K≥2⁻n1cube，原始cube中心共同半径平均逐尺度被Kr*|f|支配，故弱型常数≤∫K。恒定Lebesgue的平移卷积、自伴、Fourier、dx-Fubini、质量守恒、K*σ≤κ∫K均明确保留；不适用于坏密度背景或独立坐标多尺度。

## 正部势：所有新增线、界面与无穷尾均在全分布中

G=∫ptdt的短时部分质量T，长时sup≤CT^(1−n)/(n−1)，给weak-L^(n/(n−1))及G locallyL1、a.e.有限。任意c>0的H=(G−c)+因此全空间L1，不要求紧支撑或在线上有限。

LG=−δ0的测试绝对配对仍由Lτφ∈L1∩L∞及相同Ct⁻n控制；uniform测试生成元极限的Fourier条件ψτ|φhat|∈L1也成立。正Tonelli给G=PtG+Ut，PtG在t>0后有界，因而PtH≥H−Ut的正部半群论证不会作∞−∞。测试后t↓0得到LτH+δ0≥0于整个空间分布，不只在G有限光滑处检查源号。

先有正Radon μ，再用H L1和大cutoff MCT给μ全质量1。对h=min(G,c)有界，μ=−Lτh；Lτχε的L1缩放ε^(n−1)趋0，给μ({0})=0。新增方向线、非光滑level界面及无限尾均包含于该正源，不删除任何分布项。此源μ是本具体正部核的有限概率；上一节任意K的Lévy源ν仍可无限，两者没有混同。

## g、J和边界essential trace

严格正连续pt在紧cube、t∈[1,2]上一致正，给g=essinf_cubeG>0；局部L1给g有限。G(rx)=r^(1−n)G(x)。weak-Lq先保证V1=|{G>1}|有限，再以cube face射线Jacobian r^(n−1)和Tonelli给J=nV1<∞；不能仅从齐次推J有限。

文本的“g可归边界”可直接按essential face trace严谨理解，不需要先证明新的全方向线奇性分类或把lsc的point inf等同essinf。令g_face为所有faces的面测度essential inf。x=ry，0<r<1，给G(x)≥G(y)，故g_cube≥g_face。反过来若face正面积集上G(y)<a<g_cube，选a<a'<g_cube及r足够接近1，使r^(1−n)a<a'；这给cube内正体积的G<a'，矛盾。因此g_cube=g_face。有限J保证相应face值a.e.有限，径向论证无隐藏正则假设。未来逐点g下界认证可以比所需更强，本审计未据lsc单独断言不存在孤立低值。

齐次层积分给∫H=(n−1)J c^(−1/(n−1))/n。若真正认证0<g_low≤g，取c=g_low/n，注意这是本新模型选择，未沿用U10固定c。于是H≥g_low−c于cube，K=H/[2^n(g_low−c)]为合法majorant，成本满足

\[
C_n\le\frac{n^{1/(n-1)}J_{up}}{2^n g_{low}^{n/(n-1)}}.
\]

cube群保权的独立代数核使J=n2^n I，I为单正象限face积分。全部量正，抬至n−1幂不改变严格次序，准确突破判据为n^n I_up^(n−1)<U^(n−1)g_low^n。n3的U=8给27I²<64g³；n4的U=16给I³<16g⁴。这是构核条件成本的比较，不是从有限样本minimum推连续下界。

目标g_low=13/2500、I_up=57/100000的准确slack为56653/250000000000>0，条件成本平方为2193075/35152，独立256bit sqrt包围严格位于主192bit端点之间，展示成本7.8986291251932403343181858260407。JSON状态和本报告均只称目标；没有g/I的真实证明便不能更新8/16。

## 复现、MATLAB与终止

主导入后最小复核（stdout-only；可选--save）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/031_check.py
```

审计clone运行加 `--project /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds`。脚本先核七来源hash。MATLAB只交主串行核有限Gram和目标有理代数，实际状态在 `031_matlab_receipt.json`；不会替代这些分析极限。

本批未审cb11完整锥/退化Gamma-log表示、Green数值、continuous minimum或J积分；未运行LP/grid/优化。026依赖沿用但未重复旧认证，所有新条件已逐项核。完成manifest即停031，后续实际候选另锁源审计。


主串行MATLAB已实际exit0（session40045），同SHA196f9d47fd31b90b7d8583afcbb108ed3dc2cdc036ee87af5204bd0b7a00d9af。审计已读stdout，n3/n4 Gram及精确slack56653/250000000000、cost²2193075/35152均assert通过，vpa展示7.8986291251932403343；log SHA676b02c8b25c8370c68da42dc5bab2c42f1a6d59b7c641be2d2ea63995b9262c。脚本自身强调targets only。收尾七来源未变，本批已停止。

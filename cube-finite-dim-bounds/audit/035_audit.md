# 035 四维完整候选最终独审：接受 C4<14.753

2026-10-25。对固定四维purebody、Lebesgue、centered axis-parallel cube、所有共同正半径及一般L1输入，本次完整链审计接受

\[
\boxed{C_4<14753/1000=14.753<16.}
\]

实际组装成本不超过约14.7525411723433123415。此结论来自真实1152锥Green表示、全域g下界、2440个连续积分cell及唯一真实角帽、合法正部核和准确幂比较，不来自作者/Astra的PASS标签或任何成本采样。只接受本固定构造，未审后续11Q sharper minimum，也未借此接受三维I证书。总账由主更新，审计者没有修改。

## 锁源与独立执行

先锁 `035_source_snapshot.json` 的43个主输入，含U11完整表示、g manifest `ce84e157…`、global data `43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f`、CI4B manifest `81da109495409d65a585ba3595dbb487c8a07683027d4de308deced14488137e`、certificate `4dbb544ba34762b239adea8e4dc152987ce33c15ef855f459efb63bffa1281b0`、Astra附审、主shell/assembly及actualcost。全部实际hash列于快照，不依赖可变ledger。

两个新审计实现均首次实际通过：

- `035_representation_check.py`：1.508秒，独立生成1152锥及其rates，重算六star、三个完整rational式、全5143个g系数；三个预锁退化点的full-cone confluent residues准确核。
- `035_integral_check.py`：61.125秒，独立生成全部2440几何和固定bisection，逐cell准确重算4,932,474个非零残差系数，全部非负。256bit完整积分、帽及成本比较通过。

没有执行作者producer/checker/engine、Astra脚本或任何LP，没有读取候选数值采样来判断正性。只共享锁定的输入多项式/有理Pi、SymPy/Fraction、python-flint fmpq_mpoly/Arb算术内核；全和多项式输入先由本审计独立证明等于真实G，再用于积分。Astra报告读作有用附审，未当自动接受依据。没有数学参数修订、优化或扩网格。

## 原始1152锥与真实正常数

在a≥b≥c≥d≥0排序chamber，三个cone的逆非负参数分别为

\[
\begin{array}{ll}
I:&(a-b-c-d,b-c,c-d,d),\\
II:&(a-b,b-c,b+c-a-d,d),\\
III:&(c-d,b-c,(-a+b+c+d)/2,(a-b-c+d)/2).
\end{array}
\]

列依次ABEF、EBCD、EBDF，其det绝对1、1、2。比较a−b−c与±d给三个区域内部不交、闭包并满排序chamber；384个signed permutation chambers再覆盖R4，重复只在零体积边界。这是逆坐标/不等式证明，非体积相等推覆盖。

新实现对每锥原八个body线性forms在四ray上的符号准确核全同弱号，所以ψ在整锥真线性。rates由实际Σbody|v·ray|求和，不从作者rate表相信正常数：A..F为8,8,12,12,16,24。独立group生成元检查全cone/rate集合等变；共有240个triangulation rays。

各锥ψ球simplex体积为det/[24Πrate]，总|{ψ≤1}|=11/4608，故

\[
p_t(0)=\frac{11}{3072\pi^4t^4}.
\]

Fourier逆因子(2π)⁻⁴、det和未单位化v都保留。独立spanning Cauchy和给真正的正、保守概率半群，ψ可比|ξ|使Fourier积分绝对收敛。每锥指数积分因此给全部1152个复resolvent乘积之和；反锥恰为共轭，负实部项全部保留。

全谱球真正40个vertices亦核对：A..F上zero body-direction ranks为0,3,0,3,2,1，active normals rank为1加该rank，只有B(24点)、D(16点)是vertices。不过本证明没有仅凭vertex定理删star；下一段直接准确证四个null-star恒等式。

## time积分、π³符号及全部zero/repeated节点

对非零节点b_j=x·r_j/q_j，generic partial fractions的simple pole系数为i/Πk≠j(b_j−b_k)。取截时ε>0后积分、再取实部得到−(π/2)DD3(sign;b)。因此

\[
G(x)=-\frac1{32\pi^3}\sum_R\alpha_R DD_3(\operatorname{sign};b),
\quad\alpha_R=|\det R|/\prod q_j.
\]

独立scalar calibration：nodes(−1,−2,1,2)的DD为−1/6，直接∫dt/[(t²+1)(t²+4)]=π/12；nodes(−1,−1,1,1)的DD为−1/2，直接积分π/4。这同时核了符号、1/32及π⁴到π³的变化。

零/重复节点必须先作全锥partial fractions。若b的multiplicity为m，则pole-k系数为i^(k−4)R_(b,k)，其中R=α[u^(m−k)]Πa≠b(b−a+u)^−ma；这是把t−ib=w后各因子按w/i展开所得，不是猜sign0。新审计独立实现此Taylor系数全和。

在线外紧邻域至少两个Cauchy变量必须大；28对中任意一对可补非零body basis，最小非零det为8（本审全28对准确核）。两个受限密度各≤t/(πδ²)，另外两个basis密度各≤1/(πt)，余变量仍概率卷积，所以pt在该邻域小t一致O(1)；大t用t⁻4。于是G在线外有限且**真连续**。在body线上，一条大坐标Cauchy密度~t、其余三个~t⁻1及剩余事件固定正概率给pt≳t⁻2，故G=∞；原点直接t⁻4发散。

全锥复rational函数就是有界的真实pt，所以离active lines时b=0的所有principal parts必须准确取消。非零b的k≥2 pole时间积分实部为0；k=1的实部为−Rb/(t²+b²)，积分为−πR sign(b)/2，得到全节点公式。不能逐cone把zero sign设0；例如(0,0,1,−1)单项有1/ε发散。

预锁e1、(1,1,0,0)、(1,1,1,0)在新1152-cone residue重算中全部zero poles准确为0，Gπ³分别3/32、1/24、17/640。这些点只核实现，所有zero/repeated点的合法性来自上述全sum/rational函数与局部统一time控制，不来自三点或Astra的45点测试。

## 全star消去、整锥符号及全域g

以实际1152锥定义每ray star，独立逐项有理函数加法/gcd重算六个orbit representatives：A/C/E/F准确为0；B/D与锁定字典cross-multiply完全一致。全cone和rates在signed permutation群下等变，给S_Tr(x)=S_r(Tᵀx)和S_-r=−S_r，所以六个身份覆盖全部240ray，不漏人工projection pole。每对±ray贡献相同，将系数−1/32变为20条line reps的−1/16。

对12条B reps和8条D reps，在每个物理cone的四个ray上核projection同弱号且至少一个非零，故sign在严格cone内部恒定。case I的D符号8正，II为6正2负，III为7正1负；不是用producer的单点来假定整个region。新实现准确将二十个变换star全和与三个N17/D20 cross-multiply，全部一致。

对每case，D20的20个线性因子在物理cone代换后可定向成非负、非零线性forms，整体常数正，故D>0于u_j>0。固定g0=21/25000、πup=22/7（Arb256严格核π<πup）后，独立展开

\[
P(u)=[a^3N_{17}-g_0\pi_{up}^3D_{20}](Ru).
\]

三case非零系数数为1739、1702、1702，共5143，均非负，并逐项与源证书一致。因此generic interior内a³π³G≥g0πup³。cube中a≤1，得G≥g0。边界、D=0和zero/repeated投影点由上段的线外真连续/线上∞延拓；**没有仅凭lower semicontinuity从稠密下界反推边界值**。原点亦∞，整个闭cube闭合。

后续11Q coefficient-ratio sharper minimum未纳入035，未替换此固定门槛。

## CI4B连续majorant：全部残差准确重核

按全和真表达H=π³G=N17/D20，取齐次同阶对N20=a³N17、D20。每physical face tetra顶点x_i的首坐标均1，barycentric v_i≥0、Σv_i=1。新实现准确作N(v),D(v)的线性齐次替换。2440个D均非零且所有系数非负，故open tetra内D>0。

每cell的四个已锁Pi是非负有理数，P(v)=ΣPi v_i。新实现逐cell重算R=P D−NΣv，全部4,932,474个非零degree21系数非负；零项自然为0。于是Σv=1时0≤H≤P。真实H非负来自概率过程，不依赖N系数逐项符号。

凸性给P^(4/3)≤Σv_i Pi^(4/3)，且∫tet v_i=volume/4。因此

\[
\int_{tet}G^{4/3}\le\frac{volume}{4\pi^4}\sum_iP_i^{4/3}.
\]

π⁴来自(π³)^(4/3)，不是π³或π^(16/3)。点值、LP容差和作者修复过程均不承担证明职责；本审只核最终准确Pi和完整residual，从未导SciPy/优化。D=0的有限cell边界零3D测度可忽略，未删其邻域；全g问题已另行逐点处理，二者逻辑不同。

## 三physical tetra、shell与无叠覆盖

每physical cone取a=1截面，非负u_j变成barycentric λ_j=u_j a_j、Σλ=1，得到CI4B三根tet。其体积独立核为1/36、1/12、1/18，总1/6；覆盖来自逆cone map而非体积和。

case I做8轮最长边中点二分。新实现以精确平方边长、lex小edge处理ties，生成全部256叶子；每次中点的两个孩子在公共平面两侧并满父tet，独立核体积相加。

II/III各以q=(1,1,1)为corner的参考λ坐标，t=λ1+λ2+λ3。每shell s/2≤t≤s分为T1:t+λ1≥s、T2:t+λ1≤s≤t+λ1+λ2、T3:t+λ1+λ2≤s。三组明确barycentric系数分别是

\[
(\lambda_1+t-1,\lambda_2,\lambda_3,2(1-t)),
\quad(t+\lambda_1+\lambda_2-1,\lambda_3,2\lambda_1,2(1-t-\lambda_1)),
\quad(2t-1,2\lambda_1,2\lambda_2,2(1-t-\lambda_1-\lambda_2))
\]

（此处显示s=1参考壳）。新SymPy准确核和为1及坐标重构；各区域的非负条件与这些系数相互推出，且阈值有序，故是真正无重叠覆盖。缩放s=2⁻k、k=0…27覆盖2⁻28≤t≤1；各shell tetra再固定max(0,8−2k)轮二分。新实现生成的2440个有序顶点与证书全部准确一致，没有把证书顶点直接当作覆盖证明。

覆盖ordered体积准确为1/6−5/(36·2⁸⁴)，剩余仅II/III的两corner tetra。其非corner顶点到q距离≤√2，故t<2⁻28时r≤√2t<2⁻27。所有permutation images的剩余域均在同一个positive face的q球帽内。

## 真角帽和6/8因子

030已审n4 purebody给D=1/(48π⁴)、B=180224/π³及近q的G≤D+B/r。q与其他body-line face交点距离至少2，ε=2⁻27<1使该帽内q为唯一最近交点，故原r确实是|y−q|，不是误用非最近中心。

径向幂包络在全3D球内积分为E4(ε)。positive face的upper corner只占该球一个octant，所以其包络积分为E4/8；这里不要求G本身在q周围反射对称，只用径向上界。排序域cells乘6后，剩余所有images的并仍只由这个octant帽支付**一次**，不是6次，也不是64倍。

新256bit完整预算为

\[
I_4\le\frac6{4\pi^4}\sum_{2440}volume\sum_iP_i^{4/3}+E_4(2^{-27})/8.
\]

帽表达式约3.5270801566565500745e−9；完整上界表达式约0.0001841447010739849311214820756。准确有理上端与作者最终256bit上端相同，见 `035_integral_results.json`，并严格小于向外有理数I0=92072351/500000000000。此表达式的包围不是实际I的双侧区间。

## 合法核、一般L1/全部共同尺度及准确最终成本

固定g0=21/25000，c=g0/4=21/100000，定义

\[
K=(G-c)_+/[16(g_0-c)]=(G-c)_+/(63/6250).
\]

全cube下界给K≥1cube/16。026/031已审真实finite-frame Markov比较的每个假设在这里对应原四维八body frame：正自伴保守卷积、ψ可比|ξ|、H1 operator域、加权cutoff及weak-L4/3势。半群正部不等式给整个空间L(G−c)+=μ−δ0，μ正概率且无原点atom；包括全部body线、非光滑level界面和无限尾，非逐点源采样。weak-L势保证正部L1，不假称其有紧支撑。

对非负偶L1核、原点外正源，CND/Lévy闭包允许无限source，接触正性在同一u/σ/D下控制全部共同r。先有理尺度、再L1缩放连续及一般输入单调逼近，给一般L1弱型常数≤∫K，不是仅光滑/有界输入；未改成独立坐标多尺度。Lebesgue/Fourier/质量/Fubini及K*σ≤κ∫K均保留，不能移植坏背景。

由cube对称和face射线公式J=64I4，质量给(∫K)^3≤256I0³/g0⁴。准确值为

\[
256I_0^3/g_0^4=2275587701862626480057/708750000000000000.
\]

与候选的严格差

\[
(14753/1000)^3-256I_0^3/g_0^4
=212329869322269943/708750000000000000>0.
\]

所以C4≤∫K<14.753，低于旧16。新审计独立检查该有理比较和256bit cube root；没有用浮点14.75显示判突破，也没有沿用已纠正的14.37非正式猜测。

## 复现与最终范围

主导入后完整固定候选最小命令（stdout-only，不改旧结果；显式--save才保存新运行）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/035_check.py
```

两阶段也可分别运行 `035_representation_check.py`、`035_integral_check.py`；审计clone加 `--project /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds`。全部来源开头/结束复核，没有外部作者clone路径、优化或旧缓存依赖。典型合计约63秒。

MATLAB有限cross脚本仅交主串行核4D体积/峰值、physical tet体积、scalarπ校准和准确成本；实际状态见 `035_matlab_receipt.json`，不冒称第二套Matlab全2440系数执行。全证书g/I/source/all-L1/all-common-r链在本批已闭合，可由主据最终manifest更新U4；本审计不改账本、不提交推送、不扩招。后续11Q及三维最终验收另审。


主串行MATLAB已实际exit0（session5944），同SHA441ab8e3357714f553a8ac6eca8a39ff871a1c60fbc32ac6237abfd7519ea3d8。审计已读stdout，谱体积11/4608、peak11/3072、三tet体积、scalarπ校准及完整有理成本/差额assert均通过；log SHA99a3dbaab096a180cbdc5658f2702f07447f5de5112fd27483a92cea210c71ee。两Python阶段均首次实际PASS；便携wrapper的CLI另实际检查。收尾43来源和两checker/result关联hash均未变。035固定候选最终接受，本有限批已停止。

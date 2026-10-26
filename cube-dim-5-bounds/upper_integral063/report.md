# upper_integral063：四类区域的直接连续积分误差 pilot

本轮在固定 body16、alpha=1、Lebesgue、Fourier 逆因子 `(2*pi)^-5` 模型，以 `eps=1/2048`、`q=5/4`、三阶 jet、384-bit Arb，完成四个中心各一次全部3504锥的截时势生产。**040绝对 resolvent 尾机制能认证这些局部盒，但当前误差规模不支持宣布全域连续积分计算可行。** 粗盒虽可在 guard 通过时得到有限上界，误差可比中心值大数千至数百万倍；实用相对误差需要很细的四维尺度。

只写独立 clone `Math-Research-cube5-upper-audit-cont/cube-dim-5-bounds/upper_integral063/`。旧040、059、052、061及所有其他目录只读；无MATLAB、mesh、提交或推送。工作 clone 基准 HEAD `3c73f4e7cd65d2f7671001be41de1f228f155b1c`。任务限制禁止扩展六维，本轮未扩展。

## 1. 预登记、对象与资格

先登记 `registration.json`，再启动生产 PID59003。四个 nominal anchors 完全按任务给定；实际积分区域始终是

\[
 B_h=\{1\}\times\prod_{j=1}^4[\max(0,p_j-h),\min(1,p_j+h)],
 \qquad h\in\{2^{-6},\ldots,2^{-16}\}.
\]

这里 p 是 nominal anchor，不是 jitter 中心。中心只有提取一次普通DD jet的作用，盒包含所有原 anchor，包括原 zero/repeated 点。每个实际中心均在最小登记盒内，因此在全部11个盒内。没有删除任何零投影或投影重合位置。

有限候选依次为 `p + j/(104729*32)`，尾坐标 j 依次是 `(0,0,0,0),(-1,-2,-4,-8),(-1,-3,-7,-11),(-2,-5,-11,-19),(-3,-7,-17,-29),(-2,-11,-23,-43),(-5,-13,-31,-47),(-7,-19,-37,-49)`。候选及全部准确资格计数在登记中冻结；只取第一个280互异且落在最小盒内的中心。生产前的候选准备曾发现最初较短名单无 near-S1 合格点，随后在尚未产生任何jet时扩充为上述有限名单并登记；没有生产后的候选搜寻。

|类|nominal互异投影/零投影数|实际采用j|实际D|中心 H_eps 约值|
|---|---:|---|---:|---:|
|ordinary|161/4|`[-1, -2, -4, -8]`|123999136|0.0002610118414|
|low|280/0|`[0, 0, 0, 0]`|104729|0.0001223831996|
|nearS1|161/6|`[-2, -11, -23, -43]`|26810624|0.08774066555|
|crossing|95/6|`[-2, -5, -11, -19]`|10053984|0.00011971835|

原 crossing anchor `(1,1,1,2/3,1/16)` 有6个零投影。例如真整数顶点 `(0,-9,9,0,0)` 与其相反顶点同时投影为零；原 anchor 在每个真实盒中。这是准确的零/重复节点见证。实际中心的280个投影互异只是普通中心DD的资格；**盒内投影重合的有效性由下节完整时间 resolvent 余项支付**，不靠普通DD在盒内继续有效。

生产前的 PID58937 因 `Fraction(fmpq)` 类型转换在 guard 阶段报错，零个锥、零个jet；修复为准确字符串有理转换后，PID59003才进行唯一实际生产。失败receipt保留于外置工作目录，摘要如实记录。源码修改发生在任何jet之前。

## 2. 完整时间的机制与证明

只消费017已锁的280整数顶点 `z=144r`、3504锥及正行列式分子 d_C。040的准确DD helper仅被导入，不运行其main、旧点或旧截时参数。对实际generic中心 c，令 `k_v=D c dot z_v`、`l_v=D z_v,tail`、`a=144D eps`，有

\[
 H_\epsilon(c+(0,u))=F\Re\sum_Cd_C\int_a^\infty
       \prod_{v\in C}(s-i k_v-i l_v\cdot u)^{-1}ds,
 \quad F=\frac{D^4}{144(2\pi)^5}.
\]

由中心互异节点的 `-DD_4(.5 log(a²+k²))` 提取35个总阶数≤3的 Taylor 系数。DD递推、合并log权重均用准确 fmpq；只在全3504锥合并之后以384-bit Arb求log和F。没有按h重算jet。

设真实矩形盒相对c的端点为 u_j^-、u_j^+。每个顶点的准确线性范围为

\[
 d_v^- =D\sum_j z_{vj}(u_j^-\mathbf1_{z_{vj}\ge0}+u_j^+\mathbf1_{z_{vj}<0}),
 \quad d_v^+=D\sum_j z_{vj}(u_j^+\mathbf1_{z_{vj}\ge0}+u_j^-\mathbf1_{z_{vj}<0}),
 \quad n_v=\max(|d_v^-|,|d_v^+|).
\]

**guard直接准确比较** `n_v² < a²+k_v²`，不是用展示小数或保守浮点失败代替真实排除。最大精确比率、失败顶点及其坐标全部存于各results。若任一点失败，本批几何级数余项不适用该盒；这不说明真实H不光滑或没有其他上界。

成功时取准确有向下界 `a <= b_v <= sqrt(a²+k_v²)`，令 `q_v=n_v/b_v<1`、`q_C=max_(v in C) q_v`。对所有 `s>=a` 和整个盒，扰动resolvent比率不超过q_v。五因子乘积的三阶以上尾由

\[
 T_3(q_C)=(1-q_C)^{-5}-\sum_{n=0}^3\binom{n+4}{4}q_C^n
\]

控制。设

\[
 W_C\ge\int_a^\infty\prod_{v\in C}\min(b_v^{-1},s^{-1})ds.
\]

040 helper按排序b的断点准确积分有理幂与log，并取有向有理上端。因为 `|s-ik|>=max(s,sqrt(a²+k²))`，这是完整时间绝对积分上包。于是

\[
 |H_\epsilon-P_3|\le R_h:=F\sum_Cd_CW_CT_3(q_C).
\]

s>=a>0，所有真实空间投影包括零值和重合值的resolvent分母均非零，且无穷尾可积；因此本证明在整个盒、全部 `[eps,infinity)` 成立。不以盒内投影间隔作为额外前提。

范围采用 `V=sum_(alpha!=0)|A_alpha| prod_j max(|u_j^-|,|u_j^+|)^alpha_j`，得到 `L=max(0,A0-V-R)`、`U=A0+V+R` 的向外端点。正性来自固定正Cauchy混合。对非对称盒，准确体积及归一化矩是

\[
 |B_h|=\prod_j(u_j^+-u_j^-),\qquad
 \mathbb E_{B_h}u^\alpha=
 \prod_j\frac{(u_j^+)^{\alpha_j+1}-(u_j^-)^{\alpha_j+1}}
 {(\alpha_j+1)(u_j^+-u_j^-)}.
\]

这包括jitter导致的非零一次矩。不是把旧040对称均值硬套过来。由 `m=mean(P3)+R` 及统一区间[L,U]，凸弦给

\[
 \int_{B_h}H_\epsilon^{5/4}\le |B_h|\left[L^{5/4}
 +\frac{U^{5/4}-L^{5/4}}{U-L}(m-L)\right].
\]

这里使用向外均值上端，并与U截取；所有 powers 用Arb包围。所得为 **H_eps的局部积分上端**，不是完整G、不是I、不是范数。

几何核对例：crossing每个登记h的前两维长度h、后两维2h，故准确体积 `4h^4`，不是旧对称盒的 `16h^4`。near-S1在h=1/64时的四条边长分别为 `5/256,6/256,7/256,8/256`，体积 `105/268435456`。其余准确裁切与全部35个矩已存档。

## 3. 有限11级的结果

“最大”仅指预登记的11级、固定中心、固定三阶、固定余项测试。分三种门：guard通过即可给L=0时的有限上界；正L门要求有非平凡下界；1%门另要求 `R_h <= H_eps(c)_lower/100`。1%是生产前登记的可用性指标，不是全域积分预算。

|类|最大guard h|最大正L h|最大1%余项h|前一更粗guard尺度的精确q²约值（准确值见JSON）|
|---|---:|---:|---:|---:|
|ordinary|1/512|1/4096|1/16384|1.234718789>1|
|low|1/256|1/4096|1/16384|2.144682256>1|
|nearS1|1/512|1/8192|1/32768|1.600403217>1|
|crossing|1/256|1/4096|1/16384|1.000229163>1|

最大guard尺度（全部数值仅展示；可用端点为JSON准确有理数）：

|类 / h|L|U|mean H上端|R_h上端|凸弦积分上端|
|---|---:|---:|---:|---:|---:|
|ordinary / 1/512|0|2.427679305|2.427678021|2.427417009|7.055520438e-10|
|low / 1/256|0|209.601807|209.6018056|209.6016832|2.971008515e-06|
|nearS1 / 1/512|0|16545.81588|16545.77897|16545.69125|4.369174637e-05|
|crossing / 1/256|0|10.79303876|10.79303795|10.79291762|1.821921951e-08|

最大正L尺度（全部数值仅展示；可用端点为JSON准确有理数）：

|类 / h|L|U|mean H上端|R_h上端|凸弦积分上端|
|---|---:|---:|---:|---:|---:|
|ordinary / 1/4096|0.0001330404593|0.0003889832236|0.0003888226449|0.0001278109979|3.103791393e-18|
|low / 1/4096|1.785922598e-05|0.0002269071732|0.000226820367|0.0001044371717|1.582404894e-18|
|nearS1 / 1/8192|0.04397640058|0.1315049305|0.1298449094|0.04199300277|2.77365138e-16|
|crossing / 1/4096|8.20192431e-05|0.0001574174569|0.0001573683964|3.761228958e-05|2.504833081e-19|

最大1%余项尺度（全部数值仅展示；可用端点为JSON准确有理数）：

|类 / h|L|U|mean H上端|R_h上端|凸弦积分上端|
|---|---:|---:|---:|---:|---:|
|ordinary / 1/16384|0.0002605103517|0.0002615133312|0.0002614727653|4.611195983e-07|7.382838004e-21|
|low / 1/16384|0.0001219785535|0.0001227878457|0.000122766163|3.82963659e-07|2.869383888e-21|
|nearS1 / 1/32768|0.08704016325|0.08844116784|0.08802730946|0.0001749170804|6.654196342e-19|
|crossing / 1/16384|0.0001195540958|0.0001198826042|0.000119870261|1.426622696e-07|6.962569866e-22|

三个表同时呈现有限但很松的证书和较细的有效误差；没有只展示一个微盒成功。各results还保留每个通过尺度的mean下端、完整range、变化V、常数范围积分、余项及相对余项。

## 4. 细化趋势、瓶颈与算法局限

四类所有guard通过尺度上的R_h都严格随h减小，读回按准确有理上端检查。到最小h=1/65536时，相对余项分别约
`ordinary: 7.126939014e-06`, `low: 1.203814831e-05`, `nearS1: 0.0001436593256`, `crossing: 4.92372744e-06`.

在接近中心且q小的范围，几何级数尾从四阶起，有限结果接近四阶收紧。jitter固定后，盒持续缩向nominal而非generic中心；不能据本批有限11级声称继续任意细分R必趋零。更深细分应重新登记合法generic中心或换合流DD实现；本轮没有新增中心。

瓶颈逐类相同但程度不同：ordinary在最大guard尺度R/H约9300，low约1.71e6，near-S1约1.89e5，crossing约9.02e4。即使进入1%门，余项仍远大于384-bit算术半径；升精度不能解决这个损失。near-S1的中心本身H约0.08774，尖峰与较大的绝对cone余项共同要求更细尺度。crossing按完整时间证明通过，但其边界裁片不能当成完整对称盒。

guard本身没有把四类所有非微尺度全排除：仍到1/256或1/512；然而这些最大guard证书的范围/均值被R支配，并不表明这样的空间步长能有效积累预算。有效的正L/1%门落在1/4096～1/32768。本pilot只有四组局部嵌套盒，没有区域覆盖、没有空隙总账、没有已认证的全域U_H；不把局部成本外推为全域可接受。

下一需更强误差的数学对象是：

1. 保留锥间抵消的局部四阶积分余项/区间log组合，或直接对H_eps^(5/4)给片上余项。当前 `sum_C d_C |time remainder_C|` 丢失极强抵消；需要证明抵消在完整盒/片上持续有效，而非只在中心合并。
2. 可供空间积分消费的 **积分型余项** `int_B |H-P3|`（或有正权的高阶导数平均），及其与power上包的合法连接。不能仅把它直接替代现有点态凸弦所需的range；需一并改善range或直接power误差。
3. near-S1/边界裁片的正混合局部时间权矩或分层解析尖峰包络，减轻统一eps导致的eps负幂损失。043/047的全空间log曲率是合法备用，但在本eps和粗h上二次指数很大，并未在本轮冒充有效均值积分界。

下一步应先建立并用现有四组jet检验这些更强对象，再决定有限区域分解。没有授权或启动全域mesh。

## 5. 同eps真实短时预算与条件成本

034完整短时公式为

\[
 \|G-H_\epsilon\|_{L^{5/4}(F)}\le E_\epsilon
 =\frac56\left(\frac{14455}{31104\pi^5}\right)^{1/5}
 \pi^{-4/5}\epsilon^{6/5}.
\]

本轮实际Arb计算 `E_(1/2048)=9.6781736666585742...e-6`，这是 **Lq范数**。若连续覆盖认证 `int_F H_eps^q <= U_H`，则 `I <= (U_H^(4/5)+E)^ (5/4)`。以目标 `I<=1/25000`：

- 无另付管预算时，正规截时积分门为 `U_H <= 3.84101186903678435...e-5`；正式充分门取budget.json该gate的准确下端，不能使用其上端或四舍五入小数。
- 另付完整G管积分 `1e-6` 时，门变为 `U_H <= 3.74182806274712524...e-5`，同样消费准确下端。这只是条件组合，不是本轮取得管/覆盖。
- 若误把目标4e-5直接作为U_H，则真实合成仅给 `I<=4.16026236522835679...e-5`，超过目标。

目标 `(g,I)=(1/10000,1/25000)` 的条件成本满足 `U^4=800000<30^4=810000`，故 `U≈29.9069756244<30`。这只是准确的条件成本检查：**g≥1e-4与I≤1/25000均未在本轮取得；既有C5≤32不变。** 没有把U_H当I、把norm当积分、把正部质量的128因子放入I槽，也没有更改029的一般L1与共同尺度传递模型。

## 6. 实际运行与读回

实际唯一生产PID59003，四中心各3504锥、35个jet系数，总计14016锥，墙钟50.137486375秒；python-flint0.9.0，384 bits。每128个完整锥保存准确stencil前缀，源/登记/输入fingerprint锁定；120秒检查点未触发，未续算、未重算旧040或旧H1/64点。每中心准确stencil/cache在仓库外，归档中保留完整35系数包围与所有11级必要结果。

`verify.py` 的实际读回PID与源码hash见readback.json：重新准确核对候选资格、全部盒含原anchor及中心、裁切体积、280顶点guard、35个非对称矩、范围/均值顺序、有限最大h门及严格下降的余项序列；不运行任何H/jet/DD/log。它从系数有理端点重建保守范围/均值，再仅用384-bit Arb重做power弦，存储独立读回上端。独立外向端点可能比生产内部球多扩几个末位；两份弦值差的准确读回门是2^-280，并非新的积分误差门。端点形式保留，不用展示小数重造证书。

证明状态：040解析机制、准确裁片公式和有限数值证书已在本轮落实，读回通过；没有另一个独立数学审计代理，没有全域结论。源与库是共享依赖，hash不是数学证明的替代。

外置cache：

|类|字节|SHA256|
|---|---:|---|
|ordinary|4291635|`11ad4e0ddd6322297316b0b09f72a9167e314389faab19d7c88a7d8fd7d25708`|
|low|2700662|`49f56b1a23a7814b5f22257672d345b763f64d80048cd65caaaa4caf03ede9b5`|
|nearS1|3873680|`dfe208754f66f0e75bf95794dd49e5e3d67fe167ee80140d6c49ff5e0bbb207e`|
|crossing|3841006|`f11a3b3fe7ddce57fded688c44c0dc399b4f3e3d884ec0b71a3771478ddc1a7b`|

外置路径前缀 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_integral063/`；每文件名为 `<class>_checkpoint.json`。

## 7. 来源、复现与白名单

|只读源|SHA256|
|---|---|
|`/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz`|`8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9`|
|`/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/upper_cutoff040/cutoff_box.py`|`faa5c389976dd21d66026bbcabf8903a8e004a9342987e3dac97dffc7124aaf2`|
|`/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/upper_cutoff040/report.md`|`122cf1dd000b219e3d389d03cf8c196bec0374e62ae44b600cec017a5e3e3196`|
|`/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/upper_alpha1_034/report.md`|`d7360aaf14e935fe242a52cc7f006197925ddbd2693996c3d58c205849fab4a9`|
|`/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/upper_harnack043/report.md`|`a298ca56ba980fce4fe2542929e03befb97b528afd3d390ed9d15651da68322a`|
|`/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/audit_harnack047/report.md`|`0b680c1bd7e4aea30c027fe32cc3e51436f7dda6dd8b74ddf521c2d74bb8ca16`|
|`/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/upper_alpha1_029/report.md`|`298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81`|

生产源码SHA256 `8d93b7a02898369a5b15517441aabffef322c8b4a6eb0bfa63e57feb7882b950`；登记SHA256 `6c60629dec02b7612441b7b15d8a4430ae9a7987dffe29311de204e9a48ba15e`；生产fingerprint `3d9b0257ad1be653f0d6f3f4ab55a484c61f9eda739e0a86a0fc857c0e04148b`。

只读检查（不重算jet）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_integral063/verify.py
```

确需重建时，**先复制本目录到新的回放目录**，使输出不会覆盖本轮冻结证书；保留pilot.py字节。在复制品中移除已有results与run_receipt，运行如下命令，cache必须是新的仓库外目录：

```sh
python -B /new/replay/upper_integral063/pilot.py --work-dir /new/outside-git/cache063
```

依赖python-flint0.9.0。输入登记保存原始绝对路径，移机可仅在回放登记中把sources.source映射到相同SHA字节（040helper、017表及其余只读证明），不可换同名不同字节文件；这会自然形成新的回放登记hash/fingerprint。并不要求公开读者具有作者私人路径，也不复制全部旧输入。只有需重建时运行生产器，本轮没有为验证而重播。

归档白名单：`report.md`、`registration.json`、`pilot.py`、`verify.py`、`run_receipt.json`、四个 `<class>_results.json`、`budget.json`、`readback.json`、`actual_summary.json`、`manifest.json`。预jet guards已并入results，删除本轮重复guards输出与预登记生成草稿；外置可重建stencils、失败receipt和临时日志不入白名单。未删除/改动任何旧目录文件。未提交。

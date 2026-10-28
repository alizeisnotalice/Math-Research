# 第36轮：模类别joint-cell来源的有限weak筛选与宏观混合残差

2026-10-29。已完成真实source-cell模分类、全格求和/循环DP/有理控制及12个冻结小n模型的NumPy筛选。没有追加seed、N或改变标签。**有足够ESS的best-weak筛选没有显示随R增长的收益；joint与实际prior通常接近，均低于最佳固定表。** 同时，旧宏观窗口上的joint/prior残差随n显著减小，n15已接近机器舍入层。这是本轮更有信息的机制读数，严格全域估计由主线另行处理。

没有给出fullweak总体排序、一般 \(C_n\) 或log改善结论。一个高阈值低ESS候选也保留，未因筛选阴性而删除。

## 1. 真正joint source与质量

K固定12，R=2、3、4、6，n=3、7、15。每个source粗格 \(\mathbf j\) 的质量为 \(K^{-n}\)，选表
\[
 b(\mathbf j)=\sum_i(j_i\bmod R)\pmod R.
\]
K是R倍数，uniform粗格的每个单坐标residue均匀，因此无论其他坐标如何固定，b均匀于R类。**实际prior精确为 \(\pi_b=1/R\)**，不是用样本频数估出来的。

seed202610290360只用于事前抽一次并冻结permutation

    [11,6,4,5,1,7,2,3,9,0,8,10].

取 \(s_j=\mathrm{permutation}_j\bmod R\)，第b表
\[
 u_b(j)=\frac{(s_j+b)\bmod R+1/2}{R},\qquad
 a_b(j)=\frac{j+u_b(j)}K.
\]
每个joint格只放一个点 \((a_{b(\mathbf j)}(j_i))_i\)，再以乘积box宽 \(\delta=1/(64KR)\) 平滑。点距格边至少1/(2KR)，故box留在本joint格；各格质量、总质量严格正确，密度非负且属于L1。峰值为 \((64R)^n\)。

比较对象包括joint、全部R个fixed表的乘积源、实际prior正混合。joint抽样必须先抽source格再算b；prior先独立抽b，再uniform抽格，两者没有混淆。

全空间有结构帽
\[
 F_{\rm joint}\le\sum_bF_b=R F_{\rm prior}.
\]
它不构成增益定理，也不能由R增长本身推出log。本轮不把固定R2比较延伸成发散结论。

## 2. 真实共同cube与residue循环DP

令每坐标、每表、每source residue的权重为
\[
 z_{ibr}(x,L)=\frac1{KL}
   \sum_{j:\,j\bmod R=r}
       \int_{x_i-L/2}^{x_i+L/2}\beta_\delta(t-a_b(j))\,dt.
\]
对每张表b，循环卷积n个R维权重，得到
\[
 c_{bq}=\sum_{\sum_i r_i=q\pmod R}\prod_i z_{ibr_i}.
\]
真实joint cube平均是 \(\sum_b c_{bb}\)；fixed b平均是
\(\prod_i\sum_rz_{ibr}\)，prior取fixed平均的1/R。所有坐标和分量共用同一个L，没有query依赖的选表。

这是未归一化的非负循环DP，不再乘一次1/R，也不把各classifier分支重归一。向量化实现成本约 \(O(NnRK+NnR^3)\)，没有枚举 \(K^n\) 的大n成本。

## 3. 先运行的小例和有理控制

n3/K12、R2首先全枚举1728格与DP核对；随后同样核验R3/4/6。每R检查3个中心、5个长度，包括source点处δ短尺度、1/K、L=1与7/6。最大相对误差：

|R|全部枚举格数|每类格数|最大相对差|
|---:|---:|---:|---:|
|2|1728|864|6.63e−16|
|3|1728|576|8.29e−16|
|4|1728|432|2.22e−16|
|6|1728|288|8.34e−16|

另对n2/K12、全部四个R，使用Fraction在3/8、7/6、1三个有理L及中心(1/7,4/7)做全144格有理求和和有理循环DP，全部精确相等。例如R6、L=7/6，joint为195/196，prior为1175/1176；L=1两者均为1。此类有限点差异只验证真实输入不同，不证明weak增益。

## 4. 冻结成本、query与importance sampling

所有模型grid均62点：短尺度从δ/16几何延伸至1/(2K)，再线性至1.25；另含预定cell尺度、L=1以及10/9、7/6、11/9。L1提供finite-query maximal至少1的确定锚。冻结后未因收益修改grid。

最坏n15/R6、batch128全grid定标0.118秒，按PLAN规则冻结每case N4096，12点保守预计59秒；实际MC与Fourier诊断共14.91秒。主要数组只有batch×n×R×K及循环DP，远低于2GB；不生成或保存 \(12^{15}\) joint点表。

proposal为50%uniform，加50%uniform选择全部R+2个来源及全部62个L。先按真实source抽box点T，再抽独立cube核U，query \(X=T-U\pmod1\)。实际RN由DP计算全部来源/全部L平均，W≤2，且任意尾系数贡献至多 \(2(R+2)J\)，全部样本实际检查。

seed固定202610290361至202610290372，顺序n3/7/15，每n内R2/3/4/6。质量权重均值范围0.99131至1.00635，描述性SE约0.00461至0.00560；没有对偏低的n15/R6重新抽样。source质量的严格理由仍是每格K^−n和box质量1。

## 5. best weak筛选：没有可靠R增长信号

下表为尾ESS≥64的经验best weak；fixed列取该R下所有固定表的最佳值：

|n|R|joint|实际prior|最佳fixed|
|---:|---:|---:|---:|---:|
|3|2|1.080|1.072|1.180|
|3|3|1.041|1.037|1.234|
|3|4|1.023|1.058|1.238|
|3|6|1.029|1.023|1.261|
|7|2|1.081|1.083|1.246|
|7|3|1.057|1.059|1.291|
|7|4|1.029|1.030|1.305|
|7|6|1.014|1.016|1.313|
|15|2|1.104|1.104|1.316|
|15|3|1.069|1.070|1.409|
|15|4|1.034|1.034|1.384|
|15|6|1.016|1.019|1.427|

多数joint最佳层在阈值约1.02–1.56、ESS2600–3700，不是单样本造成的值。先用前半样本的ESS≥32扫描选阈值，再用其严格左侧浮点阈值在后半样本验证；heldout joint系数在1.009至1.094，描述性SE约0.0073至0.0110，结果没有提供R增长优势。

### 必须保留的稀疏高阈值候选

n15/R6的**未过滤**经验扫描有joint=1.42794，阈值约 \(7.494\times10^8\)，尾ESS仅13.93；最佳fixed1为1.42656、ESS3164。两者数值相近，但joint候选远不稳，注册的分半支持层未验证它。另有少数case的prior/fixed也出现低ESS高阈值峰。

因此不能根据ESS筛选表声称fullweak阴性或总体排序。本轮如实保存未过滤扫描，停止而不追seed。只对选定heldout阈值给有界Bernoulli-KL区间；未计算全阈值统一上CI，不将局部筛选当完整weak排序。

## 6. 旧宏观窗口：实际残差迅速变小

每个(n,R)事前固定64个uniform中心，检查10/9、7/6、11/9三个L。下表是在这些样本上的最大
\(\lvert A_LF_{\rm joint}/A_LF_{\rm prior}-1\rvert\)：

|n|R2|R3|R4|R6|
|---:|---:|---:|---:|---:|
|3|4.55e−4|8.45e−4|1.16e−3|2.95e−3|
|7|7.58e−9|2.23e−8|1.33e−7|8.08e−7|
|15|8.88e−16|2.00e−15|1.33e−15|5.44e−14|

这不是全x/L上界。n15的前几列已接近直接相减的浮点floor，所以另计算非零residue Fourier乘积，而不是将零差当恒等。

令
\[
 \widehat z_{ib\ell}=\sum_rz_{ibr}e^{-2\pi i\ell r/R},\quad
 T_{ib}=\widehat z_{ib0},\quad A_b=\prod_iT_{ib},\quad
 \phi_{ib\ell}=\widehat z_{ib\ell}/T_{ib}.
\]
在L>1时T>0，循环DP给以下可直接核验的残余公式：
\[
 \frac{A_LF_{\rm joint}}{A_LF_{\rm prior}}-1
 =\frac{\sum_bA_b\sum_{\ell=1}^{R-1}
       e^{2\pi i\ell b/R}\prod_i\phi_{ib\ell}}
        {\sum_bA_b}.
\]
代码同时记录右侧及将各非零乘积取绝对值后的加权帽。与直接DP比值的最大差≤2.00e−15。

n15/R6，稳定Fourier相对残差最大5.40e−14，非零乘积加权帽1.09e−13；n15/R4对应6.74e−16、1.19e−15。n15/R2的这些固定中心里乘积残差全为0，而直接比值约9e−16；这里只登记样本中的零/舍入floor，不称来源恒等。

完整周期内，各source residue的质量恰1/R，其非零characters相消；L>1时非零部分来自剩余段。这解释为什么该诊断比单纯小n weak读数更值得研究。一般全域收缩估计由主线负责；本报告不把64个中心的最大值升级为严格障碍。

## 7. 结论边界和停止

本轮验证了R增长的真实mod-cell输入和循环DP，实际prior、质量、source抽样及共同cube尺度均正确。它没有显示可靠的best weak增长；短尺度还存在少量低ESS候选，宏观旧窗口则表现出很强的residue混合。

不追加R、n、seed或改变标签追求阳性，不把固定R结构帽或R增长当log增益，不做大n正态近似。DP可计算不意味着任意joint classifier可计算：本模型只需R个residue状态；一般source格函数仍可能需要指数状态。本次结果不推渐近阶，也不排除其他分类/短尺度机制。

所有空间box积分、有限求和和循环DP有准确公式；除有理小控制外为普通浮点，未做区间舍入认证。periodic来源可用整数胞截断精确归一化转为 \(\mathbb R^n\) L1输入，内部中心损失仍可控制；本轮不以质量损失换取收益。

## 8. 复现

先运行 python3 probe.py --calibrate 完成直接/Fraction控制、成本门槛和注册，再运行 python3 probe.py 做冻结MC与Fourier诊断。实际Python3.15.0、NumPy2.5.4；脚本哈希在registration中锁定，拒绝覆盖已有结果。定标只按成本选择N4096，不能按收益调样本。

白名单为REPORT.md、PLAN.md、probe.py、registration.json、calibration.json、RESULTS.json。RESULTS约66KB，保存全部fixed表比较、支持/未过滤扫描、heldout、质量、贡献帽及36组宏观Fourier诊断。无空间大数组、cache或重复初始小例。

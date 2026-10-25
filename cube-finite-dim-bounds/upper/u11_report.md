# U11报告：四维真实body-frame Green接口

2026-10-25。base d5f647caec30f0006d3dee7c95855b604dbdd05b；只upper/u11_*，不改bridge/旧源/总账，不提交推送扩招。本有限批只表示/校准，不做成本、g/J或权重搜索；U4=16、U3=8未由本分支更新。

## 准确成品

完整推导u11_representation.md SHA63f759e850f3864390cf2af6bf39c57ead78130421389a2cfc5e853335063ea5。排序a≥b≥c≥d≥0下body谱为max(8a,4(a+b+c),2(3a+b+c+d))。增量y=(a−b,b−c,c−d,d)使三个分区本身simplicial；六rays A..F及cones AB EF、EBCD、E B D F的det为1、1、2。解析非负参数唯一表示给覆盖/不交，不以体积或点测替代。

384 signed permutations给1152 cones、240不同triangulation rays。排序Σξ=1的三个截面tetrahedron比例1/2、1/4、1/4；全ψ单位球真正极点仅B/8 orbit24及D/12 orbit16，共40。准确体积11/4608，故

p_t(0)=11/(3072pi⁴t⁴)。

正密度完整公式为(2pi)^−4 Σ_cones |detR|Π(tψ(r_j)−i x·r_j)^−1，负实cone项不能删。

Green时间积分从partial fractions推导：generic非零节点b_j=x·r_j/ψ(r_j)有

G=−(32pi³)^−1Σ_cones [|detR|/Πψ(r_j)]DD3(sign;b_j)。

π/12测试b=(-1,-2,1,2)与π/4重复测试(-1,-1,1,1)校准符号/正常数，非移植三维公式。

## 全部零/重复节点的接口

稳定源u11_green.py SHAfdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d。

- `green_float_generic(x)`：generic float discovery；近zero/repeat抛错，不返回0，不是严格数值区间。
- `green_pi3_coefficient(x)`：非零nodes允许confluent重复，返回准确有理C，G=C/pi³。
- `green_allnodes_pi3_coefficient(x)`：全部rational点的zero/repeated通路；先完整1152锥合并confluent partial-fraction coefficients i^(k−4)R_(b,k)。离active frame lines时p_t小time有界，强制所有zero pole准确消失；非零高阶pole的时间实积分为0，只剩simple residue sign和。任何未消zero pole报错，不能硬塞sign0。
- `green_cutoff_mp`：正time ε的全锥复log备用公式，任何nodes合法。遗漏0..ε时间有严格解析界3670016ε/(pi⁴dist(x,S)^4)，但mp本身roundoff未封闭，不把该tail误当总数值误差。

原点和八条body方向线G=∞，其外连续有限；node-zero hyperplanes通常不是物理奇性。唯一补登记退化符号点e1得到C=3/32，全部zero-poles准确消失；fullconfluent time积分相对差约1e-61。未重复主后续triple-edge点或成本pilot。

## 固定有限运行

12个预登记rational点为x_i=a_i+1/p_i，primes151,157,163,167及base表存prereg/results。prime模算术保证任何ray-zero或normalized-node coincidence不可能（整数系数绝对值≤144）；并非只靠有限测试宣称generic。

Python Fraction/SymPy实际核几何/volumes/ranks、1152 cones唯一性、全部12点exact generic/allnode等价、复log ε=2^-40、density与residue交叉。首两点50digit完整time积分相对差约3.93e-52及6.76e-52。12点float快式相对exact误差最大约1.29e-11；这不是全域误差界。全部固定诊断PASS。

最小检查源u11_check.py包含可重生成的原始固定点和退化测试。实际经历initial same12、epsilon Fraction输入兼容修补same12重跑、allnode实现后same12等价+预登记e1重跑；每次约2–3秒，没有新增generic点/参数/成本。重复是实现修补和新增必要节点通路的验证，不是扩大模型。一个零节点标量例(0,0,1,-1)显示逐cone积分发散，用于拒绝错误sign0。

完整geometry SHA6dc064dcdb05fc034bb9cfae177cf8a97d8fc48e941768a05dabc6f9900c7c24。详细exact fractions、所有点结果、运行源码hash在u11_results.json。

MATLAB最小校准u11_matlab.m SHA5b072aeda506cc73ef8ac119ebe980c85081390d99a91fc0238dd38b4e2e8b62已交主串行：独立重建1152 cones、base det、谱体积/p0、首genericG、π12/π4与zero-cutoff标量积分。upU11前缀，无clear/path；主串行MATLAB实际exit0/session88833，脚本hash匹配：1152 cones、volume11/4608、p0 coefficient11/3072通过；π12/π4 defect0，zero-cutoff defect8.89e-16，generic首点relative defect2.46e-12。完整stdout及log hash已入results，不把容差当区间。

## 正face权假设及适用范围

tau≥0的face附加符号在排序chamber为2tau(3a+2b+c)，不改变本三cone剖分。A..F rate分别8+6tau、8+10tau、12+12tau、12+12tau、16+18tau、24+24tau。body basis det16维持正保守密度、t^-4上界/ψ~|ξ|、H1域及cutoff；真实正部/finite-frame比较的一般L1和共同全部r假设依法保持。tau>0的active singular lines还包含face方向，不得继续只查body lines。未执行任何tau≠0点/成本；40谱球极点表仅tau0。

固定Lebesgue实际用于Fourier逆正常数、概率密度/普通卷积、谱体积与物理体积、026的能量/质量/接收体积及cube归一。此表示不适用于坏背景或独立矩形半径；不将二维3.879乘起来。

|类别|本批状态|
|---|---|
|准确分析/代数|完整谱/截面/全锥覆盖、真实G及全部节点源稿；exact rational校准PASS。新表示独审排队，不借026替代|
|实际数值|固定12 generic+e1诊断，两个完整time quad；主串行MATLAB实际exit0，generic相对差2.46e-12|
|未认证|全4D cube g下界、J上界与成本；任意泛型float误差全域界|
|未实现|本分支无U4<16/U3<8；不排除/认证任何成本族|

主后续候选/连续认证另登记、另负证明责任，本分支不重复。首批到稳定表示接口停止。

## 白名单

仅9文件：u11_inputs.json、u11_preregistration.md、u11_green.py、u11_check.py、u11_results.json、u11_representation.md、u11_matlab.m、本报告、u11_manifest.json。manifest列其余8文件hash/大小，不自记自身hash。__pycache__不归档。待正式独审/主MATLAB回执只更新report/results/manifest状态，锁定源有修补另报hash。

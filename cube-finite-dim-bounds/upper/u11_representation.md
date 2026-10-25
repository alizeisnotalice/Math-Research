# U11：四维纯body Cauchy的准确谱锥与真实Green表示

2026-10-25。固定V={(1,s2,s3,s4):sj=±1}八方向，L=−Σ_v|v·D|，ψ=Σ|v·ξ|。方向未单位化，所有时间相同。Fourier约定fhat=∫e^-iξx f，逆因子(2pi)^-4。本文件证明表示，不做g/J或成本。独立正Cauchy变量和给真实正保守半群；026/U09的一般L1/全尺度/正部势核输入不由有限点校准代替。

## 1. 排序谱公式及三个区域的准确覆盖

ψ在全部signed permutations下不变。取a≥b≥c≥d≥0。八绝对值中六项非负，总和6a+2b+2c，余项为|a−b−c+d|+|a−b−c−d|=2max(|a−b−c|,d)。因此

ψ=max(8a,4(a+b+c),2(3a+b+c+d))。

第一支在a≥b+c+d，第二支在b+c≥a+d，第三支在|a−b−c|≤d；交界值一致，三个区域内部不交、并为整个排序chamber。

令y=(a−b,b−c,c−d,d)≥0。三个区域依次是y1≥y3+2y4、y3≥y1、0≤y1−y3≤2y4。它们本身已经simplicial，不需要猜额外三角剖分。写

A=(1,0,0,0)，B=(1,1,0,0)，C=(1,1,1,0)，D=(1,1,1,1)，E=(2,1,1,0)，F=(3,1,1,1)。

|区域|非负参数表示y|x-space四rays|absdet|
|---|---|---|---|
|I|(u1+u3+2u4,u2,u3,u4)|A,B,E,F|1|
|II|(u1,u2,u1+u3,u4)|E,B,C,D|1|
|III|(u1+2u4,u2,u1,u3+u4)|E,B,D,F|2|

每个u唯一且非负，故这是覆盖和不交的解析证书，非体积/采样猜测。384个signed permutation chambers复制，共1152个full-dimensional cones，内部不交；重叠仅排序/谱边界的零4D体积。

## 2. 截面极点、全谱球极点及准确正常数

在排序chamber以Σξ=1截面，各cone的四极点是其r/Σr。Σr在A..F上为1,2,3,4,4,6；三组normalised-column det依次1/48、1/96、1/96，总和1/24，与完整排序截面相同，体积比例1/2,1/4,1/4。这些是各截面tetrahedron的全部极点；不把全部三角化ray误称全谱球极点。

body rates q=ψ(r)在A..F上为8,8,12,12,16,24。单位谱球P={ψ≤1}被每cone的simplex conv(0,r1/q1,...,r4/q4)剖分，各体积absdet/[4!Πqj]。三基础体积1/589824、1/442368、1/442368，乘384得到

|P|=11/4608，∫e^-ψdξ=4!|P|=11/192，p_t(0)=11/(3072pi⁴t⁴)。

这是独立正常数，不迁移三维p_t(0)。

全P真正极点只有40个：B/8 orbit的24个及D/12 orbit的16个。证明：P已由上述simplex并成，任何全局极点必是其某ray/q。对于z=r/ψ(r)，active halfspace normals的linear rank是1+rank{v:v·r=0}：非零项的固定sign gradient与z配对为1，zero项的正负自由改变span恰为这些v，且各v与z正交。故z为顶点当且仅当zero-direction rank3。A..F的rank依次0,3,0,3,2,1，只剩B,D。240个不同triangulation rays不全是顶点；完整几何重建算法及hash在u11_results.json。

## 3. 真实密度的1152锥公式

对任一上述cone R=[r1..r4]，ξ=Σu_jr_j，u_j≥0，ψ=Σq_ju_j。Fourier积分绝对收敛（ψ~|ξ|），换元得

p_t(x)=(2pi)^-4 Σ_cones |detR| Π_j(tq_j−i x·r_j)^-1，t>0。

全部和为实数、非负密度，conjugate cones保证虚部消失。不能把负实项截成0。各cone公式的ψ线性及Jacobian由第1节准确给出。独立spanning Cauchy和证明总质量1/保守，非有限点正值推断。

## 4. 非零节点的真实Green时间积分：符号与pi校准

G=∫_0∞p_tdt。定义b_j=(x·r_j)/q_j，α_R=|detR|/Πq_j。若各cone所有b_j非零，则其时间积分绝对收敛，可逐cone积分；相异b先有

Π(t−ib_j)^-1=Σ_j [i/Π_(k≠j)(b_j−b_k)]/(t−ib_j)。

对ε>0积分至∞，高端log项因系数和0消失，得到

I_ε=−i DD3(log(ε−i·);b1,b2,b3,b4)。

当b_j≠0，log(−ib)=log|b|−i(pi/2)sign b，故Re I_0=−(pi/2)DD3(sign;b)。因此

G(x)=−1/(32pi³) Σ_cones α_R DD3(sign;b1,b2,b3,b4)。

这一系数/符号是推导结论，不是design前提。校准b=(−1,−2,1,2)：DD3 sign=−1/6，Re I=pi/12，直接积分1/[(t²+1)(t²+4)]一致。重复非零节点按confluent DD，sign的高阶导数皆0；例(−1,−1,1,1)给pi/4。

## 5. 零/重复节点：先合并全锥的准确residues

绝不能逐cone给重复0猜sign值：例如b=(0,0,1,−1)，真实Re I_ε=1/ε−pi/2+atan ε发散。线外整体G仍可有限，因为不同cone的零pole必须先合并。

下面给不含数值扰动的全节点算法。每cone内对不同节点b按multiplicity m_b分组。其pole k=1..m_b的系数写成i^(k−4)R_(b,k)。R是真实rational expression，准确计算为

R_(b,k)=α_R [u^(m_b−k)] Π_(a≠b)(b−a+u)^(-m_a)。

注意此式中的[u^j]展开在实变量u；每因子系数为(−1)^ell binom(m_a+ell−1,ell)/(b−a)^(m_a+ell)。总i次方恰为i^(k−4)，因为4−m_b+(m_b−k)=4−k。若无其他节点，只有k4项。把所有1152cones的同(b,k) R准确相加，然后才积分。

对离所有active frame lines的x，026 two-large-jump论证给p_t(x)=O(1) as t↓0（n4）。完整p_t是上述有限rational function，因此合并后所有b=0的principal parts都必须准确为0。程序逐个检查此取消；不把非零项悄悄丢掉。

对于b≠0，k≥2的时间积分实部为0，因为

i^(k−4)(−ib)^(1−k)/(k−1) 是纯虚数。

只剩simple poles。若R_(b,1)现指全锥合并值，则对每个line外有理x（包括任意零/重复节点）有准确式

G(x)=−1/(32pi³) Σ_(b≠0) R_(b,1) sign b。

纯body的rate均整数，故G(x)=C(x)/pi³，C(x)为可准确计算的有理数。u11_green.py的green_allnodes_pi3_coefficient实现此式；先检测active frame lines，在其上返回infinite错误而不伪造有限值。generic版本和allnodes版本是同一真G的不同实现。

退化符号点e1已补登记：全锥zero poles全消，准确C(e1)=3/32；完整confluent时间积分与3/(32pi³)交叉。该单点不作为全cube最小值或成本证书。

## 6. 正时间cutoff的普遍备用式与解析误差

对任意ε>0，包括所有零/重复节点，Re I_ε及第3节全锥和永远有意义。重复节点的DD导数项是

f^(k)(b)/k!=−i^k/[k(ε−ib)^k]，k≥1，f=log(ε−ib)。

ε>0使log位于右半平面，无branch猜测。合并后的cutoff和就是∫_ε∞p_tdt。若d=dist(x,S)>0、S为八body lines，固定δ=d/32；至少two-large-jump及任意pair补basis的absdet≥8给

p_t(x)≤28/(8pi⁴δ⁴)=3670016/(pi⁴d⁴)，

从而0≤G−∫_ε∞p_tdt≤3670016 ε/(pi⁴d⁴)。det下界来自4x4±1 matrix第一行全1，减第一列后3列均为2倍整数；spanning给可补非零basis。该粗误差只支付遗漏短time，不是mp浮点roundoff界。

原点及每body方向线G=∞；其外G有限连续由026已审line分类给出。锥节点为0的hyperplanes大多不是实际奇性，必须通过上述取消或cutoff处理，不能宣布所有node-zero点G=∞。

## 7. body+face正权的可用假设（未跑tau实验）

若额外正权tau≥0乘十二face directions e_i±e_j，排序chamber中face符号为2tau(3a+2b+c)。它在整个排序chamber线性，故仍用同三个cones和1152复制；A..F的rate变为

8+6tau，8+10tau，12+12tau，12+12tau，16+18tau，24+24tau。

只要tau≥0，所有q正，真实半群由body Cauchy和加face Cauchy和独立卷积给出。body basis det16仍给Ct^-4、ψ~|ξ|、H1域、finite-direction cutoff及共同尺度，故U09完整比较/正部源接口的假设仍齐。tau>0时active singular lines还包括所有face方向；allnodes代码检测这些line。零pole取消理由对离全部active lines的点同样成立，常数可依固定tau而变。此处不计算任何tau的G点/成本，不把纯body全谱球40顶点表用于tau>0。

## 8. 有限诊断、实现范围与保留问题

12个预登记有理x_i=a_i+1/p_i，p_i=(151,157,163,167)，准确非退化不是概率：任何ray-zero或normalized-node coincidence会给整数系数d_i（|d_i|≤144）与Σd_i x_i=0。乘四prime积并逐prime取模，强制各p_i整除d_i，故d_i=0；不同normalized rays不允许此结果。所有12点exact generic系数正，与全锥复log cutoff及time/residue公式交叉。前两点50digit完整时间求积，另e1退化路径；实际误差与hash见results，不把mp精度冒充区间。

float generic快接口会拒绝近zero/repeat，拒绝点不能当积分0。exact allnodes可以给有理点严格系数；任意浮点G结果及cutoff_mp结果没有自动roundoff包络。尚无4D全face g下界或J上界，本批无新U4。

固定Lebesgue、共同正r、一般L1/all-scale输入继续由已审026/U09传递；本谱表示不改变背景或算子定义。只有后续真实cube支配和准确成本低于16才能更新上界。本批不做成本pilot，不移植三维常数。

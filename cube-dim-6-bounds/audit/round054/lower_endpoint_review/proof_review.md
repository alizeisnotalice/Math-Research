# F054：31 个低端管的全时间 Green 预算独立证明审查

本证明仅对原 body32 等时间 Cauchy frame 和正face F={(1,y):y∈[0,1]^5}成立，face measure为dy。固定q=6/5、rho=1/64。对1∈A、1≤k=|A|≤5，m=6−k、B=A补集，U_A=span(1_A,1_B)，E_A由dist(x,U_A)<rho和mean_B x∈[0,1/4]定义。

唯一已有解析输入是已审 C037 的完整、正、全时间包络

    G(x) ≤ Σ_(ordered v≠w) M_class(v,w) / (|x|² d_v(x)² d_vw(x)).

三个M的完整原dyadic区间直接读取；其定义积分未重算。992有序项分三类(192,480,320)。本批不会把它的奇异系数当真实G渐近式。

## 1. 几何扰动及所有直线的正距离

写x_A=alpha1_A+e，sum_A e=0；x_B=u1_B+z，sum_B z=0。x1=1给e1=1−alpha。令x0(u)=1_A+u1_B，s²=|e|²+|z|²=dist(x,U_A)²，则

    h=x−x0=(e−e1 1_A,z),  |h|²=s²+k e1²。

在A零和空间中，e1=(e·(e_1−1_A/k))，因此e1²≤(k−1)|e|²/k。于是|h|²≤k|e|²+|z|²≤k s²，得到|h|<sqrt(k)rho=eta_k。k=1时e=0，此论证仍成立。

|x0|²=k+m u²≥k，因此|x|≥sqrt(k)(1−rho)=r0,k。对非负x0，每个body line的投影增益(Σv_i x0_i)²/6不超过all-ones增益，故最近直线距离准确为sqrt(km/6)(1−u)。取u≤1/4并用距离函数1-Lipschitz，所有body lines同时满足

    d_v(x) ≥ d0,k = (3/4)sqrt(km/6)−sqrt(k)rho。

其正性等价于rho²<3m/32；本批五个m都严格成立。零坐标引起nearest-line ties不影响下界。

## 2. 端点1_A的所有incident pair planes

任意两个不同body代表v,w（第一坐标均1），令C为agreement block、D为difference block，则1∈C，D非空，二者平面由不交支持的v1_C、v1_D张成。平面包含1_A时，一般表示lambda v1_C+mu v1_D的每个非零系数会让对应整个支持均非零。1_A的支持是非空真子集，所以恰一个系数非零，那个支持必须等于A。由于1∈A∩C，只有C=A可能。A上的符号必须恒定；B上的符号任意。

因此恰有

    W_sigma=span(1_A,sigma_B), sigma_B∈{±1}^m / {sigma∼−sigma}

这2^(m−1)个incident planes。不同sigma modulo全负给不同B line；没有遗漏混号B。

每个W_sigma内的body代表a1_A+b sigma_B，因为第一坐标为1，必有a=1；B坐标绝对值1又强制b=±1。因此它恰含两个body方向，恰给两个ordered pairs。incident ordered总数c_k=2^m；其两代表在B的m个坐标不同，故全属class j=min(k,m)。此分类不需要枚举496planes/992pairs。

对x≥0，|Σ_B sigma_i x_i|≤Σ_B x_i，所以

    dist(x,W_sigma)²
      = |x|²−(Σ_A x_i)²/k−(Σ_B sigma_i x_i)²/m
      ≥ dist(x,U_A)²=s²。

这是真正逐plane距离下界，未把混号incident planes丢弃。

## 3. 全部nonincident planes的统一gap

若其无符号支持partition不等于A/B，至少一个块同时含一个值1的A坐标和一个值u的B坐标。对同一无符号常数lambda，这两个坐标平方误差之和至少(1−u)²/2；加入其余坐标只增大误差。非负点的signed projection gain不超过同partition的unsigned gain，故任意符号该plane在x0(u)处距离≥(1−u)/sqrt2≥3/(4sqrt2)。

若partition等于A/B而plane不incident，A上的符号必须不恒定，因为全部constant-A与任意B符号已在上一节包含。选两个相反A符号，对共同signed常数lambda，两值1的误差至少(1−lambda)²+(1+lambda)²≥2。因此距离≥sqrt2，更强。k=1时此子情况为空，不产生例外。

以上两类穷尽nonincident signed planes。减去|h|，得

    d_vw(x) ≥ g_k=3/(4sqrt2)−sqrt(k)rho>0。

正性门槛为k rho²<9/32，对k=1,…,5全成立。

## 4. 完整992项的正拆分

设N=(192,480,320)、j=min(k,m)、c_k=2^m。定义

    I_k=c_k M_j/(r0,k² d0,k²),
    B_k=Σ_l (N_l−c_k 1_(l=j))M_l/(r0,k² d0,k² g_k)。

每个background系数非负；五种向量依次为
(160,480,320)、(192,464,320)、(192,480,312)、(192,476,320)、(190,480,320)。各自加incident数后均准确992。

由前三节及C037包络逐正项支配，E_A上有

    G≤Phi≤I_k/s+B_k。

没有通过相减相关区间获得background，也未舍弃任何其他plane。可能已不奇异的incident项用1/s再放大仍合法。由nonincident gap>0，E_A内部所有S2交点只可能来自incident planes；每个plane在face中的交集至多1维、5维face测度为0。s=0时可统一以+∞解释包络，积分论证不受影响。

## 5. 真face Jacobian、maxcone和endpoint clipping

用A自由坐标deficits u_i=1−x_i≥0，i≠1。投影到零和空间e的线性映射Gram为I_(k−1)−11^T/k，det=1/k，因此deficit Lebesgue体积du=sqrt(k)de。由u_i=e1−ei，x_i≤1恰对应maxcone C_A={e1≥ei}；x_i≥0还有u_i≤1限制。

B正交分解x_B=u1_B+z有dx_B=sqrt(m)du dz。故真正face测度为sqrt(km)du de dz，而非直接平面参数测度。正规空间dim=(k−1)+(m−1)=4，距离s=|(e,z)|。

A坐标置换在其零和空间上正交，保持径向测度；k个“指定坐标为最大”的锥除tie超平面外分割全空间，所以径向球积分内C_A占1/k。此对称性在正规空间使用，不是声称物理face可任意交换固定坐标。k=1时零维A残差且fraction1，m=1时B残差零维；两者总正规维数仍4。

保留u∈[0,1/4]、maxcone及s<rho，放大其余box约束。u=0处B正性确实可能强烈clip正规球；忽略clip只会增大非负积分。绝不把以下值称为实际体积等号：

    |E_A| ≤ V_k=(1/4)sqrt(m/k)(pi²/2)rho^4,
    ∫_(E_A) s^(−q) ≤ J_k=(1/4)sqrt(m/k)[2pi²/(4−q)]rho^(4−q)。

后一式来自S3面积2pi²与径向∫r^(3−q)dr；q=6/5<4可积。设f_k=(1/4)sqrt(m/k)，则准确代数式J_k^5=(f_k*5pi²/7)^5 rho^14。

Minkowski给

    ||G||_(Lq(E_A)) ≤ U_k=I_k J_k^(5/6)+B_k V_k^(5/6)。

## 6. 全31域及交叠

保持第一坐标的S5置换是原body32 frame真实对称，故同k的区域具有相同预算，数量n_k=C(5,k−1)=(1,5,10,10,5)，总31。对同一个非负函数G，indicator_union≤sum indicator_EA，因此

    ||G||_(Lq(union E_A)) ≤ (Σ_k n_k U_k^(6/5))^(5/6)。

这一步允许任意交叠，不把signed plane个数当额外物理域，也没有把k与m交换后删去真实域。该结论是all-positive-time G的低端tube联合预算，不是全face。

## 7. 独立准确数值证明：无Arb、无旧M积分

本批使用纯Fraction及整数运算。原M每个完整dyadic区间都保留，所有正表达式取其准确upper endpoint。five-type结构、独立pi夹逼、五个固定k、最后q-union共8配置，实际一次运行PASS。

pi夹逼使用固定Machin公式pi=16atan(1/5)−4atan(1/239)。该恒等式可直接核：若alpha=atan(1/5)，则tan(2alpha)=5/12、tan(4alpha)=120/119，tan(4alpha−atan(1/239))=1；角位于(0,pi/2)，故等于pi/4。用1/5的100项与1/239的30项交错级数（均偶项数），partial sum为lower、加下一正项为upper，严格夹逼误差传播。所得pi宽度<10^-135。

根外包统一用S=10^60。对正有理X与n≥2，先算integer=floor(X S^n)，整数二分得到a=floor(integer^(1/n))；则a/S≤X^(1/n)<(a+1)/S。每次实际核l^n≤X<u^n，保存36个根外包。这里没有用浮点作为证据。

每k取eta上界、line/plane基准平方根下界，形成正d_L、g_L；r0²=k(1−rho)²准确。f_k取上界，得到I_U、B_U、V_U及J_U^5。随后分别以六次根上包I_U^6 J_U^5和B_U^6 V_U^5，构造独立U_k上界。最后对U_k^6作五次根、正权求和，再对该和的五次方作六次根。所有方向单调正确。

五个U_k独立上界均严格低于作者各自规范有理upper。最终准确比较给

    independent_upper < 10293038485441/500000000000000000
                      < (251/1000)*(4/48750)。

display约2.058607697088027e−5，仅作阅读显示。证明性fraction、strict margin、各根power证据及原M完整区间保存在results.json；旧M定义积分、kernel/native、grid、旧方向全枚举、作者writer、MATLAB均0。

## 8. 裁决及边界

PASS A054这31个E_A的全时间Green预算及所报严格有理上界。源文关于体积仅upper、混号incident完整保留、其余992项background、overlap-safe union都正确。作者环境失败发生在默认Python导入SymPy，科学配置为0；本审计没有重跑作者代码，没有失败或恢复。

高端u>3/4、中段长时间尾、管外区域及global g均未付。无全I证书、无新C6界、无新frame继承。主代理保留最终证明及后续预算组合责任。

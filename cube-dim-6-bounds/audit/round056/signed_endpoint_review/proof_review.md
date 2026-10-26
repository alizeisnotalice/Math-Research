# F056：正 B 混号端点积分独立证明审查

范围固定为原 body32 frame、正 face F={(1,y):y∈[0,1]^5}，rho=1/64、q=6/5、L=1/4；E_A、k=|A|、m=6-k、u=mean_B x、r=dist(x,U_A)均继承已完成的F054。C037完整正时间Green包络及F054的几何、完整标量区间作为冻结输入。本轮不重新证明原C037积分证书，不重算任何旧M或完整054预算。

## 1. 正坐标的混号距离

对P_sigma=span(1_A,sigma_B)，在全负等价类中选h个负号、m-h个正号的代表，1≤h≤m/2。设b=x_B≥0，P、N分别为两组坐标之和。由非负性，|P-N|≤max(P,N)。每组Cauchy–Schwarz分别给P²≤(m-h)||b||²、N²≤h||b||²≤(m-h)||b||²。因此

    dist(b,span sigma)²
      = ||b||²-(P-N)²/m
      ≥ (h/m)||b||² ≥ h u²。

全plane的A残差再加一个非负平方项，所以dist(x,P_sigma)≥sqrt(h)u。另一方面，非负b给|sigma·b|≤sum b；A投影相同，故dist(x,P_sigma)≥dist(x,U_A)=r。合起来得到

    dist(x,P_sigma)≥max(r,sqrt(h)u)。

零坐标和两组和相等均被包含。r=u=0处可以将非负包络定义为+∞，不改变积分。

## 2. 先逐点替换，再放大区域

F054的真实face Jacobian为sqrt(km)，正规坐标(e,z)维数4。A上“第一坐标最大”锥对径向测度占1/k，包括k=1零维A残差情形。故径向积分系数为sqrt(m/k)2pi²。

上一节的不等式只在实际E_A的b≥0上使用。先把原distance^(-q)替换为max(r,sqrt(h)u)^(-q)，再把其它box/positivity clipping去掉。替换后的函数在放大的正规球上仍是定义良好的非负径向函数；集合包含关系允许积分上放大，无需在球外继续断言原几何不等式，也无需假定实际clip有置换对称性。因此

    ∫_(E_A) dist(x,P_sigma)^(-q)
      ≤ J_h,k
      = sqrt(m/k)2pi² ∫_0^rho r³ ∫_0^L max(r,cu)^(-q) du dr，
    c=sqrt(h)。

这些是放大域的上积分，不是实际clip域积分等号。m=1无混号plane，此时整个k=5输入原样保留。

## 3. mean与radial积分、准确比值

因rho<cL，0<r≤rho时分点r/c位于(0,L)。直接在[0,r/c]和[r/c,L]积分，得到

    ∫_0^L max(r,cu)^(-q)du
      = q r^(1-q)/[(q-1)c] - L^(1-q)/[(q-1)c^q]。

再乘r³积分，得到

    J_h,k/[sqrt(m/k)2pi²]
      = q rho^(5-q)/[(q-1)(5-q)c]
        - L^(1-q)rho⁴/[4(q-1)c^q]。

旧unsigned放大域公式为J_0,k=sqrt(m/k)2pi²L rho^(4-q)/(4-q)。设alpha=rho/(cL)，准确相除得

    J_h,k/J_0,k
      = q(4-q)/[(q-1)(5-q)] alpha
        - (4-q)/[4(q-1)] alpha^q
      = (84/19)alpha-(7/2)alpha^(6/5)=R_h。

系数和微分正性门槛在独立配置1用Fraction核对。对0<alpha≤1，alpha^(6/5)≤alpha，所以R≥(35/38)alpha>0；R'≥84/19-21/5=21/95>0，且R(1)=35/38<1。h=1,2分别alpha=1/16、1/(16sqrt2)，故0<R_h<1。这里的比值是两个放大域上公式的比值，不是实际G或clip积分的比值。

## 4. 全部符号类别与有序项

全负不改变plane。h<m/2时，每个等价类恰有一个少数负号代表，因此N_mh=C(m,h)。m偶且h=m/2时，两互补符号向量不同，成对归为一个plane，故N=C(m,m/2)/2。另有一个unsigned plane。m=5,4,3,2,1的混号数依次为(5,10)、(4,3)、(3)、(1)、空；每行满足1+sum N=2^(m-1)。

F054已证明每个incident plane恰含两个body代表，给两个同class min(k,m)的有序pair。本轮每行再核2(1+sum N)=2^m，故混号与unsigned合计保留全部incident有序项。这里没有重新枚举方向、pair或plane。

## 5. Minkowski与冻结054完整输入

设C_k=M_j/(r0,k²d0,k²)，同一k的全部incident有序项使用相同C_k。原singular norm公式为S_k=2^m C_k J_0,k^(5/6)。对一unsigned plane与全部mixed planes分别Minkowski，得到

    new singular norm
      ≤ 2 C_k[J_0,k^(5/6)+sum_h N_mh J_h,k^(5/6)]
      = F_m S_k，
    F_m=[1+sum_h N_mh R_h^(5/6)]/2^(m-1)。

因此U_new,k=F_m S_k+B_k。无需plane间或其邻域不相交。m≥2时0<F_m<1，理论上严格改进预算公式；B_k保留原054所有nonincident正项。实际计算读取S_k、B_k的完整dyadic区间，两端均保留，未把输入当点值。

作者056中五类background完整JSON对象与054逐字一致；全部source_054_intervals的singular、background、total三对象也逐字一致。k=5的三个完整对象原样复制，且union中已保存的k5加权q幂对象也原样复制；独立批没有重新计算k5的norm、幂或根。旧054完整union区间直接读取而不重算。

## 6. 同一31域的正union

第一坐标固定，S5对称给权重(1,5,10,10,5)，总31。对共同G有1_union≤sum1_EA，故

    ||G||_(Lq(union E_A))
      ≤ [sum_k n_k U_new,k^(6/5)]^(5/6)。

重叠允许，不会遗漏k与m不同的真实区域。整条链从已积分全部t>0的C037包络出发，因此仍是低端管的全正时间预算。

## 7. 独立准确算术证书

先写theory_before_results.md，再打开作者056结果。独立代码只用Python标准库Fraction与整数幂、二分根；没有导入或运行作者程序，也没有使用Arb或浮点作证据。共10个预登记配置，一次执行。

根界尺度S=10^70：对正有理x，取整数floor(xS^n)，找其整数n次根a，再给[a/S,(a+1)/S]。每次实际核l^n≤x<u^n；19项证书保存n、准确input、lower、upper，足以直接复核幂不等式而不重做根搜索。

对alpha先向外夹逼，再用五次根外包alpha^6。R下界为A alpha_L-B(alpha_U^(6/5))_U，上界为A alpha_U-B(alpha_L^(6/5))_L，避免差式方向错误。正R的5/6次幂通过六次根外包R^5。F_m正线性求和，随后与完整S_k区间作正乘积，再加完整B_k。前四类q幂用五次根外包U_k^6，k5直接加入已保存的加权q幂区间；最后对总和五次方取六次根。

所得全部对应author alpha、R、R^(5/6)、F_m、new singular、new total及union区间，都被独立有理外包包含；前四类及union的独立upper又严格低于作者规范有理upper。最终准确比较为

    independent_union_upper
      < 13404268290817/1000000000000000000
      < old054_full_formula_interval_lower。

独立upper展示约1.3404268290815612e-5，新/旧公式比值upper展示约0.6511327199337893。完整有理端点和严格margin在results.json中；小数仅帮助阅读。这里的old lower是旧上预算公式的区间下端，不是实际G的下界。

## 8. 来源、执行与裁决

作者056 manifest固定为9305a1c06d472f5a712f4c31f580be5b3e653dd7edb01b751cef0e3413f48e71。21个来源pin涵盖其全部冻结文件、所用054原件及最终F054 manifest/report，执行前后与封存时一致。gate_release绑定最终F054 PASS和父代理release，原等待批0科学配置，continuation只做固定9 producer+1 saved check；本审计没有为其补跑任何工作。

独立批PID60298、tool chunk53815e、exit0，计算0.009299333秒；登记后139.142091秒结束，低于900秒。准备登记时曾把一个相对pin路径接到clone根而非project根，FileNotFoundError在0科学配置且登记前发生；只修正路径，未发生科学失败或重试。科学结果、receipt与代码哈希链完整。

PASS，限原31个低端E_A的全正时间Green预算改进。高端、中段长时间尾、管外、global g以及全域I仍未支付，无新C6界、无新frame继承。父代理保留最终证明及预算整合责任。

# F046 首父锥局部 fan 的独立证明审查

## 1. 真实 frame、共同时间与归一化

B044 列由32个 v=(1,±1,…,±1) 与120个 beta w 组成，其中 w 支撑大小4、首个非零坐标为+1，beta=1/80。每个四元支撑有8个 unoriented 符号向量。首符号规范只是对称 Cauchy 方向的代表选择，不能在 symbol 中补乘2。

固定坐标出现在10个支撑、每支撑8个方向中。不同坐标的乘积经其中一个坐标翻号配对抵消，即使需要重新选择全局符号，乘积仍不变。因此 Σww^T=80I，Σvv^T=32I。真实列 Gram 为(32+80beta²)I=(2561/80)I；边际 Cauchy rate 是Σ|a_i|=32+80beta=33。两式分别依赖 beta² 与 beta，不可混同；Gram 不是不存在的 Cauchy 协方差。

每列乘独立标准 Cauchy(s) 参数、全部同一 s，正概率推前定义 X_s。特征函数为 exp(−s psi)，psi=Σ|a·xi|。body32已张成R6，故psi coercive，Fourier逆变换有连续光滑密度。Cauchy加法给正保守卷积半群，积分为1；psi一次齐次给p_s(x)=s^−6p_1(x/s)。这些性质属于新152列frame，不继承旧body32点值、奇异管估计或数值结果。

独立 exact sum 验证152列互异、32个body代表、15个四元支撑各8符号、Gram及六个边际率；没有重跑B044旧2160射线检查。

## 2. 原父锥和必须细分的原因

只读旧fan第一父锥的六个整数列R，准确detR=1。正系数u生成其锥；每条非零射线唯一遇到Σu=1的标准5-simplex P。首父锥原本位于有序非负 Fourier chamber，故其子锥也在该chamber。

对120个新forms只做这个父锥的720个整数dot。准确crossing indices是2、10、26、50、82，对应

    u1−u3，u1−u3−u4，u1−u4，u1+u2−u4，u2−u4。

其它115个新forms在父锥上有固定弱符号。取系数(2,1,1,1,1,1)与(1,1,2,1,1,1)，得到两个严格chamber及父锥内部点；首cut取值+1与−1。该项的绝对值在中点产生至少beta=1/80的严格凸性差，其它正绝对值项不可能抵消。于是新psi在父锥上不线性，旧q-only更新不是合法Fourier积分恒等式。后续只核这一父锥，不接受其余20个父锥或全fan。

## 3. 顶点完整性，不是候选顶点抽查

每个闭符号区域由Σu=1和11个线性不等式定义：六个坐标非负和五个带符号cut。所有32个区域的并是P；两个不同full区域的交位于至少一个cut等式，因此full区域相对内部不交。cut为零时两种符号均保留。

任何有界多面体在5维仿射空间中的顶点，必有5个线性独立的活动约束（包括低维区域的隐含仿射等式）。若活动法向不能张成切空间的对偶，存在非零可行双向扰动，便不是顶点。反之独立的5个活动约束与Σu=1确定一个可行点，则该点是顶点。因此穷尽C(11,5)=462个系统足够且必要地覆盖全部区域顶点。

本独立实现先消去u6=1−Σ_(i≤5)ui，将5个活动等式改写成5×5系统；用标准库Fraction的独立Gauss–Jordan消元，而非作者6×6 FLINT determinant/solve。只做一次完整462枚举。237个奇异、6个唯一解在P外，其余去重为12个点。每个点、每个系统状态完整保存。对每个符号区域再用全部五个准确不等式过滤这12点，得到准确完整顶点集：12个维数5区域、20个非空低维区域，无空区域。每个低维区域的完整顶点集都包含在一个记录的full区域闭包内，故凸包也包含；所有tie边界得以覆盖。

## 4. Facets、pulling覆盖与共享面

对由原11个半空间定义的任意记录face，原边界等式与该face的适当交集都是face；其中按包含关系极大的proper交集恰为所有facets。因为任一facet必有至少一个非恒等的活动原不等式，故不会遗漏nonsimplicial或多重支撑facet。

独立reader不用作者的维数筛选来寻找facet，而先取全部非空proper zero vertex sets，再按集合包含关系保留极大项，随后验证其准确维数为d−1。这些集合与66个保存节点的facets逐一相等。

统一全局字典序最小顶点p作为apex。任意face内部点沿p出发的射线离开face时，落在一个不包含p的facet上；除低维交线外退出facet唯一。因此face被conv(p,E)覆盖，各pyramid内部不交。按维数归纳，递归triangulate E并cone到p得到完整无内部重叠的triangulation，闭包保留所有边界。独立代码按维数自底向上重建全部66个递归关系及12个根，得到恰14个full simplices。完整性由此归纳和半空间partition证明，不由体积和证明。

不同符号cell的交是双方公共face，其顶点就是完整顶点集的交。统一pulling顺序限制到该face具有相同triangulation。除此证明外，独立实算全部66对full cells的共享face，将两边所有诱导simplex交集按包含关系取极大项（不只比较预设维数的项），确认其完全相同且覆盖同维公共face。低维tie没有删除。

## 5. Primitive rays、真实symbol与物理权重

每个u顶点先精确乘R，清分母再除坐标gcd，唯一得到x=s r，s>0有理、r为primitive integer。12个r、s与保存对象一致，均留在有序chamber。各rate从全部152个真实列重新计算q(r)=Σ|a·r|，绝不读取旧q。

每个full cell的完整顶点集对每个form都有同一个弱符号，因此该form在整个cell及其锥上保持符号；新psi成为准确线性函数。将这些符号乘真实列求和，独立得到保存的完整symbol slope，并在每个顶点核对slope·r=q(r)。所有14个子锥另核152forms符号。

对六顶点子simplex，signed det(R[u1…u6])=detR·det[u1…u6]。由x_i=s_i r_i得

    det[r1…r6] · ∏s_i = detR · det[u1…u6]。

独立代码核的是带符号恒等式，不仅绝对值；各primitive determinant为±1。Fourier变量替换的测度权重使用其绝对值，不把orientation负号误当负贡献。

frame symbol为偶函数。整体负射线−r有相同rate；R6中det(−r1,…,−r6)=det(r1,…,r6)。若未来全fan将相反锥合并，两个Fourier积分为共轭，只有在正确完整覆盖后才可取2Re。这个局部正chamber对象不含额外全局负号因子2，也没有验证全部signed-permutation multiplicities。

## 6. 物理截面与参数截面明确区分

参数P的Euclidean5-volume是√6/120。物理截面R(P)的方程为ell·xi=1，ell=R^−T1=(0,−1,1,1,0,0)，||ell||²=3。线性R在该超平面上的面积因子为|detR|·||R^−T1||/||1||=1/√2。因此物理5-volume为√3/120。

独立计算还有不依赖该面积公式的验证：取五条物理边r_i−r6，它们的5×5 Gram determinant准确为3，故simplex体积直接是√3/5!。子simplex与父simplex的相对体积在两种截面相同，准确为|det[u1…u6]|。14项总和1只是覆盖证明之外的冗余校验。

这是ell·xi=1的仿射截面，不是psi(xi)=1截面。A046原记录的体积标签曾修正；本审计直接核最终数据和Jacobian，不依靠历史hash声明来推断体积。原始与最终hash分开保存，不回写作者冻结源。

## 7. 验收边界

本对象是beta=1/80新frame下首父cone的完整准确局部fan。其余20父cone、跨父cone整合、global符号/置换权重、zero/repeated Fourier node算法和任何热核/Green值均未认证。无旧tube继承，无全I/g/Cn新界。

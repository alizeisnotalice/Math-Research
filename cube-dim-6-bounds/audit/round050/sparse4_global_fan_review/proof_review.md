# F050 sparse4 全fan的独立覆盖证明审查

## 1. 输入与证明范围

对象为B044的准确152列frame：body32加beta=1/80的120条unoriented sparse4方向。其同时间正Cauchy半群、Gram和边际归一在F046已审核。旧21个parent cones完整覆盖有序正chamber的证据与原fan输入保持锁定。首parent的A046/F046完整局部证书仅继承，不重新求解462系统。

本轮独立核20个新parent的全部保存active ledger和全局拼接。接受的是全fan数学入口，尚不接受任何数值引擎、zero/repeated Fourier node处理、新S2几何、tube预算或Cn上界。

## 2. 保存边界与实际frame绑定

对每个新parent，从原fan读取准确R并核非零det及与保存parent相同。对120个新form计算wR；有正有负者恰为必须切开的超平面，其余form在整个parent上有固定弱符号。除正gcd后重建全部cuts，与保存H逐项比较。六个坐标边界也准确一致。

不同sparse4方向没有重复unoriented line，R可逆，所以不同cuts不定义同一超平面。程序还直接核normalized cuts没有相同或相反项。由此后续ledger所证明的多面体确实属于真实新frame，不能靠自行保存一套错误H来通过内部检查。

## 3. 每个顶点的完整性由全部active ledger证明

在P={u≥0,Σu=1}的仿射5空间内，顶点必须有5个独立活动边界。若活动法向未张成切空间的对偶，就有非零双向可行小扰动，与极端性矛盾。反之5个独立边界加Σu=1的可行唯一交点是顶点。低维sign cell也适用：其仿射等式作为两侧活动约束参与该秩论证。

因此穷尽C(6+m,5)个活动组合，足以列出所有sign cell的全部顶点，不需要先列出2^m个sign vectors。任何非负活动系统交点至少可按每个cut的实际符号放入一个closed cell；它在每个可行cell中仍由同样5个独立边界支撑，故不是误加入的非极端点。

独立reader只核保存记录，不求解系统。与作者的6×6系统不同，本轮消去u6=1−Σ_(i≤5)u_i，核5×5整数矩阵(H_(k,j)−H_(k,6))的det。原6×6 determinant与它相差一个固定负号，所以非奇异性严格等价。

对全部783994条lexicographic combinations逐条核：

- S记录：det准确为0。
- O记录：det非0，完整保存有理解满足Σu=1及五个活动等式，至少一坐标为负。
- I记录：det非0，所指保存顶点非负、和1并满足五个等式。

同时核完整组合序列耗尽、每个vertex的first-appearance enumeration ID和活动见证、所有ID都出现、分类数及未压缩digest。非奇异性保证保存解是唯一解，故不需要重求其值。最终247176个S、319333个O、217485个I，计783994；没有遗漏系统、候选搜索或新vertex。

## 4. Closed cells的过滤与BFS为何完整

固定sign code，将完整pool按所有带符号cuts过滤，所得恰是该bounded cell的完整顶点集。独立reader重做这一准确过滤并核rank=6，即仿射维数5；它不创建新的sign code或cell候选。

对一个full cell，cut j的zero顶点集仿射rank4当且仅当它支持一个facet。它不可能是父simplex的边界facet：否则两个仿射4维超平面重合，齐次形式必须成比例，而能cross parent的cut不可能与ui=0成比例。两个不同cuts也不能共用这一facet。因而facet相对内部存在所有父坐标和其它cut均严格的点；穿过它只翻转sign j，得到唯一邻接full cell。

独立检查对每个保存full cell的每个cut都计算该zero集合的准确rank。每个rank4 facet都必须有已保存的单bit邻居，其完整顶点集交必须恰为该facet；邻居图与记录逐项一致。记录的generic seed由正整数(1,z,…,z^5)组成，z≤5m+1，避开全部cut零点，并对应保存seed。每个cut是次数≤5的非零多项式，所以该有限选择有理论保证。

所有真实full cells的facet邻接图连通：任取两个严格内部点，可在P内部选一个generic折线中间点，使两段避开所有cut交的codimension≥2集合及退化平行情形。由于P凸，路径全在内部；它只会一次穿过一个超平面，通过其facet相对内部。因此存在一列合法单cut邻接连接两cell。

独立检查的保存图非空、连通，且对每个真实facet强制闭合。结合上述真实全图连通性，它不能只是全图的真子集。于是保存full cells确实穷尽。这个推论并不依赖BFS输出旗标或volume sum。

剩余低维sign cells也未丢弃：P中避开全部cut的点稠密，已穷尽的有限个closed full cells之并是闭集，故其闭包覆盖P的全部边界、零和tie。这避免了不必要的2^m枚举而保留完整覆盖。

## 5. Pulling覆盖不靠体积推断

完整cell的任何facet、及递归face的facet，必由至少一个原半空间等式切出。独立reader取所有非空proper zero subsets，按包含关系选极大项，再核维数d−1，复原准确facet集合。对20888个新face DAG节点，自底向上核least-vertex apex、全部不含apex的facet与cone递归，所有simplex均有正确维数且无重复。

数学归纳是：从apex向face内部点引射线，退出点落在不含apex的facet；除低维tie外退出facet唯一。因此这些pyramids覆盖face、内部不交。每个base的归纳triangulation通过cone到apex给完整无内部重叠triangulation；闭包保留边界。每个full cell的saved child列表与该递归结果完全一致。

这一facet/pyramid归纳配上上一节全cell覆盖，才证明每个parent的完整性。相对volume和1只是另一个准确冗余校验，不能弥补不完整ledger或未闭合邻接图。

## 6. 真rays、rates、sign和determinants

全局478条primitive rays逐一核坐标gcd、ordered chamber归属及真实152列的全部dot，q(r)=Σ|a·r|准确为保存有理数。每个local vertex核Ru=s r，s>0准确有理。统一排序为旧A046的12ray优先、其余primitive整数lexicographic；它限制到首parent恰为已审顺序。

每个新full cell对全部152forms独立核固定弱符号，按其真实符号重求完整symbol slope并在所有cell rays上核slope·r=q(r)。每个child的所有顶点又逐一核属于该完整cell。因此全部child cones继承全部152forms同号及真实symbol线性；没有必要再做13215×152个完全冗余符号循环。该继承逻辑明确计入操作账，不把未逐child重复实算说成实算。

所有新child都独立核primitive signed determinant及u-slice determinant，验证带符号恒等式

    det R_child ·∏s_i = det R_parent ·det U_child。

每个child非退化、体积字段准确，每parent相对体积准确和1。物理Fourier测度取|det R_child|，不把orientation符号当作积分负权。

## 7. 跨parent和跨chamber的边界

统一primitive-ray顺序使pulling在共同face上取相同组合递归。不同parent的Σu=1截面之间只是同一射线的正比例变化，保持face/vertex incidence和全局顺序；所以不会因截面不同而产生不同pulling。

此外独立重建全部126个旧parent facets：28个内部facet各有两个owner，它们的诱导4维simplex集合准确相同；70个仅有一个owner的facet各有准确ordered-chamber墙见证。计数2·28+70=126，全部边界账完整。低维公共face由pulling限制一致性和这些facet匹配保留。首parent在此只读取原certificate作为拼接数据，不重新active枚举或生成新局部fan。

新frame的两个方向orbit都在signed permutations下不变（unoriented代表的全局负号不改变Cauchy项），故psi具hyperoctahedral对称性。严格ordered正chamber内部有唯一群代表，群大小2^6·6!=46080。相邻chambers的公共墙由对应simple reflection逐点固定：相等坐标交换或最后零坐标翻号。因此墙上的诱导triangulation也一致；所有群边界均低维，积分不重复。

全群fan没有额外乘global-negative factor2。若未来使用相反锥配对，只能把两项相加重写成2Re，不能在已经列全46080群元素后再次加倍。

## 8. 全局切片及正时间Fourier入口

全局汇总由21个原parent的已核children准确组装，逐记录匹配global_fan中的parent/local ID、ray IDs、det和q，共13229cones、478rays、4266full cells；无未覆盖parent。

所有primitive rays的第一坐标正，因此可取全局xi1=1截面。每个child截面体积为|det R_child|/[120∏r_(j,1)]。独立准确求和得1/120，与ordered chamber的标准5-simplex相符。这不是各parent的Σu=1截面，也不是psi=1截面；它仍仅为覆盖证明之外的冗余校验。

对t>0，psi coercive，exp(−t psi)在R6绝对可积。每个锥上线性率q_j>0，故其Fourier积分准确为|detR|∏(tq_j−i r_j·x)^−1。完整有限cone/group求和再乘(2π)^−6给真实p_t。未代入任何t或x；没有计算热核值。

未来直接全群逻辑项数为13229·46080=609592320。该数不是运行时间或新Cn证据。rational rate/int64范围、zero/repeated Fourier nodes、数值精度、time0端点和Green积分交换均须另审。尤其正时间恒等式不能自动许可逐cone积分到0，也不让旧21cone引擎自动兼容。

## 9. 实际独立范围

唯一完整批核783994条ledger，无新solve、vertex、sign-cell搜索、首parent枚举。所有matrix/结构操作1695425，按登记低于2000000。每ledger的一个准确det与已有witness闭合计为一条验证；每cache-miss rank、其它det、ray/form dot、cell/form sign以及列明结构对象另计。普通整数乘加不是额外matrix操作。完整分类及逐parent收据保存在results.json。

接受完整数学fan；未审数值引擎、不继承旧tube/新S2、不证明全I/g或Cn改进。

# 非中心强极大算子端点常数

全部资料入口见 [CONTENTS.md](CONTENTS.md)，原始附件及历史工作档案均在本项目目录内。
研究对象为 Lebesgue 测度下轴平行非中心强矩形极大算子，统一使用
\(\Phi_n(t)=t(1+\log^+t)^{n-1}\)、自然对数及阈值 1。
本项目研究最优比值 \(C_n=\sup_f |\{M_Sf>1\}|/\int\Phi_n(f)\)。

## 当前状态（2026-11-08，第 79 轮）

- **第79轮**：固定step源的halo可用源网格顶点表示，但不能据此统一给所有G≤F定价；薄层模型的连续角点单列最优已精确求解，原输入网格可造成任意固定倍log N的对数价格损失。另证明cap2输入的有理step稠密性及固定容量下有理角点覆盖的无额外预算近似。40组精确几何、66组定价区间均通过；一般锐系数未变。见 [本轮报告](rounds/r079/report.md)。
- **第78轮**：Orlicz坐标封顶合并证明同角点真实分数覆盖可无损合成一源，适用可数列和全接收集合。618个完整源列及12个有理覆盖LP独立认证；三角固定角点族有严格分数优势，增加自然角点后差距精确消失。跨角点全局预算仍缺，一般锐系数未变。见 [本轮报告](rounds/r078/report.md)。
- **第77轮**：固定有限角点/接收小格的任意可测源分配精确归约到凸规划；32实例205组有理证书通过，12组由Mathematica精确LP复算。两个不同角点给裸分数松弛的维数损失，perspective恰修复单小格例。仍缺整个halo的维数一致源预算，一般锐系数未变。见 [本轮报告](rounds/r077/report.md)。
- **第76轮**：精确优化cap2近角点阈值成本，log B_N(τ)=N+√(2Nτ)+O(τ+1)（τ=o(N)）；临界向量动作在单格共享角点。完整高效见证迫使指数级角点熵或源复用；分离副本说明角点数可任意大、归一化复用却不随副本数增加。下一步控制带成功条件的源预算分配，一般锐上界未闭合。见 [本轮报告](rounds/r076/report.md)。
- **第75轮**：固定负角点的局部合格坏换元策略可用张量Hardy全局结算为exp(N+O(1))成本，指数低于参考尺度。重锚恒等式本身仍有源依赖、径向坍缩和向量指数失真；真实最大化折叠例可由成本筛选剔除。同一原点图还闭合了指数阈值的任意向量策略高成本筛选；多项式阈值和连续角点预算仍未闭合。见 [本轮报告](rounds/r075/report.md)。
- **第74轮**：同一个cap2紧支源可在实际动作严格成功且局部成本J_ref N^o，但固定同点参考平均及捕获成本被嵌套矩形体积比压至超多项式小。排除单点低成本直接转参考的桥梁，仍需动作/锚点的全局预算分配；条件后验尾引理亦已证明。见 [本轮报告](rounds/r074/report.md)。
- **第73轮**：一般输入与上限2输入的最优比值满足C_N/T_(N,2)=1+O(N exp(−N/100))，一般振幅归约已闭合。完整标签紧支源反例排除无条件逐点精细参考转移；局部效率候选仍待证明。一般锐上界仍需动作/锚点预算分配。见 [本轮报告](rounds/r073/report.md)。
- **第72轮**：宽径向窗口把有界输入的成本筛选降为N^(1/2+o(1))，接通禁带捕获，不再需要全局部成本弱型定理；有限c完整标签核的精确正卷积半群势表示被排除。一般上界仍缺变动作参考比较、锚点与无界余集。见 [本轮报告](rounds/r072/report.md)。
- **第71轮**：证明固定维数紧支源弱尾恒等式、临界单列极限1/2及一阶修正；给出全接收小格有效的连续下界公式。12核Python与12内核Mathematica核验通过；全输入弱型与一般锐上界仍未闭合。见 [本轮报告](rounds/r071/report.md)。
- **第70轮**：局部成本算子的物理Lᵖ范数精确为p/(p−1)，单列弱准范数为1；单次Lp筛选的最优损失仍为N量级。预算桥只需实际弱型指数γ≤(9−μ)/12≈0.64346，未证明该弱型界；限时网格未找到增长证书。见 [本轮报告](rounds/r070/report.md)。
- **第69轮**：固定参考动作下，禁用宽2K log N的临界带使最优成本增加N^(2K+o(1))；由封顶修复得到局部效率限定的窗口预算捕获。物理Jensen接口合法，仍缺变动作共用源成本与接收主体分配。见 [本轮报告](rounds/r069/report.md)。
- **第68轮**：固定动作条件总和窗口的锐L¹包络精确为1+log(H(d)/H(a))；先固定源层得到允许任意接收深度选择的缩放协变成本界。另证明固定块的混合径向—动作导数在高频无界，排除直接混合平滑路线；一般上界仍未闭合。见 [本轮报告](rounds/r068/report.md)。
- **第67轮**：条件角核保留联合符号后，一阶动作响应维数无关、径向响应为O(1/b)，并给精确有限径向L²接口。缩放与廉价填充反例排除直接全局角平均的预算捷径；一般上界仍缺兼容局部基准和接收锥分配。见 [本轮报告](rounds/r067/report.md)。
- **第66轮**：临界实际预算把任意成功块动作定位回近μ窗口；宽角网格使固定块数、有界径向模型的完整动作 halo≤(1+o(1))Q_NJ。另证明未条件化角向一阶 L² 响应维数无关，并给二阶响应反例。一般角向归约仍未闭合。见 [本轮报告](rounds/r066/report.md)。
- **第65轮**：固定块数、有界块径向模型的整个变参数 halo≤(1+o(1))Q_NJ。统一参考动作、微小空间修复及随机角网格预算分配排除该模型的主阶突破；一般真实角向输入仍未覆盖。见 [本轮报告](rounds/r065/report.md)。
- **第64轮**：角向实际模数的Jacobian精确抵消；任意角位置相关的共同块参数得到全径向预算包络。临界条件密度包络锐为 √(b log b)，明确了不加效率条件的半对数损失。全局径向积分仍有锥重数；一般 C_N 上界未闭合。见 [本轮报告](rounds/r064/report.md)。
- **第63轮**：实际预算的概率包络移除了终端单调性；固定块数的位置响应精确到 +1/N，并给共同曲率的常数项。真实角核及共享方差假设的范围已单独核查。12内核Wolfram与12线程MATLAB交叉核验完成；一般 C_N 系数区间未变。见 [本轮报告](rounds/r063/report.md)。
- **第62轮**：固定近临界单调终端的共同窗口负曲率，并扩展至至多 C log N 个平衡块的向量 Hessian；条件共享损失已推导，条件径向化存在动作依赖反例。百万维 MATLAB 与80/160位 Wolfram 核验完成。一般上界的真实矩形表示仍缺失，系数区间未变。见 [本轮报告](rounds/r062/report.md)。
- **第61轮**：临界实际模数预算下，任意非负终端输入对所有可预测反馈策略均满足 L1 振幅稳定性 O(√(logN/N))；固定连续抛物线的指定锥见证评分认证低于71.666，停止盲目加密该曲线。见 [本轮报告](rounds/r061/report.md)。一般 C_N 系数与认证下界未变。
- **第60轮**：非均匀外宽度关闭几何网格后缀逃逸，使真实点值共享缺口达到 (log N)²/N、实际 Orlicz 成本仍趋于1，秩可为 o(N/log⁴N)；但本族端点评分保持有界。33点源重优化未超过既有认证下界。见 [本轮报告](rounds/r060/report.md)。
- **第59轮**：正则冻结源在固定自然角窗内的共同动作优化仅改变径向边界 O((log N)^3/N)；几何网格给低成本条件高秩反例，但全部外区间的后缀逃逸及阈值归一化仍阻碍该族。见 [本轮报告](rounds/r059/report.md)。一般常数与认证下界未改善。
- 已复核：指数率 \(C_n^{1/n}\to B\)，以及 \(\limsup C_n/B^n\le\sqrt B\) 的证明路线；B 满足 \(\log B=1+B^{-1}\)。详细推导见 rounds/r001/report.md。
- 已实现：等参数单层见证的 Gamma–二项混合、数值优化和 Wolfram 高精度交叉复算。
- 待补证明：附件所称整个双曲示性函数族 H_n 的全局三项展开，特别是全空间 halo 与有限 L 的统一比较。
- 已完成：有限多层轮廓的严格极限归约与固定网格凸振幅优化；n=4…256 的候选均优于同维等参数单层数值对照。
- **第3轮严格下界**：n≥1024 时，log C_n ≥ n mu − sigma sqrt(n(log n−2loglog n)) − log(2n)。新增正的 sigma loglog n sqrt(n/log n) 项；首个负项系数仍为 sigma。完整证明见 rounds/r003/report.md。
- 第4轮：n=4…256 的七个多层候选已通过定向算术认证；完整端点和输入hash见 rounds/r004/certificates.json。
- **第5轮机制障碍**：任意非减轮廓、任意坐标比例向量的单一双曲接收见证，其最优 Q_n 满足 (n mu−log Q_n)/sqrt(n log n)→sigma。这是受限机制定理，不能用作 C_n 的一般上界。
- 第6轮：非负加权深层单调轮廓逐项被 Q_n 支配；分块张量示性函数的产品见证逐项被普通单层示性见证支配。均有明确固定维数深层极限限制。
- **第7轮几何障碍**：共同轴向多中心/对角缩放双曲并集满足 V(Kt)≤K V(t)。统一振幅下，仅合并各分量原有双曲见证不能超过单分量深层基准；不控制新增跨分量矩形。
- **第8轮完整机制展开**：log Q_n = n mu − sigma sqrt(n log n) + sigma[loglog n+log(sigma sqrt(2pi))]sqrt(n/log n) + o(sqrt(n/log n))，并作为一般 C_n 的严格下界。
- **第9轮重排障碍**：允许非单调标量轮廓、任意接收层及逐点选择比例后，锚定族的模块弱常数 W_n 仍精确等于 Q_n。此结论不覆盖非锚定真实矩形。
- 第10轮：闭合一维非负缺口及端点四项账本；证明仅振幅上界不能强迫纯Lp稳定性损失，必须保留Markov和模数效率。另给真实强矩形/坐标复合、真实区间/锚定区间的精确严格反例。
- 第11轮：一般输入严格归约到可预测反馈的一维上界模型 mathcalR_n，满足 Q_n≤C_n≤mathcalR_n。Bellman基准已计算至256步，反馈严格优于任意固定两步比例。
- **第12轮一般定理**：log C_n = n mu − sigma sqrt(n log n) + O(sqrt n)。证明一般首个缺口系数为sigma，并推出 C_n/B^n→0。完整一般上界见 rounds/r012/report.md。
- **第13轮一般三项定理**：第8轮公式现对一般C_n及反馈上界mathcalR_n成立；证明截停累计均值损失，再逐项核验FGL条件鞅相对中偏差。详见 rounds/r013/report.md。
- **第14轮隐式根定理**：令barPhi(z_n)=sigma z_n/sqrt n、d_n=sigma sqrt n z_n，则nmu−d_n−O(log n)≤log C_n≤nmu−d_n+O(n^(3/8)(log n)^(1/8))；因而确定任意固定阶逆对数展开。见 rounds/r014/report.md。
- **第15轮精确指数控制**：负Laplace目标的可预测最优值为exp(−n lambda c_lambda)，由唯一固定动作达到，并有逐步非负损失恒等式。此障碍仅限指数目标，不覆盖尾概率或端点；见 rounds/r015/report.md。
- **第16轮真实分位数定理**：固定c=mu−(sigma/2)sqrt(log n/n)相对c=mu的临界下分位数改善为(mu−1)log n/8+o(log n)。证明使用精确CDF参数导数与自包含Fourier密度估计；仅限该测试族，不是C_n的新对数级项。见 rounds/r016/report.md。
- **第17轮实际下界族改善**：在固定斜率1/n的指数软轮廓中，指定比例偏移使真实见证评分的对数提升(mu−1)log n/8+o(log n)；已结算严格接收与I_n分母。此为家族内改善，未缩小一般C_n上界。见 rounds/r017/report.md。
- **第18轮绝对下界项**：固定mu指数斜坡的软根已精确到o(1)，从而得到该实际下界族自身的乘法渐近等价；结合第17轮偏移，一般log C_n≥nmu−d_n−(mu+3)log n/24+o(log n)。这仍不是一般C_n的乘法等价。见 rounds/r018/report.md。
- **第19轮窄平台改善**：宽w=log n/2+loglog n的高度1平台使实际分母因子趋于1，只将软根降低1/2+o(1)；一般log C_n≥nmu−d_n+(9−mu)log n/24+o(log n)，对数级系数转为正。见 rounds/r019/report.md。
- **第20轮一般上界**：有限时域Gaussian离散屏障给log C_n≤nmu−d_n+O((log n)²)。一般上下界对数差距现为O((log n)²)，真正log n系数及一般乘法等价仍未闭合。完整证明见 rounds/r020/report.md。
- **第21轮一般对数级上界**：精确非线性Gaussian移位、C³截断与热半群积分给log C_n≤nmu−d_n+O(log n)。结合第19轮下界，一般上下界对数差距现为O(log n)，最优系数及乘法等价仍待研究。见 rounds/r021/report.md。
- **第22轮Q_n精确对数项**：任意确定性独立比例向量和任意单调轮廓的log Q_n=nmu−d_n+(9−mu)log n/24+o(log n)；机制上界匹配第19轮下界。此结论未扩展至一般C_n或反馈模型。见 rounds/r022/report.md。
- **第23轮实际相关输入反例**：一个缺角两层输入严格优于相同振幅上限下的全部常幅矩形；振幅上限1.01时认证比值≥1.020472，而常幅矩形≤1.019754。固定二维近单位收益至少gamma epsilon+O(epsilon²)，gamma=(6log(4/3)−1)/8>0；不改变一般C_n渐近结论。见 rounds/r023/report.md。
- 第24轮：高度1+1/n的满立方体近常值输入类最优比值趋于I_0(2sqrt2)，并证明三个小振幅尺度上的满立方体支撑障碍。见 rounds/r024/report.md。
- 第25轮：宏观角分块的完整锚定weak-halo仍被Q_n控制；辅助反馈热方程的精确线性响应含g log t，但未据此声称离散反馈渐近。见 rounds/r025/report.md。
- 第26轮：固定块数、次指数高度及任意确定比例的单锥接收评分≤C_m sqrt(n/∏b_j) Q_n；平衡块损失(m−1)log n/2。见 rounds/r026/report.md。
- 第27轮：共同参数的增量卷积密度在四阶起对数凹，三阶有精确反例；完成曲线接收体积及尾截断，给出精确反馈核参数导数和局部telescope。见 rounds/r027/report.md。
- **第28轮全输入固定比例障碍**：双射线全局重排使任意非负输入的固定比例反射锚定族端点精确等于对应标量Q（含有限高度）；K组固定锚点/比例的逐点最大至多KQ_n。另严格量化中等下尾带的真实一步贪婪收益p x² f/(8k²)，累计占据量仍待证明。见 rounds/r028/report.md。
- **第29轮反馈模型分离**：精确Poisson校正子、全状态比较及有界软轮廓证明liminf(log mathcalR_n−log Q_n)/log n≥(mu−1)/8；从而log mathcalR_n≥nmu−d_n+(3+mu)log n/12+o(log n)。这仅是上界模型的下界和机制障碍，未提高真实C_n下界。见 rounds/r029/report.md。
- **第30轮连续预算与熵障碍**：有限真实接收点预算有达到的精确强对偶及盒/尾认证接口；单步核移动跳点使KL随动作偏移线性，排除二次Fisher路径比较。独立比例交替仅在离散族内改善，仍低于单块参考；固定旧输入的跳跃分段积分把原网格平均1.00061修正为约0.99809，未登记连续下界。见 rounds/r030/report.md。
- **第31轮块正则性**：b≥2块核具有二次局部KL；有效界为O(bh²)。严格证明I_b(mu)→p_mu²/2及精确矩下界，同时给出b h⁴漂移障碍，禁止误用全块长统一O(h²)界。旧接收输入连续数值修复后平均约1.0000304，但仍低于单块参考。见 rounds/r031/report.md。
- **第32轮任意输入接口**：对所有可预测策略证明实际Orlicz预算的Fenchel高层尾吸收，将任意单调轮廓归约到指数/线性单终端，误差O(log n/n)。精确独立桥公式给出g log n源质量及全局对数包络反例；有限正比例参数无一阶随机序。局部接收比例小收益被统一比例对照超过，未新增C_n认证下界。见 rounds/r032/report.md。
- **第33轮正确终端结构**：证明payoff Bellman全局唯一动作根与前沿恒等式、非负running-cost等价、带中心余项的局部全动作收益；给出固定短时Gamma平移初始化和正指数外侧屏障。完整MAX递推至128仅为插值数值，global-c校准后一个独立坐标正向探针未增益，sharp全局校正仍未闭合。见 rounds/r033/report.md。
- **第34轮全控制外侧接口**：构造精确仿射/Gamma终端修补，给出自然时间管两侧的全策略边界及固定短时起点的参考桥外带源求和；纯Exp占据反例与加权大跳相消排除两种拼接捷径，结合增长屏障证明有缓冲内部的全动作跨界界。小坐标扰动获益仍被统一比例细化超过，未新增C_n认证下界。见 rounds/r034/report.md。
- **第35轮局部源与全局反例**：匹配参考局部全动作收益保留bQ/n并有可求和误差；证明固定参考多项式校正直接与两侧屏障取min在正大偏差区必失败。给出真实单分量halo外的跨分量接收方盒，以及n≥8同步二维积见证的有限分母障碍。一般锐系数仍未闭合。见 rounds/r035/report.md。
- **第36轮非线性局部闭合**：参考payoff对数凹性给出全动作Doob指数尾，控制指数化校正的Jensen缺陷；精确矩消除闭合自然时间管内MAX势，临界附加误差仅O(f_n)。全状态延拓仍未闭合。曲线/source联合变形按有限锥体积验收后无收益，三种接收点依赖方向也未胜过同输入对照。更正第34轮接收试验维数为n=64（256为求积节点）。见 rounds/r036/report.md。
- **第37轮晚时段全空间比较**：匹配仿射余量使局部非线性势在 t≥N^gamma 可与全局外屏障严格拼接；初始正缺口面积为 O(log t)，但常值偏移给出任意面积误差传播的反例，尚缺带控制损失的初始化比较。Mathematica 30组双侧变换诊断及5组真实 mixed source/receiver 均未登记新认证下界。见 rounds/r037/report.md。
- **第38轮统一面积与控制损失**：Gaussian W1比较将匹配初始化O(log T)面积扩展至所有固定gamma<1；裸全局参考移位必须支付T log(1+e^-a)，不能用O(log T)代替。精确regret递推与相对粗屏缺口给出保留反馈的传播条件，整个晚时段左外管均有mean-loss/N杀死界。15组真实mixed输入扫描中最佳全33组数值通过，仍低于同维标量且未认证。见 rounds/r038/report.md。
- **第39轮固定短时初始化及中间时段反例**：证明任意固定k≥4的全空间真值可由有限密度源+非线性势初始化；但对alpha充分接近1，现有二次势与粗屏障直接取min在t约sqrt(logN)必失败。固定正倾斜Fourier prefactor补引理闭合一步比值，失败不能靠有限密度常数修复。真实权重/比例联合优化最佳评分71.28080，仍低于同维标量且未认证。见 rounds/r039/report.md。
- **第40轮正侧局部剪切与交接尺度**：在固定正倾斜带利用inner≥Oplus强制右侧保留粗屏障，严格排除r039再进入；替代率修补Istar−epsilonD具有负Hamiltonian且保留固定k初始化。统一Mills主项给固定mu无限终端交接常数(6epsilonA3)^(-1/3)，中间单步仍未闭合。32组单步诊断未显示率修补普遍改进；几何变形最佳pruned评分71.29434，仍低同维标量且未认证。见 rounds/r040/report.md。
- **第41轮正侧中间带闭合**：有限N终端的带导数Mills引理、精确中心化payoff桥及全动作Doob尾界，将一步误差降至O(eta/t+eta⁴)，在A sqrt(log t/t)≤eta≤d epsilon严格闭合全动作比较。早期左侧剪切仍未闭合，故一般锐系数不变。Mathematica完成12组导数、120动作点及非对称分块加密复算；所有候选仍低于标量基准且未认证。见 [第41轮报告](rounds/r041/report.md)。
- **第42轮左侧参考全动作闭合**：在a=mu、无限终端、无额外源/floor范围内，真实剪切代价的接触导数缓冲、锐插值及归一化左尾控制闭合lambda介于A sqrt(log t/t)与C t^(-1/3)的全动作比较。50/80位Mathematica复算支持；匹配有限N和全空间拼接仍待证明，真实C_n账本不变。见 [第42轮报告](rounds/r042/report.md)。
- **第43轮全空间Bellman势**：有限N匹配左带保留真实仿射缺口、密度源与floor，并以对数凹score闭合早期左带外；共用常数完成自然管、正侧、晚左区和右强制clip拼接。整个屏障固定平移完成初始化，临界代价O(f_N)。端点锐系数最终换算尚待结算。20组双盒/3025矩形精确扫描及独立整数证书确认跨分量有限族面积增益，无高维下界改善。见 [第43轮报告](rounds/r043/report.md)。
- **第44轮反馈模型锐对数项**：log R_n=nmu−d_n+(3+mu)log n/12+o(log n)。一般C_n的归一化log n系数liminf/limsup夹在c_Q=(9−mu)/24与c_R=(3+mu)/12之间，宽度g=(mu−1)/8；不声称真实系数极限存在。现有20个二维见证的所有异构积接收提升在k≥2严格小于1。38点、16根及最高1476位复核支持端点计算。见 [第44轮报告](rounds/r044/report.md)。
- **第45轮全状态单向块核比较**：精确inf_r ∂c log f_(b,c)=bp_c−D_b/c，且D_b=b(mu−1)+(mu−1)/2+O(1/b)，得到驻点附近−C[1+b(c−mu)²]全状态下界及近mu窗口单向覆盖。另给真实受限矩形族的精确2epsilon²共享侧边gap；没有将kernel的1+o(1)误差误作锐模数预算误差，一般C_n区间未变。见 [第45轮报告](rounds/r045/report.md)。
- **第46轮实际预算外带吸收**：固定块数、任意非负轮廓的逐块外带贡献可压至任意N^(-K)；近mu参数下半径A sqrt(b logN)的严格阈值为A>sqrt(2muK)，低于阈值有精确预算象限反例，另给非积条带及稳定高维Fourier核。中心选参和共享侧边几何仍未闭合，真实C_n区间不变。见 [第46轮报告](rounds/r046/report.md)。
- **第47轮近最优平台修补障碍**：同一达到Q_n对数项的平台输入，在比例窗口内有精确D(t)(logN)²/N平均亏损；整体振幅修补花费D(t)(logN)²对数预算，而阈值重定位只需D(t)logN。另给内部Bregman-score控制及端点带线性障碍、真实有限矩形的有理见证。一般C_n区间未变。见 [第47轮报告](rounds/r047/report.md)。
- **第48轮下端带平方根稳定性**：实际Fenchel缺口给局部尖锐sqrt控制；精细比例bar c=mu−sqrtmu z_N/(2sqrtN)使临界score趋于−p/2，局部系数从(logN)³/N降至logN/N。修复固定Fourier截断不足导致的2^4096计算失败。未闭合上端/soft区及多块角向积分，一般C_n区间未变。见 [第48轮报告](rounds/r048/report.md)。
- **第49轮双块角向响应矩阵**：共同参数下以导数Fourier曲率和凸切线尾闭合全角向二阶矩；双侧下端位移的共同/对比方向系数分别为(3p²/2)logN/N与(p²/2)(logN)²/N。独立块动作仍有额外logN因子。同步输入/标签的整体平移保持完整V/J不变，截断则另行结算。一般C_n区间未变。见 [第49轮报告](rounds/r049/report.md)。
- **第50轮完整预算Hessian**：合并下端、soft、核曲率及实际平均约束的Schur校正，给出固定活跃面多接收公式；同一Phi₂模型以有理区间认证删核曲率或soft会误判预算曲率符号。固定规范后的相对几何最佳数值评分71.305889仍低标量且未全节点可行。一般C_n区间未变。见 [第50轮报告](rounds/r050/report.md)。
- **第51轮实际块核曲率**：证明密度/CDF二阶参数恒等式，实际Phi_N预算完整曲率通过高精度核公式及固定窗口重优化差分复核；正曲率仅为有限维数值现象。统一参数内层两步权重修复使65锚点全部数值可行，评分71.34730744仍低标量、对偶间隙3.03%，未新增认证下界。见 [第51轮报告](rounds/r051/report.md)。
- **第52轮Ward分解与soft障碍**：Fenchel密度分部积分严格闭合正象限边界，将完整曲率分解为均值平移与明确形状/响应余项；证明实际soft支持内局部系数对任意乘子均有双符号。共同参数数值驻点复核至1024维，匹配变参数公平候选评分71.35120790，仍低标量、未认证。一般C_n区间未变。见 [第52轮报告](rounds/r052/report.md)。
- **第53轮完整响应闭合**：漂移Fisher、cap乘子及推前密度卷积界严格给等双块单接收模型的(λ/P)B=O(√N(logN)²)=o(N)，完整曲率只剩Ward形状余项。修正尺度括号后，平坦几何65锚点候选评分71.35856006仍低标量、未认证；一般C_n区间未变。见 [第53轮报告](rounds/r053/report.md)。
- **第54轮严格认证与完整曲率**：显式有限L=10^8输入给log C64>71.79；旧65点固定折线见证评分≤71.668105203，形成受限机制障碍。单接收等双块的Ward形状余项闭合，完整对比曲率∼Np_mu；一般C_N对数系数区间未变。见 [第54轮报告](rounds/r054/report.md)。
- **第55轮实际双接收与规模障碍**：对称双约束精确化为平均核，完整预算响应非负；局部双锥体积带来严格一阶增益，但平衡双块有限锥要达Q_N规模，其点数和横向跨度均至少为√N量级。单接收对偶归一化权重总变差≥cN^(1/4)，排除正概率替代。双点浮点评分69.5112未改善C64认证下界；一般区间未变。见 [第55轮报告](rounds/r055/report.md)。
- **第56轮共同参数驻点闭合**：目标1、等双块正则共同驻点曲率为Np_mu+o(N)，并给出角向最优共同动作的二次偏移及四次收益公式。实际步函数有共享regret=1/4却逐高度覆盖缺口恒零，排除层覆盖替代完整regret。一般区间及C64认证下界未变。见 [第56轮报告](rounds/r056/report.md)。
- **第57轮全窗口共同曲率**：Ward响应归一化rB00=1+O(N^(-1/4))闭合非驻点Schur抵消，整个正则共同窗口曲率为Np_mu+O(N^(3/4)(logN)^3)。冻结源计算显示径向补偿大幅削弱混合收益；有界实际阶梯反例排除任意固定阶的普适共享证书。一般C_N区间未变。见 [第57轮报告](rounds/r057/report.md)。
- **第58轮校正系数与实际高维分母**：完整中心化响应给|K|=O(log³N/N²)、局部径向四阶系数O(log⁶N/N³)，尚无自然角窗的统一余项控制。512/1024维数值K变负；固定阶反例族达到反馈点值尺度需付出log²N级归一化Orlicz价格。一般C_N区间未变。见 [第58轮报告](rounds/r058/report.md)。
- 主攻：移位源曲面的统一响应尾/角余项，以及危险输入的一般秩—预算关系。
- 一般前三项现已闭合；乘法渐近等价、后续锐余项和极值结构仍待研究。

## 复现

计算优先使用WolframScript。MATLAB按需短批运行，保存结果后立即退出，不保留各子代理的常驻引擎。

MATLAB R2026a（Optimization Toolbox）：
```matlab
addpath('code');
sm_checks;
rows=sm_baseline;
```
Wolfram 15：`wolframscript -file code/crosscheck.wl`。
基准卷积、优化和Bellman递推均为确定性；第7、9轮几何诊断使用明确记录的固定种子。
结果摘要保存在 `rounds/r001/`；数值精度检查不等于区间认证。

## 研究约定

每轮分别记录严格定理、受限模型结论、数值候选与证明缺口。
每轮结束仅归档必要报告、代码和摘要；一轮一个提交。
下一轮接续点见最新轮次报告，不以低维改进宣称高维次项突破。

第2轮复现：MATLAB `rows=sm_multilayer;`，Wolfram `wolframscript -file code/multicheck.wl`。

第3轮复现：MATLAB `rows=sm_exponential;`；Wolfram `code/softproof_checks.wl`、`code/exponential_check.wl`；Python `code/soft_constants.py`。

第4轮复现：`python3 code/certify_staircase.py`；MATLAB `rows=sm_dense;`。

第5轮：完整首系数障碍及有界核心归约见 rounds/r005/report.md；单盒数值核对 `rows=sm_box;`。

第6轮：加权与分块精确支配见 rounds/r006/report.md；复现 `rows=sm_weighted;`、`code/weighted_checks.wl`。

第7轮：精确多中心膨胀障碍见 rounds/r007/report.md；MATLAB `result=sm_union_checks;`。

第8轮：完整证明与一般上界的显式维数障碍见 rounds/r008/report.md；MATLAB `rows=sm_large_soft;`，Wolfram `code/large_soft_check.wl`。

第9轮：非单调轮廓重排见 rounds/r009/report.md；MATLAB `result=sm_rearrangement_checks;`。

第10轮：缺口恒等式、缩放障碍与算子反例见 rounds/r010/report.md；Wolfram `code/deficit_checks.wl`。

第11轮：一般反馈上界归约见 rounds/r011/report.md；MATLAB `sm_bellman(h,256,outputPath)`，Wolfram `code/bellman_checks.wl`。

第12轮：一般首个次指数项与似然比鞅证明见 rounds/r012/report.md；MATLAB `result=sm_tilt_checks;`，Wolfram `code/tilt_checks.wl`。

第13轮：一般完整三项展开见 rounds/r013/report.md；MATLAB `sm_controlled_gaussian_checks`，Wolfram `code/stopping_checks.wl`。

第14轮：隐式正态根和任意固定阶逆对数展开见 rounds/r014/report.md；MATLAB `sm_root_checks`，Wolfram `code/root_checks.wl`。

第15轮：精确指数控制和局部动作偏移见 rounds/r015/report.md；MATLAB `sm_laplace_checks`，Wolfram `code/laplace_checks.wl`。

第16轮：真实下分位数的对数级改善见 rounds/r016/report.md；MATLAB `sm_quantile_derivative_checks`，Wolfram `code/quantile_derivative_checks.wl`。

第17轮：实际指数软轮廓下界的对数级改善见 rounds/r017/report.md；MATLAB `sm_soft_shift_checks`，Wolfram `code/soft_shift_checks.wl`。

第18轮：固定比例斜坡的乘法渐近及一般显式下界见 rounds/r018/report.md；MATLAB `sm_edgeworth_ramp_checks`，Wolfram `code/edgeworth_ramp_checks.wl`。

第19轮：窄平台实际下界见 rounds/r019/report.md；MATLAB `sm_plateau_checks`，Wolfram `code/plateau_checks.wl`。

第20轮：一般上界的对数平方余项见 rounds/r020/report.md；MATLAB `sm_gaussian_barrier_checks`，Wolfram `code/gaussian_barrier_checks.wl`。

第21轮：一般对数级上界见 rounds/r021/report.md；MATLAB `sm_shift_barrier_checks`，Wolfram `code/shift_barrier_checks.wl`。

第22轮：独立比例机制的精确对数项及初始化障碍见 rounds/r022/report.md；MATLAB `sm_independent_log_checks`，Wolfram `code/independent_log_checks.wl`。

第23轮：实际缺角输入与近单位振幅见 rounds/r023/report.md；Python `code/cert_capped_corner.py`，MATLAB `sm_capped_corner_checks`，Wolfram `code/capped_corner_checks.wl`。

第24轮：近常值满立方体类的精确极限与三尺度支撑障碍见 rounds/r024/report.md；MATLAB `sm_nearflat_ceiling_checks`，Wolfram `code/nearflat_ceiling_checks.wl`。

第25轮：宏观角分块的锚定障碍与反馈辅助热方程见 rounds/r025/report.md；MATLAB `sm_feedback_probe`，Wolfram `code/angular_feedback_checks.wl`。

第26轮：微观块凸预算与单锥接收的对数损失见 rounds/r026/report.md；MATLAB `sm_microblock_dual_checks`，Wolfram `code/microblock_checks.wl`。

第27轮：共同参数卷积的对数凹证书、曲线接收边界及精确反馈核见 rounds/r027/report.md；Python `code/cert_block_logconcavity.py --verify`，MATLAB `sm_microblock_halo`，Wolfram `code/block_logconcavity_checks.wl`。

第28轮：固定比例反射锚定族的全输入精确端点与带内反馈收益见 rounds/r028/report.md；MATLAB `sm_greedy_band_checks`，Wolfram `code/fixed_anchor_checks.wl`。

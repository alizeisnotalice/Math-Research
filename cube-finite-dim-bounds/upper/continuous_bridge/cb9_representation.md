# cb9：体对角 Cauchy 的20锥、正密度与真实 Green 接口

2026-10-25，第9轮固定frame首批。n=3，V={(1,±1,±1)}，L=−Σ_v|v·D|，Fourier inverse采用(2pi)^−3∫e^(ix·xi)fhat(xi)dxi。所有方向未单位化，共享同一Cauchy时间t。本分支不做g/J或域形状优化；以下不是新U3/U4。

## 1. 输入锁定及20锥的准确覆盖

主body_diagonal_frame_interface的输入SHA为8cdcd66d983404174224e8759820865938d0b54fec1c2e8a8dfe4965fde32bea，准确符号是

    psi(xi)=max(4max_j|xi_j|, 2Σ|xi_j|).

其单位球为P={max|xi|≤1/4,Σ|xi|≤1/2}。令Q=4P，则Q={max|y|≤1,Σ|y|≤2}，12个顶点是(±1,±1,0)的置换。

完整facet描述：六个y_i=±1面上，另外两坐标满足|y_j|+|y_k|≤1，是四顶点diamond；八个有符号和s·y=2面上，s_j y_j≥0、每坐标≤1，是三角形，顶点(s1,s2,0)、(s1,0,s3)、(0,s2,s3)。这些是Q的6+8条定义不等式的非冗余facet；没有其他facet。

对固定i与sign=s，j<k为其余坐标，写A=s e_i+e_j、B=s e_i+e_k、C=s e_i−e_j、D=s e_i−e_k。方形面用(A,B,C)、(A,C,D)两三角分割；各三角面保持一块。得到12+8=20个三生成元锥，完整整数列表在cb9_cones.json。

**覆盖与不重叠：** 原点在Q严格内部。每条非零射线与∂Q恰一交点，除facet边/三角剖分边外，恰落一个面三角的相对内部。因此其正锥覆盖R³，锥内部互不交叠；多重覆盖只在有限低维面，Lebesgue零测。该结论来自完整facet/radial描述，不以体积相等或有限点测试替代。

每锥ray矩阵R的|detR|准确为2；每个ray r满足psi(r)=4。若xi=Σ a_j r_j,a_j≥0，归一到面凸组合可得psi(xi)=4Σa_j。因此每个锥上都有准确线性相位，且P的体积为20·2/(6·4³)=5/48。cb9_derive.py实际SymPy逐一核20个det及合并多项式。

## 2. 完整复锥密度及全实正有理式

t>0时e^(−tpsi)可积。按上述真实锥换元并作三个绝对收敛的一维Laplace积分，得到

    p_t(x)=1/(4pi³) Σ_(20 cones) Π_(j=1)^3(4t−i x·r_j)^−1.

这里已经包括|detR|=2与(2pi)^−3。相反锥给共轭项，总和为实，但单项实部可以为负，不能删负项或把单项取绝对值。原点给20/(4pi³(4t)³)=5/(64pi³t³)，与谱单位球layer-cake校准一致。

设a=4t,u=a²，六个无向ray线性形式是x1±x2、x1±x3、x2±x3。记

    D(u,x)=Π_(i<j)[u+(xi−xj)²][u+(xi+xj)²],
    X=x1²,Y=x2²,Z=x3²,
    S1=X+Y+Z, S2=XY+YZ+ZX,
    C=X(Y−Z)²+Y(Z−X)²+Z(X−Y)²=S1S2−9XYZ,
    N(u,x)=20u4+24S1u³+4(S1²+5S2)u²+4Cu.

则完整合并后的准确密度为

    p_t(x)=a N(a²,x)/[4pi³ D(a²,x)].                  (1)

证明可直接逐锥取实部：对于ell_j=x·r_j，

    Re Π(a−iell_j)^−1 = a(a²−Σ_(i<j)ell_iell_j)/Π_j(a²+ell_j²).

每锥有三个不同无向ray，分母乘其余三个二次因子，20项相加后的N正是上式；cb9_derive.py实际展开核到四个u系数，零次系数准确为0。此为多项式恒等式，不依赖x泛型或节点互异，故重复/零线性形式同样成立，且没有丢弃任何复项。

C是显式平方和，u>0时N>0、D>0，因此(1)从代数上为正。代码直接用平方和计算C，不以abs或max(0,.)修正。t>0不存在实分母零，即使线性形式重复或为0也无异常。

## 3. 独立的正概率表示与正常数

取四个独立标准一维Cauchy变量C_t^j。以v1=(1,1,1)、v2=(1,1,−1)、v3=(1,−1,1)为列的B满足|detB|=4，v4=−v1+v2+v3。其和X_t的特征函数正是e^(−tpsi)，所以p_t是真实保守Markov密度，而非某个只在锥上的形式函数。

令

    b=(0,−(x2+x3)/2,(x1−x3)/2,(x1−x2)/2).

对第四个Cauchy变量条件化，得到另一准确正表示

    p_t(x)=t4/(4pi4) ∫_R Π_(j=1)^4[(s−b_j)²+t²]^−1 ds.   (2)

其正常数同时支付四个一维pi与basis det4。数值交叉用(2)独立积分，不导入锥积分或合并N。

若alpha_j≥0、Σalpha=1，V(alpha)=Σalpha_j b_j²−(Σalpha_j b_j)²，则Feynman参数法给一个也允许+infinity的正Green表示

    G(x)=3/(32pi³) ∫_(simplex3) 1/V(alpha) dalpha.       (3)

simplex积分用dalpha1 dalpha2 dalpha3、alpha4=1−sum，体积1/6，不是归一概率测度。证明：ΠA_j^−1=6∫(Σalpha A)^−4；积分s给5pi/[16(t²+V)^(7/2)]，再积分t4/(t²+V)^(7/2)给1/(5V)。所有项非负，Tonelli允许包含无穷极限；合并常数是3/(32pi³)。两积分正常数已实际SymPy核，见cb9_normalization.json。式(3)未用作大批成本数值，不重复root的g/J探索。

## 4. Green正有理积分：先合并，再处理退化

由(1)、a=4t,u=a²，有真正扩展值势核

    G(x)=∫_0^infinity p_t(x)dt
        =1/(32pi³) ∫_0^infinity N(u,x)/D(u,x) du.       (4)

因p_t≥0，这个换元不存在带符号级数互换问题。它与20个锥积分分别收敛与否无关。G齐次degree−2，继承全部坐标符号/置换对称。

令m为六个平方线性形式A_j=(xi±xj)²中准确等于0的个数。对非零x：

- m=0：被积在u=0为O(u)，有限；
- m=1：N与D各含一个u，约去后有有限非负极限；
- m=2：只能是coordinate axis，此时C=0，N、D各含u²，约去后仍有限；
- m=3：恰是|x1|=|x2|=|x3|>0的四条无向body-diagonal线，C=0但N只含u²，D含u³，被积正比1/u，G=+infinity。

非零点没有其他m；原点m=6，G也无穷。u→infinity时N/D~20/u²，所以尾部收敛。这个分类直接付清零form及重复平方因子的极限，不能把所有0/0点改成有限常数。

特别x=(r,±r,±r),r≠0，(1)给

    lim_(t↓0) t p_t(x)=3/(32pi³r²)>0,

故固定line点有真实对数时间发散。代码返回+infinity，不删除这些线或把数值大值截断。这里没有声称已证明完整的近line空间渐近系数。

**off-line连续极限。** 式(2)中，三重b碰撞恰对应上述四条frame线。离它们的紧集上，所有三元b簇的直径有统一正下界delta。对任意s，至多两个b_j落在s的delta/3邻域，故至少两因子≤(3/delta)²，余下两Lorentz因子用Cauchy–Schwarz积分，∫[(s−b)²+t²]^−2 ds=pi/(2t³)，得到统一p_t(x)≤C_K t（t≤1）。大时间用p_t(x)≤p_t(0)=5/(64pi³t³)。由支配收敛，G在全部非line点连续，包括zero/repeated forms。因此(4)在这些点给的值就是其他方向趋来的真实有限极限，不是任意补值。line/原点处用非负p_t和Fatou给G(x_k)→+infinity。

另∫_E G≤∫_0^1∫p_t+|E|∫_1^infinityp_t(0)<infinity对有界E成立，故G∈L1_loc，包括line邻域。由各p_t dx的Tonelli，势测度本身有该Lebesgue密度；截断后生成元的正源/极限仍由upper_theory另证，不能用此局部可积性替代。

## 5. log divided differences的准确使用边界

若一个锥的三个ell均非零且互异，部分分式及principal log给

    ∫_0^infinity Π(4t−iell_j)^−1dt
      =1/4 Σ_j log(−iell_j)/Π_(k≠j)(ell_j−ell_k).

其中log(−iell)=log|ell|−i(pi/2)sign(ell)，相反锥的虚部准确抵消。若全部锥均无zero form，可取整体实部得到

    G(x)=1/(16pi³) Σ_cone DD2(log|.|;ell1,ell2,ell3).   (5)

重复非零ell使用Hermite极限：DD2(f;a,a,a)=−1/(2a²)，DD2(f;a,a,b)=[log|b|−log|a|−(b−a)/a]/(b−a)²。该过程在重复值不为0处合法；不以“小分母容差”假装独立节点。

但zero form时，有些单锥积分出现±log发散，即使真实G有限也不能先逐项积分。x=(1,1,0)就是有限G但有zero form的例子。式(4)先完整合并后才积分，规避了这种未付取消。不能把(5)直接代入log0，也不能按各项符号删除发散项。

## 6. 两个准确退化锚点

坐标轴x=e1由(4)或(2)给

    G(1,0,0)=3/(8pi³).

edge x=(1,1,0)上，R(u)=N/D准确约为

    4(5u+2)/[(u+1)²(u+4)]
      =8/(u+1)−4/(u+1)²−8/(u+4).

因此

    G(1,1,0)=(4log2−1)/(8pi³).

cb9_normalization.json实际SymPy核该分式及积分。此点只能给cube下确界g≤G(edge)，不是完整face最小值证明。沿body diagonal的cube角点反而G=+infinity，不能沿用其他核族的角点最小。

## 7. 可调用数值模块与有限实跑

模块cb9_frame_green.py SHA256为`d23501eeaf417a736730a6e61398e27b9364c231d35ec6ee8bf573e357058bf5`，root已锁定使用。接口：

    green(x, method='auto', return_error=False)
    green(x, method='mp', dps=60)
    green_face(y,z)           # G(1,y,z)，标量/数组
    green_array(X)           # 最后一维为3
    heat_kernel(t,x,method='rational'或'cones')

默认先按max|x|齐次归一，再令u=w/(1−w)把(4)化为[0,1]正有理积分。约去准确zero form的共同u因子，body lines直接返回inf；按A_j/(1+A_j)分段。Scipy quad返回的error只是估计，不是连续区间证书；默认双精度定位用于通常尺度的face发现，高精度用mp分支，极端浮点尺度/近line仍须谨慎。本模块没有给g或J的严格误差包围。

cb9_check.py按预登记固定12泛型点+6退化类实际运行，合计18个不同空间点，未追加网格/随机点。每点在t=3/7比较全部复锥和、正有理p、真实Cauchy正积分(2)；有限G比较默认quad与60位正积分，允许式(5)的点另与完整DD/Hermite计算比较。六退化类为origin、body diagonal、axis、coordinate-plane edge、zero form及repeated nonzero forms。

实际结果：p最大相对差7.70×10^−16，有限G双精度对60位的最大相对差3.70×10^−16；DD/Hermite相对差均<10^−45；两个准确G锚点与60位积分一致。origin/bodyline的G准确分类为inf，p_t仍正常。该批约0.58秒，不是成本优化。共享frame定义、SymPy/mpmath/Python/Scipy数学库，计算自检不代独审。

主随后实际串行MATLAB_R2026a独立batch执行便携cb9_matlab.m，脚本SHA `eb5b162b0756299f8e1d97978af69b0fd1499421c3f7f9cf608820683700f181`，exit0；20个det2锥、同18点及axis/edge/infinite-line全部PASS，p最大相对差7.33×10^−16、G最大相对差3.28×10^−16。准确执行hash、原始log和parent receipt来源见cb9_matlab_receipt.json；没有clear/path/reset共享会话。该执行仍是共享模型的浮点交叉，不能替代完整g/J认证。

## 8. 来源、状态与未闭合事项

第9轮snapshot的ledger历史SHA与当前工作树不同，已用不可变base commit 4ea9d04f...实际git读取验证历史SHA，并逐项比较tasks完全一致；只有current_round、current_round_status、last_remote_confirmed_commit变化。历史hash未被刷新抹掉，详见cb9_snapshot_receipt.json。其他公式源逐文件hash均匹配，来源记录cb9_sources.json；外部design只提供猜想，不被当证明。

严格解析：完整20锥覆盖、真实p_t与G表示、全部复项/正常数、zero/repeated forms与线无限值及off-line连续极限。
实际计算：准确多项式/积分锚点和预登记18点交叉；数值接口已交root，无本分支域/成本搜索。
候选/依赖：G的正部截断、有限frame比较、全L1/尺度及真实g/J认证由upper_theory/root/独审分别接续。本式不自动证明K_c可用或成本<8。
未取得：全cube g、J严格区间、任何新U3/U4；root反馈32/64成本探路尚未低于8，仅为它的数值发现，不能当固定frame或其他核族障碍。

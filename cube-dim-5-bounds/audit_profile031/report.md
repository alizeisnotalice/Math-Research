# profile031：全rho形状族的正测度patch预算障碍独立审计

## 结论及禁止误读

**通过。** 主 round009/profile_barrier 的点预算及完整 patch 链成立，
覆盖所有有限 rho<1、所有 a0>=0，允许参数随空间 x 选择，但不允许
它们依赖尚未积分的随机形状 M。整个指定族在明确正四维面片上满足

    sup_(rho<1,a0>=0) Young证书值 < 3/25000 = 0.00012
                                           < 15198/10^8.

同一球邻域论证还直接覆盖下述正五维体积小盒，因而可以排除仅靠忽略
边界零体积集合来绕过 essential cube 方法障碍的解读。

这只是该明确证书族的预算上限，不是 G 的上界、不是 C5 的下界，也
没有新的 C5 界。真实 G(x0)=infinity 与有限预算完全兼容；本审计没有
在奇异集上对 G 作连续延拓。

## 最终锁与共同模型

只消费最终成功版本，生产源码和登记分别锁为

    e648f8e4a8c277c0b9e5bf4fc21ca8ebbd5fb5f61cf6398599448cc1517b4627
    53e887eab22f588525a4eb568b3d1b19bb9af0ecf60289e27e67c47bd3e4b06f.

结果和证明文档分别锁为

    045c87233755809c36833d8f8c7e75b32dbab7025de3febe6ac6143e4459e9b7
    69f0a857b86fe1872be09b3f12b8a91788295a28d969d0ad5e00c877f99a76ad.

六项来源见 registration.json。未运行生产脚本、其5x5 SymPy矩阵或
030的两个固定rho基线；此前结构相等/系统Python失败均不作为结果输入。

共享模型来自已审 upper_alpha1_029：alpha=1、beta=1/2，16条未单位化
方向 v=(1,+/-1,+/-1,+/-1,+/-1)，原始 S_i 独立且
`E exp(-z S_i)=exp(-sqrt(z))`。A=ΣS_i v_i v_i^T，T=ΣS_i，M=A/T。
形状概率 R 是以 T^(-1/2)/Ztilt 倾斜再推送的概率，

    Ztilt=1/(8sqrt(pi)), c1=1/(32pi^3),
    G(x)=c1 E_R[d^(3/2) B_x^(-2)],
    d=det M, B_x=x^T adj(M)x.

本审计使用029已接受的同模型正表示及S2奇异集结论，没有重建其完整
Cauchy半群/有限源传递。030只消费通用Laplace/Mellin解析框架，不消费
两固定rho的数值结论或MC诊断。

## 1. 族、正lambda及完整矩的有限性

令 W_x=Σ_(j<k)M_jk x_j x_k，u_x=(W_x+2)/12。对每个方向 v，
W_v(x)=((v dot x)^2-|x|^2)/2 在cube上多仿射，极值在顶点，顶点值为
-2、2或10。M是vv^T的凸组合，故0<=u_x<=1。于是所有rho<1都使

    lambda_i(x)=1-rho*((v_i dot x)^2-|x|^2+4)/24 > 0,
    L_x=Σlambda_i(x)P_i=1-rho u_x > 0.

rho<0时lambda>=1；0<=rho<1时lambda>=1-rho>0。后一个全cube下端
虽可随rho趋1而退化，后面的patch比较不使用这个退化下端。

族为Z=a0 L_x^(-3/2)。Young恒等式逐点给

    d^(3/2)B_x^(-2) >= 3sqrt(d)Z^2-2B_x Z^3,

故固定rho的证书函数为`c1(3a0^2 U_x-2a0^3 V_x)`，其中

    U_x=E_R sqrt(d)L_x^(-3),
    V_x=E_R B_x L_x^(-9/2).

对非零x，U_x、V_x正且有限，最大化非负a0给

    B_rho(x)=c1 U_x^3/V_x^2.

这里B_rho是能产生的最佳下函数值，不是G本身。

## 2. Mellin常数与尺度核对

置H=Σlambda_i S_i，D=Σlambda_i^(-1/2)v_i v_i^T，C=Σsqrt(lambda_i)。
R的倾斜恰使两个T指数分别为

    -1/2-5/2+3=0,  -1/2-4+9/2=0,

因此U=Ztilt^(-1)E[sqrt(det A)H^(-3)]、
V=Ztilt^(-1)E[b_A H^(-9/2)]。没有假定T与P独立。

对独立S的一次Laplace导数为
`E[S_i exp(-z lambda_i S_i)]=exp(-sqrt(z lambda_i))/(2sqrt(z lambda_i))`。
det A和b_A分别是五次/四次的多线性Cauchy–Binet多项式，因此

    E[det A exp(-zH)]=det D exp(-C sqrt(z))/(32 z^(5/2)),
    E[b_A exp(-zH)]=b_D exp(-C sqrt(z))/(16 z^2),
    b_D=x^T adj(D)x.

正Mellin积分、y=C sqrt(z)换元给

    E[det A H^(-5)]=det D/(16 C^5),
    E[H^(-1)]=2/C^2,
    V=(128/35)b_D/C^5.

常数只需`2Gamma(5)=48`、`Gamma(9/2)/sqrt(pi)=105/16`。
所有积分在零端及无穷端都收敛，lambda严格正、D正定。Cauchy作用于
`sqrt(det A H^(-5))*sqrt(H^(-1))`，得到

    U^2 <= 8pi det D/C^7.

这同时证明U有限；V的显式式证明V有限且x非零时严格正。
即使真实Green期望发散，该组U、V仍然有限。

全部lambda同乘k>0时，U乘k^(-3)、V乘k^(-9/2)，故预算U^3/V^2
不变。该尺度核对允许下一节点预算作共同比例归一。**patch比较则始终
使用原始lambda**，没有将不同x的归一化混入相对估计。

## 3. 点x0的独立谱分解与全rho范围

x0=(1,1,1,0,0)。只看首三个坐标的四个符号模式，
`(1+epsilon2+epsilon3)^2`有一次9、三次1；尾坐标四种自由选择给总数
4和12。因此点上lambda_s=1-5rho/12、lambda_l=1-rho/12。

归一lambda_l=1，令`t=sqrt((12-rho)/(12-5rho))`。不需要复建5x5
符号矩阵：平方为9的四方向首三个坐标都是(1,1,1)，尾改号交项消失，
其外积和的head块是4*11^T、tail块是4I。由frame总和16I，D的谱为

    head平行x0：4(3t+1)，
    head垂直x0的两维：16、16，
    tail两维：4(t+3)、4(t+3).

故独立得到

    det D=16384(3t+1)(t+3)^2,
    x0^T adj(D)x0=3 det D/[4(3t+1)]=12288(t+3)^2,
    C=4/t+12.

对rho<1，12-5rho>7。清分母恒等式为

    5(12-rho)-(12-5rho)=48,
    11(12-5rho)-7(12-rho)=48(1-rho)>0.

所以实际完整范围是`1/5<t^2<11/7`。rho趋负无穷时t趋1/sqrt(5)，
rho趋1^-时t趋sqrt(11/7)，均无参数遗漏。作者使用的更松范围0<t<9/7
仍正确。

代入Mellin界，预算上端为

    B_rho(x0) <= [1225sqrt(2)/(32768 pi^(3/2))]
                        [(3t+1)sqrt(t)/(144(t+3))].

本审计以谱乘积和有理平方常数检验该化简，未调用主SymPy矩阵。
函数(3t+1)/(t+3)导数为8/(t+3)^2>0。由t<9/7、sqrt(t)<8/7，
其形状因子小于17/1890；再用pi>3与sqrt(2/3)<5/6，得到

    B_rho(x0)<2975/31850496<1/10000.

这对全部rho<1统一，不是有限rho试验。

## 4. 负rho无界与趋1时的相对lambda控制

对cube中的线段，

    grad W_i(x)=(v_i v_i^T-I)x,
    ||grad W_i||_2<=4sqrt(5)<9.

因此||x-x0||_2<=eps=1/1000时，`|Delta W_i|<9eps`。
点上u_i0=k/12，k=1或5，并且

    |lambda_i(x)/lambda_i(x0)-1|
      = |rho| |Delta W_i|/(12-k rho).

两种rho符号必须分别处理：

- 0<=rho<1时，`rho/(12-k rho)<1/7`，因为k<=5且
  `12-(k+7)rho>0`。
- rho=-R<0时，`R/(12+kR)<1/k<=1`，对任意大R成立。

所以统一有`|relative difference|<=|Delta W_i|<eta=9/1000`。
负端没有使用有界rho假设；正端也没有使用趋零的1-rho下端。点上
最小lambda本来就至少7/12，正端稳定性不会在rho趋1时失效。

对每个形状P，逐分量相对比较推出

    L_x >=(1-eta)L_0,
    U_x <=(1-eta)^(-3)U_0.

## 5. adjugate Loewner及x0特征向量步骤

由lambda_x<=(1+eta)lambda_0，得到

    D_x >=(1+eta)^(-1/2)D_0,
    C_x <=(1+eta)^(1/2)C_0.

若D>=cD0>0，对非零x取x垂直空间的正交压缩，则

    x^T adj(D)x=|x|^2 det(D restricted to x-perp).

压缩保留Loewner序，四维SPD行列式单调，故
`x^T adj(D)x>=c^4 x^T adj(D0)x`。这一步并未声称所有维数下任意
矩阵adjugate映射具有未证明的序性质，而是直接证明所需二次型。

原始D0只是点归一矩阵的正数倍，故x0始终是D0及adj(D0)的特征向量。
令e=x0/sqrt(3)、A0=adj(D0)>0，A0 e=kappa e。写x=alpha e+y、
y垂直e，则

    x^T A0 x=kappa alpha^2+y^T A0 y>=kappa alpha^2,
    alpha=e dot x>=sqrt(3)-eps.

从而

    x^T A0 x >=(1-eps/sqrt(3))^2 x0^T A0 x0
                >=(1-eps)^2 x0^T A0 x0.

这填明了作者的特征向量推论，不需要随rho可能恶化的条件数或最小
特征值统一下界。每个有限rho对应SPD；rho趋负无穷时原矩阵可整体变小，
但上述相对式和预算尺度不变性仍然成立。

结合闭式V的C分母，得到

    V_x >=(1+eta)^(-9/2)(1-eps)^2 V_0.

因此整个cube球邻域上

    B_rho(x) <=[(1+eta)/(1-eta)]^9 (1-eps)^(-4) B_rho(x0).

指数9来自U三次幂与V平方的相对尺度，不遗漏x变化的四次投影损失。
本审计准确Fraction核验

    (1009/991)^9 (1000/999)^4 < 6/5.

故球邻域内所有rho与a0的预算均严格小于3/25000。

## 6. 正面patch、正体积patch与S2

作者的明确面片

    x1=1; x2,x3 in [3999/4000,1]; x4,x5 in [0,1/4000]

距x0至多1/2000<eps，四维面测度为(1/4000)^4>0。统一预算上端在
该面片逐点成立，任意空间自适应rho(x),a0(x)也不能突破它。

更进一步，完全相同的球邻域估计直接包括

    x1,x2,x3 in [3999/4000,1]; x4,x5 in [0,1/4000].

该小盒距x0至多sqrt(5)/4000<eps，五维体积为(1/4000)^5>0。
因此对essential cube方法证书也不能靠只忽略边界面来绕过障碍。
这里不需要将预算函数错误地当成齐次函数；我们在这个五维小盒上
直接应用已证的球邻域估计。

另外，

    x0=1/2(1,1,1,1,1)+1/2(1,1,1,-1,-1)

属于S2且不在S1。029已审的alpha1奇异集结论给G(x0)=infinity。
这不与本审计矛盾：U、V是受限profile族的有限矩；最优受限Young
下函数远不足以等于真实Green。面片中S2至多占零面测度；证明也没有
利用G在x0处的连续性。我们控制的是预算而非G。

不排除任意正lambda模式、M自适应参数、Z求和/条件机制、不同I预算、
或真实G的更强下界；也不能从该族失败反推任何C5下界。

## 实际有限核验、算术库及复现

唯一预登记执行PID41408，Python3.15标准库Fraction、整数多项式和
阶乘，耗时约0.000252秒，退出码0。没有SymPy/Arb依赖，也没有运行
5x5符号矩阵、rho样本、030旧参数、概率/空间样本或MATLAB。

实际核验的是：4个head符号模式计数、谱乘积的det/adj单变量多项式、
rho范围的清分母恒等式、Mellin半Gamma有理因子、预算平方常数、统一
有理点cap、patch因子，以及四维/五维小盒距离与正测度。连续全参数
结论来自上面的解析证明，有限计算不替代该证明。

白名单为本报告、registration.json、check_profile.py、results.json、
manifest.json。源码按显式--project-root读项目相对来源，默认只核hash；
registered_origin仅作出处记录，运行不依赖原外部clone路径：

    python3 -B cube-dim-5-bounds/audit_profile031/check_profile.py --project-root /checkout/cube-dim-5-bounds
    python3 -B cube-dim-5-bounds/audit_profile031/check_profile.py --project-root /checkout/cube-dim-5-bounds --run

只写独立clone的audit_profile031/，未提交、未推送、未修改主或他人文件。

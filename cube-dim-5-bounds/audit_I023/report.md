# audit_I023：I021两层stable分组的完整I资格审计

2026-10-26。**通过。** 同一16-body、alpha=3/2独立stable模型的完整Green、全部正时间和整个正unit face满足

\[
 I=\int_{[0,1]^4}G(1,y)^{10/7}dy
 <6675729/500000000000.
\]

其显式公式值约1.3351457946355969e-5，严格改善I020旧包围。该公式仍未认证I<=11/10^6，不能由上包围大于目标反推实际I大于目标。没有新g或C5突破，后续022不在本批范围。

运行句柄 `audit_I023_fixed_001` 预登记后一次成功，未重跑生产560次Gram/adj循环、旧g、旧I基线或type权重优化。只审021锁定manifest `16eeeb3e07fc15813c7ece29d415ff79c219fa116b6a5add2f4926860c1ad42b`；未提交、未改其他目录。

## 1. 共同G、尺度抵消与完整面接口

固定16条未单位化v=(1,±1,±1,±1,±1)，权1独立对称alpha=3/2 stable时间、Lebesgue测度、逆Fourier因子(2pi)^-5。令beta=3/4，各S_i独立正stable且Laplace变换为exp(-lambda^beta)，A=sum S_i v_i v_i^T，T=sum S_i，M=A/T，d=det M。

正Gaussian混合的协方差是2t^(4/3)A，密度因子(4pi)^(-5/2)。完整时间积分给 `kappa=3 Gamma(7/4)/(2^(7/2)pi^(5/2))`。矩阵函数在A=TM下为T^(-3/4)次，且

\[
 Z=ET^{-3/4}=1/[16\Gamma(7/4)],\qquad
 G=cE_RH,\quad c=3/(2^{15/2}\pi^{5/2}).
\]

这里R用T^-beta/Z倾斜原概率后推送，未假定T/P独立。对b=3/20，5b=beta准确抵消总体尺度，故

\[
 J=E_Rd^{-b}=Z^{-1}E(\det A)^{-b}.
\]

已审轨道平均/Minkowski接口亦可直接重核：H^q，q=10/7，是空间负5次齐次函数。球面积分给(8pi²/3)d^(-3/14)，全cube的160个等价正面小片给pi²/60。轨道平均后取Lq范数将det负幂变为3/20，Minkowski于是给

\[
 I\le c^{10/7}(\pi^2/60)J^{10/7}.
\]

负矩有限性由后文完整分组证书提供，避免先假设退化形状有统一条件数。非负Tonelli允许Green在某些方向取无穷；界覆盖完整面，不需要再添加时间截断、角点管或采样误差。

## 2. Cauchy–Binet系数1/20与二维商常数

对三子集J，方向行矩阵的Gram为G_J，D_J=det G_J>0。三个不同超立方体顶点不共线，首坐标为1又把线性相关转成仿射相关，因此所有三子集秩3。商空间恰为二维。

设q_i为正交补投影。Gram/Schur行列式恒等式给

\[
 w_{ik}=\det(V_J,v_i,v_k)^2=D_J\operatorname{area}(q_i,q_k)^2.
\]

本批不用作者adj3，而以有理nullspace的两行N作坐标z_i=Nv_i，H=NN^T。因为这些坐标一般不正交，正确独立核式是

\[
 w_{ik}=\frac{D_J}{\det H}\det(z_i,z_k)^2.
\]

遗漏det H会引入错误常数；本批三个代表都准确核正D_J和det H，并各用一个5阶整数Leibniz行列式平方复核该常数。

写B_J=sum_(i<k)w_ik S_iS_k。每个五方向基包含C(5,3)=10个三子集，故CB严格给

\[
 \det A=\tfrac1{10}\sum_J(\prod_{j\in J}S_j)B_J.
\]

零商向量q_i=0恰位于span J，不进入完成因子；平行/反平行商向量的面积为0，完成系数也为0。将非零q分成无向线类ell，取代表r，有q_i=t_i q_r、a_i=t_i²>0，从而w_ik=a_i w_rk。定义

\[
 X_\ell=\sum_{i\in\ell}a_iS_i,\qquad
 L_\ell=\sum_{k\notin\ell}w_{rk}S_k.
\]

每个非零完成对在其两条线类中各出现一次，所以B_J=(1/2)sum_ell X_ell L_ell。结合前一个1/10，准确得到**1/20**；不是1/10或1/40。

令x=sum a_i、l=sum w_rk，C=x*l>0，F=(prod_J S)(X/x)(L/l)。有det A=W sum theta F，theta=C/(20W)，W=16^5。tight frame VV^T=16I给每个商的sum_(i<k)w_ik=D_J16²；并有sum_J D_J=10*16³。因此sum C=20W，theta和为1。退化完成基以零权重保留，不存在除以零完成系数；至少两条非共线商线由原方向张满五维保证。

## 3. 不同组合路径的有限准确重建

只对J=(0,1,2)、(0,1,6)、(0,3,5)三个代表做有理行消元/nullspace及商面积计算。然后用五坐标改号与相邻置换生成元在16个无向body方向上的作用，准确生成三子集及(J,line) flag轨道；没有对全部560个三子集重新计算商Gram。

三子集轨道大小160、240、160，分别具有(D,非零商数,线类数)=(48,12,3)、(64,12,4)、(112,13,5)，覆盖全部560个三子集。flag有六个轨道，大小480、480、480、160、320、320，合2240组；其中两个轨道具有同一概率型，合并为五类。

独立结果与021全部概率分数/重数和组质量逐项一致：

|scope|组数|theta总权|R_X R_L的闭式|
|---:|---:|---:|---|
|15|480|3/16|2^(1/4)+2^(-1/4)|
|15|960|3/8|2^(5/4)|
|16|320|1/8|(7+8sqrt2+3sqrt3)/2^(15/4)|
|16|320|3/16|(1+1/sqrt2)3^(1/4)|
|16|160|1/8|1+sqrt2|

这些闭式由保存的正归一系数直接作3/4次幂得到，不含抽样或未认证的形状量。每个坐标的实际scope质量都为20234240，分母20W=20971520，故rho=247/256。不是仅从平均scope大小猜测每个坐标相同；程序对16个坐标分别累加检查。

## 4. 五因子独立性、scope与分数Holder

在每个组内，三个S_j的指标、线类ell的指标和其余非零商线类的指标两两不交；零商方向的系数为0，删去不改变函数。因此F的三个单变量和两个stable线性和这**五个因子**在原始产品概率下独立。不同组一般不独立，证明没有要求它们独立。

真实scope是J并非零商方向，不随ell改变。theta/rho满足对每个原始坐标的指数占用恰为1。令h_g=theta_g/rho、f_g=F_g^(-b rho)，逐坐标广义Holder给

\[
 E\prod_g f_g^{h_g}\le\prod_g(Ef_g)^{h_g}.
\]

每一步积分一个独立坐标，对依赖它的函数应用Holder；不依赖该坐标的函数提出，余下函数保持原scope条件。可先对所有非负f截断，再单调取极限。这个fractional matching步骤不在倾斜R下实施，故没有误把R的权重当作独立。

对s>0，Mellin积分准确给

\[
 m(s)=ES^{-s}=\Gamma(s/\beta)/(\beta\Gamma(s))
      =\Gamma(1+s/\beta)/\Gamma(1+s).
\]

归一正线性和的Laplace尺度是R_X=sum p_X^beta及R_L=sum p_L^beta，故

\[
 EF_g^{-s}=m(s)^5(R_XR_L)^{-s/\beta}<\infty.
\]

于是全负矩有限并有

\[
 J\le2\Gamma(7/4)m(b\rho)^{5/\rho}
 \prod_g(R_XR_L)^{-\theta_g/5}.
\]

prefactor来自Z^-1 W^-b=16 Gamma(7/4)*16^(-3/4)=2 Gamma(7/4)。指数-1/5来自b/beta=1/5。这里s=b rho=741/5120、s/beta=247/1280、5/rho=1280/247，全部准确。

## 5. Gamma与I的另一条数值资格路径

本批只算新公式，不重跑旧基线。将上式代入完整面接口，取7次幂后得到

\[
 I_{21}^{7}=\frac{3^{10}}{2^{65}60^7\pi^{11}}
 \Gamma(7/4)^{10}m(s)^{50/\rho}
 \prod_{t=1}^{5}R_t^{-2\theta_t}.
\]

独立实现使用上述五个闭代数R_t、未移位Mellin的Gamma(s/beta)/(beta Gamma(s))，以及Gamma(7/4)=(3/4)Gamma(3/4)，区别作者的移位Gamma/直接I分数幂路径。

512-bit Arb给I21^7的向外包围；随后对两个有理端点作整数7次根算法，并逐项核整数幂方向，得到95位小数精度的有理I21夹。它完全包含于021保存的384-bit输出球，且严格低于6675729/500000000000与020旧显式球下端。该有理夹仅包围**显式公式I21**，不是实际I的下界。

共享库限制：Gamma和代数函数的最终误差包围仍使用同一个python-flint/Arb底层，不能称跨库独立认证。独立性来自有理nullspace/群轨道几何、不同Gamma函数参数/代数式、I七次幂和Fraction整数根路径。没有用两个浮点中心接近替代包含关系。

## 6. 一般L1、全部尺度传递与门槛

锁定的alternative_alpha/transfer.md SHA为 `3d2acf640e1be0ee11d10fd6e29719ddc4291117975407e0175c8142573d5465`。同一G、Lebesgue正unit face、q=10/7和中心共同半径cube完全匹配，所以此完整I可直接输入

\[
 C_5\le(10/3)^{10/7}I/g^{10/7}
\]

的一般L1/全正尺度接口，不需要额外截断/边界误差。审计也以纯有理7次幂核验：配当前已审g=2821/50000000与本显式I cap，所得传递上包围仍大于32。因此新辅助I改善不能升级为C5突破。

## 7. 实际运行与精简白名单

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/audit_I023/check.py --source-root /path/to/audit-cont-repository --model-root /path/to/model-repository --run
```

默认不加--run只锁输入。固定句柄 `audit_I023_fixed_001` 一次成功约0.0193秒，Arb512位，未调用MATLAB、后台生产、旧g或022。精简白名单仅registration.md、check.py、result.json、本报告，未提交。结果保存输入hash、轨道与占用准确账本、I七次幂球及95位有理夹。

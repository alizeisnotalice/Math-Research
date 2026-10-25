# audit_I026：022完整I及024全flag权重预算审计

2026-10-26。两项分别**通过**：

- I022的三线性flag机制给完整积分 `I<13022147/10^12`。
- 在固定022分解、熵AM–GM和普通加权Holder预算族内，全部4560个任意确定权重的唯一最优预算为 `J*=2 Gamma(7/4)m(3/20)^5 Z^(-3/20)`，并给更紧完整界 `I<1620961/125000000000`。

最优I预算约1.2967687959491083e-5，其严格下端高于1.18e-5和1.1e-5，只排除**该上界预算族**达到这两个门槛；不是实际I或J的下界。没有新的C5突破。022/024各五件锁定原件、主clone及其他审计均未修改，未提交。

锁022 manifest `668f4b36fa225c865da727a1bf92f211eccfb06a35aad0ebc44cfa7071531602`，024 manifest `96f65eaa792afc290d3b83aa2284f64170386de045e54f04070b96335f93b764`。固定句柄 `audit_I026_fixed_001` 在运行前登记，输入及源码hash见result.json。

## 1. 共同模型与组合系数1/60

G仍是16条未单位化body方向、权1独立alpha=3/2 stable、Lebesgue测度的完整Green核。原始S_i独立正beta=3/4 stable，A=sum S_i v_i v_i^T，W=16^5，b=3/20。5b=beta的尺度抵消与已审轨道Minkowski接口为

\[
 J=E_R(\det M)^{-b}=16\Gamma(7/4)E(\det A)^{-b},
 \qquad I\le c^{10/7}\frac{\pi^2}{60}J^{10/7}.
\]

对无序pair J，CB中每个五方向基含10个pair，故det A=(1/10)sum_J(prod_J S)B_J。rank3商的非退化完成三元组有3个第一线类选择；继续rank2商的完成对有2个第二线类选择。因此

\[
 \det A=\frac1{60}\sum_{J,X,Y}(\prod_{j\in J}S_j)X_SY_SL_S.
\]

这个1/60=1/(10*3*2)是组合重复数，与面角因子pi²/60是两个不同来源，不能额外消去或重复加入。零完成行列式贡献为0，不妨保留于恒等式；flag只选正完成组。

## 2. 不同nullspace/Gram路径及轨道资格

只对pair(0,1)、(0,3)两个代表计算商几何，不重演作者120个pair的adj/Schur循环。rank3 nullspace的三行记N3，H3=N3N3^T，z_i=N3v_i。两个pair的Gram行列式分别16、24，第一商中非零方向恰为其余14个方向。

任取完成三方向i,j,k，正确非正交坐标公式是

\[
 \det(V_J,v_i,v_j,v_k)^2
 =\frac{\det\operatorname{Gram}(V_J)}{\det H_3}
   \det(z_i,z_j,z_k)^2.
\]

第一线类X内z_i=t_i z_r，因而第一合并系数a_i=t_i²。再用原空间rows(J,r)的rank2 nullspace N2，H2=N2N2^T，对其商坐标z'_i有

\[
 w_{sk}=\frac{\det\operatorname{Gram}(V_J,v_r)}{\det H_2}
 \det(z'_s,z'_k)^2.
\]

第二线类Y内相同比例给b_j。审计在两个代表的所有局部flag内精确验证三维商行列式平方等于a_i b_j w_sk，再各用一个5阶整数Leibniz完成系数核验。Gram分母均严格正，未假定nullspace坐标正交。

用五坐标改号及相邻置换生成元恢复轨道，得到pair轨道40、80，分别有7个双元素第一线类，以及4个双元素加6个单元素线类。10个flag轨道合4560个flag；聚合后的五种概率type、组数和每项C质量与锁022逐项准确相同。总质量62914560=60W。

这条路径只对两个pair代表及其内部第一线类求商，完整组合覆盖来自群轨道，而不是再作120次Gram生产或扫描五方向基。

## 3. 五因子独立性、全scope及退化

每个flag的两个S_j、X线类、Y线类、最终完成支持两两不交，且覆盖全部16个指标。任何第三条原body方向都不在两条方向的span内：否则三条首坐标为1的不同cube顶点仿射共线，与cube顶点的极端性矛盾。因此第一商的零方向恰为pair本身。投影到第二商后，零方向恰是pair并第一线类；再除第二线类后，剩余正完成方向就是末线性因子的支持。

所以F_f是两个单变量及三个正stable线性和的五个独立因子；不同flag一般相关，证明不要求它们独立。真实scope每组都是全16，rho=1，不得套用021的rho<1折扣。对任意flag权的active集合，scope占用仍为sum q=1。

非零商类的Gram与归一化质量均正。平行、反平行和零商方向以零面积处理，不对零完成系数除法。S_i>0几乎处处、负矩m(s)对所有s>0有限，故每个EF_f^-s有限；积分Holder可先截断再单调极限。证明同时控制所有接近退化形状的贡献，不假设统一条件数，也不删去奇异方向。

## 4. 022完整负矩及固定常数

记theta_f=C_f/(60W)，归一化因子F_f在S_i=1时为1。有det A/W=sum theta_f F_f、theta和1。普通AM–GM与Holder给

\[
 J\le K\prod_fR_f^{-\theta_f/5},\qquad
 K=2\Gamma(7/4)m(b)^5,
\]

其中R_f是三个归一stable线性和的尺度积，m(s)=Gamma(s/beta)/(beta Gamma(s))。负矩公式使用原S产品概率，不在倾斜R下强行使用独立性。指数-1/5来自b/beta=1/5，prefactor2来自16 Gamma(7/4)*W^-b。

五类归一概率有两个相同尺度，合并为四类：

|R|总theta|
|---|---:|
|(7+8sqrt2+3sqrt3)/2^(15/4)|1/8|
|3^(1/4)(1+1/sqrt2)|3/16|
|1+sqrt2|5/16|
|2sqrt2|3/8|

本批独立轨道结果恢复这些完整权重，没有删除坏type。

## 5. 024覆盖全部权重的严格最优性

对任意确定q_f>=0、sum q=1，允许零权边界。忽略q=0项仅减小算术和，AM–GM准确给

\[
 \sum\theta_fF_f\ge
 \exp[-\mathrm{KL}(q\Vert\theta)]\prod_fF_f^{q_f}.
\]

取负b次幂、期望，并用指数固定为1/q_f的普通加权Holder，得到预算

\[
 B(q)=K\exp\{b\mathrm{KL}(q\Vert\theta)
                  -(b/\beta)\sum q_f\log R_f\}.
\]

设Z=sum theta_f R_f^(1/beta)，q*_f=theta_f R_f^(1/beta)/Z。所有theta、R正，q*在整个4560维单纯形内部。逐项展开log q*，得到精确恒等式

\[
 \log B(q)=\log K-b\log Z+b\mathrm{KL}(q\Vert q^*).
\]

KL非负，且只在q=q*时为0；边界按0 log0=0连续延拓。因此q*为全单纯形唯一最优权，不仅是几个对称type的最优。type内最优条件权按原theta比例，不能无代价改为均匀flag权。

这项最优性只覆盖固定022 factorization、确定q及上述普通加权Holder指数选择。它不排除随机/形状依赖q、保留若干flag联合负矩、不同Holder架构或新的分解。对全部q都有实际J<=B(q)，所以min B仍是实际J的**上界**；min B的严格下端不是实际J/I的下界。

## 6. k0/k1层为何不新增预算信息

固定一方向后，rank4商的其余15投影若有平行类，则三条原方向线性相关；上节的cube顶点论证排除了这一情况。因此k1新增第一线类全为singleton，实际恒等式为

\[
 \det A=\tfrac15\sum_iS_iB_i,\qquad
 B_i=\tfrac14\sum_{j\ne i}S_jB_{ij}.
\]

继续022的1/3和1/2分组，得到分母120的有序pair。每个原无序pair flag恰复制两次，F、R相同，base theta各减半。uniform预算不变，而且
`sum_copies (theta/2)R^(1/beta)=Z`，全权熵最优也不变；复制权不平衡只会增加conditional KL。

原16方向没有平行类，k0第一次拆分同样全singleton，恰回到k1形式。因而两层“空商合并”不会给新信息。这不限制包含非平行类联合det的其他机制。

## 7. 短代数Z与另一数值路径

由四尺度和准确权重直接求和并化简，得到

\[
 Z=\frac32+\frac{5+3(3/4)^{1/3}}{16}(1+\sqrt2)^{4/3}
 +\frac{(7+8\sqrt2+3\sqrt3)^{4/3}}{256}.
\]

本批从独立重建的五type逐项生成sum theta R^(4/3)，其512位包围含于024的有理Z区间；短代数式与之准确对应。

数值复核采用未移位Mellin
`m(b)=Gamma(1/5)/[(3/4)Gamma(3/20)]`和 `Gamma(7/4)=(3/4)Gamma(3/4)`。令N=Gamma(7/4)^10 m(b)^50，C7=3^10/(2^65*60^7*pi^11)。则

\[
 I_{22}^7=C_7N\prod R_f^{-2\theta_f},\qquad
 I_*^7=C_7N Z^{-3/2},\qquad J_*^{10}=2^{10}N Z^{-3/2}.
\]

这区别作者I^14/J^20路径。512-bit Arb向外包围后，以Fraction整数7次/10次根生成110位有理夹，全部含于022打印球及024的100位有理区间。精确比较确认

\[
 I_{22}<13022147/10^{12},\qquad
 I_*<1620961/125000000000,
\]

且I_*预算的下端严格大于1.18e-5与1.1e-5。最优预算严格改善uniform预算也有解析解释：R不全相同，故theta!=q*，其KL严格正。

Gamma/代数函数的底层Arb库仍共享，不能称跨库独立Gamma认证。不同路径在于nullspace/轨道组合、未移位Gamma参数、Z直和、I^7/J^10代数式和Fraction最终根夹。窄区间只包围显式预算，不是实际I的双侧夹。

## 8. 同G传递、运行及白名单

已锁一般L1/全尺度transfer（SHA `3d2acf640e1be0ee11d10fd6e29719ddc4291117975407e0175c8142573d5465`）使用相同G、Lebesgue正unit face、q=10/7及中心共同半径cube。本批结果可直接输入 `C5<=(10/3)^(10/7)I/g^(10/7)`，不额外加短时/长时或管误差。当前没有新g，也没有认证C5<32。

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/audit_I026/check.py --source-root /path/to/audit-cont-repository --model-root /path/to/model-repository --run
```

默认只锁hash。固定句柄一次成功约0.0726秒，未重放作者120循环或021基线，无优化器、采样、MATLAB及后台生产。精简白名单仅registration.md、check.py、result.json、本报告，未提交。

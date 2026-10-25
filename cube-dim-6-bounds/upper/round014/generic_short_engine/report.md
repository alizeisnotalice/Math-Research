# C014：generic短势K_T、六分量梯度与p_T直接编译引擎

结论：一次编译、6组标量锚及唯一旧中心T=.032全群验证均通过。K与旧exact G−A007 H32证书严格吻合；x₂梯度与C013 exact gradG−B010已存H32 AD梯度严格吻合；p的独立有理DD5/复乘积锚通过；Euler身份在native及存档后核均含0。

K、全部6梯度和p_T一起外墙3.882001459秒，内核3.86484秒，**高于预告的约1.5秒**。本批没有继续性能优化、必要bugfix或第二次full运行。没有新空间点、其他T、MC、网格、MATLAB、扩招或push，没有新的数学C6界。

C013报告/manifest已先冻结才开始本目录。014只写本目录，010/013以及其它代理证据树均只读；当前没有活跃计算。以下是新实现的软件验证，不能登记为全域g证书。

## 1. 模型、函数及合法域

模型为六维body32 Cauchy frame，准确21-cell fan及完整符号/置换群沿用旧证据。令b_j=x·r_j/q_j、α=det/∏q_j，T>0。对于每个非零实b，使用正确branch

\[
 f_T(b)=\arctan(T/b)={\pi\over2}\operatorname{sgn}b-\arctan(b/T).
\]

因此

\[
 K_T(x)=(2\pi)^{-6}\sum_R\alpha_R D_5(f_T;b_0,\ldots,b_5).
\]

引擎直接计算此有限和；没有分别调用G、H后再相减。负b直接采用principal atan(T/b)，没有将π/2一律取正。

合法域要求**全群全部Fourier节点准确非零**。任何一个节点整数分子为0，整个point立即以exit3与unsupported_zero_node终止，不输出部分K/梯度/p。它不skip零节点、不设置sign(0)=0、不处理零节点principal parts。即使p_T本身在零节点处有意义，本generic组合引擎也一律拒绝整点。

非零重复节点保留Hermite扩展。节点符号在某个共同开邻域固定，f_T在各非零节点局部光滑；正阶jet为负的atan(b/T) jet。完整fan为有限和，因此在此合法域可逐项求节点及物理坐标导数。

输入范围与C013一致并显式检查：正int64 T分子/分母与点分母，点分子为非LLONG_MIN的int64；fan q≤10^6、ray分量绝对值≤10^6、det≤10^9。点积和den·q用__int128计算后检查int64可表示；超界exit2，不round。新T48超int64分母点未触碰。本批唯一时间为4/125。

## 2. barycentric偏导中的必需对角项

六节点互异时记

\[
 w_i=\left[\prod_{j\ne i}(b_i-b_j)\right]^{-1},\qquad c_i=f_T(b_i)w_i.
\]

值为Σc_i，节点偏导为

\[
 \partial_{b_l}D_5
 =f'_T(b_l)w_l-\sum_{i\ne l}{c_l+c_i\over b_l-b_i},\qquad
 f'_T(b)=-{T\over T^2+b^2}.
\]

这是在C013的局部sign常数情形上新增的关键项f'_T(b_l)w_l；本实现没有误用sign导数为0的简化。即使六节点全同号，f_T也不局部恒定，**没有全同号直接返回零的分支**。

有非零重复时，值用Hermite Newton表，partial_j使用copied-node D6。缓存real jets由T/(T²+(b+h)²)的级数递推，再对atan jets取负；保留到6阶，覆盖七节点DD6的全重合情形。标量锚包含六重正节点，检查了这一最高阶路径。

原坐标链保持C013已验证的映射：y_j=s_j x_(perm[j])，∂_(x_perm[j])b_m=s_j r_(m,j)/q_m。梯度对物理x求导，不再除点分母。节点偏导保留原ray标签，排序只用于Hermite表。

不采用canonical整体反号。全局负号配对时f_T奇，D5为偶；节点偏导反号，chain Jacobian也反号，原坐标梯度不变。因此720×32×21次计算用精确权2代表967680个逻辑群项。p的复乘积在全局反号下取共轭，实部同样不变。

## 3. p_T的独立复乘积与代数恒等式

\[
 \partial_T f_T(b)={b\over T^2+b^2},\qquad
 D_5\left({b\over T^2+b^2};b_0,\ldots,b_5\right)
 =\Re\prod_{j=0}^5(T-ib_j)^{-1}.
\]

后一式可由有理函数的DD乘积公式得到：D5[(T−ib)^(−1)]=i5∏(T−ib_j)^(−1)，再取虚部，Im(i·product)=Re(product)。重复节点由同一有理恒等式的Hermite延拓包含。

native只用双Arb实/虚数乘法：每个因子为(T+i b)/(T²+b²)，逐因子更新实部、虚部，最后取实部。p_T按同一α和(2π)^−6归一化。未新增acb头，未把binary64复数当作证书。

独立标量参考对T=4/125和准确有理节点，以fmpq计算函数b/(T²+b²)的Hermite DD5，并另作准确复有理乘积。五个非零案例中两个有理结果逐项相等，且包含在native p-cone的完整dyadic区间内。第六个零案例只测试拒绝，未用它绕过generic域。

## 4. 六组标量锚及独立参考

锚与013同样包括：未排序的混号互异节点、混号重复节点、六重正节点、混号三重/二重、重复有理节点，以及一组含零节点的拒绝测试；全批时间固定4/125。

φ的640bit参考采用π/2 sign−atan(b/T)表达及独立Acb正阶jet公式，和native直接atan(T/b)+real jet递推不同。每个非零案例的值及全部6个copied-node DD6偏导，由native512完整端点包含对应640bit完整向外dyadic区间。p另用上一节的准确有理恒等式复核。所有记录在anchors_native.json、anchor_validation.json；没有短radius反造区间。

所有节点的非零检查在任何缓存、同号判断或相减之前执行；本轮不存在sign/endpoint不一致，也未扩锚或放宽精度后重试。

## 5. 唯一oldc、T=.032的全群交叉

旧中心共同分母1000000007，分子

`(1000000007,998836847,892738545,807464540,467291939,96968636)`。

运行前准确比对四处中心/分母：旧exact G、A007、C013 exact gradG、B010 H32梯度；并核fan、旧G证书SHA等来源关系。所有来源一致才启动native唯一full执行。

值参考直接读取旧fmpq G的准确Q/π5和A007 H32完整球，构成1024bit外包G−H32。梯度参考读取C013新独立fmpq gradG_x2准确Q₂/π5，减B010**已存**的旧c H32 AD完整球。没有重跑G、H或H梯度。

native512的K区间及x₂梯度区间分别**严格位于**这两个旧证书差区间内部。逻辑群项967680、配对计算483840，本旧generic点zero/repeated群项为0；特殊非零重合由标量锚验证。

便于阅读的结果为

\[
 K_{.032}(c)\approx9.4191264414535010642244291\cdot10^{-6},
 \qquad p_{.032}(c)\approx0.0004132669765341514992898195764.
\]

六梯度展示为

\[
 (-2.27796954,-2.26334181,-1.18114900,-0.49999370,-0.33163748,+1.23562749)\,10^{-5}.
\]

这些只作阅读展示。K、6梯度、p与Euler均以整数端点乘2幂的完整格式保存在native_result.json。

联合缩放K_(λT)(λx)=λ^(−5)K_T(x)给

\[
 x\cdot\nabla K_T+5K_T+Tp_T=0.
\]

native累积区间含0；Python只读K、六梯度和p这八个native输出及准确x、T重建该区间，仍含0。p_T完整区间的下端严格为正。Euler只是额外身份检查，未替代K/x₂的独立旧证书交叉。

## 6. 本机依赖、构建和实际性能

复用冻结010的13个匹配头archive和现有FLINT3.6.0、GMP6.3.0、MPFR4.2.2 wheel runtime，python-flint0.9.0；没有重装、改PATH、改旧目录。只在本014 scratch解包头并编译。C++编译期ulong/fmpz=8字节、arb_struct=48字节，运行期flint_version必须3.6.0。

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/upper/round014/generic_short_engine/build_engine.py
```

build_record.json保存实际clang++ -O3 -std=c++20命令、仅修改本轮binary的@rpath修补、匹配头/runtime和binary SHA。一次编译版本，没有必要bugfix。

现成CLI：

```sh
cube-dim-6-bounds/upper/round014/generic_short_engine/scratch/generic_short full \
  /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel/fan_integer_input.txt \
  4 125 1000000007 \
  1000000007 998836847 892738545 807464540 467291939 96968636
```

这是复现命令，不表示增加第二次full执行。本批实际唯一full exit0，工具chunk457e49；外墙3.882001458987361秒、内核3.86484秒。耗时高于预估1.5秒，已如实停止，没有增加性能测试或新点。verify_saved.py实际PASS（chunk400055），它只重建存档/旧来源区间，不调用kernel。

新工具只支持准确有理point与正有理T的generic非零域，不提供连续箱证书、零节点主部或全域最小值认证。其它T、点或退化条件下的吞吐没有在本批测试。功能验收通过并不等于可直接声明global g已认证。

## 7. 有限批结账与冻结

完整运行数1≤2，唯一点1，T只有4/125；6 scalar cases，0新点、0MC/网格/MATLAB；全部进程退出，无活跃句柄。达到全部登记功能验收后停止，没有使用预留bugfix pass。

来源包括：main旧exact G与准确fan；theory clone A007 H32；冻结C013 exact gradG；upper-astra clone B010 H32 AD；010匹配headers/runtime。输入SHA在validation_results.json，所有必要文件及旧目录冻结检查在manifest.json。native源码SHA为247a1f6c6967ad52f34320b999818214488700651b46cb8e5feb3eac7064953b。

本批是新实现验证，零新数学C6界。013和更早目录保持冻结，014随本报告与manifest交付冻结。

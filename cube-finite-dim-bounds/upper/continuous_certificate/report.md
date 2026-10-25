# 11E 连续 minimum：有限批未闭合，已定位表示障碍

固定 (n=3,\tau=1/3)。本批没有证明 (G\ge13/2500) 于整个 cube，也没有构造反例；当前有效 (U_3=8,U_4=16) 不变。以下局部严格结果以共享的真实 96 锥表示正确为前提，不能替代其独立审核。本目录没有提交或推送。

## 输入与范围

已读独立 clone 的 AGENTS.md、主 011e 登记、minimum snapshot 及 weighted-frame transfer。实际 hash 在 input_hashes.json；登记为 0f5498cbfef2de1e22b0b240ca86b3e18d438f6be9bffc0174ce1f140d74031b；锁源为 ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849。

仅写本目录。未优化 τ、未计算 J、未调用共享 MATLAB、未改主状态/账本。正齐次和 cube 群将全域问题归约到 (x=(1,y,z),0\le z\le y\le1)；此处不假定 edge minimum。一般 L1、共同半径、全部尺度及最终成本接口由主负责。

## 首选表示与实际 4095 区域诊断

令每个正负配对锥的三个独立 ray 为 (r_j)，(k_j=\psi(r_j)>0)，(b_j=r_j\cdot(1,y,z)/k_j)，(W=2/(k_0k_1k_2))。generic 有限表达为

\[
8\pi^3G=\sum_C W\sum_j\frac{\log|b_j|}{\prod_{l\ne j}(b_j-b_l)}.
\]

SymPy 精确归并得到 25 个不同绝对投影，形式为

\[
8\pi^3G=\sum_{a=1}^{25}R_a(y,z)b_a(y,z)^2\log|b_a(y,z)|.
\]

grouped_formula.json 给出完整 25 项及全部人工分母因子，不省略 signed 项。每个零投影的系数含 (b_a^2)，因此在其余分母非零时按 (0^2\log0=0) 连续延拓。程序对跨零区间使用函数全局下界 (-1/(2e))，不删除负项。此事实不处理分母同时为零的交点，也不宣称所有方向奇点可去。

cc11_diagnostic.py 用 160-bit Arb、精确有理常数打印器及 dyadic 矩形，评估 4095 个区域；闭矩形与三角域相交，保留边界。结果：

|结果|评估区域数|
|---|---:|
|人工重投影分母区间包含零|2154|
|区间乘积再反演产生非有限 enclosure|499|
|有限 enclosure 但下端未过目标|1442|
|认证区域|0|

第二类属于乘积区间依赖性，不能理解为真实函数发散。第三类不能理解为反例。所有 3934 个最终未评估 leaf 矩形的准确有理端点保留在 diagnostic.json；这些 leaf 的状态没有被伪称为已评估。矩形可能重叠于边界，但没有遗漏三角域。零投影分支的保守下界修正后，在完全相同的确定性划分上重跑，计数未变；没有扩大分区预算。

一个明确障碍实例为 (y\in[3/64,1/16],z\in[1/32,3/64])：各单独线性因子可能避零，但相乘后 Arb 球扩张导致逆值非有限。更根本的是真正跨越 (8y\pm8z\pm1=0) 等人工重投影线的区域，单项有理系数本身不能直接求区间。增加精度不能消除这类参数宽度。

## 唯一表示修订：正时间截断

具体失败原因为 log 项之间的重投影相消不能被普通区间除法保存。改用

\[
H_\epsilon(x)=\int_\epsilon^\infty p_t(x)\,dt\le G(x),\qquad\epsilon=1/50.
\]

此不等式来自真实热核非负性，无须估计被丢弃短时间质量。设 (u=1/t)，得到有限有理积分

\[
H_\epsilon=\frac1{8\pi^3}\sum_CW\int_0^{1/\epsilon}\frac u2\left[\prod_j(1-ib_ju)^{-1}+\prod_j(1+ib_ju)^{-1}\right]du.
\]

对实 (u,b_j)，每个因子模平方为 (1+b_j^2u^2\ge1)。因而所有 zero/repeated 投影及原 Green 的方向奇线都可直接评价。程序不把对复数取 real 当作 holomorphic 操作：积分回调保留正负两项的有理和；它是 (u) 的亚纯函数，遇极点 Arb 自动非有限，符合 acb.integral 的解析回调要求。

cc11_regularized.py 实际用 Arb 完成 5 个点的有限积分，涵盖 edge 候选、face 方向线、body 方向线、重复投影和轴点；regularized_checks.json 记录 enclosure。点值不等于全域证明。

在候选附近的闭矩形

\[
|y-1|\le1/16384,\quad |z-2272387/5000000|\le1/16384
\]

实际得到

\[
0.00521448368971308 < H_{1/50}(1,y,z)
 <0.00527176675962013.
\]

因此这个矩形与目标三角域的交集已有局部严格下界。半宽 1/8192 的相同中心矩形下端只有 0.00518574937745，未过目标；半宽 1/4096 下端约 0.00512799848648。这是局部宽度诊断，不继续机械扩网格。详细 Arb 端点及调用次数在 selfcheck.json。

同一 (H_\epsilon) 的快速点值公式为

\[
H_\epsilon=\frac1{8\pi^3}\sum_C W f_\epsilon[b_0,b_1,b_2],\qquad
f_\epsilon(b)=\tfrac12\log(\epsilon^2+b^2).
\]

重复投影用准确 Hermite：(f'=b/(\epsilon^2+b^2))，(f''=(\epsilon^2-b^2)/(\epsilon^2+b^2)^2)，三重节点取 (f''/2)，双节点取 ([f(b)-f(a)-f'(a)(b-a)]/(b-a)^2)。rational-point 实现已对上述 5 个积分 enclosure 全部重叠；最危险样本的值为

\[
H_{1/50}(1,1,0.4544774)
=0.00524312521682739973953047695795320557\ldots.
\]

dd_benchmark.json 记录同一 generic 有理点重复 100 次实际耗时约 0.073 秒。此是点值成本，不是全域认证成本；目前点分类实现只保证精确有理输入，不承诺已解决 Chebyshev 无理节点间的精确重合识别。

## 可交下一批的严格解析误差接口

对每个锥 ray 精确检查 (k_r\ge3(|r_y|+|r_z|))，最小比值恰为 3，等号在 ((0,\pm1,\pm1))。因为 \(\psi\) 在锥上线性，且绝对值满足三角不等式，这推出全空间

\[
\psi(\xi)\ge3(|\xi_y|+|\xi_z|).
\]

于是 (F(y,z)=H_\epsilon(1,y,z)) 在 

\[
\max(|\Im y|,|\Im z|)<3\epsilon
\]

全纯。若两个虚部绝对值均不超过 (s<3\epsilon)，Fourier 绝对收敛直接给出

\[
|F(y,z)|\le H_{\epsilon-s/3}(0)
=\frac{15}{1088\pi^3(\epsilon-s/3)^2}=:M.
\]

这里 (p_1(0)=15/(544\pi^3))，且真实部分任意，不受 patch 中心位置限制。这也严格覆盖原 (G) 奇线附近的截断函数。

推荐下一批只登记固定 16×16 方形划分与三角域相交的 136 patches，每个半宽 (h=1/32)，每坐标 Chebyshev roots 插值 degree 20，ellipse 参数 \(\rho=3\)。ellipse 虚部高 (s=h(\rho-\rho^{-1})/2=1/24)，因此

\[
\epsilon-s/3=11/1800>0,\qquad M<11.906178395825.
\]

从 Chebyshev 系数界 \(|a_k|\le2M\rho^{-k}\) 和 roots 插值的 aliasing，单变量误差不超过 (4M\rho^{-n}/(\rho-1))。简单 DCT cardinal 展开给 Lebesgue 常数粗界 (\Lambda_n\le2n+1)，无需调用精细常数。两个变量的总误差不超过

\[
E_n\le(2n+2)\frac{4M\rho^{-n}}{\rho-1}.
\]

实际 Arb 在 (n=20) 得 (E_{20}<2.868313237155\times10^{-7})。若每节点 enclosure 半径 ≤(10^{-12})，节点误差放大最多 (41^2\cdot10^{-12}=1.681\times10^{-9})。二者之和远小于已见候选点对目标的 (4.31\times10^{-5}) 余量，但后者只作可行性判断，绝不代替 polynomial 全域下界。

共 (136\times21^2=59976) 节点；按上述有理点 benchmark 线性外推点值约 44 秒，未计无理节点包围、confluent fallback、系数转换和 polynomial 认证。下一批必须先解决这些接口并测真实成本，再跑固定 patch 批。推荐以区间 Bernstein 系数或严格多项式分支界验证每个插值多项式下界；仅比较节点仍然无效。对 diagonal patches 可以保守验证整个方形。

这只是建议的新有限批，尚未执行、未登记为已授权的新研究。不得把此处可证明的余项界和正数值 margin 拼接成已证 minimum。

## 最小复核与缺口

使用现有环境 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python`，依次运行本目录 `cc11_diagnostic.py`、`cc11_regularized.py`、`cc11_selfcheck.py`。需要 SymPy 与 python-flint；锁源路径与 SHA 在代码中明确检查。

自检包含 4 个非退化 rational points 的聚合表达和锁源准确点表达 Arb overlap，以及上述局部正时间盒。自检没有独立验证 96 锥到真实 Fourier Green 的全部理论；这是共享表示依赖。真正剩余缺口为：完整三角域的截断核多项式下界证书，以及相应独立审核。未产生全域 checker/certificate，因此本批不请求原 independent_audit 执行全域验收。

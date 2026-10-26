# upper_curvature128：锚点时间质量曲率的全空间 Harnack 接口

结论：新公式成立。固定五维 Lebesgue、16 条未单位化方向
\(V=\{(1,\pm1,\pm1,\pm1,\pm1)\}\)、等权共同 Cauchy 时间、逆 Fourier 因子 \((2\pi)^{-5}\)。对所有 \(\epsilon>0\)、\(c,u\in\mathbb R^5\)，有
\[
\boxed{\log H_\epsilon(c+u)\ge\log H_\epsilon(c)
 +\nabla\log H_\epsilon(c)\cdot u
 -\frac{|u|^2}{16}\frac{K_\epsilon(c)}{H_\epsilon(c)},\qquad
K_\epsilon(c)=\int_\epsilon^\infty t^{-2}p_t(c)\,dt.}
\tag{128.1}
\]
这里的系数只依赖锚点；它是直接由固定时间 Harnack 与 Jensen 得到的全空间位移下界。不能把锚点 \(K/H\) 直接当成整条线段上的统一局部 Hessian 常数。

只读来源资格已实际通过：105 的全部 4973 个残余均匹配主已接纳的 096/098 packet 中 \(\epsilon=1/64\)、目标 \(21/200000\) 的完整闭叶。107 的 4096 个与 123 的 6237 个 H128 锚点全部精确位于这些旧 H64 闭叶内。因此 123 的全部 5665 个直接失败父锚点、11330 个待证 child 均可使用
\[
L_{64}(c)=21/200000,\qquad
\boxed{q_{128,U}(c)=1024-768\frac{21/200000}{U_{128}(c)}}.
\tag{128.2}
\]
同一父 H128 的原 Delta 与四个物理 log-gradient 区间原样保留。5659 个候选的父锚点位于 child 外；全空间公式允许这种情形。

本报告未运行这些真实端点的 q 或角点资格计算，未更新任何全域 g、I、C5 或 S。新 H128 raw、保存 Delta/梯度的独立数值审计与主接纳仍是共享依赖；本分支只核其来源、类型和精确归属。主代理保留最终证明与发布责任。

## 1. 075 实际表示及固定时间常数

实际读取的 075 原件 SHA 为 `3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225`。写
\[
B=[v_1\ \cdots\ v_{16}],\quad BB^T=16I_5,\quad
A=B^T/16,\quad K=\ker B.
\]
075 已证明每条纤维均有限且连续，其处处规范正密度为
\[
p_t(x)=\frac1{1024}\int_Kf_t(Ax+w)\,dw,
\qquad f_t(z)=\prod_{j=1}^{16}s_t(z_j),\quad
s_t(z)=\frac{t}{\pi(t^2+z^2)}.
\tag{128.3}
\]
其中 K 上为 Euclidean Hausdorff/Lebesgue 测度。常数 \(1024=\sqrt{\det(BB^T)}\)，且 \(|Au|^2=|u|^2/16\)。没有把 body 方向单位化，没有改变物理坐标度量。

令 \(\ell_t=\log s_t\)。直接恒等式给
\[
\ell_t'(z)=-\frac{2z}{t^2+z^2},\quad
\ell_t''(z)=\frac{2(z^2-t^2)}{(t^2+z^2)^2},\quad
\ell_t''(z)+\frac2{t^2}
=\frac{2z^2(z^2+3t^2)}{t^2(t^2+z^2)^2}\ge0.
\tag{128.4}
\]
固定 \(t,c,u\)，取 \(a=Au\)，并在同一 K 上平移。以密度 \(f_t(Ac+w)/(1024p_t(c))\) 定义纤维概率，设
\[
\sigma=\sum_ja_j\ell_t'(z_j),\qquad
\tau=\sum_ja_j^2\ell_t''(z_j).
\]
075 的有界 score 支配允许两次微分，故
\[
D_u^2\log p_t(c)=\mathbb E\tau+\operatorname{Var}(\sigma)
\ge-\frac{2|a|^2}{t^2}
=-\frac{|u|^2}{8t^2}.
\tag{128.5}
\]
方差非负，不需要任何 Cauchy 参数矩，也不要求纤维条件坐标独立。这是对每个固定 t、所有实 c 的物理坐标 Hessian 下界。

对 \(s\mapsto\log p_t(c+su)\) 在 \([0,1]\) 积分二阶导，得到
\[
p_t(c+u)\ge p_t(c)
\exp\left(\nabla\log p_t(c)\cdot u-\frac{|u|^2}{16t^2}\right).
\tag{128.6}
\]
Hessian 常数 \(1/(8t^2)\) 与 Taylor 罚项 \(1/(16t^2)\) 的因子 2 保持正确。

## 2. 完整时间积分、导数与 Jensen

075 选取 det-16 的五方向基，证明对所有实 x
\[
0<p_t(x)\le(16\pi^5t^5)^{-1}.
\]
因 \(\epsilon>0\)，所以 \(0<H_\epsilon(c)<\infty\)，并且
\[
0<K_\epsilon(c)\le\epsilon^{-2}H_\epsilon(c)<\infty.
\tag{128.7}
\]
075 的同纤维有界 score 支配亦已覆盖完整 \([\epsilon,\infty)\) 与任意有限空间邻域。特别
\[
\nabla H_\epsilon(c)=\int_\epsilon^\infty\nabla p_t(c)\,dt.
\tag{128.8}
\]
例如 075 还给 \(\|\nabla\log p_t\|_1\le1/t\)，所以一阶导绝对值可被常数乘 \(t^{-6}\) 控制；\(t\ge\epsilon\) 排除了短时端点问题。没有移动积分端点。

定义时间概率
\[
d\nu_c(t)=\frac{p_t(c)}{H_\epsilon(c)}\,dt,\qquad t\in[\epsilon,\infty).
\]
对固定 u，随机量
\[
X_t=\nabla\log p_t(c)\cdot u-|u|^2/(16t^2)
\]
绝对值有界，故 Jensen 合法。由 (128.6)、正积分及凸函数 exp 的 Jensen 方向，
\[
\frac{H_\epsilon(c+u)}{H_\epsilon(c)}
\ge\int e^{X_t}\,d\nu_c(t)
\ge\exp\left(\int X_t\,d\nu_c(t)\right).
\]
同时由 (128.8)
\[
\int\nabla\log p_t(c)\,d\nu_c(t)
=\frac{\nabla H_\epsilon(c)}{H_\epsilon(c)}.
\]
取 log 即得 (128.1)。右侧的指数积分也有限，因为它被有限的 \(H_\epsilon(c+u)/H_\epsilon(c)\) 控制。整个推导无须限制 u 长度，跨越 S1、S2、符号面或重复投影都有效；这些集合只影响某些锥核表达式与完整 Green 的短时无限值，正时间密度与截时势均处处有限、正、光滑。

## 3. 有限截时质量给 q 的准确上界

给任意有限严格递增正有理截时
\(\epsilon_0<\epsilon_1<\cdots<\epsilon_m\)，置
\(d_j=\epsilon_{j-1}^{-2}-\epsilon_j^{-2}>0\)。对每个 \(t\ge\epsilon_0\)
\[
t^{-2}\le\epsilon_0^{-2}-\sum_{j=1}^md_j1_{[\epsilon_j,\infty)}(t).
\]
右侧在每个时间带上是左端点的倒平方，最后无穷尾仍为正 \(\epsilon_m^{-2}\)。乘非负 \(p_t(c)\) 后积分得到
\[
K_{\epsilon_0}(c)\le\epsilon_0^{-2}H_{\epsilon_0}(c)
-\sum_{j=1}^md_jH_{\epsilon_j}(c).
\tag{128.9}
\]
因此若实际可靠证书给同一锚点、同模型的
\(H_{\epsilon_0}(c)\le U_0\)、\(H_{\epsilon_j}(c)\ge L_j\ge0\)、\(U_0>0\)，则
\[
q(c):=\frac{K_{\epsilon_0}(c)}{16H_{\epsilon_0}(c)}
\le q_U(c):=\frac1{16}\left(\epsilon_0^{-2}-\frac{\sum_jd_jL_j}{U_0}\right).
\tag{128.10}
\]
方向来自 \(H_{\epsilon_j}/H_{\epsilon_0}\ge L_j/U_0\)，不能用 \(L_j/L_0\) 替代。可靠输入还自动满足 \(0\le L_j\le H_{\epsilon_j}\le H_{\epsilon_0}\le U_0\)，故
\(q_U\ge\epsilon_m^{-2}/16>0\)。若出现 \(L_j>U_0\)、非正 U、错点、错模型或非严格正截时，应拒绝该输入并核来源；不应把负 q 用进凹角点判据。可将非正的保守质量下界仅替换成 0，得到旧弱界。

单个 \(\epsilon_1=2\epsilon_0\) 给
\[
q_U=\frac1{16\epsilon_0^2}\left(1-\frac34\frac{L_{2\epsilon_0}}{U_{\epsilon_0}}\right).
\]
取 \(\epsilon_0=1/128\) 就是 (128.2)，并且 \(256\le q_U<1024\)（当 \(L_{64}>0\)）。H64 的旧闭叶结论是在新锚点 c 处的质量下界，完全不消费 H64 的梯度。新线性项与 Delta 均来自同一 H128 jet。

若将来需要叶外质量下界，可用任一旧 H64 jet 的 075 全空间 Harnack 先认证 \(H_{64}(c)\) 的下界。来源锚点不必等于 c；必须有向包围该 Harnack 右侧及 exp，并绑定同一物理点 c。本轮不计算这种新 exp，也不需要它，因为旧完整闭叶已给统一 g105 质量下界。

## 4. 与保存 H128 jet 的最小接口

给实际锚点 c 与原 H128 可靠区间 \(H(c)\in[L,U]\)、四个物理导数，保留原证书中的
\[
\Delta\le\log(L/g),\qquad
\gamma_j\in[\gamma_{j,L},\gamma_{j,U}].
\]
只用新的同点质量证据构造非负 \(q_U\)。对真实闭盒的物理偏移 \(u_1=0\)，安全下界为
\[
M(u)=\Delta+\sum_{j=2}^5\min(\gamma_{j,L}u_j,\gamma_{j,U}u_j)-q_U\sum_{j=2}^5u_j^2.
\tag{128.11}
\]
它全局凹，所以同一锚点全部 16 个闭角点的 \(M>0\) 足以证明整盒 \(H_{128}>g\)。锚点在盒外仍成立。不可逐角点切换锚点；不可把别的 epsilon 的梯度拼进来。旧 Delta 可以为负，未通过只表示该证书未闭合。

相对原 Q1024，任意偏移的安全下界准确增加
\[
(1024-q_U)|u|^2=768\frac{g105}{U_{128}(c)}|u|^2\ge0.
\]
这里只作公式说明，没有对真实候选运行旧或新余量。

`evaluate_zeroH.py` 是待主审的具体有限 consumer：从已冻结 123 final state 读取保存父证书和 raw，核 typed hash/epsilon/point，原样消费 Delta 与梯度，计算 (128.2)，检查每个 child 的 16 闭角点。总对象恰为 11330 child、最多 181280 个角点，单进程、无 pool、无新 split、不导入核/Arb/log/exp。它要求显式传入已审准备文件 SHA 及一个全新输出目录。本分支尚未运行该脚本。

## 5. 实际只读资格及依赖

`qualify_sources.py` 实际 PID **36446**、exit **0**、约 **1.996 秒**。输入先逐文件 SHA 校验，后检查：

- 105 的 4973 pending 与 source mapper 的 source leaf/certificate/record hash、精确盒边界逐字段一致；每个节点实际存在于锁定 packet 的 H64 叶，epsilon=1/64、Q256、g105 packet/已有 checker/主接纳链一致。未重做旧 log 或旧 H64 角点。
- 全 4096 个 107 与全 6237 个 123 raw manifest 锚点均在自己的精确盒内，且其 node 路径恰有一个上述旧 H64 g105 叶祖先；锚点也经 Fraction 逐坐标精确包含测试。
- 123 的 6237 任务、两个 raw manifest、final state 的全部同 epsilon cache record 逐对象一致。新 H128 的模型、kernel SHA、target、完整点、复合 cache key 匹配。
- 5665 个失败 parent 的保存 H128 证书与同点 raw 绑定。每次实际闭二分轴与 path 相符，两个 child 恰为 parent path 加 0/1；合为最终全部 11330 queue，不漏不重。全部 child 继承相同旧 H64 g105 祖先。

这些检查证明“可使用哪条旧质量下界”和“准备消费哪条保存 H128 jet”，不重新证明 raw 数值精度或保存 Delta 的有向 log 资格。原 096/098 全 g105 主接纳与 123 独立数值审计分别保留责任。

外置准备包：

`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_curvature128/prepared_candidates.jsonl`

SHA `bef5f9d45108920b7385bfca31a25a97a3ba30aa5ea4ce266cffdd2b47c31112`。只存 11330 个候选的节点、直接父 split index、typed point 与 raw/certificate hash、旧 H64 叶引用；不复制 761 MB seed，不复制任何 raw。资格结果及完整源锁在同目录 `source_qualification.json` 和 `source_qualification_registration.json`。

关键不可变输入：

|来源|SHA-256|
|---|---|
|075 report|`3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225`|
|099 report|`96e874223e5fbc78f9688f7eec574cdcdde184344b1fa809114bc87553d48e85`|
|029 report|`298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81`|
|096/098 packet manifest|`c53cef7f2626e68042b73d4bb5e80387841eb810a8fbcc1c511f9ce1d8b4780b`|
|096/098 existing checker|`5ad237758ac62c860e43e1977056b24e0ed9b366854dbd275dfb5aa3139c701d`|
|098 full state（继承，不本轮重载）|`73bc873f23f7b68f6f70132e590d28ebbe4cd619012ac9c81980f9fa500314af`|
|105 pending|`2e4639a521a88aec1b49074f1e38bf5eecb9b853d2321a4b8a2685420b8285d5`|
|105 mapper|`c1f95e6cdb69f096fcbc03c5dcceca15795f5c96a225c8b34d66498cb23b63c0`|
|123 planned tasks|`0bf33d3fd8a4ab228252062baa2d3469cb80abb3ddb46c591392a4c63d82db72`|
|123 batch1 raw manifest|`40342fbe1d0bf596dc3a31aac741237ad926b17539c8f96a5988f1e3e331b656`|
|123 batch2 raw manifest|`faf4de812b216875bebbddc802843ca8f28a9523b0f5e93a108a563374ff04e2`|
|123 final state|`6c1d4c6e075053c427d797798b5ccfc20e62e46a46151511e01a330b38213301`|
|123 seed（继承，不本轮重载）|`52fe7f9f4ecd8d4c16c1c4d48423f78bc09821770fb8e3c5d8899a9dfa4f691b`|
|H128 numerical kernel（只核引用，不导入）|`128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0`|

## 6. 本机新符号恒等式与有限有理 toy

按用户本机数学工具偏好，实际运行已有 Python/SymPy；不使用共享 MATLAB。`symbolic_check.py` PID **37223**、exit **0**，只核三个新恒等式：有限时间带望远镜、epsilon128 的 q 表达式、角点安全余量增加量。

另做一个纯人工接口 toy：声明 \(L64=105\mu\)、\(U128=140\mu\)、保存 \(\Delta=1/8\)、零梯度、偏移盒 \([-1/128,1/128]^4\)。准确计算 q=448，旧 Q1024 余量 \(-1/8\)，新余量 \(1/64\)。假设 \(H128/g\ge8/7\) 可支持该 Delta，因为 \(\log(8/7)\ge2/15>1/8\)；后者由 \(\log x\ge2(x-1)/(x+1)\)（x≥1）的导数差非负直接证明。没有调用数值 log。

此 toy 只测试接口与方向，未声称存在这样的真实 H 锚点，不能算实际 g 证据。所有真实候选角点测试仍待主审授权。

## 7. 原一般 L1、全尺度接口保持

新证书若逐叶闭合，仍仅先给相应 \(H_{\epsilon_i}\ge g\)，再以正积分 \(G\ge H_{\epsilon_i}\) 得完整 G 的面下界。包含所有奇异点、闭边界；不对分片 epsilon 求导，不对 H 错用 −4 齐次性。之后仅对 G 使用 \(G(rx)=r^{-4}G(x)\)，再进入已审 029 的同模型 essential cube 下界与完整 I 接口。

所以一般 L1、全部正共同半径的中心轴平行方体、同一接收集与成本
\(C_5\le5^{5/4}I/g^{5/4}\) 均未改变。本分支没有新增完整 g，正式 C5 上界仍依主账本为 32。

## 8. 状态与复现

独立 clone `/Users/zhengzhihao/Desktop/T/Math-Research-cube5-curvature128-sol`，基于 main `dd413e1cf61a0e376d27ba45fdbaf7ac27fc4a0c`。仅新增自身 `upper_curvature128/` 和外置 scratch；未提交、未推送，未改任何他人文件。0 新 H/梯度/核/log/exp，0 真实候选 q/角点求值，未续派任何代理。

本轮资格及符号结果已冻结，重放请复制脚本到新输出位置，避免覆盖现存结果。待主读完整证明和 consumer 后，唯一具体下一对象是已锁准备包的 **11330 个直接父 jet 的单锚点整盒新曲率资格**；不新增 H，不扩 split，不自动追加任何网格。

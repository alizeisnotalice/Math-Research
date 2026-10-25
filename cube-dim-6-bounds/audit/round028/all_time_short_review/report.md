# F028：body32 全时间短势及正时间混合独立审计

**结论：PASS。** 对固定原 body32、固定点

    x0=(50,50,50,45,27,0)/50，c=1/4100，
    ell=3141592653/10^9，u=3141592654/10^9，

独立审计证明所有 T>0 都有 K_T(x0)/T<c(ell/π)^6<c。任意正时间概率混合的全 cube essential 下端不超过 c，故这一核支配方法的成本至少 1025/16=64.0625。这是模型家族障碍；不是 C6 的下界或新上界，不排除其他 frame、其他核或其他方法。原 C6 上界 64 不变。

本目录仅为 F 独立审计。B023/B024/B028、C010/C026 和旧 F023 输入均只读；没有更改作者原件。

## 1. 确切执行范围和修复记录

实际数学审计于 UTC 2026-10-25 20:05:30.713841 开始，20 分钟总预算保持原起点。独立符号 saved-data 检查最初在首根顺序守卫失败：oldz 初始化为 0，错误守卫写成 0<oldz<a。失败前已通过 417 局部主部及 P 恒等式，尚未启动任何 native。该错误是自写 checker 的初始化条件，不是数据或数学身份失败。

原 `verify_roots.py`、`locked_inputs.json` 和 `initial_failure_receipt.json` 均保留。父任务明确允许修复；`verify_roots_repaired.py` 仅将此守卫改为 0<=oldz<a，另读 `repair_lock.json` 并继承原开始时间。第二次 saved-data 检查完整通过，PID40844；该修复不增加任何 H/PF/G 核调用。没有重置 20 分钟预算。

随后 `run_h_endpoints.py` 一次顺序运行恰两次 C010 native H512，各仅在锁定 b 上执行。两次 exit0，用时 1.1380031109 与 1.1358542442 秒；全部数学检查与两 H 结束时距初始审计 60.056947 秒。没有 native 失败、重试或更换参数。最终独立纯读回 PID41035 校验 98 个源哈希及两份准确付款，未增加核调用。

总新工作是两次完整 H native 计算；不是零核审计。新 PF full sum=0，新 G full=0，Sturm=0，新根隔离=0，MATLAB=0。未重跑 C026 的 T=1/25 或旧基线。

## 2. 有理函数、417 次分母及完整根合同

独立 `verify_roots_repaired.py` 从 B023 的 417 个实主部槽和 B024 保存的 N 出发，在每个 z=−s 处直接乘其余线性因子的短 Taylor jet，再作局部除法。全部 417 系数准确一致，359 个最高局部主部非零。N 为 414 次，D 为 417 次，D 由正 s 的线性因子组成，故 D(z)>0 (z>0)，且 N/D 没有可消极点。这里没有重合成全局大分子，也没有展开完整 D。

对每个极点及其重数，另核 N−aZ 的导数全为零；其 417 次首项准确为 64ell^6 c。因此

    N−aZ=64ell^6 cD，a>0，P=N−64ell^6cD=aZ。

这独立确定比较函数导数的符号，不依赖作者给出的 sign 摘要。保存的 squarefree 因子准确重构 Z，gcd(S,S′)=1；S 的非零系数符号变化准确为 2。两个互不交叠正有理 z 框端点异号，给至少两个根；Descartes 给至多两个根，故它们就是全部正根且均简单。准确平方不等式核对 T² 转换，得到：

|根|锁定 T 下端 a|锁定 T 上端 b|
|---|---|---|
|1|62145027/4294967296|15536257/1073741824|
|2|92903691/2147483648|185807383/4294967296|

这只核验已锁框，没有搜索、隔离、Sturm 或近似根数判断。

## 3. 两份独立 G−H 上端及全框付款

G 使用 F023 独立全 tuple Laurent 审计留下的准确 Γ：G=Γ/(128π^5)。本轮仅以新的 Arb768 π 外包把该准确系数转成有理外包，不做完整 G 计算。新的 π 外包严格处于已锁 ell、u 之间，并包含于旧 π512 外包。

H 使用 C010 已锁源码、binary、21-cone fan 和 FLINT3.6 运行库。每次保存完整 dyadic 端点；源／binary／fan／运行库哈希均在执行前核验。两次结果均为：

- 967680 个逻辑 full-group 项、483840 次实际 DD；
- 88320 个含零节点项、201408 个重复节点项；
- precision512、FLINT3.6.0、exact_global_sign_pairing=true。

没有丢弃零节点，非零重复节点走 Hermite/confluent；全局负号对的权重 2 已在 engine 中计入，不能再乘。x0 的坐标稳定子由完整 sign/permutation 枚举自动体现，不能另外乘 12。f_T(b)=atan(b/T) 在 b=0 仍合法，导数 jet 由 T/(T²+b²) 递推。源计算 H 的总 Fourier 因子为 (2π)^−6；它与 F023 的 Γ/(128π^5) 是同一归一化。

此 ABI 审计仅针对锁定 FLINT3.6 构建和本次准确输入。fan 所有 sign/permutation 下原始节点分子上界 600，分母上界 8000，int64 及 __int128 交叉乘积安全；两个 T 的准确整数参数也在 int64 内。没有把此局部输入安全性外推成任意参数引擎保障。

用完整准确端点构造 K_upper(b)=G_upper−H_lower。两次

    c*a−(u/ell)^6 K_upper(b)

严格为正，向下取至 10^−18 格点仍有以下安全付款：

|根|独立安全付款下界|K(b) 仅展示近似|
|---|---|---|
|1|130791170401/62500000000000000|1.436430694326986×10^−6|
|2|232263757011/500000000000000000|1.0087094990298506×10^−5|

精确付款、H/K 的全部分子分母在 `native_endpoint_results.json`；展示小数不参与比较。B028 两份 PF640 球均严格位于独立 G−H 外包内部，作为交叉证据；付款没有使用 PF 球。B028 原 verify_saved.py 曾有 exponent 字符串解析错误，作者保留失败原件并以 check.py 加 int 转换修复；本审计只读已冻结结果，并自行把 exponent 准确转为整数，没有执行作者 writer/checker。

## 4. 全正时间与任意概率混合

严格固定 F(T)=cT−(π/ell)^6 K_T(x0)，所以

    F′(T)=−P(T²)/(64ell^6D(T²))。

P 的初始及末端符号为负，两个正根简单且完整，故 F′ 的三个符号为 +,−,+。K≥0 且随 T 增加，因而上述付款事实上在各自整个 [a,b] 上成立。F(0+)=0，G(x0)<∞ 给 F(+∞)=+∞。完整临界点、两框严格付款及这些端点条件共同推出 F(T)>0 对全部 T>0 成立。没有把两个样本值直接推广为全时间，也没有混用相反势函数的导数／无穷端符号。

`theory.md` 给完整接受论证。`uniform_mixture_proof.md` 独立补足 essential 与任意 μ 的关键：dist∞(x0,S2)≥1/20；在 U=B2(x0,1/100) 上，正参数推前及大参数秩至少 3 给所有 t>0 的共同 p_t 上界 4960(1600√6)^6/π^6。只对独立 triple 扩 basis，以正测度先作变量替换，再由 fixed-t 连续性把 a.e. 上界提升到逐点，不在零概率 fiber 条件化，也不假定任意三列独立。

由此所有 q_T=K_T/T 在 U 连续且受同一常数支配。任意 Borel 概率 μ，即使支持趋零、无界且无逆时间矩，Qμ=∫q_T dμ 仍连续。x0 的每个开邻域与 cube 内部交集均有正六维测度，故点上全时间界确实迫使 essential 下端 ≤1/4100。该模型成本因此 ≥1025/16（保留的 π 缩放实际给严格更强的不等式）。

一般 L1／全正尺度沿用已独审 F009 的正 Markov maximal 常数 1，并用 Cauchy 缩放 r^−6q_T(x/r)=q_(rT)(x)。没有把 sup_T 放入卷积后视为同一个核。分布源分裂若要求有限 δ0 系数需 ∫T^−1dμ<∞；本审计不作该分裂，故没有这一额外假设。

## 5. 冻结证据与范围

主要入口：`root_identity_results.json`、`native_input_manifest.json`、`H_root1.json`、`H_root2.json`、`native_endpoint_results.json`、`final_readback.json`、`audit_receipt.json`。输入锁及两份独立源码均保留。`manifest.json` 列本目录全部文件、字节数和 SHA256；完整外部源哈希保存在 locked/repair/root 记录中。

关键已锁输入：B024 report 662af6b462cc8fb9a8a72a07cc9b1aa4ae243208c726f1773cc9eb3325a36eaf；B028 report 9d6c88ea66a3f69f728abc44a332d819af0eee1ec950158a3afe27f01b3ded60；B028 manifest bce8ac93a01df53f24e4c22745f0427425bb5dacc0814008b35d987b82cdfbc4；F023 exact Γ results a7135d22c023d1dab9251b24f282810a3e67c785e0aa6709326db7ef2933deb9。

C010 binary SHA256 429cdccca45a69c1de6c04725231297e00ee8d3951ccbb1a6141feb934f13bb8；source e857e480efae68ec99fd82f97c64f4756b6bff465476e3a788d9c283f5b751b7；fan 7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812。最终独立 K/payment results SHA256 1cfa61e6e4bd70b35f34a2d69395105060da99dec6cf5a036a6bb8c049e6ff48。

没有遗留运行进程或下一批候选。无新 frame 数值、下界搜索、C6 晋级、commit/push；父任务保留最终证明与集成责任。

# F052：零节点 Green 有限部分独立审计

**PASS，仅限标量理论接口及下述精确参考。** 一般正、有限、张成 R6 的 Cauchy frame，在真实 S2 外 Green 有限且连续，S2 上 Green 无穷。对已证明完整、归一正确的有限 fan，正 cutoff Hε 的完整 Laurent 和中，负幂 −1、−3、−5 在 offS2 点准确抵消，常数 functional 给出

    G(x) = Γ(x)/(128π^5),  x ∉ true S2.

零节点的全部 functional 数据设0，是固定点常数系数提取，不是 sign 在0可微，也不能据此获得梯度或 Hessian。本轮没有新 C6 界、全域 Green floor 或 production engine 放行。

## 独立数学核验

[proof_review.md](proof_review.md) 自包含证明一般 frame 的 rank-greedy 估计。离 S2 至少 δ 的参数和，全部“大参数”列秩必须至少3；选一个独立 triple 扩为六列 basis，用正推前和真实 determinant，三个 t 尾因子与另外三个 t^-1 峰因子相消，得到一致 p_t=O(1)。它无需“任意三列独立”，且不在零概率 fiber 上条件化。固定 t 连续性和开邻域处理 a.e.→pointwise，结合大时间 t^-6 给 offS2 Green 连续。S2 点以包含该点的 pair 扩 basis，在剩余参数正概率事件上给 p_t≥c_x t^-2，因此必须拒绝有限值。

独立核源文的商密度、Gram Jacobian、coercivity常数、rational all-pairs 规范、完整 Fourier fan 的相位和权重、正 cutoff 交换、Hermite jet 展开及归一。所有零/重复节点保留，不能单锥先积分到0。完整 fan 和真 S2 gate 是前提，不能用有限部分数值替代。

B052 的 check_scalars.py 已只读审查，未执行作者12个 scalar，也未把其正式六节点配置当成真实完整 frame。

## 唯一实际独立批

先登记并核定25个输入哈希，等待/准备338.578079秒后启动实际审计时钟。唯一独立批 PID **53463**、tool chunk **c9675b**、exit **0**，计算耗时 **0.230572250秒**，无重试或子进程。配置共4项（上限8）：

1. 一族 symbolic normalized atan jets，k=0,…,5，正非零节点常数、高阶零常数、负节点奇性及1/(128π5)归一。
2. 新六节点配置 [0,0,0,0,0,3]，Fraction负二项卷积与直接Taylor remainder一致；完整通道为 (1/243,−1/81,1/27,0)，符号时间导数与 resolvent 相位一致；附全节点取负不变性。
3. 独立 axis6 frame 的 x=(1,2,3,0,0,0)，完整64个orthants逐项准确求和，Γ=16/15，L1=L3=L5=0。其全部15个真pair残差最小1，直接密度积分给 G=1/(120π5)，与有限部分严格相等。
4. 同一 tiny frame 的 x=(1,1,0,0,0,0)，完整64项给 Γ=−96、L1=64、L3=L5=0。真pair最小残差0；直接 Hε=π^-6[1/ε−3π/4+(3/2)atanε+ε/(2(1+ε²))] 发散。常数有限部分 −3/(4π5) 只能是诊断，明确标 **S2_REJECT_INFINITE_GREEN**。

这两项 tiny analytic reference 是 **2次真实完整参考、合计128个orthant项**，不是“0所有核工作”。生产 body32/152frame fullkernel、G/H/native 引擎、MATLAB、空间网格、旧baseline均0。全部128项及两个15-pair检查保存在 [results.json](results.json)。最小独立实现 [verify_reference.py](verify_reference.py) 不导入或执行作者源；采用 Fraction 负二项卷积，区别于作者的 Hermite 矩阵逆、SymPy导数残数及三角DD运行。文件有一次执行守卫，避免意外重跑。

## 来源与局限

已锁 B052 manifest 175171c3be5757699e33cd3895b29a2a5e6ef4e5edcdb5f4c4cff9f6bc8223d4；theorem f4ead1ef2785595fa71603f6a2ec651f5a31d16dd19296b963bfc4cd584233b3；report 998c14afeb6c94509dcc0f9a019707c00c2fd13b7a483abc59eacdbc60255966。完整绝对路径与25项哈希见 [input_manifest.json](input_manifest.json)，执行前后全部匹配。

唯一前期读取事件是 theorem 尚未创建，后续准备时读到排字未修的瞬态；作者最终封存前已修复，最终源按上述哈希重新绑定。无数学失败、无计算恢复或重试。等待与最终实际耗时分列在 manifest 中。

未来production标量实现仍须逐点绑定真实frame、rank6、完整fan、multiplicity和全部pair严格offS2证书；零/重复节点相等性与符号必须准确。不得复用未经证明的新frame几何，不得跳项或再乘已合并global负号2。**gradient/Hessian推广=false，production新引擎独审=false，053放行=false，新C6界=false。** 主代理保留最终证明及放行责任。

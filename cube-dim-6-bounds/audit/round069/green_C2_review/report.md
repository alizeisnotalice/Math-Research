# B069：真 S3 外 Green C² 与导数有限部契约独立审计

**结论：PASS，限于真 S3 补集的局部完整 Green C²、统一 cutoff 导数收敛，以及由此得到的完整 fan 导数常数项契约。** 理论未发现阻断。7 项新预登记配置实际一次通过；不重跑作者的 384 个 axis tuple-field 项，不重跑 052 baseline，不运行 body32 Green 核或 body32 全 triples。

真 S3 在这里指所有至多三列的**实际线性张成空间**的并，不是“恰好所有不可微点”的判定。作者没有证明 S3 上 Green 的实际 gradient/Hessian 存在或发散，本审计也不作此结论。没有全正面 I、全局 g、C6 或 Cn 新证书，主代理保留最终证明和生产接口放行责任。

## 1. 源文件冻结与实际验收对象

作者目录为 `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/upper/round069/green_C2_off_S3/`。在本轮登记和科学运行之前，实际验证其最终 manifest 的 10 个产物 pin、input manifest 的 12 个预执行 pin 全部一致；科学程序末尾再次核对不变。

|原件|SHA-256|
|---|---|
|manifest.json|527cb3097e55d346e2867c7f58b0ed66f8a8a590ffe26b07979e9f6e111f8ce8|
|theorem.md|a1689203f2c6f36765535247847f6ab29e591e1182c0c3e5d7d1b691fb6b65c7|
|derivative_contract.py|0fa8f5b49b6acb3cd95453dcca6b0f8b2c4a37026dff8fa49c65420dd08e06e8|
|report.md|ae9427c045e28e13cc9aa25ed500429562233339e6a94ffee387bf1207cd979c|
|results.json|ba36ada1067ac2fbd3f4039fcc8ab719a73a579593c1fd22ff2cf2756c261727|

读取了完整 theorem、contract、check_fixed、最终 report/manifest/receipts，以及被引用的 052 scalar theorem；核对了 040 的 S3 正 moment 障碍与 048 的 labelled band-Hessian 公式。本轮不把已有 scalar 有限部直接升级成导数结论。准确源边和作者保存结果检查见 `source_binding.json`。

作者旧结果实际是 10 个新配置、6 个 axis 字段各 64 orthants，共 384 tuple-field 项，PID 73307、保存 tool receipt exit 0；本审计只是读取并绑定这些保存记录，没有把它们算作本轮独立重跑。

## 2. 理论审查：四大参数与二阶正混合控制

令非零列 a_j 张成 R⁶，正谱权重已吸收到列中，ψ(ξ)=Σ|a_j·ξ|。取固定可逆六列基 N，C=||N⁻ᵀ||²_F，c=(1+C)⁻¹，S=Σ||a_j||₁。由 ||ξ||≤√C||Nᵀξ||₂≤√Cψ(ξ) 和 1+C≥√C，可用 ψ≥c||ξ||。

在 dist(x,S3)≥δ 的 fiber 中，以 R=δ/(2S) 划分大参数。若所有大参数列的秩≤3，其贡献落在真 S3；剩余贡献范数≤RS=δ/2，矛盾。因此存在四个独立的大参数。必须按秩选取，不能声称任意四个大参数独立；四元组中的依赖以及重复列均未忽略。

固定独立四元组 J，沿 U_J=span(A_J) 消去这四变量，Jacobian 为 √det(A_JᵀA_J)。余变量在二维正交商上的 exponent 恰是原 ψ 的限制，因为被选列在商空间中消失。它仍满足同一 coercivity c，故 quotient 峰为

\[
 (2\pi)^{-2}2\pi\int_0^\infty r e^{-ctr}\,dr
 =\frac1{2\pi c^2t^2}.
\]

每个被选 Cauchy 密度在 |z|>R 上≤t/(πR²)，四者相乘给 t⁴；正推前/coarea 加独立四元组的 union bound 得 p_t(x)≤M_δt²。无需给同一 fiber 赋唯一 selector，也无需假设事件互斥。邻域 δ′<δ、固定 t 连续性、δ′↑δ 的步骤正确补齐 a.e. 到逐点的升级。

对有理列，共同分母 d_A 使 d_A⁸det(A_JᵀA_J) 为正整数；故每项逆平方根≤d_A⁴。用 π>2 得粗有理包

\[
\bar M_\delta=\frac{\binom m4 d_A^4}{64c^2R^8}.
\]

此处分母 64、分母幂 d_A⁴ 和二维 quotient 维度都匹配。它只保证定性可用，未给实用覆盖规模。

导数部分另取完整正混合：固定 N，把余变量写成 tMw，w 的概率律不依赖 x。对每个六因子 product kernel k，Cauchy log score 满足 |ℓ′|≤1/t、|ℓ″|≤2/t²。因此 b=N⁻¹u、d=N⁻¹v 给

\[
 |D_uk|\le\|b\|_1k/t,\quad
 |D_{uv}k|\le\bigl(\|b\|_1\|d\|_1+2\sum_i|b_id_i|\bigr)k/t^2.
\]

固定 t 的六因子全局峰给概率积分下的导数 dominator，积分后得到相同相对密度界。这一步使用**完整**混合，未对大参数事件、阈值、rank selector 或移动 fiber 求导。因而将它与前面的 p_t≤M_δt² 结合没有遗漏 indicator 导数。

## 3. 局部 C² 与 ε³、ε²、ε 余项

在正距离闭邻域上，小时间的 0、1、2 阶导数分别由常数乘 t²、t、1 支配；大时间由固定六因子峰及 score 分别给 t⁻⁶、t⁻⁷、t⁻⁸。两个时间端点都可积。

先定义这些导数积分的局部一致极限，再由函数及一、二阶导数的一致收敛识别极限为 G 的真实导数，避免预先假设 G 可微。对 H_ε=∫_ε^∞p_tdt 得

\[
0\le G-H_\epsilon\le\bar M_\delta\epsilon^3/3,
\quad |D_uG-D_uH_\epsilon|\le C_1(u)\bar M_\delta\epsilon^2/2,
\quad |D_{uv}G-D_{uv}H_\epsilon|\le C_2(u,v)\bar M_\delta\epsilon.
\]

坐标或混合导数用显式 C₁/C₂，未遗漏换基代价。由此严格闭合 G∈C²(R⁶\S3)。052 的 off-S2 scalar 连续性不足以替代这里的 off-S3 资格。

040 已证原 body32 的正 all-time J₂ 在 S3 上为无穷。这个结论没有被本定理消去；它不等于有符号 Hessian 必发散。本接口对 S3 返回 unsupported 的含义只是其假设失败。

## 4. 先正 ε 求导，再完整 Laurent 消去

对已证明完整且 multiplicity 正确的有限 fan，正 ε 的 H_ε 采用真实 f_ε(b)=atan(b/ε) 的 DD5。Cauchy contour 公式对每个 label 求导：插入该 label 的一个 copy；二次求同一 label 才出现 factor 2。两个不同 label 即使数值相同，也不因此得到 factor 2。线性节点 Jacobian 的物理二阶链式组合在作者公式 (8) 中完整保留了 diagonal 与 mixed 项。

固定节点 tuple 的 Hermite jet 展开中，非零节点值常数为 (π/2)sign(b)，正阶 jet 为 O(ε)；零节点奇数阶 normalized jet 是 (−1)^((j−1)/2)/(jε^j)。因此 DD6 只有 ε⁻¹、ε⁻³、ε⁻⁵ 负渠道，DD7 还可有 ε⁻⁷。常数 functional 的人工零高阶 jets 只是此展开的系数，不能先视为随点可微的函数。

第 3 节已证明完整物理导数有有限极限。对有限完整 fan 逐项展开后，如果某个负渠道非零，最高负幂就阻止该极限存在。因此**完整物理字段**的负渠道必须精确消去，随后才可取 DD5/DD6/DD7 的常数。归一化为 (2π)⁻⁶(π/2)=1/(128π⁵)。这里不需要点的邻域上逐锥一致 Laurent 展开，也未声称每锥各自消去。

据此作者补上了 052 明确缺失的导数论证。单个 formal 常数、小数近零 residues、off-S2 资格或不完整 fan 都不能替代上述顺序与前提。

## 5. 七项实际独立配置

固定 9 列有理 frame 由新的六列有理加权/剪切基，加 2a₀、a₀+a₁、a₁+a₂ 构成，完整原值见 `registration.json`。它与作者 axis fixture 不同；所有计算保持权重，不把非单位列归一成 unit directions。

|配置|独立计算与实际结果|
|---|---|
|1|Cauchy score 的全实变量界化为 (y²−1)²、2y⁴+2y²+4、2y⁴+6y² 三个非负多项式；精确符号通过。|
|2|新加权点 (7/6,32/15,12/5,32/7,80/11,192/13) 的全部 84 triples，以独立 Fraction 消元和增广 Gram 行列式比值求残差，与作者 contract 每行一致；最小距离平方 2101/180。依赖 triples 未跳过。|
|3|点 (1/2,0,0,4/7,5/11,0) 的全部 36 pairs 距离平方最小 25/121，但 84 triples 最小为 0；确为 off-S2/on-S3，tail 与 derivative request 均正确拒绝。|
|4|大参数列顺序 [0,1,6,7,3,4] 的前四列仅 rank 2，greedy 选 [0,1,3,4] 得 rank 4。其实际和距 S3 的平方为 25/121。全部 126 个四元组中 87 个独立；以共同分母 30030 缩放后的 Gram 行列式均满足所需整数下界。|
|5|独立 Fraction 基逆矩阵重构六个 C₁、36 个 C₂ 和精确余项，与 contract 一致；ε 减半分别给 8、4、2 的比。半径到边界、负半径、非正 ε、错 frame hash、缺 fan/multiplicity/ref、非法字段、零列、缺秩和浮点输入均正确拒绝。|
|6|在 labelled 变量中先微分 DD5(z⁸+3z⁷−2z⁵) 的 complete-homogeneous 多项式，再代入 [0,0,1/2,1/2,−2,3]。同 label 二阶为 10，不同但等值 label 二阶为 5；另一等值 label pair 为 6，与 copied DD 的正确因子严格一致。|
|7|新 DD7 tuple [−2,−2,0,0,0,0,1/2,1/2] 用独立 normalized-Hermite 矩阵和真实 atan jets 得 −5639π/500+35/(4ε)−1/(3ε³)+O(ε)。与作者 formal function 一致，非零主部被明确保留；未报告为 Green 导数。|

这些测试检查新 theorem/contract 的错误敏感边界；解析证明不由这些有限例子代替。配置 5 的非空 fan proof reference 只是 API 边界 fixture，返回值明确标注外部 fan 证明未由该函数认证，并不是本轮新 fan 证书。

本轮实际运行命令是 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_independent.py`。PID **74270**，tool chunk **3c56ef**，exit **0**，程序报告 **0.461138292 秒**；180 秒硬限制未触及。7/7 通过、无重试。完整准确输出、源 hash 和脚本 hash 在 `results.json`。脚本拒绝覆盖已有收据；本批已经停止，不应在封存目录重跑。

## 6. 接受范围及未来 oracle 必须保留的前提

可接受的是局部 full-G C² 定理、显式 ε 误差阶、严格 off-真-S3 资格算法和条件导数常数公式。未来生产 oracle 仍必须绑定全部实际列及权重、完整 fan proof/multiplicity、点、所有 triples 的实际秩和精确残差、请求的物理字段，保留零/重复节点及准确符号。

`bind_derivative_request` 只检查这些声明的结构与 frame hash，并明确不认证外部 fan 证明。通过它的布尔参数不是 fan 证明本身。若计算的是有限 band [ε,T] 而不是 H_ε=[ε,∞)，还必须另外支付上端时间尾，不可直接使用本轮小时间差额冒充完整 G jet。

凡 ball/rectangle 跨 S3，不能使用本定理给完整 Green 的 C² 余项；应保留既有正 band 的匹配导数/覆盖路线或另证适用结果。粗 M 常数不保证便宜的全域 cover，本轮没有 body32 full-G Hessian 数值、全正面连续 I 或 C6/Cn 端点。

写入仅限本审计目录，067 未改。0 body32 kernel、0 body32 triples、0 全 fan/64 orthants 求和、0 旧 384/052 重跑、0 grid、0 MATLAB、0 递归代理、0 stage/commit/push。最终判定：**PASS_LOCAL_C2_AND_QUALIFIED_DERIVATIVE_CONTRACT；全域数值目标仍未支付。**

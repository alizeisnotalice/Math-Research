# Round002 Hadamard 矩与完整 I 预算独立审计

2026-10-27。独立 clone：`/Users/zhengzhihao/Desktop/T/Math-Research-cube7-moment-audit-r002`，起点 `0478c861dcbbfbbc8be472f7617ef4aa261f86eb`。本审仅写 `cube-dim-7-bounds/audit/round002/hadamard_moment/`。作者锁定源为独立 moment clone 中 `upper/round002/hadamard_moment/` 的九文件；completeproof.md SHA256 为 `dde15c7d6ece9d6b3f5cab2280afcf3745aebd42eedc95d2fc4a8204b601aa74`。

**结论：接受作者 Hadamard 分块矩上界和完整 I 上预算，未发现实质修复项。** 8 个块确切分割 64 个原方向，块独立性、determinant Minkowski、自负幂 AM–GM 方向、Gamma/Dirichlet 矩和 full-boundary896 因子均成立。独立 768-bit 复评被作者512-bit四个主表达式区间严格包含；精确有理比较确认

`I_H768upper < author512dyadicupper < R_I=178444228577/20000000000000000 < I_Gamma768lower`，

以及 `M_H768upper<M_Gamma768lower`。新的细 dyadic 端只是同一解析表达式的算术锐化。没有全域 g 证书，故没有新 C7 数值上界。

## 1. 范围、独立性与锁源

作者 completeproof.md、mincheck.py、registration、results、sources、receipt、manifest 与 stdout 已完整读取。作者 PASS 不作证明前提；下文重新说明各解析接口。九源在本次检查前后哈希相同，manifest 的每一项字节数/哈希和 receipt 的四个运行哈希均一致。sources 指向的 AGENTS 与 R001 两个解析证明哈希也吻合。本 clone 的全部 83 个 R001 文件在实际检查前后哈希相同，未改动封存成果。

已按根 AGENTS 使用独立 clone；本轮明确 no-commit/no-push 指令优先于一般归档规则。本审未导入或执行作者 mincheck，未读/重跑 Monte Carlo、未运行旧 baseline science、无 Green 节点、指数变体、优化、MATLAB、子代理或提交推送。旧 Gamma 解析表达式仅按规定固定重评，用于上界表达式间的严格比较。

独立有限路径不同：作者由 character-code 生成 coset，本审由三条 dual parity 的 syndrome fibres 重建同一分割；64 个 determinant 由本审手写整数 fraction-free Bareiss 计算，没有调用作者或 SymPy determinant 实现。固定 Gamma 表达式用 768-bit 的 collapsed log-Gamma 路径，与作者512-bit block-product路径不同；I 系数从原始 C、H、896 表达式计算。共享 SymPy、FLINT/Arb 基础库，因此不声称独立算术内核。

## 2. 原过程、S 的尺度与 Gaussian 表示

原方向是 `(1,±1,…,±1)` 的64列、每列 weight1，所有一维 Cauchy 的时间共同为 t，未单位化或删换方向。Y_i iid Gamma(1/2,rate1)、S_i=1/(4Y_i)。作者在 §1 用 J(c)=∫exp(−u²−c²/u²)du 推出 Laplace 变换：c>0 时可在正紧 c 区间内微分，代换 w=c/u 给 J'(c)=−2J(c)，c↓0 用 e^(−u²) 支配，得到 J(c)=sqrtπ exp(−2c)/2。令 c=sqrtλ/2，即 `Eexp(−λS)=exp(−sqrtλ)`。所有 4 和时间尺度都吻合。

给定 Σ=Σ_iS_i v_i v_i^T，Gaussian 协方差2t²Σ，特征函数 exp(−t²ξ^TΣξ) 对 S 平均后为 exp(−tΣ|v_i·ξ|)。它仍是原共同时间的同一半群。条件密度的 `(4πt²)^(-7/2)`、Tonelli 后 t 积分因子 Γ3/(4π^(7/2))=`1/(2π^(7/2))` 以及 Q^−3 均正确。允许 G 在 S3 上取∞不影响非负表示，也不导致删掉奇集邻域。

## 3. 八 cosets、tight-frame 与全部 Cauchy–Binet 系数

作者 code 来自 r=1,…,6 的六个非零 F2³ characters，r=1,2,4 使映射 injective，因此 code 大小8；它在 F2^6 有8 cosets，覆盖所有64列且无重叠。每个 coset 是 canonical retained-Hadamard matrix 的行符号变换与列置换，故每块7×8的 Gram=8I7。

独立重建方式：对六个方向 sign bits t1,…,t6，按

`(t1+t2+t3, t1+t4+t5, t2+t4+t6) mod2`

分组。这三条 parity 正好消去原六 characters，并线性独立，kernel 大小8。实际8个 fibres 每个8项，排序后与作者 coset_masks 完全相同。每块重新建立矩阵，检查 Gram=8I7，所有64方向仅出现一次。不能随便选八列就调用下述矩公式；tightness 与所有七阶 minors 是本审实际核过的条件。

每块删去任一列的七阶 minor 模512，既由 Hadamard8 inverse/cofactor 的解析公式成立，也由本审仅8×8=64次 Bareiss 准确重算通过，signed值逐项与作者相同。平方 Cauchy–Binet 系数因此全部2^18，且 `8·512²=8^7`。每个随机块 Σ_g=Σ_(i∈g)S_i v_i v_i^T 因 tight frame spanning 和每个 S_i>0 而 SPD。八块使用不相交的原独立变量，故独立；行符号变换及列置换保持 determinant law 相同。这里没有把重叠块当独立，也没有声称随机加权块仍等于标量 I。

## 4. 自足 SPD determinant Minkowski 与 AM–GM 方向

取 SPD A、B。令 λ_i>0 为 A^−1/2BA^−1/2 的特征值，ρ=(Πλ_i)^(1/n)。k阶 elementary symmetric sum 的 C(n,k) 个正项中，每个 λ_i 出现 C(n−1,k−1) 次；普通 AM–GM 给每个该和≥C(n,k)ρ^k。k=0,n 为等式，求和得到 det(I+A^−1/2BA^−1/2)≥(1+ρ)^n。乘 detA 再取正 n 次根，正好是

`det(A+B)^(1/n) ≥ detA^(1/n)+detB^(1/n)`。

作者的自足证明无需额外来源或未核 theorem。对八 SPD 块迭代合法。置 W_g=detΣ_g^(1/7)>0，有 detΣ^(1/7)≥Σ_gW_g。函数 x^−1/2 递减，所以

`detΣ^(-1/14) ≤ (Σ_gW_g)^(-1/2) ≤ 8^(-1/2) Π_gW_g^(-1/16)`，

后一式由 `ΣW≥8(ΠW)^(1/8)` 反向负幂得到。转换回 determinant 后指数恰为 −1/112，而非 −1/16 或 −1/56。对不相交八组独立变量取期望，非负 Tonelli/独立性与后面的有限单块矩给

`M≤8^(-1/2)[E detΣ_g^(-1/112)]^8`。

两次不等式的方向都正确；这给实际未知 M 的上界，不能把表达式差 M_Gamma−M_H 当实际 M 的误差或下界。

## 5. 单块负矩 exact Gamma 参数

全部八个 Cauchy–Binet 项的系数为2^18，代 S_i=1/(4Y_i)，七次乘积每项带4^−7，故

`detΣ_g=(2^18/4^7) Σ_j Π_(i≠j)Y_i^−1 =16(Σ_iY_i)/(Π_iY_i)`。

本审精确核 `512²/4^7=16`。这里无漏一个 S 或错误把八变量之和当七变量之和。

T=Σ_(i=1)^8Y_i、P_i=Y_i/T，positive-simplex变换 Jacobian T^7。联合 Gamma 密度分解给 T~Gamma4、P~Dirichlet(1/2)^8 且独立。于是对任意 a>0，

`detΣ_g^(-a)=2^(-4a) T^(7a) Π_iP_i^a`。

7a=8a−a，分母 sumY 被支付一次。T矩与 Dirichlet beta 积分相乘，Γ4 消去，得到

`E detΣ_g^(-a)=2^(-4a) Γ(4+7a)/Γ(4+8a) [Γ(1/2+a)/Γ(1/2)]^8`。

所有形状参数正，矩有限；原 proof 从 Gamma product change of variables 推 beta integral，足够自足。固定 a=1/112 给参数65/16、57/14、57/112。代八次幂及8^−1/2后，等价 collapsed式为

`M_H=2^(-25/14)[Γ(65/16)/Γ(57/14)]^8[Γ(57/112)/Γ(1/2)]^64`。

本审768-bit从这个 log-Gamma式计算，未调用作者单块 product计算路径，亦未启动其他 a 变体。

## 6. 全 boundary896 与完整 I，非单 shape 对称

任意 SPD Σ 的固定 shape 不必有坐标对称性。沿 cube 全 boundary 的 ray坐标 x→rx，Jacobian r^6，先积分 exp(−r²Q(x))，得到

`∫∂K Q^(-7/2)dσ = [2/Γ(7/2)]∫R7exp(−Q)dz = (16π³/15)sqrt(detΣ)`。

q=7/6、3q=7/2 时，shape kernel detΣ^−1/2Q^−3 的 boundary Lq norm 是 `(16π³/15)^(6/7)detΣ^(-1/14)`。先在完整 boundary 上做 Minkowski；限制到 `{N^−1I≤Σ≤NI}` 后可处理有限核，再用 MCT及有限矩，避免循环假定 G∈Lq。非负 Gaussian表示包括全部 t>0，既不截短时间也不删 singular caps。

只有在 expectation 后才使用原64等权 body law 的 signed-permutation 对称：14个face、每face64个orthants给 `∫∂K G^q=896I`。作者顺序正确，没有对不对称的单 shape 强加face symmetry。最终 `I≤C^qH M^q/896=M^(7/6)/(105·2^(25/6)π^(13/12))`。本审符号等式与另一 C,H,896 路径确认了所有2、π幂次。

旧 M_Gamma 的 §5 reconstruction 同样是解析 Cauchy–Binet/weighted AM–GM矩，不是MC baseline：边际7/64与 `ES^(-1/128)=4^(1/128)Γ(65/128)/Γ1/2` 给 `(1/4)[Γ(65/128)/Γ1/2]^64`。它与 M_H 是同一未知 M 的两个上界；比较表达式大小合法，未改过程。

## 7. 实际有限运行与严格端点比较

登记先于数学执行。一次独立 mincheck exit0，无断言失败或实现修复，脚本计时0.035464667016640306秒；Python3.15.0、SymPy1.14.0、python-flint0.9.0，Arb768bit。精确有限范围仅8组、64 minors、固定Gamma与boundary表达式，无其它matrix sweep。全部数据在results.json；代码和原源哈希在receipt/manifest。

作者512-bit区间的 M_group、M_Gamma、I_group、I_Gamma 均严格包含本审768-bit exact rational interval。上下端使用 Arb `.lower().fmpq()` 与 `.upper().fmpq()` 序列化，随后 Fraction 精确比较，不依赖打印小数。M_H约0.08730273429780108714944687916，I_H约8.92221142884962178793463878578×10^−6，只作可读展示。实际严谨结论是

`M_H768upper<M_Gamma768lower`，

`I_H768upper<author512dyadicupper<R_I<I_Gamma768lower`。

results.json 保存了本审完整 rational upper 与 R_I减它的正有理差。R_I=178444228577/20000000000000000 本来已经有效；若主集成采用 finer512或768端，只是同一固定表达式的算术锐化，无新增模型、矩不等式或科学采样。

## 8. 结论和未解决义务

本审接受 completeproof.md 的解析声明及 registered finite scalar/partition证书。无实质数学反例或待修复缺口。标准 SPD谱分解、普通AM–GM、Gamma/Dirichlet变换、Gaussian积分、Minkowski、Tonelli/MCT 是明确解析工具，有限代码仅核整数结构和标量区间。

得到完整 I≤I_H<R_I，尚未证明有用的全局 `g0≤essinf_K G`。因此不能由本审或其更细 rational endpoint 声称任何新 C7 上界；条件成本仍需独立 g 证书。既有R001成果和阈值均未改动。

复现命令，从本审独立 clone 根执行：

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-7-bounds/audit/round002/hadamard_moment/mincheck.py
```

此轮在注册固定范围成功后停止；不另开权重或指数路线。

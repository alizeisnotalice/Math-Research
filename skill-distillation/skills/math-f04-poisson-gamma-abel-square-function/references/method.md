# Poisson–Gamma 与 Abel 平方函数接口：方法与验收边界

1. 非负自伴L在Hilbert空间上定义P_t=e^{-t√L}；谱参数λ≥0的乘子e^{-t√λ}，欧氏L=−Δ且Fourier约定匹配时为e^{-t|ξ|}。

2. 明确g积分是∫₀∞t||∂_tP_tf||²dt。谱参数λ>0时∫₀∞tλe^{-2t√λ}dt=1/4，λ=0时为0；Tonelli/谱定理得(1/4)||(I−Π_kerL)f||²。

3. 单频plane wave仅是乘子核验，不是L²(Rⁿ)函数；欧氏L²版本零频集合零测，可由Plancherel给(1/4)||f||²。去t权会变成(1/2)||L^{1/4}f||²并要求此量有限。

4. Gamma密度用shape α>0与rate r>0明确给r^αs^{α−1}e^{-rs}/Γ(α)；对Laplace参数u≥0积分e^{-su}得(r/(r+u))^α。扩展至复u须Re u>−r并声明复幂分支，一般混合按实际参数核验。

5. 有限Fourier/谱和先核Abel极限，再凭适用的支配或强收敛推广；逐谱收敛不自动控制极大值。

6. 新增邻近P-1cd2c10e2d8842c2限半空间、半有限von Neumann代数的列/行BMO–HMO与标准Poisson核，平方范数含完整梯度及球Carleson积分；BMO→HMO方向与部分对偶步骤引Mei且省细节。该表示不含Gamma shape/rate、Abel极限或任意中心立方体比较，不能替换本Hilbert谱恒等式。

7. P-ba2136b0727e6a09为粗糙域zero-boundary −Lu=H−divΞ，corkscrew+n-AR、n≥2、real uniformly elliptic L，H/Ξ紧支撑有界；modifiedPR_p已可解且1<p≤2时作者报告q∈(1−ε,p)外推，但其p18/p21打印atom normalization与λ预算不相合。Th5.7是两侧平均积分，proof继承Th5.1；外推原子估计细节留给读者。先用已修source卡与显式PDE入口，不将打印atom theorem当已认证，也不混同Gamma/Abel半群。

8. ROOT已检查热核正原子→L¹弱型扩展，仅限正连续热核；立方体指标边界、Poisson-Gamma私有接口和与平均尺度族的比较仍需另证。

9. 粗糙域PDE显式迁移入口（P-ba2136，非Gamma/Abel）：输入Ω corkscrew+n-AR、实强椭圆L、已可解modifiedPR_p且1<p≤2；对zero-boundary −Lu=H−divΞ先计算δH/Ξ的T²_q area norm（δ^(−n−1)dx），再在data远离2B时用Th5.7的边界截断N梯度均值≲域内annulus梯度均值；两侧都是fint。q∈(1−ε,p)外推需原子与annulus控制，作者细节留读者；其p18/p21原子定义/构造/λ^p预算指数不相合，暂不调用打印atom theorem作已认证入口。energyσ^(1−2/q)的标准修正仅独立代数建议，覆盖和收敛另证，详见TEAM_EF source amendment。

### 正连续热核：有限正原子到 L¹ 弱型

ROOT 的独立推导见 [局部核验记录](audit-root.md)。调用时先把其前提作为假设明确写出：固定维数 n，h_t≥0；对任意紧集 K、固定 x 与 t>0，y↦h_t(x−y) 在 K 连续；对固定紧支撑 f_R，t↦h_t*f_R(x) 在 (0,∞) 连续；有限正原子 ν 已满足 α|{sup_t h_t*ν>α}|≤Cν(Rⁿ)。有限原子界本身不由以下约化证明。

取 f≥0 于 L¹，令 f_R=f1_{[-R,R]^n}。用网格划分支撑立方体，并把每格的 f_R 质量放到格内一点，得到总质量恰为 ||f_R||₁ 的 ν_j。固定 x,t 时，紧集上的空间连续性使 h_t*ν_j(x)→h_t*f_R(x)；先对可数个正有理 t 同时成立。时间连续性保证 sup_{t>0} 等于 sup_{t∈Q_{>0}}。于是 {sup_t h_t*f_R>α} 包含于 liminf_j {sup_t h_t*ν_j>α}；Fatou 和有限原子假设给出测度≤C||f_R||₁/α。再让 R↑∞，正核单调收敛使层集递增且 ||f_R||₁↑||f||₁，得到同一常数 C 的 L¹ 弱型界。一般实/复 f 使用 |h_t*f|≤h_t*|f|。

这里的关键是正核、紧集空间连续性、正时间连续性和已知原子估计；不能把点态网格逼近直接移植到球或立方体指标核，亦不能对不连续或有符号核照搬。

## 不可省略的限制

生成元幂次、Gamma 参数化不同会给出不同核。逐频率 Abel 极限不自动推出算子极限。 粗糙域source打印atom归一化与系数预算有gap；局部化两侧均为平均积分。标准指数修正只局部代数，不作为原文已证明结论。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。

局部证据补充见 [ROOT记录](audit-root.md)；其审读者与范围以记录为准。

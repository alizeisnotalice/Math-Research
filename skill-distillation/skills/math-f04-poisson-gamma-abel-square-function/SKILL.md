---
name: math-f04-poisson-gamma-abel-square-function
description: "用于Poisson–Gamma 与 Abel 平方函数接口的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **交接准备状态（结构性）**：本副本只整理技能结构与引用可移植性；未进行全文深读或数学验证。原评级、可用状态、假设、证明限制和未决条件均沿用原件。结构检查不构成科学验证。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# F04 · Poisson–Gamma 与 Abel 平方函数接口

当前证据状态：**待证据验收：可执行方法工作流，不是已证明的课题结论**。全文转换、论文阅读和证明核验是三个独立状态。（此行声明已被下方【修订】标注替代：该 skill 经外部审核提升为「部分可用」，详见对应修订行与 references/verification-log-20261006.md）
  【修订 2026-10-06（外部审核 R21）】三条核心恒等式经 sympy 符号精确验证（∫tλe^{−2t√λ}=1/4、去 t 权得 (1/2)‖L^{1/4}f‖²、Gamma-Laplace=(r/(r+u))^α），P-1cd2c10e/P-ba2136 守卫与源卡一致；证据状态提升为：部分可用：基础推导已检查，专题证据仍待完整验收。粗糙域迁移细节与私有接口比较待做。

## 输入与产出

半群生成元及谱/Fourier 约定、Poisson 参数 t、Gamma 形状/尺度、Abel 正则化和平方函数归一化。 若选邻近粗糙域PDE入口，另给n≥2、corkscrew Ω及n-Ahlfors regular边界、实强椭圆L、H/Ξ紧支撑有界、zero-boundary variational solution和modifiedPR_p可解证书（1<p≤2）；该输入与Hilbert semigroup/Gamma输入分模型。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 非负自伴L在Hilbert空间上定义P_t=e^{-t√L}；谱参数λ≥0的乘子e^{-t√λ}，欧氏L=−Δ且Fourier约定匹配时为e^{-t|ξ|}。
2. 明确g积分是∫₀∞t||∂_tP_tf||²dt。谱参数λ>0时∫₀∞tλe^{-2t√λ}dt=1/4，λ=0时为0；Tonelli/谱定理得(1/4)||(I−Π_kerL)f||²。
3. 单频plane wave仅是乘子核验，不是L²(Rⁿ)函数；欧氏L²版本零频集合零测，可由Plancherel给(1/4)||f||²。去t权会变成(1/2)||L^{1/4}f||²并要求此量有限。
4. Gamma密度用shape α>0与rate r>0明确给r^αs^{α−1}e^{-rs}/Γ(α)；对Laplace参数u≥0积分e^{-su}得(r/(r+u))^α。扩展至复u须Re u>−r并声明复幂分支，一般混合按实际参数核验。
5. 有限Fourier/谱和先核Abel极限，再凭适用的支配或强收敛推广；逐谱收敛不自动控制极大值。
6. 新增邻近P-1cd2c10e2d8842c2限半空间、半有限von Neumann代数的列/行BMO–HMO与标准Poisson核，平方范数含完整梯度及球Carleson积分；BMO→HMO方向与部分对偶步骤引Mei且省细节。该表示不含Gamma shape/rate、Abel极限或任意中心立方体比较，不能替换本Hilbert谱恒等式。
7. P-ba2136b0727e6a09为粗糙域zero-boundary −Lu=H−divΞ，corkscrew+n-AR、n≥2、real uniformly elliptic L，H/Ξ紧支撑有界；modifiedPR_p已可解且1<p≤2时作者报告q∈(1−ε,p)外推，但其p18/p21打印atom normalization与λ预算不相合。Th5.7是两侧平均积分，proof继承Th5.1；外推原子估计细节留给读者。先用已修source卡与显式PDE入口，不将打印atom theorem当已认证，也不混同Gamma/Abel半群。
8. ROOT已检查热核正原子→L¹弱型扩展，仅限正连续热核；立方体指标边界、Poisson-Gamma私有接口和与平均尺度族的比较仍需另证。
9. 粗糙域PDE显式迁移入口（P-ba2136，非Gamma/Abel）：输入Ω corkscrew+n-AR、实强椭圆L、已可解modifiedPR_p且1<p≤2；对zero-boundary −Lu=H−divΞ先计算δH/Ξ的T²_q area norm（δ^(−n−1)dx），再在data远离2B时用Th5.7的边界截断N梯度均值≲域内annulus梯度均值；两侧都是fint。q∈(1−ε,p)外推需原子与annulus控制，作者细节留读者；其p18/p21原子定义/构造/λ^p预算指数不相合，暂不调用打印atom theorem作已认证入口。energyσ^(1−2/q)的标准修正仅独立代数建议，覆盖和收敛另证，详见TEAM_EF source amendment。

## 证据与失败处理

生成元幂次、Gamma 参数化不同会给出不同核。逐频率 Abel 极限不自动推出算子极限。 粗糙域source打印atom归一化与系数预算有gap；局部化两侧均为平均积分。标准指数修正只局部代数，不作为原文已证明结论。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

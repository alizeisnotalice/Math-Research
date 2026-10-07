# P318 固定端点能量：首批冻结来源合同

本卡从首批 H01 冻结资产逐字节复用，无新增全文阅读、无 future 来源导入，也不增加主张或案例。

- 论文：He，arXiv:2311.05409v1（2023-11-09），9 页；身份及原文问题以[阅读卡](papers/TEAM_HI-P-318df19635ebf8e5.json)为准。
- 原冻结卡：`reviews/TEAM_N/integration/stage16-batch2/generated-v3/skills/math-h01-uniform-cramer-moderate-deviations/references/papers/TEAM_HI-P-318df19635ebf8e5.json`；复制卡 SHA-256：`459be296f9eb3537d28b717a5d20a313cbc674ab898465d883b77238fe021320`。
- 原 PDF：`corpus/papers/P-318df19635ebf8e5/paper.pdf`，SHA-256 `318df19635ebf8e52c12d13ad29ed476599b17b9dffd9abe13932f1f452a1db6`。
- 转换正文：`corpus/papers/P-318df19635ebf8e5/paper.md`，SHA-256 `bc6b9ef2001c2474b449058d7a2135630a563d3c8f35fecc64ba8c8851c324fb`；行号按物理换行计。
- 能量定义：PDF p.3，正文 128–143；Proposition 2.2 及完整局部证明：PDF pp.3–4，147–174。
- 概率主定理 Theorem 1.4：PDF pp.2–3，91–120；假设：PDF p.1，28–50。原文引用的 Hu–Lee 路径 MDP 的原始证明没有在此独立认证。

确定性合同为 T>0、σ²>0，I_T(f)=∫₀ᵀ(f′)²/(2σ²)，f 绝对连续且 f(0)=0；其余函数定义能量为 +∞。任意实数 a 下，固定 f(T)=a：
a²=(∫₀ᵀf′)²≤T∫₀ᵀ(f′)²=2σ²T I_T(f)，而 f(t)=at/T 达到等号，故最小值为 a²/(2σ²T)。
这是确定性变分结果，不需要把它当作某一随机过程的 LDP。

若用于该文的概率分支，先核 iid X_i、μ=EX₁>0、σ²=Var X₁>0；Λ 在全体实参数有限，sup_a[−Λ(a)]<∞；
存在 θ∈(0,1]、v>1、b>0 使 Eexp(θ|X₁|^v)≤e^b。固定 r∈(0,μ)、t>0，正尺度 a_n 满足 a_n/n→0 与 √n/a_n→0。
居中插值路径在 C([0,T],R) 的引用 MDP 使用速度 a_n²/n、上述 I_T，并分别按闭集上界/开集下界使用。
Theorem 1.4 的两个严格命中尾只有对数极限 −μ³t²/(2σ²r)，没有有限阈值概率、正态尾比或鞍点前因子。

阅读卡中保留的缺口仍适用：早命中证明原文仅写 analogous；原 [0,1] 路径的未命中 τ=∞ 需扩展路径或另核指数小的未命中事件；
原文 Cramér 共轭在支持端点的有限最大值断言有反例，须使用 sup；Exp(1) 与 Poisson(1) 示例不能直接满足打印的全部假设。
本合同不认证全部路径 LDP 引文、不提供统一方差分层或中心立方体迁移结论。

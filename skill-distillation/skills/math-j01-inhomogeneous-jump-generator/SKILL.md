---
name: math-j01-inhomogeneous-jump-generator
description: "用于非齐次跳跃过程真实生成元的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# J01 · 非齐次跳跃过程真实生成元

> 审核边界：本流程整理给定对象的定义、前提、推导与证据；它本身不构成专题定理证明。数值或玩具例只核对对应实例。来源结论须在账本列出的条件内使用，未独立核验的迁移必须标注“未证迁移”。本轮逐条断言、原文定位与实际案例记录见工作包 `audit_current/jk/`。
> 文献定位见 `references/handoff-evidence-index.csv`；完整原文应以工作包 `evidence/papers/<paper_id>/paper.pdf` 及其 SHA-256 为准，自动转换文本只作检索辅助。

## 输入与产出

时间-状态依赖转移率/跳核、漂移扩散项、补偿截断约定、测试函数和初始状态。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 从小时间转移概率定义真实生成元 L_t f(x)=lim_{h↓0}(E[f(X_{t+h})|X_t=x]−f(x))/h；先确认极限及测试函数域。
2. 按模型写漂移/扩散与 ∫[f(x+z)−f(x)−χ(z)∇f(x)·z]ν(t,x,dz)，若有跳出域或杀死，明确其边界/杀死项。
3. 检查 Lévy 可积性、局部有界性与时变强度；空间非齐次是 rate/kernel 依赖 x，时间非齐次是依赖 t，二者分别记录。若模型是非齐次纯跳 Q-process，还须登记 conservative/stable Q 与是否爆炸；Zhang (arXiv:1511.05011v2) Theorem 3.1 给出该模型的非爆炸等价判据，Theorem 3.2(a) 给出漂移充分条件，(b) 的必要性另需局部紧可分度量状态空间、紧集上局部跳率界和嵌入链 weak-Feller。不得仅由局部有界率推出非爆炸。
4. 用 f≡1 和坐标函数（在合法域内）检查质量守恒与漂移；由生成元推出演化方程时注明时间不齐次算子次序。
5. 引用前向 Kolmogorov 方程时限定其测试集与时间：Feinberg–Shiryaev (2021), Theorem 3.5 (PDF p. 7; citing [9, Theorem 3]) assumes Assumption 2.3: each state has locally bounded total rate on finite time intervals. It covers each `(q,s)`-bounded set `B` and gives the forward equation only for almost every `t∈(u,s)`. Its minimality statement is within `\hat P` for functions satisfying the initial boundary limit, terminal-time absolute continuity, and the equation; uniqueness additionally requires a regular transition function. Example 3.2 (PDF pp. 6–7; citing [10, Example 2]) has both integrals in the right-hand side infinite for `B=X`, so that equation is undefined there. Do not extend the theorem to every measurable `B` or every `t`.

## 证据与失败处理

从独立安装位置访问文献 PDF 或共享审计记录时，先按 [便携证据访问指南](references/portable-audit-access.md) 定位证据包并调用共享解析器。

冻结系数只是一种近似，不是原生成元；跳到定义域外必须用实际落点/墓地状态处理。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

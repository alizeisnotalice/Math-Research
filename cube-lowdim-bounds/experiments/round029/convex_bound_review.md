# 固定 R 的 beta 凸性列独立审读

结论：以下解析关系无实质阻断；本审读未执行新计算。设 `1<alpha<=beta<2`、`R>0` 固定、`p=beta/alpha`，且 `psi_gamma=(r^(-gamma)-R^(-gamma))_+`。定义

`Phi(t)=(t+R^(-alpha))^p-R^(-beta), t>=0`。

则 `Phi(0)=0`、`psi_beta=Phi(psi_alpha)`，且

`Phi'(t)>=mu=p R^(alpha-beta)>0`。

故 `psi_beta>=mu psi_alpha`，不要求 `mu>=1`。逐跳凸性切线不等式在同一空间点给出

`A psi_beta(x)>=chi(x) A psi_alpha(x)`，

其中 `chi(x)=Phi'(psi_alpha(x))=p min(r,R)^(alpha-beta)>=mu`。这是 **chi 乘旧生成元** 的点态下界，不是 `A(chi psi_alpha)`；也不能把样条项的生成元整体乘 chi。

若旧候选 `K1=c1 psi_alpha+S1` 满足 `c1>=0`、`A psi_alpha>=0`，则

`K_beta=(c1/mu) psi_beta+S1>=K1`，且 `A K_beta>=A K1`。

因此这给出新族的保证可行基线。`beta=alpha` 时准确退化为原核。以上点态推导可在离轴且 `r!=R` 处使用；轴与截断边界不应强行填入有限生成元点值，转分布条件仍须使用已审的局部可积/Fubini桥，而非旧有限源测度接口。

用于发现的旧 float 缓存仍只提供数值列。严格认证时必须以旧 `A psi_alpha` 的真正严格下包围，乘 chi 的严格正区间，并与样条生成元联合包围；候选奇异系数需非负。外区 `r>=R` 有 `chi=mu`，保守列对上述基线只复现旧奇异生成元下界，不凭空增加外区裕量；样条组合仍需认证。

成本不能由可行支配推断改善。奇异核半积分为

`I_gamma=gamma R^(2-gamma)/(2-gamma)`。

保证可行替换的奇异成本比准确为

`I_beta/(mu I_alpha)=(2-alpha)/(2-beta)>1`（`beta>alpha`）。

所以保持 S1 不变的基线成本必升；新优化必须借助内区更强生成元或重配样条来抵消涨价。目标函数与最终准确成本须使用 beta 的完整质量，不能沿用 alpha 的积分权重。本关系是新族探索接口，不是新弱型上界或成本改进证书。

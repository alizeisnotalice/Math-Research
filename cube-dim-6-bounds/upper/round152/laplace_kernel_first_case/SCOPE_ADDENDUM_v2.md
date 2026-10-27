# C152 v2追加口径：group-reduced贡献，不是全R6固定pair事件

本说明追加于原PROOF、single_case及HANDOFF之后；这些原文件
保持原样，不覆盖或改写。数学公式与唯一lambda30求值不变。

130给出的完整恒等式是

  beta2=(46080/pi^3) sum_(C in reduced21fan) alpha_C
                    sum_(i<j) E_C[k_h(A_i)k_h(A_j)].

因全pair-sum在所用群作用下不变，原全空间积分可以约化到该
21-cone fan并付46080群因子。对单独pair，群作用一般会置换其
标签，不能直接保持同一个原始pair不动而把同一因子搬过去。

本首例的量

  (46080/pi^3)*(1/171992678400)*900*I30^2

准确含义是：**group-reduced 21-cone正表达式内，cone0的单个
pair(16,17)项的正贡献**。它不是原始固定pair(16,17)在整个R6
上的事件积分，也不建立那种逐pair全空间等式。

它可作为完整约化正和的一个部分项；不额外乘496，不额外复制
group factor，不自动给其他pattern重数。证书字段
partial_mu2_contribution_interval必须按本说明解释，且保留
is_full_beta2_upper=false、new_M_upper=false、new_C6_bound=false。

v2准备程序把本说明自身的准确hash加入frozen_sources，另在
frozen_sources中保存显式contribution_scope与全R6固定pair否定
标记。原single_case在求值前检查这些源hash，结果再绑定整个
frozen_sources文件hash，从而将本口径与将来唯一证书绑定。

general2D signed积分仍是完整上界的当前瓶颈。首例只认证一个
真实可分离原型，不证明完整beta2上端，不替代149 full compiler
的9961模式覆盖、signed项/尾与全域误差账本。

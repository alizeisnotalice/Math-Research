# R21 六层加权结算：代数核验

`six_layer.wl` 已用Mathematica 15.0.0实际运行，输出见 `six_layer.log`。取 τ=bv/c，w_h=τ^(6−h)，六层和的精确望远镜是 c f₆−bvτ⁵f₀；W=Σw_h=1+O(v)，前五权重之和为O(v)。若 f₀=A₀vM，则端点项除以a恰为 A₀b⁶v⁷M/(c⁵a)，在M/a≤S时为O(v⁷S)。这些恒等式没有假设f₅=o(·)。

设 A=A_v、k=K+ε、x=M/E、y=S/E，两个slack为 h₁=a(1+kv⁶y)−(A−a)x 和 h₂=1−(1−a)y。脚本核得目标比值上界减去1+x恰等于

  h₁/(A−a) + a k v⁶ h₂/[(A−a)(1−a)]。

在0<a<1、A>a、k≥0及两slack非负下，负差的违例集为False。带入a=v^n、A=(1+v)^n，n≥2时归一化前因子极限为1。几何源预算和层间关系本身仍由外部证明提供；这里不证明全输入上界。

复现：`wolframscript -file six_layer.wl`。脚本没有私有路径或自动时限，日志为实际Mathematica15 ARM输出。几何与同层接口已由 proofs/021-sharp-settlement.md 和 reviews/021-settlement.md 审定，真实锐上界由其中列明的R14/R17/R19/R20依赖组合取得；软件只核算式。

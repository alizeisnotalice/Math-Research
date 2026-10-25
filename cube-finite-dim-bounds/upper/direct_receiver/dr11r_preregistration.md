# 11R直接接收对偶预登记

2026-10-25。只upper/direct_receiver/，不改任何旧锁源、总账或其他项目，不提交推送扩招。用户要求本轮完成后等待，故本批结束不追加试验。真实main卡/账本/035输入hash已锁dr11r_inputs.json：读取时U3=8，U4严格14753/1000；14.371等候选不冒充已审。

固定最多六个配置，n3/n4各同三种。两个receiver boxes内部不交，h各box一个常数，r1=1、r2=2（不是半径1/2）。所有边界取闭/半开不影响Lebesgue体积。

每维n中V=4^n，配置如下：

1. volume_critical：B1=[−16,0]×[0,4]^(n−1)，B2=[0,4]×[0,4]^(n−1)。
2. face_adjacent：B1=[−4,0]×[0,4]^(n−1)，B2=[0,4]×[0,4]^(n−1)。
3. vertex_adjacent：B1=[−4,0]^n，B2=[0,4]^n。

先证明a_i(y)=|B_i∩Q(y,r_i)|/(2r_i)^n在每坐标breakpoint grid cell上multi-affine，全空间Th最大值落有限断点笛卡尔积。每维断点仅l_ij±r_i,u_ij±r_i；不另加数值网格、不只检查receiver中心。然后Python Fraction准确列所有coeff vectors，Pareto去冗余，解2变量LP并给准确primal h和有限input原子dual certificate。仅这六配置，允许六进程独立并行；不为CPU占用重算。

理论scope另列：完整一般L1/共同r的条件对偶定理及有限空间近似误差预算。h≡1/2给有限两尺度U≤2只是传统finite-scale envelope，不算新能力；若提到2−2^-n亦明确是固定kernel envelope基准。有限成功不替代任意E、任意r(x)的一致mass引理。需指出一般measurable E未必可由内含full boxes逼近；不给未证strong duality。

停止于六配置准确证书/反例和上述条件接口。MATLAB给相同六例断点load/质量检查脚本由主串行，独立变量前缀，无clear/path。无新Cn上界，无下一轮搜索。

实际配置计算前的设计缩小：原separated草案仅复述无重叠负载的平凡情况，未运行。按主要求不把已覆盖结论当新能力，替换为volume_critical以检验uniform5/4是否还能由per-box h改善；仍只六配置，没有第七例。配置JSON为最终权威输入。

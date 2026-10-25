# A049 二维边际短预算独立审计：PASS

接受原 body32 模型、A3 中段域 D、rho=1/64、mean_B∈[1/4,3/4]、epsilon=1/1000 的严格结论

    ||∫_0^epsilon p_s ds||_(L6/5(D))
       <66002272277/250000000000000000。

本结论覆盖 D 上全部实际参数事件。没有认证其它区域、长时间、全 I/global g 或 C6 改善。

独立证明确认以下接口：

- (X1,X4)=(U+V,U−V)，U,V 独立 Cauchy(16s)，逆 Jacobian1/2；六列 basis determinant32、retained pivot2、free Schur determinant16，固定坐标的真实边际积分也成立。
- 余4坐标的 flip 只用于全密度 p，准确 orthant 因子1/16。actual P/B 不单独声称对称，只用≤p转移 L1 上界。
- B三维块的二维零和残差上，坐标泛函范数 sqrt(2/3)，所以任取固定 j=4 的真实范围为[1/4−sqrt(2/3)/64,3/4+sqrt(2/3)/64]。无需因“某j”再付3倍。
- 保留这一坐标范围、放大其余4个正坐标至整个orthant，得到 Cjoint=(128/π²)∫(1−z²)^−2dz，以及∫_D p_s≤Cjoint s²。
- 043 的两大参数事件和 actual 不交补事件完整。为 L∞ 放大的坏事件可重叠，插值仍用 actual B≤p 的一次 L1 预算，未把扩大后事件误当互斥。
- L∞–L1 插值给 s^(4/3)、s^(5/3)，Minkowski/单调收敛给3/7、3/8及epsilon^(7/3)、epsilon^(8/3)。

数学审查未发现阻断性缺口。完整点态/事件/测度论证明见 proof_review.md。

实际最终源锁后仅一次4配置批：PID50203、exit0；本机 SymPy+Arb768，实际数学0.027817秒，源锁至完成13.913228秒。源锁实际开始 UTC2026-10-25T21:51:12.311233；之前218.311233秒等待/理论准备单独记账，不计入720秒实际上限。

前两配置为准确结构及符号；后两为固定 Cjoint 和固定新短预算。独立 Cjoint≈15.13131638048143112946，新总上端≈2.64009089107041×10⁻⁷。严格判定使用完整有理端点，未使用这些展示小数。

旧043的 Dbig、Mbad、oldbudget只从完整保存dyadic端点读取。新预算用两种单调端点的六次根表达式重算，其准确上端之和严格小于规范阈值；另与作者完整区间相容。保存旧预算只用于上包公式间的准确比较，未重算它，也不声称实际范数按该比例降低。

作者最终 pins：manifest d70a5d3684f3d2e2373e6348b55cf0cbfac51486bb766b4e68c78e01437bfba0；results de7ebaa9e2cec7d64bcd33435224daa5a205d900ead76e01129b802b68900a8e。完整 whitelist 和所有外部输入哈希吻合，计算前后源不变。

0作者writer/checker、0旧043预算重算、0frame/fan/triples全枚举、0kernel/native/G/H/K/p、MATLAB、MC、grid、优化、扩招、commit/push。全部句柄关闭。有限批至此停止，主保留最终证明和集成责任。

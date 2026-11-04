# 基础部分逐项核查

本表以原始 PDF 的 PDFKit 恢复文本为索引。每项的“已证”均指本次正文给出了独立证明；它不包含后续章节才会使用的固定模型或全体集合结论。

| 原页与命题 | 假设与量词 | 复核结论 | 新位置 | 修订理由 |
|---|---|---|---|---|
| pp.4--5，引理 1.1 | $f\in L^1_{loc}$，严格阈值 $\lambda>0$，三种基 | 已证 | `fnd:measurability` | 中心基以固定边长卷积连续性处理；非中心基用严格见证的外扩，避免对不可数参数上确界直接作错误可测性推断。 |
| p.5，引理 1.2；定义 (3)--(4) | $0<|E|<\infty$ | 已证 | `fnd:top-level`, `fnd:def-halo` | 保留严格水平集和 Lebesgue 归一化。 |
| p.6，注 1.3--1.4 | 中心分数阶模型；外层权 | 已证的定义性说明 | `fnd:center-density`, `fnd:outer-weight` | 明确外层权不等于内禀 $w\,dx$ 平均，也不等于两权范数问题。 |
| p.6，引理 2.1；命题 2.3 | $r>0$；所有可测 $E$ | 已证 | `fnd:scaling` | 变量替换显示齐次度为 $\beta$；唯一指数为 $\beta/n$。 |
| pp.6--7，定理 2.2 | $\lambda>0$；三种基 | 已证 | `fnd:raw-divergence` | 保留从单位方体的远场下界到伸缩发散的完整量纲计算。 |
| pp.7--8，定理 3.1 | 单一严格见证 $B$；$0<\alpha<1$ | 已证且最优 | `fnd:double-density` | 对称差估计以凸函数端点最大值实现，不能将两个单侧误差相加后称为最佳。 |
| p.9，推论 3.2、注 3.3 | $\alpha\uparrow1$，固定 $\theta$ | 已证 | `fnd:rho-expansion` | 只给固定参数展开；不暗示对 $\theta\downarrow0$ 一致。 |
| p.9，引理 4.1 | 两个矩形的对称差不超过 $\varepsilon|R_0|$ | 已证 | `geo:rectangles` | 使用严格正确的 $A_n=(1+\varepsilon)/(1-\varepsilon)^{n-1}$。 |
| pp.9--10，引理 4.2 | 两中心方体；$0<\varepsilon<1/2$ | 已证 | `geo:centered-cubes` | 显式说明逐坐标交长的推导。 |
| pp.10--12，定理 5.1 | 上界另需 $\varepsilon_\alpha<1/2$ | 已证 | `geo:explicit-bounds` | 修正历史重写中的括号遗漏：非中心上界为 $(1+2h_n(\varepsilon_\alpha))^n\alpha^{-1/(1-\theta)}$；旧式只写了 $1+2h_n(\varepsilon_\alpha)$。 |
| p.12，推论 5.2、注 5.3 | $\alpha\uparrow1$；三种基 | 已证 | `geo:endpoint-orders`, `geo:scope` | 将中心项限定为可验证的首项下界；高维全体 sharp 首项仍未由本节证明。 |

## 依赖与范围

本部分只使用 Lebesgue 测度、严格水平集和有限测度的测试集合。它提供后续权模型和切向模型需要的尺度窗与紧性输入，但不产生一般凸体、内禀加权极大算子、全体中心最优首项或 PDE 正则性结论。

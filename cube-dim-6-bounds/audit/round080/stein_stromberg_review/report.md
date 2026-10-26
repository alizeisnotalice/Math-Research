# A080：经典覆盖固定分配式独审

结论：**PASS_SPECIFIED_ALLOCATION_DIAGNOSTIC_ONLY**。
在 beta>1、整数 m>=1 的既定范围内，080 的完整选择、覆盖、壳多重性和一般
L1 全尺度传递成立，给出指定的六维分配式

\[
K(\beta,m)=\left(1+\left(\frac{\beta+2}{\beta}\right)^6\right)
\left(1+m\left(\frac{\beta+1}{\beta}\right)^6\beta^{6/m}\right).
\]

全参数解析证明及冻结准确事实支持
\(K(\beta,m)>1602/25=64.08>64\)。这只排除该分配表达式改善64；不是实际最优
覆盖常数或 C6 的下界，不排除更锐利几何、其他选择规则或其他覆盖证明。beta=1
未被纳入。没有新 C6/Cn 上界，主代理保留最终数学责任。

## 绑定和独立范围

作者目录：`Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round080/stein_stromberg_finite6/`。
作者 manifest SHA256 为
`5b69ee9803ff0a133d0ebf95fe0423da878710192270b3f197f817ecea1a319f`；
theorem SHA256 为
`6a44d397fce20e250e55ee85bdcb7022006e4b8ec4102d5dd85c44a9ea21113b`。
全部 manifest 文件、输入绑定与运行收据核对一致。

原始 Stein–Strömberg 1983 PDF 只读 numerics/078 冻结文件，SHA256 为
`ecb7c3bbbf804545fd68795f5628d242ba8d03f2a3ebeab67806af9c8da6c615`。
本审只提取局部 PDF 第3–5页，即印刷261–263；未重读全文或重新下载。准备中
PyMuPDF 不可用，改用已有 pypdf 完成局部读取；这不是科学配置或科学程序重试。
主五个 Fraction 事实只读和静态审核，没有调用其程序重新计算；MATLAB 收据只读，
本审没有 MATLAB 调用。公开原文提供选择结构，当前全参数64.08诊断属于内部
量化结果，本审不附加发表或穷尽新颖性声明。

## 关键审核结论

完整推理见 `proof_review.md`，其中明确补足了容易误用的环节：

1. 每次选所有仍合格球中半径最大的球。因为旧增量和权重一经生成不再改变，
   分数只增不减，不合格球不能后来变合格。若未选球半径为 r，在第一个更小
   半径入选前它已不合格；若没有更小半径，终止时即得同一结论。因此所需
   大球控制只使用半径至少 r 的先前球。积分作用于未选球的并集，不重复加其体积。
2. 归一化最小贡献球半径为1后，r_j>=beta 给出
   beta*r_j+beta<=(beta+1)r_j；而 beta>1 给出 r_j>1，保证 j 确实在 k 之前，
   可用 k 的中心条件。r_j=beta 由 outer 部分包含。beta=1 不能照抄这一步。
3. 每壳中互斥的是 I_j，而不是原球。I_j 包含在半径(beta+1)b的一个球内，
   分母至少为 |B|(beta*a)^6，故壳成本与六维体积因子准确。m个半开壳覆盖
   [1,beta)，outer覆盖[beta,infinity)，等半径和壳端点均未漏计。
4. 对每个中心点用原 witness 半径 R，再以 R/beta 为 small ball。star expansion
   恰回到原 witness；有限覆盖、权重 Tonelli 恒等式及全尺度 L1 传递均成立。
   固定 R 平均的连续性和 Lebesgue 内正则性使 compact-to-level-set 步骤完整。
   这里不需要把原文辅助 maximal operator 当成与 centered operator 相等。
5. exp(u)/u 的导数符号给出对所有正实 m 的松弛下界，因此覆盖所有正整数 m。
   三个 beta 区域无空隙；严格性来自正级数尾和有限 beta，不要求 minimizer
   本身严格。封存五事实足够支持64.08阈值，无需参数搜索。

## 新执行及限制

事前固定4个配置，仅一次科学调用，全部 PASS。一个一维五区间 fixture 的
beta=3/2，球依次为(center,radius)=(-4,4)、(4,4)、(59/10,3)、(8,1)、(10,1)。
新 fixture 只检验维数无关的选择逻辑，不能替代六维普遍证明。

|配置|实际独立核验|结果|
|---|---|---|
|1|源绑定；精确 greedy；较大球拒绝后的持续性；非平凡增量|PASS|
|2|同一 fixture 的 outer 条件、壳体积；r=beta 的符号边界|PASS|
|3|全参数 minimization 导数、首因子单调性、六维壳比例与 K 组装|PASS|
|4|同一 fixture 的原尺度 witness 和一个精确阶梯函数的加权 L1 传递|PASS|

选中 A、B、D、E，拒绝半径3的 C；准确增量长度为8、8、1、2，总计19。
E 的入选中心分数恰等于1。在 x=9，outer 权重为2/3，壳权重为1，均满足所需
分配。对 f=indicator[-10,12]、lambda=1/2，全部 star 平均准确为1，加权积分
准确为19。没有用浮点采样或空间网格替代这些有理恒等式。

科学 PID86910、guard PID86902，exit0；CPU0.186138秒，科学墙钟0.189545秒，
含清理的总墙钟0.234080秒，子进程已回收。CPU软硬限均30秒，工作170秒加
清理10秒的单程序总限180秒，只允许清理自己的新进程组。工具 chunk `b8004a`。
未超时，未自动重试，未追加参数或科学配置。

主五事实未重算，0参数网格、0旧baseline、0body32 kernel、0MATLAB、0递归代理、
0stage/commit/push。原080、078及本人的074冻结目录均未改动。审计已停止。

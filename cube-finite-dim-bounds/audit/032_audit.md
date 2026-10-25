# 032 cb11完整96锥与退化Green表达：独立审计

2026-10-25。结论：主最终导入的 cb11 报告SHA `2e8126458ee3b5b6a472be612645dabd411e8705271c8195a18f9c85cc8b34c2`、模块SHA `ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849` 的三维完整pointwise表达通过。本审包括96锥覆盖、符号/Jacobian、真实Fourier密度、谱体积/峰值/time尾、截断时间DD、全部zero/repeated实部极限、精确全Γ聚合及有限/无穷分类；不只是旧36点或两个物理锚点。

未发现需通知minimum/J暂停使用真实公式的数学缺口。不过这仍不是跨节点层的连续box区间算法，也不是g/J或新上界认证。module的float/mp值不成为严格端点，generic DD系数不能直接在含coalescence的box除小差；下游需另证连续误差。四维u11不在本批。

## 锁源、预锁例与独立路线

按011i登记，先锁 `032_source_snapshot.json` 的16项来源；作者manifest SHA为 `7b22d25442f8d13837ece005849d2d142a18fd0e33d57543663e244de05701fa`。`032_cases.json`在计算前锁定十个参数例，含六类节点、τ0/正τ的face方向切换、body奇线、原点。没有添加空间grid或用小epsilon外推有限性。

独立Python首次执行通过，内部约0.124秒，无修补/重复。新代码从真实body/face谱线性forms直接生成全部96锥、每锥全10个forms的单符号判定，并以实函数gε(b)=log(ε²+b²)/2的Hermite递推重建全96项。与作者48个反锥pair的分支表实现不同。

原module仅只读执行green_expression作重建后的准确系数比较；没有调用其green/heat、cone-checker/limit-checker或time积分。此执行会载入作者依赖NumPy，不冒称完全未执行作者代码；独立计算本身使用SymPy1.14.0、Fraction、python-flint/Arb256bit。十例原始x下的全96系数与module归一后的系数经正确尺度映射完全一致。031已审真实正权过程作为明确概率/一般L1比较背景，不重复。

## 96锥的真实覆盖和谱线性

对任意频率，将其绝对坐标排序为a≥b≥c≥0。body谱和为3a+b+c+|a−b−c|，即max(4a,2(a+b+c))；face和为4a+2b。设A=(1,0,0)、B=(1,1,0)、C=(2,1,1)、D=(1,1,1)。ABC、BCD的逆坐标分别为

\[
(a-b-c,b-c,c),\qquad(b-c,a-b,b+c-a).
\]

在排序chamber中，第一组非负恰对应a≥b+c，第二组非负恰对应a≤b+c。它们的像覆盖整个chamber，公共BC仅为边界；独立准确矩阵计算给det分别±1。严格坐标大小和非零符号唯一确定一个signed permutation chamber，48个此类chamber在内部不重叠，边界是有限平面零测集。因此96锥确实覆盖R³ a.e.，不是从数量、体积或有限采样推覆盖。

新checker另外对每锥、每个原body/face方向v，检查三个v·ray_j全同弱号。故对所有u_j≥0，绝对值谱在该锥真正线性。ray rates由真实forms重新求和，得到

\[
k_A=4+4\tau,\ k_B=4+6\tau,\ k_C=8+10\tau,\ k_D=6+6\tau.
\]

τ≥0时都严格正。主cones.json全部锥、signed det、rate整数系数与独立构造一致；共有50条不同ray、48对反锥。τ0时face过程退为零，但96个细分锥仍是body谱的合法细分，没有凭空保留face jump源。

## Fourier密度、谱体积与完整大时间尾

031正权过程给特征函数exp(−tψτ)、严格正连续概率密度pt。ψτ≥2|ξ|保证Fourier积分绝对收敛。各锥ξ=Σu_j r_j、u_j≥0、Jacobian |detR|=1，逐个正实部指数积分得到

\[
p_t(x)=\frac1{8\pi^3}\sum_{96}\prod_{j=1}^3(tk_j-i x\cdot r_j)^{-1}.
\]

Fourier inverse正常数(2π)⁻³、未单位化ray和真实τ均已保留。反锥与本锥的项严格共轭，才可写48个2Re；这不允许取abs、删负实部项或改变符号。

在每个u锥内，谱单位球是Σk_j u_j≤1的simplex，体积1/(6Πk)。48个ABC加48个BCD给

\[
\operatorname{Vol}\{\psi_\tau\le1\}
=\frac5{6(1+\tau)(2+3\tau)(4+5\tau)}.
\]

其6倍除8π³给pt(0)的t⁻³系数。因此

\[
0<p_t(x)\le p_t(0)=\frac5{8\pi^3(1+\tau)(2+3\tau)(4+5\tau)t^3},
\]

其中上界来自Fourier三角不等式，正性来自真实概率过程。t≥T的Green尾即峰值系数除2T²。τ0系数5/64，τ1/3系数15/544；独立符号核全部一致，未用36点证明峰值或尾界。

## 独立实DD证明：分支与所有confluence

设每锥b_j=x·r_j/k_j、W=1/Πk_j。ε>0时，time被截离0，全部96个单锥积分绝对有限。对互异实节点，直接有实有理恒等式

\[
\frac{d}{d\epsilon}\,DD_2(g_\epsilon;b_1,b_2,b_3)
=-\operatorname{Re}\prod_j(\epsilon-i b_j)^{-1},
\quad g_\epsilon(b)=\tfrac12\log(\epsilon^2+b^2).
\]

checker准确核了恒等式；其中未加负号的 ReΠ(t−ibj)⁻¹ 在t=ε处为

\[
\epsilon(\epsilon^2-\sum_{i<j}b_ib_j)\big/\prod_j(\epsilon^2+b_j^2).
\]

DD与积分在ε→∞均趋0，因此它们相等。对固定ε>0，time integrand的绝对值≤t⁻³，统一可积；DD也按实Hermite导数连续延拓，所以恒等式涵盖全部重复/zero节点。它等于作者右半平面principal complex log DD的实部，因为实节点下取实部与线性DD/Hermite导数可交换。故没有依赖跨负实轴log分支的隐含选择。

当b≠0固定时，gε(b)=log|b|+O(ε²)，gε'(b)=1/b+O(ε²)，gε''(b)/2=−1/(2b²)+O(ε²)。零节点处gε(0)=logε、gε'(0)=0准确成立。非零x不可能同时垂直于一个独立锥的三条ray，所以只有以下六类；独立Hermite递推得到完整规则：

| 节点 | γ | finite |
|---|---:|---|
| a,b,c互异非零 | 0 | Σ log|bj|/Πk≠j(bj−bk) |
| a,a,b非零且不等 | 0 | [log|b|−log|a|−(b−a)/a]/(b−a)² |
| a,a,a非零 | 0 | −1/(2a²) |
| 0,b,c非零且不等 | 1/(bc) | log|b|/[b(b−c)]+log|c|/[c(c−b)] |
| 0,b,b非零 | 1/b² | (1−log|b|)/b² |
| 0,0,b非零 | −1/b² | log|b|/b² |

所有b的正负在分母和导数中保留，只有log参数可取abs。这里直接计算实部，0,0,b情形不会误留下纯虚1/ε；非零x没有gε''(0)的实ε⁻²项。全部三节点为0仅在原点，那里由pt(0)直接判ε⁻²发散。

六类实部都是γ logε+finite+O(ε²)，对固定x,τ有限锥求和后无更慢的未列无穷项。此余项是pointwise的，不声称跨非零差/零节点box的一致常数。

## 全Γ聚合与真实有限/无穷分类

由真实pt≥0，Gε=∫ε∞pt随ε↓0单调趋真实扩展G。必须先对全部锥加权聚合Γ=ΣWγ和finite。若Γ<0，Γlogε→+∞；若Γ=0，G等于有限log表达。若Γ>0，会得到Gε→−∞，与正性矛盾，因此合法输入的准确表达不可能出现。Γ=0时上面的实部余项已排除更慢发散，不是仅凭数值Gamma接近0判有限。

独立十个预锁参数例覆盖全部六类。关键准确事实（Γ均指未乘1/(8π³)、按输入原x计算）：

| x,τ | 全Γ | 结论 |
|---|---:|---|
| (1,0,0),1/3 | 0 | finite；局部32锥γ正、16锥γ负，必须全合 |
| (1,1,0),0 | 0 | 精确G=(4log2−1)/(8π³) |
| (1,1,0),1/3 | −79/612 | true +∞，新增face jump线 |
| (1,1,1),1/3 | −6/17 | true +∞，body jump线 |
| (0,9,8),1/3 | 0 | finite，含0,b,b及a,a,b |
| (8,1,1),1/3 | 0 | finite，含a,a,b |
| (8,1,0),1/3 | 0 | finite，含a,a,a |

原点另有peak系数15/544的ε⁻²发散。完整十例的精确rational常数/log系数、局部符号计数及有限点Arb包围均存 `032_results.json`；原点无强塞log公式。它们是对实现的有限检验，完整表示的证明来自上面覆盖/Fourier/实DD及六类穷尽，不来自这些点的通过。

## module尺度、精确输入与连续证书边界

module先将x按s=max|x|归一，返回scale及归一后的系数。checker用全部96原x系数与返回值比较：constant、Γ、各log系数乘s²，log参数除s后逐项一致。Γ符号不受此归一影响，finite的Γ=0时Σlog系数=0，所以最终除s²正确。若要解释Γ≠0时的原始截时ε表达，须同时使用ε/s；不能把归一化Γ误标为原x的数值系数。

module在合并log之前保留signed knots的所有分母/导数；ΣA_q+Γ=0的准确检查有效。finite时可用log(q/qref)减少共同log正常数的取消，不能删除signed项。float坐标/tau解释为其十进制拼写，若要精确τ=1/3应传Fraction或字符串；它不会替调用者证明近似Chebyshev点的位移误差。

auto条件检查或80位mp都是数值稳定性措施，不能当作舍入的严格区间保证。这里独立Arb只给预锁finite rational点的真实值包围。连续盒跨zero/coalescence需新误差证明，ε截时策略也还须付小time残差；这些不因pointwise表示通过而自动成立。

## 实际工具、复现和停止

最小复核（主导入后，stdout-only；可选--save）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/032_check.py
```

审计clone执行加 `--project /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds`。入口核全部16主来源hash。MATLAB最小脚本只交主串行核实DD恒等式、带±号的8个scalar confluence极限及peak；它不是新增空间点批或作者36点的重复，也不代替解析证明。实际状态见 `032_matlab_receipt.json`。

本批无谱成本、g/J grid、time数值外推、优化或LP；不审下游minimum/J/四维u11。完整pointwise表示通过，连续候选证书另锁源。最终manifest完成后停032，不改作者/旧证据/账本/其他项目，不提交推送或扩招。


最终锁数学清单时，新032 scalar MATLAB脚本仍待主串行实际exit/log，未收到回执，不宣称已运行或通过。脚本SHA aa4f5ecfbd22cb85b902a137843563df2857d059ef2cfcba747baedd29c1df74保持；本次已有实际独立Python准确检查。此补充cross不阻塞完整解析表示结论；若主稍后实际运行，可另存后续receipt于主总报告，不把本pending记录冒作已成功。收尾16项主来源及checker/result关联hash再核未变，032数学有限批已停止。

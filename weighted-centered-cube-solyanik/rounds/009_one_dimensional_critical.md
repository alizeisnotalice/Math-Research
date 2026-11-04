# 第九轮：一维无限质量权的内禀临界发散

## 目标

证明一维内禀共同中心等指数模型在临界及次临界区间的发散，不把外层周期条带结论直接移植到内禀分母。

## 严格结论

设 \(\mu=w\,dx\) 是非零局部有限绝对连续测度且 \(\mu(\mathbb R)=\infty\)。对任意 \(m\ge2\)，
\[
C^{(m),\rm int}_{c,1,w}(\lambda)=\infty,
\qquad0<\lambda\le(1-1/m)^m.
\]

证明先以累计质量坐标把 \(\mu\) 化成无穷长度区间上的 Lebesgue 测度，再使用周期 \(m\)-分割。任意区间的第 \(i\) 个分块质量与总质量的 \(1/m\) 相差至多两个端点区间；更强地，选择总质量为整数的中心区间时，周期性使每个分块质量恰为总质量的 \(1/m\)。向补集 \(A_1^c\) 加入有界的有限正质量集合 \(F\subset A_1\)：它只严格增加第一个因子，其他因子不变。于是乘积显式为
\[
\left(1-\frac1m\right)^{m-1}\left(1-\frac1m+\frac{\mu(F)}N\right)>
\left(1-\frac1m\right)^m,
\]
其中 \(N=\mu(Q)\in\mathbb N\)。每个点均有这种共同中心见证，故光晕为全空间并有无穷 \(\mu\)-质量。

这是真正的内禀构造：每个因子的分母仍为同一 \(\mu(Q)\)，没有借用外层 Lebesgue 条带结论。

## 计算核验

`critical_infinite_certificate.py` 对 \(m=2,3,4,8\) 的有理参数实际验证精确恒等式
\[
\prod_i p_i-((m-1)/m)^m=((m-1)/m)^{m-1}\mu(F)/N>0.
\]
`critical_infinite_check.wl` 用 Wolfram Language 独立复算同一精确恒等式。两者只核验有限代数步骤；一般测度、零权区间和全空间光晕由 TeX 中的证明完成。

实际运行于 2026-11-05：Python 3 的四组精确余量分别为
\(1/238,8/2349,9/2624,823543/1017118720\)，均为正；
WolframScript 15.0.0 对同四组返回 `True`。`critical_infinite_smoke.tex`
以 XeLaTeX (TeX Live 2026) 成功编译为两页 PDF；逐页视觉检查未见缺字、公式溢出或重叠。

## 交付物

- `critical_infinite_extension.tex`：可并入主文的自足定理与证明；
- `critical_infinite_certificate.py`：Python 精确有理数复算；
- `critical_infinite_check.wl`：Wolfram Language 独立符号复算；
- `critical_infinite_smoke.tex` 及其 PDF：最小编译检查入口与输出。

## 范围与下一步

本轮只证明一维且 \(\mu(\mathbb R)=\infty\)。有限总质量时光晕的总质量不能自动发散，必须使用分母趋零或可变公共集合的独立构造；高维情形也不由本证明推出。

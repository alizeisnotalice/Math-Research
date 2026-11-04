# 课题2计算证书

本目录只核验可以明确分离为有限代数、有限几何或数值扫描的部分。运行时间为 2026-11-04；工具为 Python 3.15、Wolfram Language 15.0.0 和 MATLAB R2026a Update 4。

## 实际运行

```
python3 exact_certificate.py       exit 0
python3 scan_independent.py        exit 0
wolframscript -file symbolic_check.wl    exit 0
MATLAB: run_matlab_file(matlab_check.m)  completed; all displayed pass=1
python3 water_filling.py           exit 0
python3 star_geometry_check.py     exit 0
python3 -m py_compile *.py         exit 0
```

`exact_certificate.py` 用 `fractions.Fraction` 实际核对了

\[
\lambda_m^*=\left(\frac{m-1}{m}\right)^m,
\qquad s_m(\lambda_m^*)=1,
\qquad \Phi_m(\lambda_m^*)=0
\]

在 \(m=2,3,4,8\) 上均为精确有理数恒等式。独立 Python 和 MATLAB 浮点扫描了每个阈值两侧 \(10^{-9},10^{-6},10^{-3}\) 的点；\(\Phi_m\) 在下侧均为负、上侧均为正。Wolfram Language 还返回

\[
s_m(1-t)=t+\frac{m-1}{2m}t^2+
\frac{(m-1)(2m-1)}{6m^2}t^3+O(t^4).
\]

这些结果只支持标量参数和级数的正确性，不证明多线性常数的传递。

## 全权族下界的可实现星型构造

设 \(0<\delta<1/4\)，中心块 \(B_0=[-\delta,\delta]^n\)，并对每个 \(v\in\{-1,1\}^n\) 取卫星块 \(B_v=v+B_0\)。令 \(E=B_0\)，权在 \(B_0\) 上有总质量 1，在每个 \(B_v\) 上有总质量

\[
q=(1-\eta)\frac{1-\alpha}{\alpha},\qquad 0<\eta<1,
\]

其余处为零。可取每块上的常数密度，因此 \(w\in L^1_{\rm loc}\)。对每个 \(x\in B_v\)，以 \(x\) 为中心、半径 \(1+2\delta\) 的 \(\ell_\infty\) 方体包含 \(B_0\cup B_v\)，不接触其他卫星块；故

\[
M_{c,w}\mathbf 1_E(x)\ge\frac1{1+q}>\alpha.
\]

于是

\[
\frac{w(\{M_{c,w}\mathbf 1_E>\alpha\})}{w(E)}
\ge 1+2^nq
=1+2^n(1-\eta)\frac{1-\alpha}{\alpha}.
\]

令 \(\eta\downarrow0\) 得到全权族公式下界所需的逼近方向。`exact_certificate.py` 在 \(n=1,2,3\)、\(\alpha=1/2,3/4,9/10\) 上用有理数验证了严格性和误差，`scan_independent.py` 与 `matlab_check.m` 另外扫了 \(\alpha=0.99,0.999\)。此构造本身是严格的连续密度构造，图见 `star_geometry.svg`；它仍依赖正文另证的中心方体覆盖上界才能得到等号。

`star_geometry_check.py` 还以 \(\delta=1/10\) 精确检查了 \(n=1,2,3\)：每个卫星的见证方体包含中心和自身，且因 \(2-2\delta>1+2\delta\) 不会接触符号不同的卫星。对任意 \(0<\delta<1/4\) 同一不等式给出所需间隔余量。

## 加权指数的标量优化

`water_filling.py` 考虑

\[
\max\left\{\prod_i(1-d_i)^{a_i}:d_i\in[0,1),\ \sum_i d_i=s\right\},\qquad a_i>0.
\]

对活跃坐标，KKT 方程为 \(d_i=1-a_i/\Lambda\)，其中
\(\Lambda=\sum_{i\in A}a_i/(|A|-s)\)；不活跃坐标满足 \(a_i\ge\Lambda\)。程序用精确分数枚举活跃集，并与独立的有理网格比较。四组测试均验证网格目标值不超过 KKT 目标值。

这可用于论文中解释不等指数的标量缺陷分配；它不自动给出共同见证方体或原始多线性算子的界。

## 限制

- 数值扫描、有限原子质量和 KKT 证书不代替 Besicovitch 覆盖、Boolean 原子化平滑或一般集合的传递证明。
- 临界以下的外层发散、固定权端点系数以及外层集中模量结论均未由此计算认证。
- 任何进入主文的算子等式仍需由主证明单独检查量词、严格水平集和共同中心条件。

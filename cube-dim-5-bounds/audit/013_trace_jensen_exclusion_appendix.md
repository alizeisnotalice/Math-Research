# 013 附录 — trace-Jensen Gaussian 下界目标的严格排除

令
\[
D=(\det A)^{-1/2},\qquad B=A^{-1}.
\]
在五维，AM--GM 给
\[
\operatorname{tr}B\ge5(\det B)^{1/5}=5D^{2/5}.
\]
所以
\[
{D\over(\operatorname{tr}A^{-1})^2}\le {D^{1/5}\over25}.
\]
第 013 主审计已给出
\[
g_{\rm mix}=2(2\pi)^{-5/2}E\left[{D\over(\operatorname{tr}A^{-1})^2}\right].
\]
因 \(z^{1/5}\) 凹，Jensen 和原点峰值恒等式
\[
 E D=(2\pi)^{5/2}p_1(0),\qquad
 p_1(0)={35509\over764411904\pi^5}
\]
给出
\[
g_{\rm mix}\le {2\over25(2\pi)^2}\,p_1(0)^{1/5}. \tag{*}
\]
指数完全相消：外面的 \((2\pi)^{-5/2}\) 与
\(((2\pi)^{5/2})^{1/5}=(2\pi)^{1/2}\) 合为 \((2\pi)^{-2}\)。峰值公式来自
第 004 轮独立 body-fan 审计，而非浮点锚点。

`013_trace_jensen_exclusion_check.py` 用纯 `Fraction` 作了严格比较。它由
Machin 恒等式和交错 arctan 级数得到有理下界
\[
\pi>{1231847548\over392109375}>{333\over106}.
\]
将此下界代入 (*) 的五次幂上界，脚本精确验证
\[
\left({2\over25(2\pi)^2}p_1(0)^{1/5}\right)^5
 <\left({11\over100000}\right)^5.
\]
其有理 \(\pi\) 下界产生的展示上界约为
\(8.769664596\times10^{-5}\)，已严格小于 \(1.1\times10^{-4}\)；主
256-bit Arb 值约 \(8.769656397\times10^{-5}\) 与此一致但不是本排除的
证明依据。

因此 trace-Jensen 的该特定混合接口不可能认证
\(g_{\rm mix}\ge11/100000\)。应停止这一目标的数值认证搜索；此结论不排除
保留更多矩阵信息、不同不等式或其他 Gaussian 结构所给的更强下界。

# Carathéodory、Krein–Milman 极点：方法与验收边界

1. 核对凸性和具体拓扑；区分有限维 Carathéodory 条件与局部凸紧集的 Krein--Milman 条件。

2. 在 R^d 中对 conv(S) 内的点，用仿射相关性将有限凸组合化简到至多 d+1 个点。

3. 对 Hausdorff 局部凸空间中的紧凸集，只在 Krein--Milman 假设下推出它是极点闭凸包。

4. 检查目标究竟需要有限凸组合、闭包还是积分表示。

5. 若文献只是提到 Krein 空间或 Krein 算子，须确认存在真正的凸极点论证。

## 不可省略的限制

无限维 Krein--Milman 一般只给闭凸包，不保证有限凸组合。Krein 空间是带不定内积的空间，不能凭同名认定为 Krein--Milman 结果。

# 有限维Carathéodory消元：已检查基础推导

x=Σ_{i=1}^mλ_i a_i，a_i∈R^d，λ_i>0，Σλ_i=1，m>d+1。
(a_i,1)∈R^{d+1}线性相关，故有非零c满足Σc_i a_i=0、Σc_i=0。
c含正负项。取t=min_{c_i>0}λ_i/c_i>0，设λ_i'=λ_i−tc_i。
则所有λ_i'≥0，仍表示同一x、总和1，至少一项归零。
重复消元得到至多d+1个点。此结论不自动证明无限维紧凸集的Krein–Milman定理。

局部证据补充见 [SOL_HN记录](audit-sol_hn.md)；其审读者与范围以记录为准。

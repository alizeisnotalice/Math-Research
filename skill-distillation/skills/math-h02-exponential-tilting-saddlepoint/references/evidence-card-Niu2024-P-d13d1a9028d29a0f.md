# 独立阅读卡：Niu, Ray Choudhury & Katsevich (2024)

**Paper ID:** `P-d13d1a9028d29a0f`  
**身份/版本：** Ziang Niu, Jyotishka Ray Choudhury, Eugene Katsevich, “The saddlepoint approximation for averages of conditionally independent random variables,” arXiv:2407.08915v3; 64 页。v3 脚注称本文被纳入较长稿 arXiv:2407.08911。本轮读的是 v3。  
**附件映射：** H02 E0215 与 E0255 为同 SHA 重复附件；计作同一篇。  
**SHA-256:** `d13d1a9028d29a0f21d767d68ecde717a0300ece7938472782cbb15b7a332ea9`  
**覆盖：** PDF p.1–64、全文转换 lines 1–3575；逐页含附录 A–H。PDF 公式视觉回看 pp.5–7（定义、定理、LR 参数）、28–29（Mills/条件 Berry–Esseen 辅助式）、61–64（Appendix H 关键极限和常数）。身份核对 arXiv 题录。

## 定理的逐项范围

对每行条件于 `F_n` 独立、条件均值 0 的变量 `W_{in}`，令
\[
K_{in}(s)=\log E(e^{sW_{in}}\mid F_n),\qquad K_n(s)=n^{-1}\sum_iK_{in}(s).
\]
Theorem 1 假设以下条件全部成立：

1. **CSE 分支：**条件尾 `P(|W_{in}|≥t|F_n)≤θ_n e^(−βt)`，固定 `β>0`，`θ_n` 为 `F_n`-可测、有限 a.s. 且 `θ_n=O_P(1)`；**或 CCS 分支：**`|W_{in}|≤ν_{in}` a.s.，`ν_{in}` 为 `F_n`-可测，且 `n^(−1)Σ_iν_{in}^4=O_P(1)`。
2. 平均条件方差 `n^(−1)Σ_i E(W_{in}²|F_n)=Ω_P(1)`。
3. 阈值 `w_n` 为 `F_n`-可测且 `w_n=o_P(1)`。

此时以概率趋于 1 存在固定邻域内唯一鞍点根 `K_n′(ŝ_n)=w_n`。置总 CGF `Λ_n=nK_n`、总和阈值 `x_n=nw_n`，
\[
\lambda_n=\hat s_n\sqrt{nK_n''(\hat s_n)},\qquad
r_n=\operatorname{sgn}(\hat s_n)\sqrt{2\{\hat s_nx_n-\Lambda_n(\hat s_n)\}}.
\]
则条件尾满足
\[
P\!\left(n^{-1}\sum_iW_{in}\ge w_n\mid F_n\right)
=\left[1-\Phi(r_n)+\phi(r_n)\left(\lambda_n^{-1}-r_n^{-1}\right)\right](1+o_P(1)).
\]
这是完整 LR 近似的相对误差式，不是有限样本概率界；未给明确误差常数或最终误差率。原文的阈值条件排除固定非零 cutoff 大偏差；一般依赖、置换检验依赖及 studentized 统计量也不由此覆盖。Theorem 1 本身没有额外非格点假设。

## 简化前因子是附加推导

只考虑正上尾，要求 `w_n>0` 以概率趋于 1，因此在唯一根事件上 `ŝ_n>0,r_n>0`。Appendix H（尤其 pp.61–64）给出 `λ_n/r_n→_P1`。定义 `M(r)=rQ(r)/φ(r)`，`Q=1−Φ`。将完整 LR 表达式除以 `φ(r_n)/λ_n` 得精确恒等式
\[
\frac{Q(r_n)+\phi(r_n)(\lambda_n^{-1}-r_n^{-1})}
     {\phi(r_n)/\lambda_n}
=1+\frac{\lambda_n}{r_n}(M(r_n)-1).
\]
对 `r>0`，`r²/(r²+1)<M(r)<1`。因此 LR/`(φ/λ)` 在 1 以下，亏损不超过 `λ_n/[r_n(r_n²+1)]`。若再有 `r_n→_P∞`，上式比值趋于 1，且
\[
\frac{\phi(r_n)}{\lambda_n}
=\frac{\exp\{\Lambda_n(\hat s_n)-\hat s_nx_n\}}
 {\hat s_n\sqrt{2\pi\Lambda_n''(\hat s_n)}}.
\]
所以**在该定理的模型和 shrinking-cutoff 假设内**，并附加正尾和 `r_n→∞`，简化前因子是相对渐近等价式。不能把这一步独立提升为全范围鞍点定理，也不能称作上/下界；相对余项率未量化。

## 文中局部常数问题

Appendix H p.62 公式 (H.23) 声称对 `r≥1`，`r e^(r²/2)Q(r)≥1/2`。该常数不正确（例如 `r=1` 左边约 `0.261`）。从本卡的 Mills 下界可以得到
\[
r e^{r^2/2}Q(r)>\frac{r^2}{\sqrt{2\pi}(r^2+1)}\ge\frac1{2\sqrt{2\pi}},\qquad r\ge1.
\]
该处证明只需要一个绝对正数常数，因此以 `1/(2√(2π))` 替换后局部论证仍然成立。记录为局部笔误/修复，不据此否定整篇主定理。

## 阅读中的其它可复核信息

- Appendix H Lemma 29 控制 `λ_n/r_n→_P1`，但该结论依赖 Theorem 1 的 shrinking-cutoff 设定。
- PDF p.28 的标准 Gaussian Mills 计算给出上尾前因子方向；其对象是 Gaussian `Q(u)`，不支持一般倾斜分布的有限样本相同不等式。
- 条件 Berry–Esseen 辅助结果与 tilted product law 是主证明的组成部分；Theorem 1 最终仍只给 `o_P(1)`，没有显式速率。
- 本卡数学结论/外推界限与数值例子分开登记；实例运行不能替代上述证明。

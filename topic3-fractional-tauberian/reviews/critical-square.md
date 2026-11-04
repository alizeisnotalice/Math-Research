# 中间临界点两层一致过渡：独立复审

结论：`fix:critical-uniform` 的三项均正确。复审限于组合四阶式、末端截断面积和内层
\(Q(\mu)\) 的分段；没有重复运行已有符号或数值计算。

## 组合四阶式

令 \(s=tu\)、\(\beta=1+\lambda t\)。直接的统一代数展开为
\[
 w_B=t(1-2u)-2\lambda u t^2+2\lambda(u-u^2)t^3+O_\Lambda(t^4),
\]
\[
 w_C=t(1-2u)+2(1+\lambda)u t^2-2\lambda(u-u^2)t^3+O_\Lambda(t^4).
\]
两式平均给出 \(t(1-2u)+ut^2\)，半差给出
\(ut^2[1+2\lambda-2\lambda t(1-u)]\)。因为有限个输入取最大值是 Lipschitz 映射，
在 \(0\le u\le1/2\) 上
\[
 \max\{w_B,w_C\}=t(1-2u)+ut^2+
 ut^2|1+2\lambda-2\lambda t(1-u)|+O_\Lambda(t^4)
\]
确实一致。乘以面积 Jacobian \(8t\,du\) 后，主体的误差为 \(O_\Lambda(t^5)\)，
而主体三阶项正是 \(t^3[1+8I(\lambda,t)]\)。这里没有遗漏共同的三阶项：
\(w_B,w_C\) 的三阶项互为相反数，已进入绝对值中的 \(t\)-修正。

## 末端截断

在 \(u=1/2+O_\Lambda(t)\)，最大候选的主形为
\[
 t^{-1}\max\{w_B,w_C\}=-2(u-1/2)+m_\lambda t+O_\Lambda(t^2),
 \qquad m_\lambda=\max\{-\lambda,1+\lambda\}.
\]
其正部三角形在 \((u,t^{-1}w)\) 平面的面积为
\(m_\lambda^2t^2/4+O_\Lambda(t^3)\)。乘回 \(8t^2\) 后给出
\(2m_\lambda^2t^4+O_\Lambda(t^5)\)，与正文完全一致。特别地，
\(m_\lambda\ge1/2\) 保证在 \(u=1/2\) 处最大值为正，故主体区间没有隐藏的内侧截断。

## 内层 \(Q(\mu)\)

令 \(\lambda=-1/2+\mu t\)、\(a=1+2\mu\)。则
\[
 I(\lambda,t)=t\int_0^{1/2}u|a-u|\,du+O_M(t^2),
 \qquad 2m_\lambda^2=\frac12+O_M(t),
\]
对 \(|\mu|\le M\) 一致。因此
\[
 Q(\mu)=\frac12+8\int_0^{1/2}u|a-u|\,du.
\]
在 \(a\le0\)、\(0\le a\le1/2\)、\(a\ge1/2\) 三段分别积分，得到
\[
 Q(\mu)=
 \begin{cases}
 -\frac16-2\mu,&\mu\le-\frac12,\\
 \frac56-a+\frac83a^3,&-\frac12\le\mu\le-\frac14,\\
 \frac76+2\mu,&\mu\ge-\frac14,
 \end{cases}
\]
边界值连续：在 \(\mu=-1/2\) 均为 \(5/6\)，在 \(\mu=-1/4\) 均为 \(2/3\)。
分段与绝对值零点是否落在积分区间完全对应，无遗漏。

## 建议

数学上可通过。为使“一致”不被误读，内层公式的余项应写成 \(O_M(t^5)\) 或在
公式紧邻处明确常数依赖于预先固定的 \(|\mu|\le M\)；正文文字虽已说明“有界 \(\mu\)”，
但显示公式中的裸 \(O(t^5)\) 不够精确。

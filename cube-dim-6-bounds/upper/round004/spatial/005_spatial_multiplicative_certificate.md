# 第004轮：截时势的空间乘性盒证书

对于一维尺度 \(t\) Cauchy 密度
\(f_t(z)=t/[\pi(t^2+z^2)]\)，有
\[
 \left|{d\over dz}\log f_t(z)\right|={2|z|\over t^2+z^2}\le {1\over t}.
\]
积分给任意实 \(a,z\)
\[
 e^{-|a|/t}f_t(z)\le f_t(z+a)\le e^{|a|/t}f_t(z). \tag{1}
\]

令 \(V\) 为32条 body 方向代表，\(A(c)=\sum_{v\in V}c_vv\)。其系数 \(\ell^1\)
gauge 满足
\[
 \inf_{A(c)=h}\sum_v|c_v|=\|h\|_\infty. \tag{2}
\]
理由是允许系数取负后，单位球是 \(\operatorname{conv}(\pm V)=[-1,1]^6\)。
一边由三角不等式，另一边由 cube 的顶点凸组合得到；最小值实际达到。

在32维参数积分中作平移 \(z\mapsto z+c\)，(1)--(2) 给真实 body32 密度
的逐点比较
\[
 e^{-\|h\|_\infty/t}p_t(x)\le p_t(x+h)
 \le e^{\|h\|_\infty/t}p_t(x). \tag{3}
\]
可先对测试函数证明，再由连续密度升级为逐点式；没有 Jacobian 因子，因为参数
平移保 Lebesgue 测度。

对 \(H_\epsilon=\int_\epsilon^\infty p_tdt\)，由 \(t\ge\epsilon\) 有
\[
 e^{-\|h\|_\infty/\epsilon}H_\epsilon(x)
 \le H_\epsilon(x+h)
 \le e^{\|h\|_\infty/\epsilon}H_\epsilon(x). \tag{4}
\]
这在 \(S_2\) 及其交叠上也有效：\(H_\epsilon\) 是有限连续正函数，未取
\(\epsilon\downarrow0\)。

若一个 \(\ell_\infty\)-box 的中心有准确区间 \([L,U]\)，半径为 \(r\)，则全 box
上 \(H_\epsilon\in[e^{-r/\epsilon}L,e^{r/\epsilon}U]\)。因此：

* cube 下界可由所有 box 的 \(e^{-r/\epsilon}L\) 最小值给出，随后用
  \(G\ge H_\epsilon\)；
* 正 face 积分上界可用 \(\sum_B|B|(e^{r_B/\epsilon}U_B)^{6/5}\)，再以
  第002轮 whole-face 短时范数预算加入；
* 不需要为 \(S_2\) 另设 pointwise 帽或从其补集延拓。

这是一套严格连续证书，但均匀网格昂贵。取 \(r=\epsilon/4\) 时，乘性损失仅
为 \(e^{\pm1/4}\)。cell side 是 \(2r=\epsilon/2\)：cube 每轴需
\(\lceil4/\epsilon\rceil\) 个、正 face 每轴需 \(\lceil2/\epsilon\rceil\) 个。
实际 sizing 表显示，\(\epsilon=1/20\) 也需要约 \(2.62\cdot10^{11}\) 个 cube
cell 和 \(1.024\cdot10^8\) 个正 face cell；更小 \(\epsilon\) 更不可行。故它不应替代
全域 uniform grid。

最有前景的有限用法是：将 (4) 用于少数已认证 node 的局部 box 扩张、候选极小
点的排除或 singular-line 邻域的统一覆盖；其余平滑区域仍需 Chebyshev 或另一套
有完整误差账的自适应方法。任何自适应 tree 也必须记录 box 覆盖、中心区间、
半径和指数因子，不能仅按数值采样丢弃 leaf。

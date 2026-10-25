# 第 002 轮五维 body fan 独立审计

本审只读数值分支的七个排序 chamber 锥和锚点程序。它重新以整数 Leibniz 行列式、Fraction 截面体积和所有 16 个 body 线性型计算，不导入生产检查器。

## 覆盖结论

七个锥都包含于 \(a_0\ge a_1\ge\cdots\ge a_4\ge0\)，且每一个 body 线性型在每锥上同弱号。去掉始终正的 \((1,1,1,1)\) 型后，七个锥的严格符号签名两两不同；所以它们的内部没有正体积重叠。对截面 \(\sum_i a_i=1\)，锥的相对体积为
\[
 \frac{|\det R|}{\prod_j {\bf1}\cdot r_j}.
\]
七项和和排序 chamber 的值均为 \(1/120\)。因此不交内点的锥并在 chamber 中具有满截面体积，给出 a.e. 覆盖。边界重叠为有限超平面并，对 Fourier 积分没有影响。

该论证足以支持锥内符号线性化、Jacobian 和 a.e. Fourier 锥分解；它不认证截时积分或连续 \(g/I\) 证书。

## 峰值常数

在原点，完整 signed-permutation 轨道有 \(2^5\,5!=3840\) 个相等贡献。七个基本锥的精确和为
\[
 \sum_{C\ \mathrm{base}}\frac{|\det R_C|}{\prod_jq_{Cj}}
 =\frac{35509}{91729428480}.
\]
因此
\[
 p_1(0)=\frac{35509}{764411904\pi^5},\qquad
 p_t(0)=\frac{35509}{764411904\pi^5t^5}.
\]
这严格替代较粗的 \((48\pi^5t^5)^{-1}\) 上界，可用于后续复管常数；不把浮点锚点当作该身份的证据。

## 截时 DD4 接口

令 \(b_j=x\cdot r_j/q_j\)，\(f_\epsilon(b)=\tfrac12\log(\epsilon^2+b^2)\)。对五个锥射线，正确公式的符号是
\[
 H_\epsilon(x)=-(2\pi)^{-5}
 \sum_C\frac{|\det R_C|}{\prod_jq_{Cj}}
 f_\epsilon[b_1,b_2,b_3,b_4,b_5].
\]
这里方括号是四阶实 divided difference。因为
\[
 \partial_\epsilon f_\epsilon[b_1,\ldots,b_5]
 =\Re\prod_{j=1}^5(\epsilon-i b_j)^{-1},
\]
而四阶差商在 \(\epsilon\to\infty\) 趋零，时间尾积分带一个负号。生产实现不得误用三维 DD2 的正号。

任何重复或零节点在 \(\epsilon>0\) 时必须按 Hermite 差商连续延拓；重数 \(m\) 的单节点值为 \(f_\epsilon^{(m-1)}(b)/(m-1)!\)。本审脚本对互异有理节点和五重节点分别作了精确符号核查。\(\epsilon\downarrow0\) 前必须先合并全部锥的 principal parts，不能逐锥删去零节点。

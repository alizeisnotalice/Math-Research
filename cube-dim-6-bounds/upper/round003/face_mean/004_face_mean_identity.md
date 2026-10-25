# 第003轮：body32 face 均值的逆范数恒等式

令 \(X_1=\sum_{v\in V}vC_v\)，其中 \(C_v\) 是独立尺度一的标准 Cauchy
变量，\(p_1\) 是其密度，\(R=\|X_1\|_\infty\)。body32 Green 势齐次次数为
\(-5\)。令
\[
 M=\int_{[0,1]^5}G(1,y)dy,\qquad
 J=\int_{\partial[-1,1]^6}G(y)dS(y).
\]
则严格有
\[
 \boxed{\quad M={1\over384}\mathbb E[R^{-1}].\quad} \tag{1}
\]

证明：由 Tonelli、\(p_t(x)=t^{-6}p_1(x/t)\) 和 \(r=1/t\)，
\[
 J=\int_{\partial Q}\int_0^\infty t^{-6}p_1(y/t)dtdS(y)
  =\int_{\partial Q}\int_0^\infty r^4p_1(ry)drdS(y).
\]
cube-radial 变换 \(x=ry\) 的 Jacobian 是 \(r^5drdS(y)\)，所以右端为
\(\int_{\mathbb R^6}p_1(x)/\|x\|_\infty dx=\mathbb E[R^{-1}]\)。没有额外的
径向6因子。body frame 的坐标置换与反射对称把边界分成6个坐标方向、每个2个
符号面、每个面32个正交象限，故 \(J=6\cdot2\cdot32M=384M\)。

期望有限：原点邻域以 \(p_1\in L^\infty\) 和 cube-radial 体积元控制
\(R^{-1}dx\sim r^4drdS\)，远离原点则 \(R^{-1}\le1\)。所以 (1) 不需要把
未知发散量相减，也不忽略边界的零测奇异线。

令 \(q=6/5\)。由于 face 体积为1，Jensen 的方向为
\[
 I=\int_{[0,1]^5}G(1,y)^qdy\ge M^q. \tag{2}
\]
若另有 cube 点严格上界 \(g\le U\)，则 body32 正部成本至少为
\[
 6^q{I\over g^q}\ge6^q(M/U)^q. \tag{3}
\]
要由均值单独排除成本 \(<64\)，需要
\[
 {M\over U}>{64^{5/6}\over6}={16\over3}. \tag{4}
\]
这里没有 Jensen 上界；均值本身也不能给 \(I\) 的上界。

## 有限可行性诊断

本轮实际运行四个独立 scramble 的 Sobol 批次，每批 \(2^{15}\) 个32维 Cauchy
参数；只计算 \(R^{-1}\)。得到
\(\mathbb E[R^{-1}]\) 的 replicate 均值约 \(0.02148828\)，相应
\(M\) 诊断值约 \(5.59591\cdot10^{-5}\)。若未来 Laurent 证书给
\(U\approx3.49\cdot10^{-5}\)，则 \(M/U\) 的诊断值约1.60，低于 (4) 的
5.333。它表明“只认证全 face 均值再 Jensen”的路线在该 \(U\) 量级下看起来
不足；这不是严格反证，也不排除多盒 \(L^{6/5}\) 下界或其他机制。

最小严格批次应是：以完整 cone 截时公式认证 \(\int_{[0,1]^5}H_\epsilon\) 的
下界，先与 \((16/3)U\) 比较。若仍不通过，就停止均值路线，不扩展到大网格；
若通过，再用 whole-face 短时预算和 Laurent 点上界完成证书。

## 全时间 Minkowski 预算的独立审核

对正 face，slice 的准确 \(L^1(\mathbb R^5)\) 边际是
\(32t/[\pi(1+1024t^2)]\)。与全 slice 的峰值界插值时会出现
\(32^{5/6}\)。但 body frame 对五个自由坐标逐个反射不变，故
\(p_t(1,y)\) 在32个正交象限上的 \(L^{6/5}\) 质量相同；限制到
\([0,\infty)^5\) 的范数精确乘以 \(32^{-5/6}\)。二者完全抵消。

因此 primary 的正-face 全时间预算积分
\[
 B=\int_0^\infty\min(C_0t^{-6},C_ft^{-4})^{1/6}
 \left[{t\over\pi(1+1024t^2)}\right]^{5/6}dt
\]
常数方向正确，其中 \(C_0\) 是全局峰值常数、\(C_f\) 是正 face 峰值常数。
Minkowski 给的是 \(\|G(1,\cdot)\|_{L^{6/5}([0,1]^5)}\le B\)，不是该真实
范数的下界。若 \(g\le U\)，任何只使用这个预算来证明成本小于64的方法都
必须有 \(B<16g/3\le16U/3\)。所以 primary 的严格 \(B>16U/3\) 仅排除
这套 \(L^\infty\)--\(L^1\) 插值加全时间 Minkowski 预算；它不排除 body32
模型，也不能与上节 Jensen 下界混用。

此前我曾暂时担心该预算漏掉 \(32^{5/6}\)，经上述正象限因子复核后已纠正；
本报告以本节为最终口径。

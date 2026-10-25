# 第004轮补充：二阶截时势比较与顶点证书

令 \(A:\mathbb R^{32}\to\mathbb R^6\) 为 \(Ac=\sum_vc_vv\)。Gram 恒等式
\(AA^T=32I\) 给正交右逆 \(R=A^T/32\)。在 coarea 表示
\[
 p_t(x)=C_A\int_{\ker A}\prod_v f_t((Rx+w)_v)\,dw \tag{1}
\]
中，\(\|Rh\|_2^2=\|h\|_2^2/32\)，且 \(Rx\perp w\)。

## log-半凸：face 上界

\(\log f_t\) 的二阶导不小于 \(-2/t^2\)。所以
\(\prod_vf_t(z_v)e^{\|z\|_2^2/t^2}\) log-convex。将 (1) 乘
\(e^{\|x\|^2/(32t^2)}\)，积分内剩下固定正权 \(e^{-\|w\|^2/t^2}\)。每个
固定 \(w\) 的 integrand 对 \(x\) log-convex，Hölder 表明积分仍 log-convex。
故
\[
 \nabla^2\log p_t\succeq-{I\over16t^2}.
\]
对 \(t\ge\epsilon\)，再乘非负二次指数并作正时间混合，得到
\[
 H_\epsilon(x)e^{\|x\|^2/(32\epsilon^2)}\ \text{log-convex}. \tag{2}
\]

在五维 face box \(B=\prod_i[c_i-r_i,c_i+r_i]\)，令 \(U_v\) 为32个顶点的
严格上界。由 (2)、顶点 barycentric 权和 AM--GM，
\[
 \int_BH_\epsilon^{6/5}
 \le |B|\exp\left({(6/5)\sum_ir_i^2\over32\epsilon^2}\right)
 {1\over32}\sum_{v\in\mathcal V(B)}U_v^{6/5}. \tag{3}
\]
这只给积分上界，不能给 \(g\) 的下界。若 \(r_i=\epsilon\)，指数损失为
\(e^{3/16}\)，正 face cell side 为 \(2\epsilon\)，覆盖数
\(\lceil1/(2\epsilon)\rceil^5\)。

## 坐标二阶 ODE：cube 下界

固定坐标方向 \(e_i\)，取 \(c=R e_i=V^Te_i/32\)。每个 \(|c_v|=1/32\)，
\(\sum c_v^2=1/32\)，\(\sum|c_v|=1\)。并且
\[
 (\log f_t)''\le {1\over4t^2},\qquad |(\log f_t)'|\le {1\over t}.
\]
对 (1) 内的参数平移求二阶导，逐点有
\[
 \partial_{ii}F\le\left[{1\over t^2}+{1\over128t^2}\right]F.
\]
由正积分和 \(t\ge\epsilon\)，
\[
 \partial_{ii}H_\epsilon\le k_\epsilon^2H_\epsilon,\qquad
 k_\epsilon={\sqrt{129/128}\over\epsilon}. \tag{4}
\]
若一维截面端点为 \(a<b\)，最大值原理对 \(f''-k^2f\le0\) 给
\[
 f(x)\ge{\sinh(k(b-x))f(a)+\sinh(k(x-a))f(b)\over\sinh(k(b-a))}. \tag{5}
\]
逐坐标施用正的一维算子，给五维 box 的顶点下界。所有顶点至少为 \(L\)、半径
同为 \(r\) 时，整 box 有统一下界
\[
 H_\epsilon\ge L\cosh(k_\epsilon r)^{-5}. \tag{6}
\]
这和 (3) 可在同一 cell 分别服务 cube \(g\) 下界和 face \(I\) 上界。它不是
对 \(G\) 的上界；cube 侧只用 \(G\ge H_\epsilon\)。

## 复参数相对管：可合法但须局部使用

对实 \(z\)、\(|a|<t\)，直接平方分母可验证
\[
 |f_t(z+ia)|\le {f_t(z)\over1-(a/t)^2}. \tag{7}
\]
式 (1) 使这不是形式 delta 代换：若 \(x\mapsto x+ih\)，直接在 coarea
积分中把 \(Rx\) 改成 \(Rx+iRh\)。当每个 \(|(Rh)_v|<t\) 时，(7) 给统一
可积支配，从而由支配收敛定义真实的解析延拓并有相对界。

特别对单轴 \(h=se_i\)，\(|(Rh)_v|=|s|/32<t\)，故
\[
 |p_t(x+ise_i)|\le p_t(x)(1-s^2/(1024t^2))^{-32}. \tag{8}
\]
对 \(H_\epsilon\) 取 \(|s|<32\epsilon\) 并用 \(t\ge\epsilon\) 即得同一相对
管界。多自由坐标同时复移时应逐个检查 \(|(Rh)_v|<\epsilon\)；不能把单轴
宽度32\(\epsilon\)无条件套到五维 max-norm 管。该管可改善局部 Chebyshev 常数，
但本轮未据此启动连续网格。

用于全cube的G下界时，只需先覆盖正face：对0<rho≤1，由G(rho y)=rho^-5G(y)≥G(y)，并利用坐标反射/置换可扩至整个cube。不能将这个齐次关系误用于固定epsilon的H_epsilon。

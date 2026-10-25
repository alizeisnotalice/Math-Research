# 五维面短时常数的商空间峰值证书

固定大方向 \(v=(1,1,1,1,1)\)，并用坐标 \(u_j=x_j-x_1\)（\(j=2,\ldots,5\)）描述商空间；\((x_1,u)\mapsto x\) 的 Jacobian 为一。其余十五条 body 方向的商像为
\[
 -2\mathbf1_S,\qquad\varnothing\ne S\subseteq\{1,2,3,4\}.
\]
故商过程的谱函数为
\[
 \psi_\perp(\xi)=2\sum_{\varnothing\ne S\subseteq[4]}
 \left|\sum_{j\in S}\xi_j\right|. \tag{1}
\]

设 \(B=\{\psi_\perp\le1\}\)。因为 \(e^{-\psi_\perp}\ge0\)，商密度在零点取得全局上确界。齐次积分给
\[
 C_\perp=p^\perp_1(0)=(2\pi)^{-4}4!|B|. \tag{2}
\]

`round002_quotient_vertices.py` 对十五条 form 的全部 455 个三元组做准确秩检查。430 个秩三组给出全部 90 个顶点：顶点的充要条件是零 form 秩为三，再加非零线性 gradient 为四维 active normal span。脚本还枚举所有 1,611 个 support vector，并以精确等式与仿射秩得到 370 个 facet（350 个四顶点、20 个五顶点）。

五顶点 facet 的唯一仿射 circuit 的 Radon 两部分大小为 \(2\) 和 \(3\)。删除其中任一固定一侧的各顶点，得到标准无重叠四面体剖分；脚本选择最小局部指标更小的一侧。于是全部边界有
\[
350+8\cdot2+12\cdot3=402
\]
个精确单形。逐个取以原点为顶点的行列式体积，得到
\[
 |B|=\frac{14455}{191102976}. \tag{3}
\]
所有顶点、facet、circuit 和单形体积均使用整数或 `Fraction` 运算；脚本断言顶点满足全部 1,611 个 halfspace、仿射系数确实湮灭所有五个齐次点，并断言 (3)。

由 (2)--(3)，
\[
 C_\perp=\frac{14455}{127401984\pi^4}. \tag{4}
\]
若第一个 Cauchy 参数足够大以命中 \(x_1=1\)，其密度至多 \(256t/\pi\)。给定其余变量，沿该参数积分后只剩商密度；商密度由 (4) 上界。对十六个大变量事件并集，
\[
 p_t(1,y)\le C_Ft^{-3},\qquad
 C_F=\frac{4096}{\pi}C_\perp
 =\frac{14455}{31104\pi^5}. \tag{5}
\]
第一坐标边际仍为 \(16t/[\pi(1+256t^2)]\le16t/\pi\)。所以第 001 轮的 \(q=5/4\) 插值和 Minkowski 论证保持成立，只需把旧 \(C_F\) 换成 (5)：
\[
 \|T_\varepsilon(1,\cdot)\|_{L^{5/4}([0,\infty)^4)}
 \le16^{-4/5}\frac56
 \left(\frac{14455}{31104\pi^5}\right)^{1/5}
 \left(\frac{16}{\pi}\right)^{4/5}\varepsilon^{6/5}. \tag{6}
\]
相对第 001 轮的面上确界常数，(5) 的比值为
\(14455/2654208\)。这只改进短时误差预算；尚未给出截时 Green 的连续 \(g/I\) 证书或新的 \(C_5\) 上界。

## 最小复核

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper/round002_quotient_vertices.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper/round002_quotient_volume_audit.py
```

# N01 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

S-0c1c25fd019e7c3c JohnEllipsoid lines44–52、273–276、350–351；R^n满维compactconvex，Johnmaxinscribed与Löwnerminoutside不同；sym度量√n/generaln。

## 证明骨架与已检查推导

Affine map JohnE→unitball→contactmoment∑ci ui=0、∑ci uiuiᵀ=I（sym可简化）→support/凸体包含；variationproof同时check可行与行列式增加。

本次实际检查：cube[-1,1]^n containsB2且eachvertexnorm√n，factor√n；平面regulartriangleinradius.5/circumradius1违√2，generalfactor2紧。

## 正反例与失败处理

正例：对 K=[-1,1]^n，最大体积内接椭球为 B_2^n，且 K⊂sqrt(n)B_2^n，因为 ||x||_2<=sqrt(n)||x||_infinity<=sqrt(n)。

条件缺失例：对任意非对称凸体直接套用 K⊂sqrt(n)E；一般 John 包含因子是 n，sqrt(n) 需要中心对称。

## 中心立方体迁移推断

用椭球归一化凸几何参数，或在体积/容量估计前将凸体夹在 Euclidean 球之间。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

√n isgeo度量containment，不给cube maximal 弱型常数；affinechange须另外checkfamily与质量的 Jacobian 因子。

对称性和内接/外接约定会改变常数。John 域或 Loewner 演化的同名文献并未证明椭球结论。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。

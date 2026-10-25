# 第 011 轮：完整面内空间消元与严格 I 上界

## 结论与证据边界

本轮固定五维 16 条未单位化 body 方向、相互独立的 \(\alpha=3/2\) stable 分量，
以及 Lebesgue 测度。沿用 Astra upper001 已证明的正形状表示后，完整消去四个面内空间坐标，证明
\[
 I:=\int_{[0,1]^4}G(1,y)^{10/7}dy
 \le c^{10/7}\frac{\pi^2}{60}\mathbb E_R\det(M)^{-3/14}. \tag{1}
\]
这里的 \(\pi^2/60\) 比单独扩张到整平面的 \(\pi^2/12\) 小五倍；这个因子来自全部
立方体面的精确径向投影与坐标置换对称。它没有使用面内单调性、角点最小原则或数值采样。

一次预登记的有限 Python/256-bit Arb 计算认证
\[
 \mathbb E_R\det(M)^{-3/14}\le Q_*=2.4576305439372284\ldots,
\]
从而得到本轮实际的完整上界
\[
 \boxed{I\le1.9382520037464934\ldots\times10^{-5}<\frac1{50000}.} \tag{2}
\]
这是一条全时间、全正面的严格积分上界；无须另加短时间或晚时间误差。
它仍高于目标 \(11\cdot10^{-6}\)，且本轮没有证明所需的全域 \(g\) 或改进 \(C_5\le32\)。

## 1. 正形状表示及对称性

令 \(\beta=3/4\)，取 16 个独立正 stable 变量
\(\mathbb E e^{-\lambda S_v}=e^{-\lambda^\beta}\)。记
\[
 T=\sum_v S_v,\quad P_v=S_v/T,\quad
 M=\sum_vP_vvv^T,\quad v=(1,\pm1,\pm1,\pm1,\pm1).
\]
\(M\) 几乎处处正定，每个对角元为一，\(\operatorname{tr}M=5\)。把原概率用
\(T^{-\beta}\) 倾斜并推送到 \(P\)，所得概率记为 \(R\)。归一化与 Green 表示为
\[
 \mathbb ET^{-\beta}=\frac1{16\Gamma(7/4)},\qquad
 G(x)=c\mathbb E_R H(M,x), \tag{3}
\]
\[
 c=\frac3{2^{15/2}\pi^{5/2}},\qquad
 H(M,x)=\det(M)^{-1/2}(x^TM^{-1}x)^{-7/4}. \tag{4}
\]
这些是 Astra upper001 报告 (6)--(7) 的正表示。其原始条件 Gaussian 核包含
\(\exp[-x^TA^{-1}x/(4t^{4/3})]\)，先对 \(t\) 作 Gamma 积分，再消去 \(T\) 即得 (3)。
因此 \(x_1=1\) 的早时间离原点衰减已被精确积分，未被 \(p_t(0)\) 替换。
等式允许沿方向线取 \(+\infty\)，适合 Tonelli。

\(vv^T\) 的集合在全部五个坐标的符号变化与置换下闭合：若变换后的首坐标为负，
将整个方向改号不改变 \(vv^T\)。独立同分布的 \(S_v\) 因而被置换，\(T\) 不变，
故 \(R\) 的矩阵分布也在这些共轭变换下不变。不能只证明保持第一坐标的反射对称；
本轮的 \(1/5\) 因子还使用了交换第一坐标与其余坐标的对称。

## 2. 临界幂次的精确球面/立方体恒等式

置 \(q=10/7\)、\(d=\det M\)、\(b=3/14\)。对固定正定 \(M\)，
\[
 H(M,x)^q=d^{-5/7}(x^TM^{-1}x)^{-5/2}. \tag{5}
\]
这是 \(-5\) 次齐次函数。由五维 Gaussian 积分或椭球体积的极坐标公式，
\[
 \int_{S^4}(\theta^TM^{-1}\theta)^{-5/2}d\omega(\theta)
 =\frac{8\pi^2}{3}d^{1/2},
\]
所以
\[
 \int_{S^4}H(M,\theta)^q d\omega
 =\frac{8\pi^2}{3}d^{-3/14}. \tag{6}
\]
例如可用 \(\int_{\mathbb R^5}e^{-x^TM^{-1}x}dx=\pi^{5/2}d^{1/2}\) 核对全部常数：
径向积分是 \(\Gamma(5/2)/2\)，故球面积系数为 \(2\pi^{5/2}/\Gamma(5/2)=8\pi^2/3\)。

对立方体每个面 \(F\subset\partial[-1,1]^5\)，令 \(\theta=x/|x|\)。径向投影的
Jacobian 为
\[
 d\omega=\frac{x\cdot n_F}{|x|^5}d\sigma(x)=|x|^{-5}d\sigma(x),
\]
因为 \(x\cdot n_F=1\)。乘以 (5) 的齐次缩放因子 \(|x|^5\) 后恰好抵消。因此对每个固定 \(M\)，
\[
 \sum_{F\text{为十个面}}\int_F H(M,x)^q d\sigma(x)
 =\frac{8\pi^2}{3}d^{-3/14}. \tag{7}
\]
面交界的径向像为球面零测集，不产生额外项。

对 (7) 用 Tonelli 取 \(R\) 期望，再用第 1 节的完整 signed-permutation 对称。
十个面的每个正交象限，共 \(10\cdot16=160\) 块，其期望积分相同，故得到精确式
\[
 \int_{[0,1]^4}\mathbb E_R H(M,(1,y))^qdy
 =\frac{\pi^2}{60}\mathbb E_Rd^{-3/14}. \tag{8}
\]
最后只在形状期望上使用一次 Jensen：\((\mathbb E_RH)^q\le\mathbb E_RH^q\)。
结合 (3) 即得 (1)。下面的负矩上界同时证明 (8) 右端有限，因而无需预先假定其可积。

作为条件正积分接口，若将 \(M\) 写成
\(\begin{pmatrix}1&a^T\\a&C+aa^T\end{pmatrix}\)，则
\[
 H(M,(1,y))^q=\det(C)^{-5/7}
 [1+(y-a)^TC^{-1}(y-a)]^{-5/2}.
\]
这本来是一个四维椭圆 Cauchy 正积分；式 (7)--(8) 利用全部面和形状对称把它完整消掉。

## 3. 负行列式矩的全域可证 Gamma 上界

令 \(w_B=\det(V_B)^2\)，只对非退化的五方向子集求和。Cauchy--Binet 给
\[
 d=\sum_Bw_B\prod_{v\in B}P_v.
\]
frame 恒等式 \(\sum_vvv^T=16I\) 给 \(W:=\sum_Bw_B=16^5\)。方向反射群在 16 个方向上传递，
所以每个方向的包含权重相同；其总和为 \(5W\)，故
\(\sum_{B\ni v}w_B/W=5/16\)。对 \(w_B/W\) 作加权 AM--GM 得到全 simplex 内部的不等式
\[
 d\ge16^5\prod_{v=1}^{16}P_v^{5/16}. \tag{9}
\]
此式在 \(P_v=1/16\) 同时取等，利用了全部 Cauchy--Binet 权重；本轮不重枚举旧五基。

固定
\[
 b=\frac3{14},\quad a=\frac{5b}{16}=\frac{15}{224},\quad
 \gamma=5b-\beta=\frac9{28}\in(0,\beta).
\]
由 (9)、倾斜概率定义和 \(P_v=S_v/T\)，
\[
 \mathbb E_Rd^{-b}
 \le\frac{16^{-5b}}{\mathbb ET^{-\beta}}
 \mathbb E\left[T^\gamma\prod_{v=1}^{16}S_v^{-a}\right]. \tag{10}
\]
在原独立概率下，\(T^\gamma\) 对每个 \(S_v\) 递增，\(\prod S_v^{-a}\) 对每个坐标递减，
因此独立坐标的反向关联不等式给
\[
 \mathbb E\left[T^\gamma\prod_vS_v^{-a}\right]
 \le\mathbb ET^\gamma\left(\mathbb ES^{-a}\right)^{16}. \tag{11}
\]
这里不是假定 \(T\) 与 \(P\) 独立。一个自足证明是：单坐标有
\(\operatorname{Cov}(f(X),h(X))=\frac12\mathbb E[(f(X)-f(X'))(h(X)-h(X'))]\le0\)
当两函数单调方向相反；按独立坐标逐次条件化并归纳即可。先对两函数截断后应用，再单调取极限，
即可处理无界情形。下列显式矩有限，因而右端确为有限上界。

由 Laplace 变换的 Gamma/Mellin 积分，
\[
 \mathbb ES^{-a}=\frac{\Gamma(1+a/\beta)}{\Gamma(1+a)},\qquad
 \mathbb ES^\gamma=\frac{\Gamma(1-\gamma/\beta)}{\Gamma(1-\gamma)},\qquad
 \mathbb ET^\gamma=16^{\gamma/\beta}\mathbb ES^\gamma. \tag{12}
\]
负矩式来自 \(s^{-a}=\Gamma(a)^{-1}\int_0^\infty\lambda^{a-1}e^{-\lambda s}d\lambda\)；
正矩式来自
\(s^\gamma=\gamma\Gamma(1-\gamma)^{-1}\int_0^\infty(1-e^{-\lambda s})\lambda^{-1-\gamma}d\lambda\)。
两个积分均非负，故可用 Tonelli；\(0<\gamma<\beta\) 保证正矩收敛。
代入 (10)--(12) 并整理 16 的幂次，得本轮固定公式
\[
 \boxed{\displaystyle
 \mathbb E_Rd^{-3/14}\le Q_*=
 16^{5/14}\frac{\Gamma(7/4)\Gamma(4/7)}{\Gamma(19/28)}
 \left[\frac{\Gamma(61/56)}{\Gamma(239/224)}\right]^{16}.} \tag{13}
\]
本式直接给所有退化形状邻域的全域可积支配；没有把 \(\det M\) 有限采样值当作正下界。

## 4. 一次实际有限批次及准确缺口

本轮按 `round011_shape_integral_registration.json` 登记的唯一批次实际运行：
Python 精确有理运算与 python-flint 0.9.0、256-bit Arb，约 0.0245 秒，退出码零。
检查 tight frame 及全部 \(120\cdot32=3840\) 个 signed-permutation 方向映射；
计算固定有理指数、一个低阶精确矩和 (13) 的 Gamma 公式。空间点、求积、Monte Carlo、QMC
调用次数均为零；共享 MATLAB 未使用。

| 本轮对象 | 严格数值含义 |
|---|---:|
| \(Q_*\) | \(2.4576305439372284727\ldots<5/2\) |
| 完整 \(I\) 上界 (1)、(13) | \(1.9382520037464933788\ldots\cdot10^{-5}<1/50000\) |
| 同一负矩界配旧整平面因子 | \(9.6912600187324668944\ldots\cdot10^{-5}\) |
| 达到目标 \(I\le11\cdot10^{-6}\) 的充分负矩门槛 | \(1.3947585727271260714\ldots\) |
| 当前显式 \(I\) 上界 / 目标 | \(1.7620472761331757989\ldots\) |

全精度 Arb 球见结果 JSON，表中省略号不作为有向舍入后的有理端点使用。
简便可用严格有理上界是 \(I<1/50000\)，更精细可用 \(I<1.938253\cdot10^{-5}\)。

为说明负矩上界还有真实空间而非已被证明不可能，另以精确权重矩核对
\[
 \mathbb E_Rd=16^5\mathbb E_R(P_1P_2P_3P_4P_5)
 =\frac{\beta^5\,5!}{(\beta)_5}=\frac{648}{1463}.
\]
由凸性只有下基准
\(\mathbb E_Rd^{-3/14}\ge(648/1463)^{-3/14}=1.1906559399\ldots\)，
低于充分门槛 \(1.394758\ldots\)。因此当前 Gamma 上界超过门槛不排除目标。
形状 Jensen 本身也可能有损失；本轮未声称 (1) 是最优空间上界。

这条新接口避免了 upper_integral019 已排除的“将峰值从早时间一路积分”的机制。
它直接界完整 Green 面积分，不消费第 010 轮短时/角点管预算；第 010 轮管接口仍可用于
以后局部保留形状的求积或降低 Jensen 损失。

剩余实质缺口有两条可区分的路线：给 \(\mathbb E_Rd^{-3/14}\) 更紧的全域可证上界，
或者保留形状之间的相关性而减少 \((\mathbb EH)^q\le\mathbb EH^q\) 的损失。
本轮没有扩展次数、搜索参数或追加空间单元，也没有把这些任务声明为已完成。

## 5. 复现与来源锁定

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper/round011_shape_integral.py
```

必要新文件仅为本报告、登记、脚本和结果 JSON。读取的 Astra upper001 报告、源码、验证
JSON，以及登记和脚本 SHA-256 均锁定在结果 JSON 中。该上游路径为
`/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_upper001/`；
本轮只读，不修改、不重新运行其验证。

本轮结果 JSON SHA-256：
`5b77e6a098f314e8d3cac1fbbc1625ea846e6baeec562f839e9d9d298846211d`。
重运行包含实际耗时的 JSON 时，整份 hash 可以变化；固定公式和整数/有理证书应保持。
主 clone 只读；本轮仅新增 upper/round011 文件，未提交、未推送、未覆盖既有成果。

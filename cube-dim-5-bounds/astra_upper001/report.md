# Astra upper 001：trace 机制排除与保留形状的精确矩证书

2026-10-26。独立 clone：`Math-Research-cube5-astra-upper`；仅写本目录，未提交。

**本轮结论。** 对已固定的 16-body、独立 α=3/2 stable 模型，整个
“对符号轨道使用 Jensen，最后只留下 tr(A⁻¹)”的全域下界机制，都有严格天花板

\[
 g_{\rm tr}\le 8.38102676972744590524\ldots\,10^{-5}
 <84/10^6<1/10^4.
\]

这是方法本身的天花板，不是 Green 函数的上界，更不是 C5 的下界或上界。
因此，继续仅提高该 trace 期望或逐时倾斜的数值精度不能达到现有 g 目标。

替代路线是先精确消去时间和总体尺度，保留一个对角线恒为 1 的随机相关矩阵。
它的所有多项式矩均可精确计算。原分数幂 minorant 约束可以变成纯多项式不等式，
从而将全域 g 认证改写为有限的有理多项式证书，避免 16 维概率积分。
本轮已实际做出一个低阶全域证书，但仅给 g≈7.209e-6，弱于主线已得的倾斜证书。
新接口的完整性证明不表示所需目标证书已找到；C5 仍未改进。

## 1. 共同归一化与先积分时间

令 β=3/4、r=7/4，V={(1,±1,±1,±1,±1)}。取 16 个独立正 stable 变量

\[
 \mathbb E e^{-\lambda S_v}=e^{-\lambda^\beta},\qquad
 A=\sum_{v\in V}S_vvv^T,\quad D=\det(A)^{-1/2}.
\]

条件 Gaussian 的协方差是 2t^(4/3)A，其特征函数平均后恰为
exp(-t Σ|v·ξ|^(3/2))。正 stable 变量几乎处处正且有限，V 张成 R5，故 A 正定。
对 x≠0 置 a=xᵀA⁻¹x。Tonelli 及 z=a/(4t^(4/3)) 给

\[
\begin{aligned}
G(x)&=(4\pi)^{-5/2}\mathbb E\left[
 D\int_0^\infty t^{-10/3}e^{-a/(4t^{4/3})}\,dt\right]\\
&=\kappa\mathbb E[D a^{-7/4}],\qquad
\kappa=\frac{3\Gamma(7/4)}{2^{7/2}\pi^{5/2}}.
\end{aligned} \tag{1}
\]

等式允许 +∞，不需要事先删除奇异方向。x=0 时 G=∞。
本处 covariance 的因子 2、stable 时间及 Γ(7/4) 都不可省略。

## 2. 完整 trace-Jensen 机制的严格天花板

A 的分布在坐标符号共轭下不变，D 同时不变。对固定 A 的 32 个符号矩阵 U，
凸函数 z↦z^(-7/4) 给

\[
 2^{-5}\sum_U((Ux)^TA^{-1}(Ux))^{-7/4}
 \ge\left(\sum_i(A^{-1})_{ii}x_i^2\right)^{-7/4}
 \ge(\operatorname{tr}A^{-1})^{-7/4}
\]

对 x∈cube 成立。因此该机制的最强精确期望值是

\[
 g_{\rm tr}=\kappa\mathbb E[D(\operatorname{tr}A^{-1})^{-7/4}]. \tag{2}
\]

对 A⁻¹ 的五个正特征值作 AM–GM：tr A⁻¹≥5D^(2/5)，再用 z^(3/10) 的凹性，

\[
 g_{\rm tr}\le\kappa5^{-7/4}\mathbb E D^{3/10}
 \le\kappa5^{-7/4}(\mathbb E D)^{3/10}. \tag{3}
\]

下面给出 ED 的完全初等、有限几何上界，无需引用 Brascamp–Lieb。

取五方向基 B，其行是

\[
\begin{pmatrix}
1&1&1&1&1\\1&-1&-1&-1&1\\1&-1&-1&1&-1\\
1&-1&1&-1&-1\\1&1&-1&-1&-1
\end{pmatrix},\qquad\det B=-48.
\]

让保持第一坐标的 16 个符号矩阵作用于 B，每个所得基仍有绝对行列式 48，
且 V 中每个方向恰出现 5 次。这由符号作用的传递性直接得到，脚本亦逐项整数检查。
设 a₁=s₁(0)=2Γ(2/3)/(3π)，对 16 个基的指数函数用 Hölder，得到

\[
\begin{aligned}
p_1(0)&=(2\pi)^{-5}\int_{\mathbb R^5}e^{-\sum_v|v\cdot\xi|^{3/2}}d\xi\\
&\le\frac{a_1^5}{48}\left(\frac5{16}\right)^{10/3}
=:P_*.
\end{aligned} \tag{4}
\]

具体地，每个 Hölder 因子是 exp[-(16/5)Σ_(v∈B_U)|v·ξ|^(3/2)]；
16 个因子的 1/16 次几何平均正好是原被积函数，每个因子的积分由五个独立坐标求出。
由条件 Gaussian 在零点的公式，ED=(4π)^(5/2)p₁(0)。代入 (3)--(4)，即得

\[
 g_{\rm tr}\le\kappa5^{-7/4}[(4\pi)^{5/2}P_*]^{3/10}
 <84/10^6<1/10^4. \tag{5}
\]

实际 256-bit Arb 包围见 verification.json，约 8.3810267697274459e-5。
这个排除不限制保留整个符号轨道的下界，也不限制下面的完整相关矩阵接口。

## 3. 径向消元：精确概率与十个形状变量

令 T=ΣS_v、P_v=S_v/T、M(P)=ΣP_v vvᵀ，所以 A=TM。每个对角元 M_ii=1，
tr M=5。M 只依赖十个非对角相关量 r_ij=ΣP_v v_i v_j，落在

\[
\mathcal K=\operatorname{conv}\{vv^T:v\in V\}
\]

这个十维紧多面体中。P 几乎处处是 simplex 内点，M 几乎处处正定。

负幂的 Gamma 积分给，对任意 z_v>0，

\[
\mathbb E(\sum_v z_v S_v)^{-\beta}
=\frac1{\Gamma(1+\beta)\sum_vz_v^\beta},\qquad
\mathbb ET^{-\beta}=\frac1{16\Gamma(7/4)}. \tag{6}
\]

定义概率 R 为用 T^(-β) 对原概率倾斜后推送到 P 的分布；没有声称 T 与 P 独立。
由于 (1) 中矩阵函数对 A 是 -β 次齐次，得到精确正表示

\[
\boxed{G(x)=c\,\mathbb E_R H(M,x),\quad
c=\frac3{2^{15/2}\pi^{5/2}},\quad
H(M,x)=\det(M)^{-1/2}(x^TM^{-1}x)^{-7/4}.} \tag{7}
\]

Arb 实际包围 c≈0.00094737431379484146449。
这一步同时消去了时间和总体尺度，没有 determinant-Jensen 或 trace-Jensen。
由 (6) 直接得到精确变换

\[
\boxed{\mathbb E_R(\sum_vz_vP_v)^{-\beta}=
\frac{16}{\sum_vz_v^\beta}.} \tag{8}
\]

本轮证明只依赖 Gamma 积分和给定 Laplace 变换，不依赖 normalized-stable 文献定理。

## 4. 所有需要的多项式矩都是有限有理运算

若 L_j=Σ_v a_(vj)P_v，令

\[
F(z)=\frac{16}{\sum_v(1-\sum_j a_{vj}z_j)^\beta}.
\]

在原点邻域展开 (8)，对多重指标 k、|k|=n，有

\[
 \mathbb E_R\prod_jL_j^{k_j}
 =\frac{\prod_j k_j!}{(\beta)_n}[z^k]F(z). \tag{9}
\]

P 位于 compact simplex，因此局部解析展开与期望交换合法。β=3/4、a 有理时，
固定次数只涉及有限有理运算。十个 r_ij 的 a 列就是十个不同 Walsh 符号字符。

实际有限验证给：

\[
\begin{gathered}
\mathbb Er_{ij}=0,\quad
\mathbb Er_{ij}r_{kl}=\tfrac17\,1_{\{(i,j)=(k,l)\}},\quad
\mathbb E\|M-I\|_F^2=20/7,\\
\mathbb Er_{12}^4=3/55,\quad \mathbb Er_{12}^6=141/4807,\\
\mathbb Er_{12}r_{13}r_{23}=5/77,\quad
\mathbb Er_{12}^2r_{13}^2=\mathbb Er_{12}^2r_{34}^2=17/385.
\end{gathered} \tag{10}
\]

因此不能把相关量当作独立变量，也不能只靠“它在 I 附近”做未经控余项的 Taylor。
尤其非零三角矩记录了符号轨道中 trace 所丢掉的信息。

另用独立的 Laplace 导数展开交叉核验了登记的九个权重矩。若

\[
(-\partial_t)^k e^{-t^\beta}=e^{-t^\beta}
\sum_j a_{k,j}t^{j\beta-k},
\]

则 a_(k+1,j)=(k-jβ)a_(k,j)+βa_(k,j-1)，并由 Gamma 积分得

\[
\mathbb E_R\prod_iP_i^{k_i}=
\frac1{(\beta)_{\sum k_i}}
\sum_{j_i}\frac{(\sum j_i)!}{16^{\sum j_i}}\prod_i a_{k_i,j_i}. \tag{11}
\]

verify.py 中 (9) 和 (11) 是两个独立实现，有限输出完全相等。

## 5. 全域 g 的纯代数证书接口

令 d=det M、b=xᵀadj(M)x。在 M 正定且 x≠0 时

\[
 H(M,x)=d^{5/4}b^{-7/4}. \tag{12}
\]

选有理系数多项式 P₀(r,x)，构造 h=d²P₀。以下两项都是**纯多项式约束**：

\[
 P_0\ge0,\qquad 1-d^3P_0^4b^7\ge0
 \quad\text{在 }\mathcal K\times\partial Q. \tag{13}
\]

它们蕴含 0≤h≤H；这只是对 (12) 的正数四次幂比较。边界奇异矩阵具有零概率，
而 (13) 可在闭多面体上直接检验。利用 (9)，

\[
 L(x)=\mathbb E_R[d^2P_0(r,x)] \tag{14}
\]

是可精确输出的有理系数空间多项式。最终用 Arb 取得一个有理 c_L<c，再认证
c_L L(x)≥g₀ 于 cube 各面，即给所需全域 Green 下界。齐次性把面下界扩到内部；
原点为∞。同时对 P₀ 做联合坐标符号/置换平均，可把空间域缩到
x=(1,y)、0≤y₄≤y₃≤y₂≤y₁≤1，不引入任何 nonlinear Jensen 损失。

有限验证形式可以选：精确 Bernstein 系数、带可验有理系数的 SOS 恒等式，或
带舍入证明的区间细分。优化器结果本身不是证书。用 P 的 15 维 simplex
参数化 K 是最直接但较大的实现；利用十维相关多面体的描述可减变量。
式 (13) 的总次数较高，当前未实现目标次数的优化或全域 positivity 验证。

### 这个 minorant 族为何没有 trace 的固定天花板

以下是存在性/完整性结论，**不是计算复杂度保证**。在紧域 K×∂Q 上，d、b≥0，
且有界。正定处定义 f=H/d²=d^(-3/4)b^(-7/4)。令

\[
f_N=\min\{N,f\},
\]

在 d b=0 处延拓为 N，则 f_N 连续、有界，并随 N 增加。
第 6 节给 f≥1/35，所以 f_N（N≥1）有统一正下界。
由正定性几乎处处成立，

\[
\mathbb E_R[d^2f_N]\uparrow\mathbb E_RH=G/c.
\]

左侧为连续的空间函数。如果某目标 g₀ 严格低于整个 cube 面的**逐点** inf G，
连续函数的递增开覆盖及 compactness 保证某个有限 N 已在整个面上超过 g₀/c。
对 (1-ε)f_N 做足够精确的有理多项式一致逼近，利用 f_N≥1/35，可得到 P₀>0、
P₀≤(1-ε/2)f_N，同时保留期望的目标余量。于是 (13) 的第二残差具有统一严格
正余量 1-(1-ε/2)^4。需要逼近时，可把 max(d,0)、max(b,0) 代入上述截断式，
先连续延拓至环境盒，再用 Bernstein 逼近。

严格正的有理多项式在 simplex×box 上，经足够高的 Bernstein 次数提升后有
正系数（固定多项式的提升系数一致趋近对应网格值），从而存在有限有理证书。
这证明该族在上述严格目标假设下不受 trace-Jensen 天花板限制。
尚未证明实际 inf G>1/10000，也未找到实现该目标的 P₀；不可倒过来把完整性当成目标成立。

## 6. 一个实际完成的低阶证书

令 μ 是 M 的最小特征值，其余为 λ₂,…,λ₅。tr M=5、|x|²≤5 给
b≤5d/μ。因此

\[
d^3b^7\le5^7\mu^3\prod_{j=2}^5\lambda_j^{10}
\le 5^7(15/43)^3(50/43)^{40}<35^4. \tag{15}
\]

第二步是带权 AM–GM，权重为 3、10、10、10、10，总和 43。
最后一步为脚本实际核验的纯有理严格比较。于是 P₀=1/35 满足 (13)，

\[
 H(M,x)\ge d^2/35. \tag{16}
\]

Cauchy–Binet 展开 d=Σ_(|B|=5) det(V_B)²∏_(i∈B)P_i；枚举恰好 4368 个五子集，
其中 3008 个非退化，系数和是 16⁵。用两个基交集大小 k=0,…,5 分类，结合 (11)，
精确求得

\[
 \mathbb E_Rd^2=\frac{126404864}{474619145}.
\]

代码没有逐对抽样：令 w_B=det(V_B)²，计算
b_j=Σ_(|J|=j)(Σ_(B⊇J)w_B)²，然后二项反演获得所有有序基对的交集权重。
代入 (7)、(16)，256-bit Arb 认证

\[
 G(x)\ge\frac{c}{35}\frac{126404864}{474619145}
 =7.20895134274030898246\ldots\,10^{-6}
\quad(x\in Q).
\]

这个值弱于已有约 1.543e-5 的倾斜下界，**不替换主账本**。作用是展示新接口
从解析变换、有限有理矩到全域下界的一条实际运行闭环；没有用随机最小值。

## 7. 交付、运行与剩余缺口

实际本机运行：Python 3.15.0、SymPy 1.14.0、python-flint/Arb 256 bit。
最终登记内固定检查约 0.05 秒；没有 Monte Carlo/QMC，也没有新增空间网格。
复现：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/astra_upper001/verify.py
```

必要文件只有本报告、registration.md、verify.py、verification.json；源输入和
复现代码 SHA-256 在 JSON 中逐项保存。主线尚未审阅本轮证明，建议按 (1)、(5)、
(8)、(13)、(15) 五个接口分别审计。

剩余实质工作是寻找保留形状/空间相关性的 P₀ 并验证 (13)--(14)，同时另行完成 I
的严格上界。十维相关矩与高次 positivity 仍可能很昂贵；本轮没有声称已解决该
计算瓶颈。可优先用已有精确矩筛查低次数对称多项式，然后只认证一个预登记候选。

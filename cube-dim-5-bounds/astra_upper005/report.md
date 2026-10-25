# Astra upper 005：二次全相关辅助函数，g≥4583/10⁸

2026-10-26。本轮严格取得固定16-body、α=3/2 stable模型的全域 Green 下界

\[
\boxed{G(x)\ge4583/10^8=0.00004583\quad(x\in[-1,1]^5).}
\]

它提高了004的7/156250。本轮也严格排除固定二次行列式minorant下的整个全相关
**线性**辅助函数族达到1/20000。二次辅助函数只验证三个登记候选，未被该线性族
排除覆盖。I尚未满足联合成本，C5≤32没有改善。

## 1. 继承接口与连续线性族天花板

沿用001--004的精确表示和有理端点：

\[
G=cE_RH,\quad c=3/(2^{15/2}\pi^{5/2}),\quad
H=d^{5/4}b^{-7/4},\quad
\ell_L(d)=A_Ld-B_Ud^2\le d^{5/11},
\]

其中d=det M、b=xᵀadj(M)x，0<d≤1，且ℓ_L≥0。全部端点分数保存在checkpoint.json。
Young不等式给任意Z≥0：

\[
H\ge\frac{11}{4}\ell_L(d)Z^7-\frac74bZ^{11}. \tag{1}
\]

记W=Σ_(i<j)M_ij x_i x_j。004已证明−2≤W≤10于全cube。
首先考虑Z=a(1+ηW)、a≥0、0≤η≤1/2。在角点x=𝟙，以未圆整的ℓ=A d−Bd²定义
U(η)=Eℓ(1+ηW)^7、V(η)=Eb(1+ηW)^11，均正；解析最大化a给必要天花板
cU^(11/4)/V^(7/4)。

用A_U>A、B_L<B构造U_U≥U，并取c_U>c。checkpoint.py实际准确生成

\[
R(u)=(1/20000)^4V(u/2)^7-c_U^4U_U(u/2)^{11},\qquad0\le u\le1.
\]

该77次有理多项式的78个Bernstein系数全部严格为正，整区间无需细分。
因此固定此ℓ的所有线性Z、a≥0、0≤η≤1/2都无法达到1/20000。
这是明确候选族的上限，不是G或C5的上界，也不排除下述非线性Z。

## 2. 三个二次辅助函数与全域非负性

本轮恰好检查

\[
Z=a f_\kappa(W),\qquad
f_\kappa(W)=1+W/8-\kappa W^2,\quad
\kappa\in\{1/512,1/256,1/128\}.
\]

f_κ凹，且在W的两个端点−2、10严格正，具体为

|κ|f(−2)|f(10)|
|---|---|---|
|1/512|95/128|263/128|
|1/256|47/64|119/64|
|1/128|23/32|47/32|

因此Z≥0在全cube成立，(1)无需额外高维positivity假设。
按照登记规则，由角点U=Eℓ_L f^7、V=Eb f^11计算最佳a=(U/V)^(1/4)的Arb包围，
向下取1/1000网格的有理数；随后全部做完整空间验收，角点仅用于参数选择。

|κ|准确a|本轮严格全域g|
|---|---|---|
|1/512|671/1000|4553/10⁸|
|1/256|27/40|4583/10⁸|
|1/128|17/25|2279/50000000|

没有追加第四个κ，也没有用小数最小值构建全域结论。

## 3. 任意精度图矩与独立校验

因为f²次，所需矩增加到

\[
EdW^n,\ Ed^2W^n\quad(n\le14),\qquad EbW^n\quad(n\le22),
\]

最高相关图次数26。Euler图整数递推及D_n=∏_(j=0)^(n−1)(4j+3)与004完全相同。
本机Apple clang不支持超过128位的_BitInt，未将失败编译当作运行结果；最终改用
另一代理独立实现并检验的GMP任意精度后端
`upper/round012_graph_moments_gmp.hpp`，锁定SHA-256为
`7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8`。

本轮graph_moments.cpp仅负责空间多项式展开，调用该后端的B_k。矩、阶乘、
二项乘积及空间累计均为mpz_class，不受固定128位容量限制。边的指数使用五位
编码；本轮最大单边指数不超过24，低于31，最大总阶26也在后端的显式检查范围内。
后端逐项验证0≤B_k≤D_|k|，并用绑定源哈希的journal保存已经完成的整数矩。
kernel明确关闭隐式时间中断，在每个完整W阶数后sync并输出进度。

完整空间生产实际完成143289个规范图矩，最大B位数为131，证实不能继续使用signed128。
GMP阶段约217.2秒。持久化缓存和临时原始表不归档；中断后可以利用缓存重放空间
累加，不能把半完成的空间累计重复相加。

所有结果进入两组准确交叉检验后才作下界认证：

1. 与004的全部既有空间多项式逐项有理相等；
2. 每个所需阶数在角点的值，与独立三组Cauchy–Binet/Stieltjes算法提取至22阶的
   有理矩完全相等。三组权重仍为10、2、−2，方向数1、5、10。

任意精度后端另有其独立高阶pilot报告，包括26阶单边矩的131位整数及18阶三角图
的独立Stieltjes计算；本轮不把那三个pilot等同于全部空间生产的独立证明。

## 4. 28,561个完整面系数给出全域结论

将f^7、f^11作准确有理展开，代入(1)并取期望，得到X_i=x_i²的对称有理多项式

\[
F(X)=\frac{11}{4}a^7E_R[\ell_L(d)f_\kappa(W)^7]
-\frac74a^{11}E_R[b f_\kappa(W)^{11}]. \tag{2}
\]

总X次数最高23，每个X_i次数最高12。完整partition系数保存在verification.json；
可携精确矩保存在space_moments.json。所有数都是原始有理数，不能由展示小数重建。

利用G和模型的带符号置换对称性，只需验收x₁=1、X₂,…,X₅∈[0,1]。
verify_global.py将(2)准确转换为四维、每个变量次数12的Bernstein系数。
每个候选共13⁴=28561个系数，逐个有理比较后发现最小值均在(12,12,12,12)。
由于Bernstein基非负且和为1，这些最小系数给整个连续面的下界。
“最小恰为角点值”是完整系数核验后的结果，不是预设的极值原理。

最佳κ=1/256、a=27/40时，最小系数为

\[
L=\frac{113191546367073879688956437920125691089139674463974443421317398051467}
{2339680948874089473235836222342600438430872633344000000000000000000000}.
\]

L严格正，乘以已核有理c_L<c得到

\[
c_LL=
\frac{107234763566675158859659345654848680739883593890039193435759673167060784421984331}
{2339680948874089473235836222342600438430872633344000000000000000000000000000000000000}
>4583/10^8.
\]

G的−7/2次齐次性将面下界扩到cube内部；原点G=∞。这完成了全域证书。
Python全空间检验约4.1秒。

## 5. 兼容的空间强区接口

等待后端期间，另以一次1820个四方向子集的有限整数计算，补全002的点态Hölder下界。
令

\[
Q_4(x)=E_Rb(x)^2=A\sum_i x_i^4+B\sum_{i<j}x_i^2x_j^2,
\quad A=65024/168245,\quad B=860544/1043119.
\]

spatial_bridge.py精确得到A，再用已审角点Q₄(𝟙)核定B。由002的同一证明，

\[
G(x)\ge K Q_4(x)^{-7/8},\quad
K=c\frac{(648/1463)^{5/2}}{(126404864/474619145)^{5/8}}.
\]

Arb和有理比较给出两个严格空间强区：

- Q₄≤3280621/10⁶时，G≥1/10000；
- Q₄≤515999/250000时，G≥3/20000。

又A−B/2<0，故Q₄≤(B/2)||x||⁴，分别可使用更简单的
||x||²≤14/5、||x||²≤11/5两个充分条件。这不是整个cube的目标证书。
例如(1,1,1,1,0)在这个空间接口中仅得到约5.50e-5，说明剩余弱区不只包含极近体对角线
的点。单方向ridge补偿尚未在本轮计算或认证。

## 6. 复现、白名单与剩余问题

完整重算需要GMP6.3.0开发文件及锁定的round012头文件。以实际本机路径为例：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper005/checkpoint.py
clang++ -O2 -std=c++17 -isystem /opt/homebrew/Cellar/gmp/6.3.0/include -I /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper cube-dim-5-bounds/astra_upper005/graph_moments.cpp -L /opt/homebrew/Cellar/gmp/6.3.0/lib -lgmpxx -lgmp -o /tmp/astra_upper005_graph_moments
/tmp/astra_upper005_graph_moments /tmp/astra_upper005_graph_moments.json
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper005/verify_global.py /tmp/astra_upper005_graph_moments.json
```

已保存的有理空间矩可直接用于审阅/复核，无需先重跑约三分钟的后端：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper005/verify_global.py cube-dim-5-bounds/astra_upper005/space_moments.json
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper005/spatial_bridge.py
```

精简白名单和逐文件SHA-256见manifest.json，包括本报告、登记、三个验证脚本及
必要JSON、空间kernel源码、可携空间矩与最终有理证书。不归档二进制、GMP安装副本、
临时图矩大表或journal；原始图矩表SHA绑定在可携结果中。

本轮新g仍距离现有联合成本目标很远。二次ℓ(d)在d→0仅为O(d)，小于真实
d^(5/11)，会损失接近低秩形状的重要贡献。下一步可提高该分数幂minorant，或使用
保留单/多方向stable尺度的条件下界覆盖空间弱区。没有宣称本轮三个二次候选已经
穷尽所有非线性辅助函数，也没有新的C5上界。

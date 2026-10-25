# Astra upper 004：全部十相关的全域多项式证书，g≥7/156250

2026-10-26。固定 α=3/2、16 条未单位化 body 方向的独立 stable 模型。
本轮严格取得

\[
\boxed{G(x)\ge\frac7{156250}=0.0000448\quad(x\in[-1,1]^5).}
\]

它严格提高了 002 的 1/27000。所有空间方向由完整四维面上的2401个有理
Bernstein系数覆盖；没有把角点检查或小数最小值当成全域下界。
当前 C5 上界仍为32，未完成 g 与 I 的突破成本。

## 1. 点态下界、辅助函数的全域非负性

继承已锁定的 001--003 归一化：

\[
G(x)=cE_RH(M,x),\quad c=\frac3{2^{15/2}\pi^{5/2}},\quad
H=d^{5/4}b^{-7/4},\ d=\det M,\ b=x^T\operatorname{adj}(M)x.
\]

M 是16个 vvᵀ 的随机凸组合，正定几乎处处、对角线恒为1，0<d≤1；R 在带符号
坐标置换下不变。003 锁定的有理二次函数 ℓ_L(d)=A_Ld−B_Ud² 满足
0≤ℓ_L≤d^(5/11)，常数准确值保存在本轮 corner_checkpoint.json 和 verification.json。

Young 不等式对任意 Z≥0 给

\[
H\ge\frac{11}{4}\ell_L(d)Z^7-\frac74bZ^{11}. \tag{1}
\]

本轮保留全部十个相关项，取

\[
W(M,x)=\sum_{i<j}M_{ij}x_ix_j,\qquad Z=a(1+\eta W).
\]

对每条 v，W_v(x)=Σ_(i<j)v_iv_jx_ix_j 是每个坐标分别仿射的函数，故在cube上的
最小值在顶点取得。顶点处 v·x 是五个±1的和，因此
W_v=((v·x)²−5)/2≥−2。对 vvᵀ 的凸组合仍有 W≥−2。
所以0≤η≤1/2时 Z≥0（a≥0），(1) 在整个cube成立。

三个预登记参数及准确有理 a 为

|η|a|严格全域 g|
|---|---|---|
|1/32|137/200|2007/50000000|
|1/16|689/1000|1107/25000000|
|1/8|333/500|7/156250|

a 按登记规则由角点最佳 a 的Arb包围向下选择1/1000网格点；角点仅用于选择参数。
每个候选随后都通过了独立的完整空间有理证书。

## 2. 独立的角点三组计算

在x=𝟙处，W=Σ_vP_v w_v，w_v=((Σ_i v_i)²−5)/2 只取10、2、−2，方向数分别为
1、5、10。corner_checkpoint.py 先用三组 Cauchy–Binet 系数及基对交集分类
计算 Ed W^n、Ed²W^n，再用广义 Stieltjes/Laplace 有限导数算所有所需有理矩。
它沿用003的有限混合矩证明，分母改为

\[
\{(1-10z)^\beta+5(1-2z)^\beta+10(1+2z)^\beta\}/16,
\quad\beta=3/4.
\]

三组 Gram 矩阵的对角元为 a+5b+10c、非对角元为 a+b−2c，因此沿𝟙的特征值为
5a+9b+2c，其余四个特征值为4b+12c。故det的组生成多项式是
(5a+9b+2c)(4b+12c)^4，且𝟙ᵀadj(M)𝟙的系数由5(4b+12c)^4给出。
基对产生174个非零指数模式，全部精确处理，无随机取样。

三个候选的角点值都通过必要筛查，但 JSON 明确标为“非全域证书”。下一节
完全不同的图矩算法在构造全部空间系数后，与这些角点有理矩逐个准确相等。

## 3. 完整空间系数：有限整数 Euler 图矩

记 r_e=M_ij，边e对应K5的十条无向边。对多重指标 k，|k|=n，由001的变换

\[
E_R(1-\sum_et_er_e)^{-3/4}
=\frac{16}{\sum_{v\in V}(1-\sum_e t_ev_iv_j)^{3/4}}
\]

得到：当多重图k有奇度顶点时 Er^k=0；当全部顶点度数为偶数时，定义

\[
D_n=\prod_{j=0}^{n-1}(4j+3),\quad B_k=D_nE_Rr^k,\quad B_0=1,
\quad C_m=3\prod_{j=1}^{m-1}(4j-3),
\]

则精确递推为

\[
\boxed{B_k=\sum_{0<j\le k,\ j\text{为Euler图}}
C_{|j|}\prod_e\binom{k_e}{j_e}B_{k-j}.} \tag{2}
\]

证明：符号字符 v_iv_j 的乘积平均在每个顶点度数偶时为1，否则为0；展开分母，
非零非常数项均带负号。将分母与逆级数相乘，比较系数，再乘4^n即可得到 (2)。
递推右侧都是整数，阶数严格下降；没有浮点矩计算。

graph_moments.cpp 用S5置换归一化多重图，只需3389个不同的规范图矩，
共归一化71782个输入图。它展开det(M)的73项、det(M)²的2021项以及b的315个
余子式项，配合 W^n 的有限多项式展开，准确求出

\[
E_RdW^n,\ E_Rd^2W^n\ (0\le n\le7),\qquad
E_RbW^n\ (0\le n\le11)
\]

作为X_i=x_i²的对称多项式。代码只保存一个partition代表的单项式系数，恢复时
对所有不同的坐标排列求和；不是把聚合系数误当成单项系数。

### 整数实现与可复核性

所有大整数采用有符号128位，矩递推及系数累计的加、乘均带compiler overflow检查；
任何溢出都会立即报错，不能绕回。最高图次数17，D₁₇约2.44236e24；计算中还核验
0≤B_k≤D_n。小整数乘法最大可由11!及有限det系数的绝对和约束，远低于64位限制。
有限精度整数不是近似数值运算。

计算结果交给 Python Fraction 合并不同阶数的D_n分母。对每一个所需n，空间
多项式在𝟙的值都与第二节独立三组算法的有理结果完全相等；同时恢复
Ed=648/1463、Ed²=126404864/474619145及 Eb(x)=(216/385)Σx_i²。
这些交叉检验记录在space_moments.json。

## 4. 完整正面上的有理 Bernstein 证书

将上述精确矩代入 (1)，得到对称有理空间多项式

\[
F_{a,\eta}(X)=\frac{11}{4}a^7\sum_{n=0}^7\binom7n\eta^n
\{A_LEdW^n-B_UEd^2W^n\}
-\frac74a^{11}\sum_{n=0}^{11}\binom{11}n\eta^n EbW^n. \tag{3}
\]

它满足 E_RH(M,x)≥F_(a,η)(x₁²,…,x₅²)。最高X总次数12（b额外带一个X次数），
每个坐标最高次数6。verification.json 完整保存每个候选的有理partition系数。

带符号坐标置换对称性和G的−7/2次齐次性允许只验收x₁=1，其他x_i²∈[0,1]。
verify_global.py 将 (3) 代入X₁=1，精确转换到四维、每个变量次数6的Bernstein基。
转换公式逐坐标为

\[
\mathcal B_k=\sum_{j\le k}f_j\prod_{i=1}^4
\frac{\binom{k_i}{j_i}}{\binom6{j_i}}.
\]

每个候选的2401个系数全部通过有理比较，最小系数都准确位于k=(6,6,6,6)。
Bernstein基非负且和为1，所以这个最小值是整个连续面的严格下界；不需要网格外推。
它恰与角点值相等是**完整系数核验后的结论**。

对于最强的η=1/8、a=333/500，完整最小Bernstein系数是

\[
L=\frac{542862901331884932496974940874070007125668534129721}
{11479405938548387840000000000000000000000000000000000}.
\]

乘以已验证的有理c_L<c，实际得到

\[
c_LL=
\frac{514294368632971553675557760658888820627722847837516632821541753}
{11479405938548387840000000000000000000000000000000000000000000000000}
>\frac7{156250}.
\]

这里L>0已由有理证书确认，所以用c下端相乘的方向合法。
对非零cube内点x，令ρ=||x||∞≤1，点态齐次性给G(x)=ρ^(-7/2)G(x/ρ)≥7/156250；
原点G=∞。因此该结论覆盖整个cube及所有边界层。

## 5. 复现与精简白名单

实际本机环境：Apple clang 21.0.0，C++17；Python 3.15.0；python-flint/Arb256位。
三组必要检查约0.21秒；C++图矩阶段约0.16秒；全域Python证书约0.19秒。
这些时间只是实际运行记录，不参与数学推理。

从独立clone根运行：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper004/corner_checkpoint.py
clang++ -O2 -std=c++17 cube-dim-5-bounds/astra_upper004/graph_moments.cpp -o /tmp/astra_upper004_graph_moments
/tmp/astra_upper004_graph_moments /tmp/astra_upper004_graph_moments.json
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_upper004/verify_global.py
```

白名单：本报告、registration.md、corner_checkpoint.py/.json、graph_moments.cpp、
verify_global.py、space_moments.json、verification.json、manifest.json。
不归档可重建的二进制或/tmp图矩中间文件；该中间文件的SHA-256绑定在verification.json。
完整Bernstein系数可从保存的有理空间多项式准确重建，结果散列也已保存。

依赖锁定的astra_upper001/verify.py与astra_upper003/verification.json，输入SHA-256
保存在corner_checkpoint.json；绝不能从展示小数重新建立本证书。

## 6. 结论边界

新全域g已获本轮精确数学及有限算术证书，待主代理独立审计。004未证明任何
整个全W线性族的天花板；003的单相关族排除不能移用于这里。后续若研究全W族
连续参数上限，会另行登记与保存，不能改变当前证书。
I尚未足够小，故本轮不修改C5≤32，也不把g的改进当成C5改善。

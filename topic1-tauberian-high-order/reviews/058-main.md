# R058 main：独立完整数学审查

2026-11-05。**最终裁决：PROOF-VALID。无必要数学补丁。**

本次只审查新版主文件，未修改作者源、baseline、论文仓库或他人文件，未编译、提交或发布。审查完成后按暂停要求待命，不开启新研究。

## 1. 版本、供货与裁决范围

实际全文读取并核对的对象为：

`/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/scratch/r058-eighth-main/main.tex`

**SHA256：`413a02da01cccd08f33e25709b5ac579f19d092dbfd4aa965804e675c22f68ed`**，524 行。已实际读取与同目录 `main.baseline.tex` 的完整 `diff -u`。改动集中在精确八阶主陈述、ordinary 第八项、五阶参数及全重数 little-o 传递；既有精确约化、临界发散、有限见证选择及显式包络证明保持原样。

高阶供货为 R057：

`scratch/r057-eighth-global-assembly/global_eighth_assembly.md`

源 SHA256：`e3ed7614d8775cd57fe6792ab6235bfdb8ec0029ec94532564052ba67d78277c`。

实际读取其 §0、§§6–11 的准确结论与量词，并全文读取两份正式终审：

* `scratch/r057-review-global/audit.md`，SHA256 `76a80e950cfe11bbbaaa73137ddfe176558a5abc593572bf2606b46a3cdd9a4f`；
* `scratch/r057-review-adversarial/audit.md`，SHA256 `cd12441096e803a07674ecb7182148c43cc4ef55308a4cffb6e307abaa0737a9`。

上述三个 hash 本次均实际重核。两审均最终 PROOF-VALID，明确供应全部固定整数 (n\ge2)、任意正有限测度可测输入、全部正半径、严格水平集的精确八阶式；下界对每个充分小实参数成立。它们已经核闭新二维补集、同一 N8、一次 J/删除账本、全输入上确界以及两尺度 lower，故 main 在此基础上使用真正极限合法，不再只是局部模型或形式展开。

本签署不代替正在独立转写审查的 `sections/eighth.tex` 全文保真审。main 只需其 `eighth:theorem` 供应 R057 的上述定理；该外部标签保持。本文主文件自身的证明、系数和重数传递已独立核查，未将未审转写内容当作新的数学供货。

## 2. 精确八阶陈述及二维新费用

摘要、主定理和自然参数式统一取

\[
a=v^n=2^nt,\qquad v=2t^{1/n},\qquad
L_n(t)=\left[1-\frac a{(1+v)^n}\right]^{-1},
\]
\[
K_n=\binom n2\frac{n^2}{27},\qquad
\widehat J_n=\frac{(n-1)n^3(2609+1971n+81n^2)}{26244}
                +\frac{16}{27}\mathbf1_{\{n=2\}}.
\]

这与已审 R057 完全一致。尤其二维指示项已经是最终总系数中的 (16/27)，并未沿用旧单尺度 (8/27)。旧正常模型和真实质量分母给 (15694/6561)，新增 inside 精确费用及匹配两尺度 lower 再加 (8/27)，最终为 (17638/6561)。main 未把两笔不同来源的费用混淆，也没有把下界构造的安全可删预算当作上界必要性。

由 (av^j=2^{n+j}t^{1+j/n})，主式三个修正的 (2^{n+6},2^{n+7},2^{n+8}) 因子正确。二维明确式为

\[
C_2(1-t)=L_2(t)+\frac4{27}av^6-\frac{176}{243}av^7
                   +\frac{17638}{6561}av^8+o(av^8),
\]

其中 ((2+26/9)(4/27)=176/243)。自然归一化以准确的 (P_{7,n}) 为中心，分母为 (av^8)，与 (t) 式等价。

高阶部分始终限定固定 (n\ge2)，并未把精确八阶供货扩大到 (n=1) 或维数一致。原 (n\ge1) 包络仍独立成立。所有下界范围写成每个充分小实参数，未从离散序列 lower 升格为完整 liminf。ordinary 第八项没有额外非零或统一符号声明。

## 3. ordinary (b_{n,0},\ldots,b_{n,8}) 的全部低维碰撞

令 (v\to0)。解析函数的绝对收敛级数为

\[
L_n-1=\sum_{k\ge1}v^{kn}(1+v)^{-kn}.
\]

在 (v^{n+j}) 上令 (\ell=j-(k-1)n\ge0)，则

\[
[v^{n+j}]\{v^{kn}(1+v)^{-kn}\}
=(-1)^\ell\binom{kn+\ell-1}{\ell}
=(-1)^{j-(k-1)n}\binom{n+j-1}{j-(k-1)n}.
\]

所以 (k\le1+\lfloor j/n\rfloor)，main 的有限和正确，且已经包含所有整数幂碰撞。对 (j=8)，上限分别为：(n=2) 时 5；(n=3,4) 时 3；(5\le n\le8) 时 2；(n\ge9) 时 1。main 的通式没有把第七项的 (n\ge8) 单项分界错误搬给第八项。

固定 (n) 时，在零点的解析性使截至 (v^{n+8}) 的余项为 (O_n(v^{n+9})=o_n(v^{n+8}))。把真实修正只加入 (j=6,7,8)，再乘 (v^{n+j}=2^{n+j}t^{1+j/n})，正好得到 main 的 (b_{n,j}) 及其 little-o。无把包络形式系数当真实系数的步骤。

特别 (b_{3,5}=2^8[-\binom75+\binom72]=0)。(j=1) 仅有 (k=1)，故 (b_{n,1}=-n2^{n+1})，首个归一化极限正确。

原保留的第七项非零说明也已核查：(n\ge8) 时单项为负；(n=2,\ldots,7) 的六个有理数均正确，三维为正，其余为负。该说明只涉及 (j=7)，没有暗含 (b_{n,8}\ne0)。

## 4. 五阶参数与统一六阶余项

精确参数仍为

\[
\sigma_m(\varepsilon)=m\{1-(1-\varepsilon)^{1/m}\},\qquad m\ge1.
\]

凹性给 (\sigma_m\ge\varepsilon)，(1-e^{-x}\le x) 给 (\sigma_m\le-\log(1-\varepsilon))。对数级数给

\[
0\le\sigma_m-\varepsilon\le\frac{\varepsilon^2}{2(1-\varepsilon)}.
\]

第五项的准确系数为

\[
d_m=\frac{(m-1)(2m-1)(3m-1)(4m-1)}{120m^4};
\]

这是 (\varepsilon^5) 项，不是 bracket 自身的五次项。main 的 (A_m) 到 (\varepsilon^4)，所以 (\varepsilon A_m) 准确到五阶，记号没有混淆。

直接对全部实 (m\ge1) 求导，

\[
\sigma_m^{(6)}(z)=\prod_{r=1}^5(r-1/m)(1-z)^{1/m-6}.
\]

当 (0\le z\le\varepsilon\le1/2) 时，每个乘积因子非负，且乘积不超过 (5!=120)；最后一因子不超过 (2^6=64)。因而第六导数在 (0\) 与 7680 之间，五阶 Taylor 余项介于 (0) 与

\[
\frac{7680}{6!}\varepsilon^6=\frac{32}{3}\varepsilon^6
\]

之间，阈值和常数独立于 (m)。这也直接支持后文使用的 (\varepsilon\le\varepsilon A_m\le\sigma_m)。(m=1) 时所有高阶系数及余项准确为零，未漏端点。

因此旧四阶截断在二维留下 (O(\varepsilon^5)) 的同阶问题已经消除：本轮误差为 (O(\varepsilon^6)=o(\varepsilon^5))。

## 5. 全部重数的统一 little-o 和自然极限

这里保留原模型的全部输入类：共享一个方体/半径，以公共交集 (F) 的正有限测度归一化；各单独 (E_i) 允许无限测度。

对所有 (m\ge1) 和 (0<\varepsilon\le1/2)，同一个区间给

\[
\varepsilon\le s=\sigma_m(\varepsilon)
\le-\log(1-\varepsilon)\le\log2<1,
\qquad s\le2\varepsilon.
\]

因此全部 (m) 同时处于精确约化的有限支路，未使用依赖 (m) 的起效区间。

固定 (n\ge2)，置 (p_8=1+8/n\le5)。线性 ordinary 余项是同一个标量函数 (r_n(t)=o(t^{p_8}))。定义

\[
\omega_n(\delta)=\sup_{0<t\le\delta}\frac{|r_n(t)|}{t^{p_8}},
\]

则由 little-o 定义，(\omega_n(\delta)\to0)，不需要连续性。因此

\[
\sup_{m\ge1}\frac{|r_n(\sigma_m(\varepsilon))|}{\varepsilon^{p_8}}
\le 2^{p_8}\omega_n(2\varepsilon)\to0.
\]

这是真正统一的余项模量，未交换逐 (m) 的极限与 supremum。

再记 (z_m=\varepsilon A_m(\varepsilon))。由前节统一界，

\[
0\le s-z_m\le(32/3)\varepsilon^6,\qquad
\varepsilon\le z_m\le s\le2\varepsilon.
\]

对九个固定指数 (p=1+j/n\in[1,5])，均值定理在同一实正区间给

\[
|s^p-z_m^p|\le p(2\varepsilon)^{p-1}(32/3)\varepsilon^6
=O_p(\varepsilon^{p+5}).
\]

固定 (n) 下有限求和为 (O_n(\varepsilon^6)=o_n(\varepsilon^{p_8}))，包括最紧的二维 (p_8=5)。因此 `eq:uniform-coefficients` 的全 (m) 一致 little-o 正确，也允许 (m=m(\varepsilon)\to\infty)。

自然参数推论更直接：主定理已经给同一个标量函数

\[
g_n(t)=\frac{C_n(1-t)-P_{7,n}(t)}{2^{n+8}t^{1+8/n}}-\widehat J_n\to0.
\]

精确约化及 (s\le2\varepsilon) 给

\[
\sup_{m\ge1}|g_n(\sigma_m(\varepsilon))|
\le\sup_{0<t\le2\varepsilon}|g_n(t)|\to0.
\]

分母 (a_sv_s^8=2^{n+8}s^{1+8/n}) 正确，未错用 (\varepsilon) 参数分母。这个推论不需要 Taylor 截断，也未把逐重数极限错误提升为一致极限。

## 6. 保留的精确约化及临界发散

本次完整重读该证明，并以实际字节比较确认它与 baseline 完全相同，仍然有效。

正向：乘积严格大于 (\lambda) 的共享见证，由 AM–GM 得总缺陷严格小于 (s_m(\lambda))；并集界给公共交集的同方体缺陷不大于总缺陷。因此多线性 halo 包含于所需线性 halo，适用单个输入无限测度的原输入类。

反向：对固定正有限测度 (F) 和线性 halo 内紧集，每点固定严格半径，在中心移动时保留 ((s-d_x)/(2m)) 余量；紧性只选择有限个原见证。在其有限并 (D\setminus F) 的 Boolean 原子上等分，能同时控制每个见证。向 (F^c) 扩张可测分割后，交集准确是 (F)。每个移动方体的各缺陷严格小于 (s/m)，于是乘积严格大于 ((1-s/m)^m=\lambda)。内正则性和对 (F) 取上确界合法；不要求构造一组同时实现整个线性 halo 的输入。

临界：周期条带的每条在任意整数长度区间占 (1/m)，删除有界正有限 (F) 后，任意中心的足够大整数边长方体都包含 (F)。各因子不小于 (1-1/m>0)，至少一个严格增大，所以即使 (\lambda=\lambda_m^*) 也属于严格 halo。公共交集为 (F)，全空间 halo 给无穷。允许无限测度的单个 (E_i) 正是这里需要的定义，main 没有在后文缩窄它。

## 7. 有限全点见证、外部上界和周期孔包络

有限见证 lemma 与外部上界组成的整块已实际字节比较，保持 baseline 原样；本次也独立重核其全点逻辑。

固定半径平均值在中心连续。紧集取得有限允许半径列表后，各点选实际可用的最大值，所得函数 Borel。按半径降序选择尚未覆盖中心，闭方体使新中心距离严格大于所有先选半径。最小允许半径正且中心有界，故半径分离 packing 强制选择次数有限；不需要先把 K 离散化。有限算法确实穷尽所有层并覆盖整个 K。

对私有 (x\in P_i)，若指定半径大于 (r_i)，在更早半径层结束时它必已被另一方体覆盖，与私有性矛盾，故所需 (r(x)\le r_i) 对所有私有点成立。这个条件保留了后续 actual witness 接口。

Besicovitch 族的 (2^n) 重叠由共同交点的中心正交象限分类直接得出；闭边界零坐标归非负，不影响严格中心分离。于是严格平均余量给原外部估计。(K\subset \mathcal H\setminus E) 不必开；它是可测差集中的紧集，Lebesgue 内正则性足够，未把“差集开”当假设。密度定理使 (E\setminus\mathcal H) 零测，恢复完整 halo。

周期孔 lower 中 (\rho<\min\{(1-\alpha)^{1/n},1/2\}) 给 (2\delta<r) 和严格平均值。(r+2\delta=1/N) 同时排除其他孔正体积侵入，(r+\delta<1/N) 保证见证留在有限外盒。输入的内点几乎处处属于 halo，先 (N\to\infty) 再 (\rho\uparrow\rho_\alpha) 得包络；主定理在 (t<2^{-n}) 时的 (L_n) 正是该分支。没有把有限构造 lower 声称为全部半径 upper。

## 8. 本次实际计算与静态核查

未重复 R057 已成功的全维 WL/Fraction 批次。本次另实际运行本机 Python3 标准库 `fractions.Fraction` 和 `math.comb` 作窄范围抄录核查，所有 `assert` 通过，输出为：

```text
b7_normalized = [-31280/243, 676/9, -21424/81,
                 -84575/243, -8836/9, -169442/81]
Jhat2 = 17638/6561
K2 = 4/27
planar_seventh = -176/243
b35 = 0
baseline_n2_j8 = 256
ordinary_b28 = 1737988096/6561
sigma_remainder_bound = 32/3
duplicate_labels = 0
external_main_refs = [eighth:theorem]
reduction_block_unchanged = true
greedy_external_block_unchanged = true
```

另对 (m=1,2,3,7,1000003) 的二至五阶参数系数作有理 falling-product 抄录比较，全部准确相等；全 (m) 的结论依据第 4 节的解析导数公式，不靠有限样本外推。ordinary (b_{2,8}) 包括 baseline 的 256，不是只有 (2^{10}\widehat J_2)。这些计算仅核标量及字面一致性，不替代集合几何。

main 自身无重复标签。`eq:quarticparameter` 已内部一致替换为 `eq:quinticparameter`；实际检索保留的 `global.tex`、`sharp.tex` 无旧标签外部引用。`thm:endpoint`、`eq:sharp`、`lem:finitegreedy`、`prop:upper`、`prop:lower` 和 `eq:external` 均保留；主文件唯一引用但不自定义的标签是模块约定的 `eighth:theorem`。已核 main 没有额外依赖被移除的 `prop:nonexact`。未编译，因此不把本次静态标签核查冒充 PDF/排版验收。

## 9. 最终签署

**绑定 SHA256 `413a02da01cccd08f33e25709b5ac579f19d092dbfd4aa965804e675c22f68ed`：PROOF-VALID。**

精确八阶主陈述与已经两审通过的 R057 全输入供货吻合，ordinary 全部低维碰撞、五阶参数及其统一六阶误差、全重数 little-o 和自然参数极限均闭合；旧精确约化、有限全点 greedy 和包络证明完整保留且独立重核。无必要数学补丁。本报告只签 main 的数学与其供货接口，不替代模块转写的独立保真审或后续实际编译。审查完成后待命，不开展新研究。

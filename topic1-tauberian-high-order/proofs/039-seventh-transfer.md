# R039：七阶主渐近假设下的系数与重数一致传递

固定维数 \(n\ge2\)。本稿只证明参数传递与系数代数，不证明也不提前确认一般维数的七阶主渐近。以下结论关于 \(C_n\) 的含义全部以假设 \(H_n\) 为前提；若主假设未成立，本稿仅作为条件材料保留。不主张维数一致性。

## 1. 主假设及条件系数公式

记
\[
 v=2t^{1/n},\quad a=v^n=2^nt,\quad
 L_n(t)=\left[1-\frac{a}{(1+v)^n}\right]^{-1},\qquad
 K_n=\binom n2\frac{n^2}{27},\quad p_*=1+\frac7n\le\frac92.
\]
**假设 \(H_n\)。** 当 \(t\downarrow0\) 时，
\[
 C_n(1-t)=L_n(t)+K_nav^6-\left(n+\frac{26}{9}\right)K_nav^7
              +o_n(av^7).                                           \tag{H}
\]

**命题1。** 在 \(H_n\) 下，
\[
 C_n(1-t)=1+\sum_{j=0}^7 b_{n,j}t^{1+j/n}+o_n(t^{p_*}),              \tag{1}
\]
其中
\[
 \boxed{b_{n,j}=2^{n+j}\left[
 \sum_{q=1}^{1+\lfloor j/n\rfloor}
 (-1)^{j-(q-1)n}\binom{n+j-1}{j-(q-1)n}
 +\mathbf1_{\{j=6\}}K_n
 -\mathbf1_{\{j=7\}}\left(n+\frac{26}{9}\right)K_n
 \right],\quad 0\le j\le7.}                                         \tag{2}
\]

证明。对充分小的 \(v\)，绝对收敛的几何级数给
\[
 L_n(t)-1=\sum_{q\ge1}v^{qn}(1+v)^{-qn}.
\]
要得到 \(v^{n+j}\)，必须有 \(q\le1+\lfloor j/n\rfloor\)，所需二项式次数为
\(r=j-(q-1)n\)。该项系数是
\[
 (-1)^r\binom{qn+r-1}{r}
 =(-1)^{j-(q-1)n}\binom{n+j-1}{j-(q-1)n}.
\]
截至 \(v^{n+7}\) 的剩余项为 \(O_n(v^{n+8})\)。将 (H) 的第六、七阶修正加入对应项，再用
\(v^{n+j}=2^{n+j}t^{1+j/n}\)，即得 (1)、(2)。这也保留了全部低维幂次碰撞。证毕。

## 2. 全部多线性重数的自然参数传递

使用已经证明的共享方体、公共交集归一化的精确约化，记
\[
 \sigma_m(\varepsilon)=m[1-(1-\varepsilon)^{1/m}],\qquad m\ge1.
\]
当 \(0<\varepsilon\le1/2\) 时，
\[
 \varepsilon\le\sigma_m(\varepsilon)\le-\log(1-\varepsilon)
 \le2\varepsilon,
 \qquad C_n^{(m)}(1-\varepsilon)=C_n(1-\sigma_m(\varepsilon)).          \tag{3}
\]
这里约化同时适用于全部 \(m\)：事实上
\(\sigma_m(\varepsilon)\le\log2<1\)，均在有限约化分支。

**命题2。** 在 \(H_n\) 下，写 \(s=\sigma_m(\varepsilon)\)、
\(a_s=2^ns\)、\(v_s=2s^{1/n}\)，则关于全部整数 \(m\ge1\) 一致有
\[
 \boxed{C_n^{(m)}(1-\varepsilon)
 =L_n(s)+K_na_sv_s^6-\left(n+\frac{26}{9}\right)K_na_sv_s^7
       +o_n(\varepsilon^{p_*}).}                                    \tag{4}
\]
等价地，可在自然参数下写成统一归一化极限
\[
 \sup_{m\ge1}\left|
 \frac{C_n^{(m)}(1-\varepsilon)-L_n(s)-K_na_sv_s^6}{a_sv_s^7}
 +\left(n+\frac{26}{9}\right)K_n\right|\longrightarrow0.              \tag{5}
\]

证明。将 (H) 的标量余项记为 \(r_n(t)=o_n(t^{p_*})\)，令
\[
 \omega_n(\delta)=\sup_{0<t\le\delta}\frac{|r_n(t)|}{t^{p_*}}
 \longrightarrow0\qquad(\delta\downarrow0).
\]
由 (3)，
\[
 \sup_{m\ge1}\frac{|r_n(\sigma_m(\varepsilon))|}{\varepsilon^{p_*}}
 \le2^{p_*}\omega_n(2\varepsilon)\longrightarrow0.
\]
这证明 (4)；用 \(a_sv_s^7=2^{n+7}s^{p_*}\) 直接归一化即得 (5)。模量只来自同一个标量余项，因此没有 \(m\) 依赖，也没有把对固定 \(m\) 的 little-o 直接交换为上确界。证毕。

## 3. 四阶参数截断已足够

定义四阶参数截断的因子
\[
 A_m(\varepsilon)=1+\frac{m-1}{2m}\varepsilon
 +\frac{(m-1)(2m-1)}{6m^2}\varepsilon^2
 +\frac{(m-1)(2m-1)(3m-1)}{24m^3}\varepsilon^3.
\]
则关于全部 \(m\ge1\) 一致有
\[
 \sigma_m(\varepsilon)=\varepsilon A_m(\varepsilon)+O(\varepsilon^5).
                                                                         \tag{6}
\]
例如直接求五阶导数，
\[
 \sigma_m^{(5)}(z)=
 \prod_{q=1}^4(q-1/m)(1-z)^{1/m-5},
\]
在 \(0\le z\le\varepsilon\le1/2\) 上绝对值至多 \(24\cdot32=768\)，所以 (6) 的余项可统一界为 \(\frac{32}{5}\varepsilon^5\)。\(m=1\) 时所有高阶项和余项均为零。

**命题3。** 在 \(H_n\) 下，关于全部整数 \(m\ge1\) 一致有
\[
 \boxed{C_n^{(m)}(1-\varepsilon)
 =1+\sum_{j=0}^7 b_{n,j}
       [\varepsilon A_m(\varepsilon)]^{1+j/n}
       +o_n(\varepsilon^{1+7/n}).}                                  \tag{7}
\]

证明。先将命题1代入 (3)，其标量余项按命题2的模量论证一致为
\(o_n(\varepsilon^{p_*})\)。两底数 \(s\) 和
\(\varepsilon A_m(\varepsilon)\) 均在 \([\varepsilon,2\varepsilon]\) 内。
对 \(1\le p\le p_*\le9/2\)，均值定理及 (6) 给
\[
 |s^p-[\varepsilon A_m(\varepsilon)]^p|
 \le C_p\varepsilon^{p+4}.
\]
对八个固定指数有限求和，误差为 \(O_n(\varepsilon^5)\)。因为
\(5-p_*\ge1/2>0\)，该误差仍是 \(o_n(\varepsilon^{p_*})\)。因此已有的四阶参数截断足够，无需计算第五阶参数系数，也无需重列逐幂组合的长表。证毕。

## 4. 新第七系数在全部维数均非零

此处是 (2) 定义的系数的独立代数事实；只有在 \(H_n\) 下才能把它称为真实 \(C_n\) 的系数。
若 \(n\ge8\)，\(\lfloor7/n\rfloor=0\)，故
\[
 2^{-(n+7)}b_{n,7}
 =-\binom{n+6}{7}-\left(n+\frac{26}{9}\right)\binom n2\frac{n^2}{27}<0.
                                                                         \tag{8}
\]
第一项是 \(L_n\) 贡献的唯一负单项，另一项也严格为负。这对所有整数
\(n\ge8\) 是直接解析证明，不由有限数值范围外推。

余下六个维数只需准确有理算术：

| \(n\) | \(L_n\) 的 \(v^{n+7}\) 系数 | \(b_{n,7}/2^{n+7}\) |
|---:|---:|---:|
| 2 | \(-128\) | \(-31280/243\) |
| 3 | \(81\) | \(676/9\) |
| 4 | \(-240\) | \(-21424/81\) |
| 5 | \(-275\) | \(-84575/243\) |
| 6 | \(-804\) | \(-8836/9\) |
| 7 | \(-1715\) | \(-169442/81\) |

这些数全部非零。特别 \(b_{2,7}=-16015360/243\)，与已证二维七阶公式的基线展开及修正相符；\(n=3\) 的第七系数为正，故不能将“全维非零”加强成“全维为负”。结合 (8)，得到 \(b_{n,7}\ne0\) 对每个整数 \(n\ge2\) 成立。

## 5. 实际计算证据及其范围

已实际读取 `r039-transfer-computation/coefficients.wl`、日志及 JSON：Mathematica 15.0 的 15 个子核分别执行 \(n=2,\ldots,16\)，记录了 15 个不同 kernel ID 和 PID；逐项 `SeriesCoefficient` 与 (2) 的基线有限和全部一致，二维回归为 \(-16015360/243\)。这些只核系数代数，绝不验证主假设 (H)。

原 `coefficients.log` 的空零维集合经过 `Lookup` 后输出
`Missing[KeyAbsent, n]`，不能据此将那一行报告成成功返回空列表。已有 `summary.py` 的实际 `summary.log` 从 JSON 的逐项布尔记录独立复核，明确返回 `zero dimensions: []`。本目录的 `verify_fractions.py` 又用 Python 标准库 `Fraction` 和 `comb` 实际逐项核验了全部15个有理结果、基线和零标记，输出全匹配与空零列表；上表就是其中 \(n=2,\ldots,7\) 的必要有限证据。原日志的 `Reduce` 正性违反集合返回 `False` 表示该违反集合为空；(8) 已另给一般维数解析证明。

本稿只读取旧源及计算缓存，未改仓库；不依赖一般维数局部价格或总装的审核结果去取消 (H)。

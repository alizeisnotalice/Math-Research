# upper_shift091：五方向链提升给完整正 band 的新平移比

本轮独立证明，对固定 body16、共同时间 Cauchy 模型，任意 x∈R5、t>0、面内位移 u1=0，令 r=||u||∞，则
\[
\boxed{R_{\rm new}(r/t)^{-1}\le {p_t(x+u)\over p_t(x)}
\le R_{\rm new}(r/t),}\tag{1}
\]
\[
R_{\rm new}(z)=\exp\{2\operatorname{asinh}(z/4)+8\operatorname{asinh}(z/16)\}
=\left(\sqrt{1+z^2/16}+z/4\right)^2
 \left(\sqrt{1+z^2/256}+z/16\right)^8.
\tag{2}
\]
这是一条全空间、全部正时间的正密度比值定理，不依赖锥节点非退化、梯度求值、p 的凸性或只在盒角点成立的比值。

当 z>0 时，\(R_{\rm new}(z)<e^z\)。对同一真实盒、同一 anchor、同一时间下端，若旧 Q-asinh 参数取合法的 \(B_Q^2\ge\max_Bu^TQu\)、Q=I4/4−J4/36，新公式还严格小于该旧公式。下文给出完整解析比较；不把不同盒、不同时间或不同 anchor 的旧数值当作比较输入。

六个预登记比例及严格认证的显示值如下：

| z=h∞/a | Rnew 显示近似 | log Rnew 显示近似 | Rnew/已保存 exp(z) 球的比值上端（显示） |
|---|---:|---:|---:|
| 1/8 | 1.1331362107968384 | 0.124989196176294 | 0.999989196235 |
| 1 | 2.703664224496367 | 0.994607973138849 | 0.994622484023 |
| 2 | 7.098174408991937 | 1.959837625490749 | 0.960633444104 |
| 4 | 42.20247490740483 | 3.742478866417194 | 0.772965290615 |
| 8 | 842.9988137587104 | 6.736965550834448 | 0.282794597383 |
| 16 | 76146.45395641754 | 11.24041379067855 | 未重算 exp(16)；已证 log Rnew<16 |

z=8 时另准确验证 Rnew<843。完整 outward 端点保存在 JSON，不使用显示小数作为证书。六项都严格低于同 h/time 的旧 asinh 保存下端；前五项也严格低于已有 exp 保存下端。仅新算六个代数 R 及六个 log 资格表达式，旧 exp/asinh 无重算，零 H、stencil、梯度或 Green 调用。没有新 I、g、C5 或 Cn 结论。

## 1. 单 Cauchy 比值与固定纤维

保持未单位化方向 V={(1,±1,±1,±1,±1)}、各时间权1、α=1、Lebesgue 与逆 Fourier 因子 (2π)^−5。写 B 为方向矩阵、K=ker B、A=Bᵀ/16、\(f_t(z)=\prod_{j=1}^{16}s_t(z_j)\)。075 已证明处处正、有限的固定纤维表示
\[
p_t(x)=\frac1{1024}\int_Kf_t(Ax+w)dw.\tag{3}
\]
任意固定提升 a 满足 Ba=u，则 a−Au∈K。对 w 作固定平移，K 上 Hausdorff 测度不变，故
\[
p_t(x+u)=\frac1{1024}\int_Kf_t(Ax+w+a)dw.\tag{4}
\]
比例的分母 p_t(x)>0。这里不改变纤维 Jacobian、不取 Cauchy 变量矩，也不要求对 p 求导。

对单坐标，067 的精确比值可直接重证。设
\[
\rho_t(a)=\left(\sqrt{1+a^2/(4t^2)}+|a|/(2t)\right)^2.
\]
a=0 时比值为1；a≠0 时 ρ>1 且 \((\rho-1)^2=\rho a^2/t^2\)。二次式
\(\rho[t^2+(z+a)^2]-(t^2+z^2)\) 的首项为ρ−1>0，判别式为零，因此处处非负。反向位移给另一侧，于是
\[
\rho_t(a)^{-1}\le s_t(z+a)/s_t(z)\le\rho_t(a)
=e^{2\operatorname{asinh}(|a|/(2t))}.
\tag{5}
\]
对16个因子相乘，再作 (3)–(4) 的正积分，保留同一双侧常数。该步对任意合法固定 lift 成立，不需要 lift 为 Euclidean 最小范数右逆。

## 2. 五方向链提升及所有退化情形

将尾四坐标按绝对值递减排序，记索引排列 π，
\(w_k=|u_{\pi(k)}|\)，\(w_1=r\ge w_2\ge w_3\ge w_4\ge w_5=0\)。同值按原坐标索引排序；零坐标的 sign 固定为+1。记各尾符号为σj。

令 v0=(1,−σ1,−σ2,−σ3,−σ4)，vk 依次翻转排序中的前 k 个尾坐标，k=1,…,4。这五列都属于原 V。定义
\[
a_0=-r/2,\qquad a_k=(w_k-w_{k+1})/2\quad(1\le k\le4).
\tag{6}
\]
其它11个方向的系数为零。差值 \(d_k=w_k-w_{k+1}\ge0\) 且 Σdk=r，所以 Σak=0、Σ|ak|=r。

对排在第 j 位的尾坐标，其合成值为
\[
\sigma_{\pi(j)}\left[-a_0-\sum_{k<j}a_k+\sum_{k\ge j}a_k\right]
=\sigma_{\pi(j)}\left[r-\sum_{k<j}d_k\right]
=\sigma_{\pi(j)}w_j=u_{\pi(j)}.
\]
首坐标为Σak=0，故 Ba=u。相等绝对值只产生零系数；零符号选择不改变合成值。u=0 时五个系数全零，比值准确为1，不除以 r。lift 只依赖本次固定位移 u，和积分纤维变量 x、w、t 无关。

由 (5)，这个具体 lift 的双侧指数是
\[
2\operatorname{asinh}(r/(4t))
+2\sum_{k=1}^4\operatorname{asinh}(d_k/(4t)).\tag{7}
\]
asinh 在非负半轴递增且凹，四个非负 dk 总和为r，Jensen 恰给
\(\sum_k\operatorname{asinh}(d_k/(4t))\le4\operatorname{asinh}(r/(16t))\)。这证明 (1)–(2)。用于 Jensen 的变量是四个非负 gap，不是空间盒顶点，也不是 p_t 本身。比如绝对值 (r,3r/4,r/2,r/4) 给四等 gap，达到此一步的等号；这不表示真实混合密度的比值上界可达到。

## 3. 与 exp 及匹配的旧 Q-asinh 的全参数比较

Rnew(z) 在 z≥0 单调递增，Rnew(0)=1。因 asinh(s)<s 对 s>0，
\[
\log R_{\rm new}(z)<2z/4+8z/16=z\quad(z>0).
\tag{8}
\]
所以它严格改进075的 exp 位移公式，包含任意大小的正比例，不仅六个有限资格点。

再把一个 z/4 和四个 z/16 视为五个非负数，凹性给
\[
\log R_{\rm new}(z)
\le10\operatorname{asinh}(z/10).\tag{9}
\]
对同一真实盒及 anchor，令 h∞=max||u||∞，并令旧 BQ² 为 \(\max u^TQu\) 的合法上端。085 的准确平方完成已证明
\(u^TQu\ge(5/24)||u||∞^2\)，因此
\[
B_Q^2\ge(5/24)h_\infty^2.
\]
对 z=h∞/a>0，原五基 asinh 公式满足
\[
10\operatorname{asinh}\frac{\sqrt{B_Q^2}}{2\sqrt5a}
\ge10\operatorname{asinh}\frac{z}{4\sqrt6}
>10\operatorname{asinh}\frac{z}{10}
\ge\log R_{\rm new}(z),\tag{10}
\]
因为 \((4\sqrt6)^2=96<100\)。这是严格匹配 h∞、BQ² 和 a 后的解析比较，不依赖重新求旧 R。若旧证书用了更保守的 BQ² 上端，结论仍成立。

式 (10) 只比较指定的旧 Q-asinh 公式，不声称优于任何可能的其它合法比值证书。数值端点可能受舍入影响，因此 consumer 仍应取新向上端、旧向上端及其它合法 cap 的 min。旧 anchor、盒或时间下端不匹配时，旧数值不能用于 (10) 或新的 band。

## 4. 真正的盒、任意 anchor 和完整正时间 band

对固定 anchor c 和同首坐标的真实四尾盒 B，h∞ 按全部真实端点相对 c 计算。无论 c 是否居中、是否在 B 内，对于任意 u=x−c、x∈B，以及 s∈[0,1]，都有 ||su||∞≤h∞。因此 (1) 对整个盒及所有 anchor-to-box 线段统一成立；不需要以盒角点的 p 值推断内部。

固定 J=[a,b]，a>0、b≤∞。因为 t≥a 时 \(R_{\rm new}(h∞/t)\le R_{\rm new}(h∞/a)\)，正时间积分给
\[
\boxed{\frac{m_J}{R_{\rm new}(h∞/a)}
\le B_J(c+u)\le R_{\rm new}(h∞/a)m_J,}
\quad B_J(x)=\int_Jp_t(x)dt,\ m_J=B_J(c).
\tag{11}
\]
同一式对 su 成立。无穷尾由075的 \(p_t\le1/(16\pi^5t^5)\) 可积峰值保证有限；没有删时间尾，也没有改变 a 或原 ε。若盒退化为 anchor 本身，则 h∞=0，R准确为1；若只是部分坐标零宽度，h∞ 仍按实际端点定义，定理不变。

将 (11) 用在既有正 band，保留其原质量外端 \([m_L,m_U]\)：range 为 \(m_L/R_U\) 到 \(m_UR_U\)。085 的均值公式也可用同一个 R_U：
\[
\operatorname{mean}B_J\le m_U\min\left\{R_U,
1+\|\mu\|∞/a+(25/48)M_2R_U/a^2,\ C_J^{\rm existing}\right\}.
\tag{12}
\]
随后仍按原时间划分相加、使用合法全盒range与均值作 \(z^{5/4}\) 的凸弦。本轮不生产 band、不改变中心、不改 ε、不重算 H；(11)–(12) 是供后续另行登记消费的接口，不是新的全层积分证书。

## 5. 纯整数向上 API 与有限严格资格

`ratio.py` 只依赖标准库。核心 `ratio_interval(z,bits=192)` 对 (2) 的两个平方根使用整数 isqrt 及准确 Fraction 输入，保留向下/向上根，再作正整数幂，最后向外量化到分母2^bits。每个根都检查 \(l^2\le q<u^2\)。z=0 特判返回(1,1)。没有 exp、log、asinh 或 Arb 依赖。

公开调用：
```python
from ratio import new_ratio_upper, combine_ratio_upper
R_U = new_ratio_upper(h_infinity, time_lower, bits=192)
R_U = combine_ratio_upper(h_infinity, time_lower,
                          existing_upper_bounds=(old_R_U,), bits=192)
```
输入使用 Fraction、准确有理字符串或整数，拒绝浮点输入。time_lower 必须正，h∞ 必须非负；该接口固定用于 u1=0 的同首坐标盒。内部缓存最多4096个 (z,bits) 结果，复用常见盒/时间比例；缓存不改变准确算术和模型。`chain_lift(tail)` 可返回明确的五列和系数，供检查或更细方向比值研究。

有限计算预登记12个有理尾向量，覆盖零、并列绝对值、零符号、正负混合、单坐标、四等 gap 和实际079保存几何偏移。准确核验 Ba=(0,u)、Σa=0、l1=r、非负 gap 与确定排序。

标量批只计算 z∈{1/8,1,2,4,8,16} 的六个新代数 R 和六个 log(R) 资格表达式，384-bit Arb。每个192-bit整数 API 外区间完整包含其 Arb 球，且保存的 log 上端严格小于 z。没有为比较重算 exp 或旧 asinh：

- 前五个比例读取074 h=1/64的已存同一裁片，band 下端依次1/8、1/64、1/128、1/256、1/512；其真实 h∞ 恰为1/64，旧 exp/旧 asinh完整外端均保留。
- z=16读取072 h=1/128、band下端1/2048的旧 asinh外端；该真实裁片 h∞=1/128，保存的最大Q值为1/16384。没有保存 exp(16) 数值输入，因此只用本轮 log 上端<16核证与 exp(16) 的比较，不补算旧值。

注册中的旧结果字段与源文件逐项相同，h∞ 和 h/a 的匹配在实际程序中准确检查。六个新整数向上端全严格低于旧 asinh 下端，前五全严格低于已存 exp 下端。来源身份以 SHA 锁定；其数学有效性依赖相应旧报告和本报告的匹配检查，不靠 hash 本身代替证明。

独立 `verify.py` 不导入生产 ratio 函数、不调用 Arb，用不同缩放公式
\[
R_{\rm new}(z)=\frac{(\sqrt{16+z^2}+z)^2\cdot(\sqrt{256+z^2}+z)^8}{16^9}
\]
进行224-bit整数二分平方根包围。该路径的区间包含保存 Arb 球，并包含于192-bit API区间；所有旧端点比较相同。它另独立重构12个 lift，不重新求任何 log 或旧比值。

## 6. 一般 L1 模型接口与运行收据

比值只依赖 r/t，符合 α=1 共同尺度缩放；原 body方向、密度、时间权、正源和 Fourier归一化全部不变。029 的一般 L1、全空间有限正源、所有正尺度与共同半径中心cube传递链保持原样。后续若使用它改善合法的 I 上端，仍需独立合法的 essential cube g 才可代入 \(C_5\le5^{5/4}I/g^{5/4}\)。本批不改变这些资格，不把局部比例当作全域g，也不把正部质量的128因子移入I或 band 槽位。

生产 PID76850、session36437：先登记、armed并发句柄后运行；六公式/六log与有限lift的计算耗时0.019804625秒，实际exit0。独立整数验证PID77098，0.004171333秒，exit0。按当前用户规则没有wallcap/alarm，有限前缀逐项保存；本轮无失败或恢复。没有H、stencil、梯度、Green、Gamma、旧exp或旧asinh重算，所有程序已结束。

登记、源码、完整球、旧端点副本、独立有理复核与本报告的hash见manifest。旧085/076及其它来源保持原字节，仅新增本目录，未提交或推送、未使用MATLAB。再次只读复核入口为 `python -B verify.py`；生产器有已完成checkpoint时拒绝重放。真实全层旧band再包络应另行登记并验证其范围/均值/凸弦账本。
